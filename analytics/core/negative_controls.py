"""Matched null controls for GAPH SNV candidates.

Two controls answer different questions:

* matched callable positions test whether GAPH candidates differ from the
  technically observable background;
* unobserved alternate alleles at the same position test whether the exact ALT
  carries information beyond the site itself.

The implementation samples a bounded, deterministic candidate set per
strategy.  It never materializes every possible allele in a target gene.
"""

from __future__ import annotations

import bisect
import concurrent.futures
import csv
import gzip
import hashlib
import heapq
import json
import math
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from bin.fetch_gnomad_variants import fetch_region_variants_recursive

from .clinvar_validation import clinvar_label, info_value, path_metadata, split_strategies
from .conservation import annotate_track, parse_tracks


DNA_BASES = ("A", "C", "G", "T")
CONTEXT_PRIORITY = ("cds", "utr", "exon", "intron")
INVALID_SEGMENT_FLAGS = frozenset({"low_mapq", "non_primary", "ambiguous_event_allele"})
CONTROL_VERSION = 1
CALLABLE_BLOCK_VERSION = 1
MATCHED_POOL_SIZE = 5


@dataclass(frozen=True)
class NegativeControlAnalysis:
    matched_summary: pd.DataFrame
    matched_context_summary: pd.DataFrame
    matched_ecdf: pd.DataFrame
    same_site_summary: pd.DataFrame
    manifest: dict
    manifest_path: Path
    matched_path: Path
    same_site_path: Path
    conservation_path: Path
    permutations: int


def build_negative_controls(
    *,
    run_dir: Path,
    variant_annotations_tsv: Path,
    alignment_segments_tsv: Path,
    target_features_tsv: Path,
    genes_tsv: Path,
    target_sequences_dir: Path,
    clinvar_vcf: Path,
    clinvar_regions_bed: Path,
    strategies: list[str],
    sample_size_per_strategy: int = 25_000,
    permutations: int = 1_000,
    seed: int = 20_260_721,
) -> NegativeControlAnalysis:
    """Build or load both negative-control analyses for one completed run."""

    if sample_size_per_strategy < 1:
        raise ValueError("negative-control sample size must be >= 1")
    if permutations < 100:
        raise ValueError("negative-control permutations must be >= 100")

    outdir = run_dir / "analytics" / "negative_control"
    outdir.mkdir(parents=True, exist_ok=True)
    matched_path = outdir / "matched_callable.snv.tsv.gz"
    same_site_path = outdir / "same_position_alt.snv.tsv.gz"
    conservation_path = outdir / "matched_callable.phyloP100way.tsv.gz"
    manifest_path = outdir / "manifest.json"
    callable_blocks_path = outdir / "callable_blocks.tsv.gz"
    callable_manifest_path = outdir / "callable_blocks.manifest.json"

    expected_inputs = {
        "version": CONTROL_VERSION,
        "variant_annotations": path_metadata(variant_annotations_tsv),
        "alignment_segments": path_metadata(alignment_segments_tsv),
        "target_features": path_metadata(target_features_tsv),
        "genes": path_metadata(genes_tsv),
        "clinvar_vcf": path_metadata(clinvar_vcf),
        "clinvar_tbi": path_metadata(Path(f"{clinvar_vcf}.tbi")),
        "strategies": sorted(strategies),
        "sample_size_per_strategy": sample_size_per_strategy,
        "matched_pool_size": MATCHED_POOL_SIZE,
        "seed": seed,
        "conservation_track": "phyloP100way",
    }
    if _cache_is_valid(manifest_path, expected_inputs, [matched_path, same_site_path, conservation_path]):
        manifest = json.loads(manifest_path.read_text())
        return _load_analysis(
            matched_path,
            same_site_path,
            conservation_path,
            manifest,
            manifest_path,
            permutations,
            seed,
        )

    genes = _read_genes(genes_tsv)
    contexts = _read_disjoint_contexts(target_features_tsv, genes)
    focal = _sample_focal_snvs(
        variant_annotations_tsv,
        contexts,
        strategies,
        sample_size_per_strategy,
        seed,
    )
    if focal.empty:
        raise ValueError("No normalized GAPH SNVs were available for negative controls.")

    _build_or_load_callable_blocks(
        alignment_segments_tsv,
        callable_blocks_path,
        callable_manifest_path,
    )
    sequences = _read_target_sequences(target_sequences_dir, set(focal["gene_id"]))
    tentative_matched, same_site = _generate_control_options(
        focal=focal,
        callable_blocks_path=callable_blocks_path,
        contexts=contexts,
        genes=genes,
        sequences=sequences,
        seed=seed,
    )
    observed_controls = _collect_observed_control_keys(
        variant_annotations_tsv,
        tentative_matched,
        same_site,
    )
    matched = _finalize_matched_options(tentative_matched, observed_controls)
    same_site = _finalize_same_site_options(same_site, observed_controls)

    clinvar = _read_clinvar_snv_annotations(clinvar_vcf, clinvar_regions_bed)
    gnomad_keys, gnomad_complete, gnomad_error = _fetch_gnomad_presence(same_site)
    same_site = _annotate_same_site_controls(same_site, clinvar, gnomad_keys, gnomad_complete)
    _write_tsv(same_site_path, same_site)

    conservation_rows, conservation_manifest = _annotate_matched_conservation(
        matched,
        conservation_path,
    )
    matched = matched.merge(conservation_rows, on="variant_key", how="left", validate="many_to_one")
    _write_tsv(matched_path, matched)

    manifest = {
        "inputs": expected_inputs,
        "complete": bool(gnomad_complete and conservation_manifest.get("status") == "complete"),
        "focal_candidate_count": int(len(focal)),
        "matched_row_count": int(len(matched)),
        "matched_focal_count": int(matched.loc[matched["role"] == "observed", "focal_id"].nunique()),
        "same_site_row_count": int(len(same_site)),
        "same_site_focal_count": int(same_site.loc[same_site["role"] == "observed", "focal_id"].nunique()),
        "gnomad_complete": gnomad_complete,
        "gnomad_error": gnomad_error,
        "conservation": conservation_manifest,
        "callable_blocks": str(callable_blocks_path),
        "matched_tsv": str(matched_path),
        "same_site_tsv": str(same_site_path),
        "conservation_tsv": str(conservation_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return _summarize_analysis(
        matched,
        same_site,
        manifest,
        manifest_path,
        matched_path,
        same_site_path,
        conservation_path,
        permutations,
        seed,
    )


def _cache_is_valid(manifest_path: Path, expected_inputs: dict, outputs: list[Path]) -> bool:
    if not manifest_path.exists() or not all(path.exists() for path in outputs):
        return False
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    if manifest.get("complete") is False:
        return False
    return manifest.get("inputs") == expected_inputs


def _load_analysis(
    matched_path: Path,
    same_site_path: Path,
    conservation_path: Path,
    manifest: dict,
    manifest_path: Path,
    permutations: int,
    seed: int,
) -> NegativeControlAnalysis:
    matched = pd.read_csv(matched_path, sep="\t", compression="gzip", keep_default_na=False)
    matched["phyloP100way"] = pd.to_numeric(matched["phyloP100way"], errors="coerce")
    same_site = pd.read_csv(same_site_path, sep="\t", compression="gzip", keep_default_na=False)
    return _summarize_analysis(
        matched,
        same_site,
        manifest,
        manifest_path,
        matched_path,
        same_site_path,
        conservation_path,
        permutations,
        seed,
    )


def _write_tsv(path: Path, frame: pd.DataFrame) -> None:
    frame.to_csv(path, sep="\t", index=False, compression="gzip", lineterminator="\n")


def _read_genes(path: Path) -> dict[str, dict[str, object]]:
    frame = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
    required = {"gene_id", "chromosome", "begin", "end", "sequence_length"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Genes table missing columns: {', '.join(sorted(missing))}")
    genes = {}
    for row in frame.itertuples(index=False):
        genes[str(row.gene_id)] = {
            "chrom": str(row.chromosome).removeprefix("chr"),
            "begin": int(row.begin),
            "end": int(row.end),
            "length": int(row.sequence_length),
        }
    return genes


def _read_disjoint_contexts(
    path: Path,
    genes: dict[str, dict[str, object]],
) -> dict[str, list[tuple[int, int, str]]]:
    frame = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        keep_default_na=False,
        usecols=["gene_id", "feature_type", "target_start0", "target_end0"],
    )
    by_gene: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    for row in frame.itertuples(index=False):
        feature = str(row.feature_type).lower()
        if feature in CONTEXT_PRIORITY:
            by_gene[str(row.gene_id)][feature].append((int(row.target_start0), int(row.target_end0)))

    result = {}
    for gene_id, gene in genes.items():
        feature_intervals = by_gene.get(gene_id, {})
        boundaries = {0, int(gene["length"])}
        for intervals in feature_intervals.values():
            for start, end in intervals:
                boundaries.add(max(0, start))
                boundaries.add(min(int(gene["length"]), end))
        ordered = sorted(boundaries)
        disjoint = []
        for start, end in zip(ordered, ordered[1:]):
            if end <= start:
                continue
            context = "other"
            for candidate in CONTEXT_PRIORITY:
                if any(left < end and right > start for left, right in feature_intervals.get(candidate, [])):
                    context = "other_exon" if candidate == "exon" else candidate
                    break
            if disjoint and disjoint[-1][2] == context and disjoint[-1][1] == start:
                disjoint[-1] = (disjoint[-1][0], end, context)
            else:
                disjoint.append((start, end, context))
        result[gene_id] = disjoint
    return result


def _context_at(intervals: list[tuple[int, int, str]], position: int) -> str:
    starts = [item[0] for item in intervals]
    index = bisect.bisect_right(starts, position) - 1
    if index >= 0 and intervals[index][0] <= position < intervals[index][1]:
        return intervals[index][2]
    return "other"


def _stable_rank(seed: int, *parts: object) -> int:
    payload = "|".join([str(seed), *(str(part) for part in parts)]).encode()
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "big")


def _sample_focal_snvs(
    path: Path,
    contexts: dict[str, list[tuple[int, int, str]]],
    strategies: list[str],
    limit: int,
    seed: int,
) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", compression="gzip", nrows=0).columns.tolist()
    columns = [
        "variant_key",
        "gene_id",
        "event_type",
        "target_start0",
        "genomic_start1",
        "ref",
        "alt",
        "strategies",
    ]
    if "lookup_status" in header:
        columns.append("lookup_status")
    strategy_set = set(strategies)
    heaps: dict[str, list[tuple[int, str, dict[str, object]]]] = defaultdict(list)

    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=columns,
        keep_default_na=False,
        chunksize=200_000,
    ):
        chunk = chunk[chunk["event_type"].astype(str).eq("snv")]
        chunk = chunk[
            chunk["ref"].astype(str).str.len().eq(1)
            & chunk["alt"].astype(str).str.len().eq(1)
            & chunk["ref"].astype(str).str.upper().isin(DNA_BASES)
            & chunk["alt"].astype(str).str.upper().isin(DNA_BASES)
        ]
        if "lookup_status" in chunk.columns:
            chunk = chunk[chunk["lookup_status"].astype(str).eq("ok")]
        for row in chunk.itertuples(index=False):
            gene_id = str(row.gene_id)
            target_pos = int(row.target_start0)
            record_base = {
                "gene_id": gene_id,
                "variant_key": str(row.variant_key),
                "target_pos": target_pos,
                "genomic_pos": int(row.genomic_start1),
                "ref": str(row.ref).upper(),
                "alt": str(row.alt).upper(),
                "context": _context_at(contexts.get(gene_id, []), target_pos),
            }
            for strategy in split_strategies(str(row.strategies)):
                if strategy_set and strategy not in strategy_set:
                    continue
                record = {**record_base, "strategy": strategy}
                token = f"{gene_id}:{record['variant_key']}"
                rank = _stable_rank(seed, strategy, token)
                heap = heaps[strategy]
                item = (-rank, token, record)
                if len(heap) < limit:
                    heapq.heappush(heap, item)
                elif rank < -heap[0][0]:
                    heapq.heapreplace(heap, item)

    rows = []
    for strategy, heap in heaps.items():
        for _negative_rank, _token, record in sorted(heap, key=lambda item: (-item[0], item[1])):
            rows.append(record)
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame = frame.sort_values(["strategy", "gene_id", "variant_key"], kind="mergesort").reset_index(drop=True)
    frame.insert(0, "focal_id", [f"focal_{index:09d}" for index in range(len(frame))])
    return frame


def _read_target_sequences(directory: Path, gene_ids: set[str]) -> dict[str, str]:
    sequences = {}
    for gene_id in sorted(gene_ids):
        path = directory / f"{gene_id}.fa.gz"
        if not path.exists():
            raise FileNotFoundError(f"Missing target sequence for gene {gene_id}: {path}")
        parts = []
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if not line.startswith(">"):
                    parts.append(line.strip())
        sequences[gene_id] = "".join(parts).upper()
    return sequences


def _build_or_load_callable_blocks(
    segments_path: Path,
    blocks_path: Path,
    manifest_path: Path,
) -> None:
    expected = {
        "version": CALLABLE_BLOCK_VERSION,
        "alignment_segments": path_metadata(segments_path),
        "invalid_flags": sorted(INVALID_SEGMENT_FLAGS),
        "support_unit": "tax_id, else taxname, else ortholog_gene_id",
    }
    if _cache_is_valid(manifest_path, expected, [blocks_path]):
        return

    with tempfile.TemporaryDirectory(prefix="callable_blocks_", dir=blocks_path.parent) as tmp_name:
        sorted_path = Path(tmp_name) / "segments.sorted.tsv"
        env = dict(os.environ)
        env["LC_ALL"] = "C"
        with sorted_path.open("w") as sorted_handle:
            process = subprocess.Popen(
                ["sort", "-t", "\t", "-k1,1", "-k2,2", "-k3,3", "-k4,4n", "-k5,5n"],
                stdin=subprocess.PIPE,
                stdout=sorted_handle,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
            )
            assert process.stdin is not None
            with gzip.open(segments_path, "rt", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                required = {
                    "gene_id",
                    "strategy",
                    "ortholog_gene_id",
                    "tax_id",
                    "taxname",
                    "target_start0",
                    "target_end0",
                    "is_primary",
                    "qc_flags",
                }
                missing = required - set(reader.fieldnames or [])
                if missing:
                    process.stdin.close()
                    process.kill()
                    raise ValueError(f"Alignment segments missing columns: {', '.join(sorted(missing))}")
                for row in reader:
                    if str(row["is_primary"]).lower() != "true":
                        continue
                    flags = {item for item in str(row.get("qc_flags") or "").split(",") if item}
                    if flags & INVALID_SEGMENT_FLAGS:
                        continue
                    start = int(row["target_start0"])
                    end = int(row["target_end0"])
                    if end <= start:
                        continue
                    unit = row.get("tax_id") or row.get("taxname") or row.get("ortholog_gene_id")
                    if not unit:
                        continue
                    process.stdin.write(
                        f"{row['strategy']}\t{row['gene_id']}\t{unit}\t{start}\t{end}\n"
                    )
            process.stdin.close()
            stderr = process.stderr.read() if process.stderr is not None else ""
            return_code = process.wait()
            if return_code != 0:
                raise RuntimeError(f"sort failed while building callable blocks: {stderr.strip()}")

        block_count = _collapse_sorted_segments(sorted_path, blocks_path)

    manifest = {
        "inputs": expected,
        "block_count": block_count,
        "blocks_tsv": str(blocks_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def _collapse_sorted_segments(sorted_path: Path, blocks_path: Path) -> int:
    fields = ["strategy", "gene_id", "target_start0", "target_end0", "callable_species"]
    block_count = 0
    current_group: tuple[str, str] | None = None
    current_unit: tuple[str, str, str] | None = None
    merged_start = 0
    merged_end = 0
    events: dict[int, int] = defaultdict(int)

    def add_merged_interval() -> None:
        if current_unit is not None and merged_end > merged_start:
            events[merged_start] += 1
            events[merged_end] -= 1

    def write_group(writer: csv.DictWriter) -> int:
        if current_group is None or not events:
            return 0
        written = 0
        depth = 0
        previous = None
        strategy, gene_id = current_group
        for position in sorted(events):
            if previous is not None and position > previous and depth > 0:
                writer.writerow(
                    {
                        "strategy": strategy,
                        "gene_id": gene_id,
                        "target_start0": previous,
                        "target_end0": position,
                        "callable_species": depth,
                    }
                )
                written += 1
            depth += events[position]
            previous = position
        return written

    with sorted_path.open() as source, gzip.open(blocks_path, "wt", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line in source:
            strategy, gene_id, unit, start_text, end_text = line.rstrip("\n").split("\t")
            group = (strategy, gene_id)
            unit_key = (strategy, gene_id, unit)
            start = int(start_text)
            end = int(end_text)
            if current_group is not None and group != current_group:
                add_merged_interval()
                block_count += write_group(writer)
                events.clear()
                current_unit = None
            if current_unit != unit_key:
                if current_unit is not None:
                    add_merged_interval()
                current_unit = unit_key
                merged_start, merged_end = start, end
            elif start <= merged_end:
                merged_end = max(merged_end, end)
            else:
                add_merged_interval()
                merged_start, merged_end = start, end
            current_group = group
        add_merged_interval()
        block_count += write_group(writer)
    return block_count


def _depth_bin(depth: int) -> str:
    if depth < 5:
        return "1-4"
    if depth < 10:
        return "5-9"
    if depth < 20:
        return "10-19"
    if depth < 50:
        return "20-49"
    return "50+"


def _intersect_blocks_with_contexts(
    blocks: list[tuple[int, int, int]],
    contexts: list[tuple[int, int, str]],
) -> tuple[
    list[tuple[int, int, int]],
    dict[tuple[str, str], list[tuple[int, int]]],
]:
    searchable = []
    strata: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    context_index = 0
    for start, end, depth in blocks:
        searchable.append((start, end, depth))
        while context_index < len(contexts) and contexts[context_index][1] <= start:
            context_index += 1
        index = context_index
        while index < len(contexts) and contexts[index][0] < end:
            context_start, context_end, context = contexts[index]
            overlap_start = max(start, context_start)
            overlap_end = min(end, context_end)
            if overlap_end > overlap_start:
                strata[(context, _depth_bin(depth))].append((overlap_start, overlap_end))
            index += 1
    return searchable, strata


def _depth_at(blocks: list[tuple[int, int, int]], position: int) -> int:
    starts = [item[0] for item in blocks]
    index = bisect.bisect_right(starts, position) - 1
    if index >= 0 and blocks[index][0] <= position < blocks[index][1]:
        return blocks[index][2]
    return 0


def _weighted_position(
    intervals: list[tuple[int, int]],
    cumulative: list[int],
    rng: np.random.Generator,
) -> int:
    offset = int(rng.integers(0, cumulative[-1]))
    index = bisect.bisect_right(cumulative, offset)
    previous = cumulative[index - 1] if index else 0
    start, _end = intervals[index]
    return start + (offset - previous)


def _generate_control_options(
    *,
    focal: pd.DataFrame,
    callable_blocks_path: Path,
    contexts: dict[str, list[tuple[int, int, str]]],
    genes: dict[str, dict[str, object]],
    sequences: dict[str, str],
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    focal_groups = {
        (strategy, gene_id): group.copy()
        for (strategy, gene_id), group in focal.groupby(["strategy", "gene_id"], sort=False)
    }
    matched_rows = []
    same_site_rows = []

    for row in focal.itertuples(index=False):
        common = {
            "focal_id": row.focal_id,
            "strategy": row.strategy,
            "gene_id": row.gene_id,
            "context": row.context,
            "role": "observed",
            "option": 0,
            "variant_key": row.variant_key,
            "chrom": str(genes[row.gene_id]["chrom"]),
            "pos": int(row.genomic_pos),
            "ref": row.ref,
            "alt": row.alt,
        }
        same_site_rows.append(common)
        for option, alt in enumerate((base for base in DNA_BASES if base not in {row.ref, row.alt}), start=1):
            same_site_rows.append(
                {
                    **common,
                    "role": "control",
                    "option": option,
                    "alt": alt,
                    "variant_key": f"{common['chrom']}:{common['pos']}:{row.ref}>{alt}",
                }
            )

    current_key = None
    current_blocks: list[tuple[int, int, int]] = []

    def process_group(key: tuple[str, str] | None, blocks: list[tuple[int, int, int]]) -> None:
        if key is None or key not in focal_groups or not blocks:
            return
        strategy, gene_id = key
        searchable, strata = _intersect_blocks_with_contexts(blocks, contexts.get(gene_id, []))
        gene = genes[gene_id]
        sequence = sequences[gene_id]
        for row in focal_groups[key].itertuples(index=False):
            depth = _depth_at(searchable, int(row.target_pos))
            if depth <= 0:
                continue
            depth_group = _depth_bin(depth)
            common = {
                "focal_id": row.focal_id,
                "strategy": strategy,
                "gene_id": gene_id,
                "context": row.context,
                "depth_bin": depth_group,
                "role": "observed",
                "option": 0,
                "variant_key": row.variant_key,
                "chrom": str(gene["chrom"]),
                "pos": int(row.genomic_pos),
                "target_pos": int(row.target_pos),
                "ref": row.ref,
                "alt": row.alt,
                "callable_species": depth,
            }
            matched_rows.append(common)
            intervals = strata.get((row.context, depth_group), [])
            if not intervals:
                continue
            cumulative = np.cumsum([end - start for start, end in intervals]).astype(int).tolist()
            rng = np.random.default_rng(_stable_rank(seed, strategy, row.focal_id))
            seen = set()
            target_options = MATCHED_POOL_SIZE * 3
            for _attempt in range(target_options * 100):
                target_pos = _weighted_position(intervals, cumulative, rng)
                if target_pos >= len(sequence) or sequence[target_pos] != row.ref:
                    continue
                alt = DNA_BASES[int(rng.integers(0, len(DNA_BASES)))]
                if alt == row.ref:
                    continue
                genomic_pos = int(gene["begin"]) + target_pos
                variant_key = f"{gene['chrom']}:{genomic_pos}:{row.ref}>{alt}"
                if variant_key == row.variant_key or variant_key in seen:
                    continue
                seen.add(variant_key)
                matched_rows.append(
                    {
                        **common,
                        "role": "control",
                        "option": len(seen),
                        "variant_key": variant_key,
                        "pos": genomic_pos,
                        "target_pos": target_pos,
                        "alt": alt,
                        "callable_species": _depth_at(searchable, target_pos),
                    }
                )
                if len(seen) >= target_options:
                    break

    with gzip.open(callable_blocks_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (str(row["strategy"]), str(row["gene_id"]))
            if current_key is not None and key != current_key:
                process_group(current_key, current_blocks)
                current_blocks = []
            current_key = key
            current_blocks.append(
                (int(row["target_start0"]), int(row["target_end0"]), int(row["callable_species"]))
            )
        process_group(current_key, current_blocks)

    return pd.DataFrame(matched_rows), pd.DataFrame(same_site_rows)


def _collect_observed_control_keys(
    annotations_path: Path,
    matched: pd.DataFrame,
    same_site: pd.DataFrame,
) -> set[tuple[str, str]]:
    wanted: dict[str, set[str]] = defaultdict(set)
    for frame in (matched, same_site):
        if frame.empty:
            continue
        controls = frame[frame["role"] == "control"]
        for row in controls[["variant_key", "strategy"]].itertuples(index=False):
            wanted[str(row.variant_key)].add(str(row.strategy))
    if not wanted:
        return set()

    found: set[tuple[str, str]] = set()
    for chunk in pd.read_csv(
        annotations_path,
        sep="\t",
        compression="gzip",
        usecols=["variant_key", "strategies"],
        keep_default_na=False,
        chunksize=250_000,
    ):
        subset = chunk[chunk["variant_key"].astype(str).isin(wanted)]
        for row in subset.itertuples(index=False):
            key = str(row.variant_key)
            requested = wanted.get(key, set())
            for strategy in split_strategies(str(row.strategies)):
                if strategy in requested:
                    found.add((key, strategy))
    return found


def _finalize_matched_options(
    frame: pd.DataFrame,
    observed: set[tuple[str, str]],
) -> pd.DataFrame:
    if frame.empty:
        return frame
    keep_rows = []
    for _focal_id, group in frame.groupby("focal_id", sort=False):
        focal = group[group["role"] == "observed"]
        controls = group[group["role"] == "control"].copy()
        controls = controls[
            [
                (str(row.variant_key), str(row.strategy)) not in observed
                for row in controls.itertuples(index=False)
            ]
        ]
        controls = controls.drop_duplicates("variant_key").head(MATCHED_POOL_SIZE)
        if focal.empty or controls.empty:
            continue
        keep_rows.append(focal.iloc[[0]])
        controls.loc[:, "option"] = np.arange(1, len(controls) + 1)
        keep_rows.append(controls)
    if not keep_rows:
        return frame.iloc[0:0].copy()
    return pd.concat(keep_rows, ignore_index=True)


def _finalize_same_site_options(
    frame: pd.DataFrame,
    observed: set[tuple[str, str]],
) -> pd.DataFrame:
    if frame.empty:
        return frame
    keep_rows = []
    for _focal_id, group in frame.groupby("focal_id", sort=False):
        focal = group[group["role"] == "observed"]
        controls = group[group["role"] == "control"].copy()
        controls = controls[
            [
                (str(row.variant_key), str(row.strategy)) not in observed
                for row in controls.itertuples(index=False)
            ]
        ]
        controls = controls.drop_duplicates("variant_key")
        if focal.empty or controls.empty:
            continue
        keep_rows.append(focal.iloc[[0]])
        controls.loc[:, "option"] = np.arange(1, len(controls) + 1)
        keep_rows.append(controls)
    if not keep_rows:
        return frame.iloc[0:0].copy()
    return pd.concat(keep_rows, ignore_index=True)


def _read_clinvar_snv_annotations(
    clinvar_vcf: Path,
    regions_bed: Path,
) -> dict[str, str]:
    tabix = shutil.which("tabix")
    if tabix is None:
        raise FileNotFoundError("tabix is required for negative-control ClinVar annotation.")
    result = subprocess.run(
        [tabix, "-R", str(regions_bed), str(clinvar_vcf)],
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode not in (0, 1):
        raise RuntimeError(f"tabix ClinVar query failed: {result.stderr.strip()}")

    labels_by_key: dict[str, set[str]] = defaultdict(set)
    for line in result.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        chrom, pos, _record_id, ref, alts = fields[:5]
        ref = ref.upper()
        if len(ref) != 1 or ref not in DNA_BASES:
            continue
        label = clinvar_label(info_value(fields[7], "CLNSIG"))
        for alt in alts.split(","):
            alt = alt.upper()
            if len(alt) != 1 or alt not in DNA_BASES:
                continue
            key = f"{chrom.removeprefix('chr')}:{int(pos)}:{ref}>{alt}"
            labels_by_key[key].add(label)

    annotations = {}
    for key, labels in labels_by_key.items():
        if labels == {"benign"}:
            annotations[key] = "B/LB"
        elif labels == {"pathogenic"}:
            annotations[key] = "P/LP"
        elif labels == {"excluded_vus"}:
            annotations[key] = "VUS"
        else:
            annotations[key] = "Other"
    return annotations


def _cluster_positions(
    positions: set[int],
    max_gap: int = 20_000,
    max_span: int = 200_000,
) -> list[tuple[int, int]]:
    if not positions:
        return []
    ordered = sorted(positions)
    clusters = []
    start = previous = ordered[0]
    for position in ordered[1:]:
        if position - previous > max_gap or position - start > max_span:
            clusters.append((start, previous))
            start = position
        previous = position
    clusters.append((start, previous))
    return clusters


def _fetch_gnomad_presence(frame: pd.DataFrame) -> tuple[set[str], bool, str]:
    if frame.empty:
        return set(), True, ""
    positions_by_chrom: dict[str, set[int]] = defaultdict(set)
    for row in frame[["chrom", "pos"]].drop_duplicates().itertuples(index=False):
        positions_by_chrom[str(row.chrom)].add(int(row.pos))
    tasks = [
        (chrom, start, end)
        for chrom, positions in positions_by_chrom.items()
        for start, end in _cluster_positions(positions)
    ]
    variants = []
    errors = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(4, max(1, len(tasks)))) as executor:
        futures = {
            executor.submit(fetch_region_variants_recursive, chrom, max(1, start - 1), end + 1): (
                chrom,
                start,
                end,
            )
            for chrom, start, end in tasks
        }
        for future in concurrent.futures.as_completed(futures):
            chrom, start, end = futures[future]
            try:
                variants.extend(future.result())
            except Exception as exc:  # network/API failures must not create partial evidence
                errors.append(f"{chrom}:{start}-{end}: {exc}")
    if errors:
        return set(), False, " | ".join(errors[:5])
    keys = set()
    for variant in variants:
        ref = str(variant.get("ref") or "").upper()
        alt = str(variant.get("alt") or "").upper()
        if len(ref) == 1 and len(alt) == 1 and ref in DNA_BASES and alt in DNA_BASES:
            chrom = str(variant.get("chrom") or "").removeprefix("chr")
            keys.add(f"{chrom}:{int(variant['pos'])}:{ref}>{alt}")
    return keys, True, ""


def _annotate_same_site_controls(
    frame: pd.DataFrame,
    clinvar: dict[str, str],
    gnomad_keys: set[str],
    gnomad_complete: bool,
) -> pd.DataFrame:
    frame = frame.copy()
    frame["clinvar_class"] = frame["variant_key"].map(clinvar).fillna("Not Found")
    frame["clinvar_found"] = frame["clinvar_class"].ne("Not Found")
    if gnomad_complete:
        frame["gnomad_found"] = frame["variant_key"].isin(gnomad_keys)
    else:
        frame["gnomad_found"] = ""
    return frame


def _annotate_matched_conservation(
    matched: pd.DataFrame,
    output_path: Path,
) -> tuple[pd.DataFrame, dict]:
    unique = (
        matched[["variant_key", "chrom", "pos"]]
        .drop_duplicates("variant_key")
        .sort_values(["chrom", "pos", "variant_key"], kind="mergesort")
    )
    rows = [
        {"variant_key": str(row.variant_key), "chrom": str(row.chrom), "pos": str(int(row.pos))}
        for row in unique.itertuples(index=False)
    ]
    track = parse_tracks("phyloP100way")[0]
    try:
        summary = annotate_track(
            rows=rows,
            track=track,
            max_block_bp=250_000,
            max_gap_bp=50_000,
            remote_retries=3,
            retry_sleep_seconds=5.0,
            precision=6,
        )
    except RuntimeError as exc:
        for row in rows:
            row[track.name] = ""
        summary = {
            "track": track.name,
            "status": "failed",
            "error": str(exc),
            "unique_positions": len({(row["chrom"], row["pos"]) for row in rows}),
            "annotated_positions": 0,
        }
    frame = pd.DataFrame(rows)
    _write_tsv(output_path, frame)
    frame[track.name] = pd.to_numeric(frame[track.name], errors="coerce")
    return frame[["variant_key", track.name]], summary


def _summarize_analysis(
    matched: pd.DataFrame,
    same_site: pd.DataFrame,
    manifest: dict,
    manifest_path: Path,
    matched_path: Path,
    same_site_path: Path,
    conservation_path: Path,
    permutations: int,
    seed: int,
) -> NegativeControlAnalysis:
    matched = matched.copy()
    matched["phyloP100way"] = pd.to_numeric(matched["phyloP100way"], errors="coerce")
    matched_summary = _matched_summary(matched, ["strategy"], permutations, seed)
    matched_context_summary = _matched_summary(
        matched,
        ["strategy", "context"],
        permutations,
        seed + 1,
    )
    matched_ecdf = _matched_ecdf(matched)
    same_site_summary = _same_site_summary(
        same_site,
        permutations,
        seed + 2,
        bool(manifest.get("gnomad_complete")),
    )
    return NegativeControlAnalysis(
        matched_summary=matched_summary,
        matched_context_summary=matched_context_summary,
        matched_ecdf=matched_ecdf,
        same_site_summary=same_site_summary,
        manifest=manifest,
        manifest_path=manifest_path,
        matched_path=matched_path,
        same_site_path=same_site_path,
        conservation_path=conservation_path,
        permutations=permutations,
    )


def _group_key(values: object) -> tuple[object, ...]:
    return values if isinstance(values, tuple) else (values,)


def _paired_values(frame: pd.DataFrame, value_column: str) -> tuple[np.ndarray, list[np.ndarray]]:
    observed = []
    controls = []
    for _focal_id, group in frame.groupby("focal_id", sort=False):
        focal_values = pd.to_numeric(
            group.loc[group["role"] == "observed", value_column],
            errors="coerce",
        ).dropna()
        control_values = pd.to_numeric(
            group.loc[group["role"] == "control", value_column],
            errors="coerce",
        ).dropna()
        if focal_values.empty or control_values.empty:
            continue
        observed.append(float(focal_values.iloc[0]))
        controls.append(control_values.to_numpy(dtype=float))
    return np.asarray(observed, dtype=float), controls


def _resampled_statistics(
    observed: np.ndarray,
    controls: list[np.ndarray],
    statistic,
    permutations: int,
    seed: int,
) -> tuple[float, np.ndarray]:
    if observed.size == 0 or not controls:
        return math.nan, np.array([], dtype=float)
    observed_stat = float(statistic(observed))
    rng = np.random.default_rng(seed)
    null = np.empty(permutations, dtype=float)
    for index in range(permutations):
        draw = np.fromiter(
            (values[int(rng.integers(0, len(values)))] for values in controls),
            dtype=float,
            count=len(controls),
        )
        null[index] = float(statistic(draw))
    return observed_stat, null


def _empirical_two_sided(observed: float, null: np.ndarray) -> float:
    if not np.isfinite(observed) or null.size == 0:
        return math.nan
    lower = (1 + int(np.sum(null <= observed))) / (len(null) + 1)
    upper = (1 + int(np.sum(null >= observed))) / (len(null) + 1)
    return min(1.0, 2.0 * min(lower, upper))


def _matched_summary(
    frame: pd.DataFrame,
    group_columns: list[str],
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    columns = [
        *group_columns,
        "matched_focals",
        "observed_median",
        "null_median",
        "null_ci_low",
        "null_ci_high",
        "median_difference",
        "empirical_p",
    ]
    if frame.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    grouper = group_columns[0] if len(group_columns) == 1 else group_columns
    for raw_key, group in frame.groupby(grouper, sort=True):
        key = _group_key(raw_key)
        observed, controls = _paired_values(group, "phyloP100way")
        observed_stat, null = _resampled_statistics(
            observed,
            controls,
            np.median,
            permutations,
            _stable_rank(seed, *key),
        )
        null_median = float(np.median(null)) if null.size else math.nan
        row = dict(zip(group_columns, key))
        row.update(
            {
                "matched_focals": len(controls),
                "observed_median": observed_stat,
                "null_median": null_median,
                "null_ci_low": float(np.quantile(null, 0.025)) if null.size else math.nan,
                "null_ci_high": float(np.quantile(null, 0.975)) if null.size else math.nan,
                "median_difference": observed_stat - null_median,
                "empirical_p": _empirical_two_sided(observed_stat, null),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def _matched_ecdf(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for strategy, group in frame.groupby("strategy", sort=True):
        observed, controls = _paired_values(group, "phyloP100way")
        if observed.size == 0:
            continue
        pooled_controls = np.concatenate(controls)
        combined = np.concatenate([observed, pooled_controls])
        low, high = np.quantile(combined, [0.005, 0.995])
        if not np.isfinite(low) or not np.isfinite(high):
            continue
        grid = np.linspace(low, high, 121) if high > low else np.array([low])
        observed_ordered = np.sort(observed)
        observed_fractions = np.searchsorted(observed_ordered, grid, side="right") / len(observed)
        control_fractions = np.mean(
            [np.searchsorted(np.sort(values), grid, side="right") / len(values) for values in controls],
            axis=0,
        )
        for label, fractions in (
            ("GAPH", observed_fractions),
            ("Matched callable", control_fractions),
        ):
            rows.extend(
                {
                    "strategy": strategy,
                    "set": label,
                    "phyloP100way": float(score),
                    "fraction_leq": float(fraction),
                }
                for score, fraction in zip(grid, fractions)
            )
    return pd.DataFrame(rows)


def _as_boolean(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().map({"true": True, "false": False})


def _same_site_summary(
    frame: pd.DataFrame,
    permutations: int,
    seed: int,
    gnomad_complete: bool,
) -> pd.DataFrame:
    columns = [
        "strategy",
        "metric",
        "matched_focals",
        "observed_rate",
        "null_rate",
        "null_ci_low",
        "null_ci_high",
        "enrichment_ratio",
        "empirical_p",
    ]
    if frame.empty:
        return pd.DataFrame(columns=columns)
    frame = frame.copy()
    frame["clinvar_found"] = _as_boolean(frame["clinvar_found"])
    metrics = [("clinvar_found", "ClinVar found")]
    if gnomad_complete:
        frame["gnomad_found"] = _as_boolean(frame["gnomad_found"])
        metrics.append(("gnomad_found", "gnomAD found"))

    rows = []
    for strategy, group in frame.groupby("strategy", sort=True):
        for value_column, label in metrics:
            observed, controls = _paired_values(group, value_column)
            observed_stat, null = _resampled_statistics(
                observed,
                controls,
                np.mean,
                permutations,
                _stable_rank(seed, strategy, value_column),
            )
            null_rate = float(np.mean(null)) if null.size else math.nan
            rows.append(
                {
                    "strategy": strategy,
                    "metric": label,
                    "matched_focals": len(controls),
                    "observed_rate": observed_stat,
                    "null_rate": null_rate,
                    "null_ci_low": float(np.quantile(null, 0.025)) if null.size else math.nan,
                    "null_ci_high": float(np.quantile(null, 0.975)) if null.size else math.nan,
                    "enrichment_ratio": (
                        observed_stat / null_rate
                        if np.isfinite(observed_stat) and np.isfinite(null_rate) and null_rate > 0
                        else math.nan
                    ),
                    "empirical_p": _empirical_two_sided(observed_stat, null),
                }
            )
    return pd.DataFrame(rows, columns=columns)
