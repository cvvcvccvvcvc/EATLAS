"""Consequence-matched target-space null for GAPH SNV candidates.

For each sampled GAPH SNV, controls are possible SNVs from the same gene,
target context, genomic REF>ALT substitution, and RefSeq VEP consequence. The
control allele must not be observed by the same GAPH strategy. Conservation is
the measured outcome and is therefore not part of matching.
"""

from __future__ import annotations

import bisect
import gzip
import hashlib
import heapq
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .clinvar_validation import directory_metadata, path_metadata, split_strategies
from .conservation import annotate_track, parse_tracks
from .vep_consequences import annotate_vep_consequences


DNA_BASES = ("A", "C", "G", "T")
CONTEXT_PRIORITY = ("cds", "utr", "exon", "intron")
CONTROL_VERSION = 2
MATCHED_POOL_SIZE = 5
CANDIDATE_POOL_SIZE = MATCHED_POOL_SIZE * 3
CANDIDATE_FOCAL_CHUNK_SIZE = 2_000


@dataclass(frozen=True)
class TargetSpaceNullAnalysis:
    summary: pd.DataFrame
    consequence_summary: pd.DataFrame
    ecdf: pd.DataFrame
    manifest: dict
    manifest_path: Path
    matched_path: Path
    conservation_path: Path
    vep_cache_path: Path
    resamples: int


def build_target_space_null(
    *,
    run_dir: Path,
    variant_annotations_tsv: Path,
    target_features_tsv: Path,
    genes_tsv: Path,
    target_sequences_dir: Path,
    strategies: list[str],
    sample_size_per_strategy: int = 25_000,
    resamples: int = 1_000,
    seed: int = 20_260_721,
) -> TargetSpaceNullAnalysis:
    """Build or load the target-space null for one completed run."""

    if sample_size_per_strategy < 1:
        raise ValueError("target-space-null sample size must be >= 1")
    if resamples < 100:
        raise ValueError("target-space-null resamples must be >= 100")

    outdir = run_dir / "analytics" / "negative_control"
    outdir.mkdir(parents=True, exist_ok=True)
    matched_path = outdir / "target_space_null.snv.tsv.gz"
    conservation_path = outdir / "target_space_null.phyloP100way.tsv.gz"
    vep_cache_path = outdir / "vep_consequences.sqlite"
    manifest_path = outdir / "manifest.json"

    expected_inputs = {
        "version": CONTROL_VERSION,
        "variant_annotations": path_metadata(variant_annotations_tsv),
        "target_features": path_metadata(target_features_tsv),
        "genes": path_metadata(genes_tsv),
        "target_sequences": directory_metadata(target_sequences_dir),
        "strategies": sorted(strategies),
        "sample_size_per_strategy": sample_size_per_strategy,
        "matched_pool_size": MATCHED_POOL_SIZE,
        "candidate_pool_size": CANDIDATE_POOL_SIZE,
        "candidate_focal_chunk_size": CANDIDATE_FOCAL_CHUNK_SIZE,
        "seed": seed,
        "matching": ["gene_id", "target_context", "ref", "alt", "vep_primary_consequence"],
        "conservation_track": "phyloP100way",
    }
    if _cache_is_valid(manifest_path, expected_inputs, [matched_path, conservation_path]):
        manifest = json.loads(manifest_path.read_text())
        return _load_analysis(
            matched_path,
            conservation_path,
            vep_cache_path,
            manifest,
            manifest_path,
            resamples,
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
        raise ValueError("No normalized GAPH SNVs were available for the target-space null.")
    sampled_focal_count = len(focal)

    sequences = _read_target_sequences(target_sequences_dir, set(focal["gene_id"]))
    focal, reference_mismatch_count = _validate_focal_reference(focal, genes, sequences)
    reference_valid_focal_count = len(focal)
    if focal.empty:
        raise ValueError("No sampled GAPH SNVs matched the target reference sequence.")

    focal_annotations, focal_vep = annotate_vep_consequences(
        focal[["variant_key", "gene_id", "chrom", "pos", "ref", "alt"]],
        vep_cache_path,
    )
    focal = _merge_vep(focal, focal_annotations)
    focal = focal[focal["vep_status"].eq("ok")].reset_index(drop=True)
    if focal.empty:
        raise ValueError("VEP returned no target-gene consequences for sampled GAPH SNVs.")

    candidates, generated_candidate_count, candidate_vep = _annotate_candidate_controls(
        focal,
        contexts,
        genes,
        sequences,
        vep_cache_path,
        str(focal_vep["release"]),
        seed,
    )
    if candidates.empty:
        raise ValueError("No consequence-matched target-space control candidates were available.")

    observed_controls = _collect_observed_control_keys(variant_annotations_tsv, candidates)
    matched = _build_matched_rows(focal, candidates, observed_controls)
    if matched.empty:
        raise ValueError("No consequence-matched target-space controls were available.")
    matching_diagnostics = _matching_diagnostics(focal, matched)

    conservation_rows, conservation_manifest = _annotate_conservation(matched, conservation_path)
    matched = matched.merge(conservation_rows, on="variant_key", how="left", validate="many_to_one")
    _write_tsv(matched_path, matched)

    manifest = {
        "inputs": expected_inputs,
        "complete": conservation_manifest.get("status") == "complete",
        "sampled_focal_count": sampled_focal_count,
        "reference_valid_focal_count": reference_valid_focal_count,
        "reference_mismatch_count": reference_mismatch_count,
        "vep_annotated_focal_count": int(focal["focal_id"].nunique()),
        "matched_focal_count": int(matched.loc[matched["role"] == "observed", "focal_id"].nunique()),
        "matched_row_count": int(len(matched)),
        "matched_control_count": int((matched["role"] == "control").sum()),
        "matching_by_consequence": matching_diagnostics.to_dict(orient="records"),
        "generated_control_candidate_count": generated_candidate_count,
        "consequence_matched_candidate_count": int(len(candidates)),
        "focal_vep": focal_vep,
        "candidate_vep": candidate_vep,
        "conservation": conservation_manifest,
        "matched_tsv": str(matched_path),
        "conservation_tsv": str(conservation_path),
        "vep_cache": str(vep_cache_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return _summarize_analysis(
        matched,
        manifest,
        manifest_path,
        matched_path,
        conservation_path,
        vep_cache_path,
        resamples,
        seed,
    )


def _cache_is_valid(manifest_path: Path, expected_inputs: dict, outputs: list[Path]) -> bool:
    if not manifest_path.exists() or not all(path.exists() for path in outputs):
        return False
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    return manifest.get("complete") is not False and manifest.get("inputs") == expected_inputs


def _load_analysis(
    matched_path: Path,
    conservation_path: Path,
    vep_cache_path: Path,
    manifest: dict,
    manifest_path: Path,
    resamples: int,
    seed: int,
) -> TargetSpaceNullAnalysis:
    matched = pd.read_csv(matched_path, sep="\t", compression="gzip", keep_default_na=False)
    matched["phyloP100way"] = pd.to_numeric(matched["phyloP100way"], errors="coerce")
    return _summarize_analysis(
        matched,
        manifest,
        manifest_path,
        matched_path,
        conservation_path,
        vep_cache_path,
        resamples,
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
    return {
        str(row.gene_id): {
            "chrom": str(row.chromosome).removeprefix("chr"),
            "begin": int(row.begin),
            "end": int(row.end),
            "length": int(row.sequence_length),
        }
        for row in frame.itertuples(index=False)
    }


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
                "pos": int(row.genomic_start1),
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
        with gzip.open(path, "rt") as handle:
            sequences[gene_id] = "".join(
                line.strip() for line in handle if not line.startswith(">")
            ).upper()
    return sequences


def _validate_focal_reference(
    focal: pd.DataFrame,
    genes: dict[str, dict[str, object]],
    sequences: dict[str, str],
) -> tuple[pd.DataFrame, int]:
    rows = []
    mismatches = 0
    for row in focal.to_dict(orient="records"):
        gene_id = str(row["gene_id"])
        target_pos = int(row["target_pos"])
        sequence = sequences[gene_id]
        if target_pos < 0 or target_pos >= len(sequence) or sequence[target_pos] != str(row["ref"]):
            mismatches += 1
            continue
        rows.append({**row, "chrom": str(genes[gene_id]["chrom"])})
    return pd.DataFrame(rows), mismatches


def _merge_vep(frame: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
    renamed = annotations.rename(
        columns={
            "status": "vep_status",
            "consequence_terms": "vep_consequence_terms",
            "transcript_id": "vep_transcript_id",
            "mane_select": "vep_mane_select",
            "canonical": "vep_canonical",
            "impact": "vep_impact",
            "variant_class": "vep_variant_class",
        }
    )
    return frame.merge(renamed, on=["variant_key", "gene_id"], how="left", validate="many_to_one")


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


def _generate_candidate_controls(
    focal: pd.DataFrame,
    contexts: dict[str, list[tuple[int, int, str]]],
    genes: dict[str, dict[str, object]],
    sequences: dict[str, str],
    seed: int,
) -> pd.DataFrame:
    rows = []
    unique_focal = focal.drop_duplicates(["variant_key", "gene_id"])
    for focal_row in unique_focal.itertuples(index=False):
        intervals = [
            (start, end)
            for start, end, context in contexts.get(str(focal_row.gene_id), [])
            if context == focal_row.context and end > start
        ]
        if not intervals:
            continue
        cumulative = np.cumsum([end - start for start, end in intervals]).astype(int).tolist()
        sequence = sequences[str(focal_row.gene_id)]
        gene = genes[str(focal_row.gene_id)]
        rng = np.random.default_rng(_stable_rank(seed, focal_row.gene_id, focal_row.variant_key))
        seen = set()
        for _attempt in range(CANDIDATE_POOL_SIZE * 250):
            target_pos = _weighted_position(intervals, cumulative, rng)
            if target_pos >= len(sequence) or sequence[target_pos] != focal_row.ref:
                continue
            pos = int(gene["begin"]) + target_pos
            variant_key = f"{gene['chrom']}:{pos}:{focal_row.ref}>{focal_row.alt}"
            if variant_key == focal_row.variant_key or variant_key in seen:
                continue
            seen.add(variant_key)
            rows.append(
                {
                    "control_group": f"{focal_row.gene_id}|{focal_row.variant_key}",
                    "focal_consequence": focal_row.primary_consequence,
                    "gene_id": str(focal_row.gene_id),
                    "context": focal_row.context,
                    "variant_key": variant_key,
                    "chrom": str(gene["chrom"]),
                    "pos": pos,
                    "target_pos": target_pos,
                    "ref": focal_row.ref,
                    "alt": focal_row.alt,
                }
            )
            if len(seen) >= CANDIDATE_POOL_SIZE:
                break
    return pd.DataFrame(rows)


def _annotate_candidate_controls(
    focal: pd.DataFrame,
    contexts: dict[str, list[tuple[int, int, str]]],
    genes: dict[str, dict[str, object]],
    sequences: dict[str, str],
    vep_cache_path: Path,
    vep_release: str,
    seed: int,
) -> tuple[pd.DataFrame, int, dict[str, object]]:
    unique_focal = focal.drop_duplicates(["variant_key", "gene_id"]).reset_index(drop=True)
    matched_parts = []
    generated_count = 0
    summaries = []
    for start in range(0, len(unique_focal), CANDIDATE_FOCAL_CHUNK_SIZE):
        focal_chunk = unique_focal.iloc[start : start + CANDIDATE_FOCAL_CHUNK_SIZE]
        candidates = _generate_candidate_controls(focal_chunk, contexts, genes, sequences, seed)
        if candidates.empty:
            continue
        generated_count += len(candidates)
        annotations, summary = annotate_vep_consequences(
            candidates[["variant_key", "gene_id", "chrom", "pos", "ref", "alt"]],
            vep_cache_path,
            release=vep_release,
        )
        summaries.append(summary)
        candidates = _merge_vep(candidates, annotations)
        candidates = candidates[
            candidates["vep_status"].eq("ok")
            & candidates["primary_consequence"].eq(candidates["focal_consequence"])
        ]
        if not candidates.empty:
            matched_parts.append(candidates)

    summary = _merge_vep_summaries(summaries, vep_release, vep_cache_path)
    if not matched_parts:
        return pd.DataFrame(), generated_count, summary
    return pd.concat(matched_parts, ignore_index=True), generated_count, summary


def _merge_vep_summaries(
    summaries: list[dict[str, object]],
    release: str,
    cache_path: Path,
) -> dict[str, object]:
    status_counts: dict[str, int] = defaultdict(int)
    for summary in summaries:
        for status, count in dict(summary.get("status_counts", {})).items():
            status_counts[str(status)] += int(count)
    first = summaries[0] if summaries else {}
    return {
        "status": "complete",
        "base_url": first.get("base_url", ""),
        "release": release,
        "options": first.get("options", {}),
        "requested": sum(int(item.get("requested", 0)) for item in summaries),
        "cached": sum(int(item.get("cached", 0)) for item in summaries),
        "queried": sum(int(item.get("queried", 0)) for item in summaries),
        "batch_count": sum(int(item.get("batch_count", 0)) for item in summaries),
        "status_counts": dict(sorted(status_counts.items())),
        "cache_path": str(cache_path),
    }


def _collect_observed_control_keys(
    annotations_path: Path,
    candidates: pd.DataFrame,
) -> set[tuple[str, str]]:
    wanted_keys = set(candidates["variant_key"].astype(str))
    if not wanted_keys:
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
        subset = chunk[chunk["variant_key"].astype(str).isin(wanted_keys)]
        for row in subset.itertuples(index=False):
            for strategy in split_strategies(str(row.strategies)):
                found.add((str(row.variant_key), strategy))
    return found


def _build_matched_rows(
    focal: pd.DataFrame,
    candidates: pd.DataFrame,
    observed_controls: set[tuple[str, str]],
) -> pd.DataFrame:
    candidates_by_group = {
        group_id: group.drop_duplicates("variant_key")
        for group_id, group in candidates.groupby("control_group", sort=False)
    }
    rows = []
    for focal_row in focal.itertuples(index=False):
        group_id = f"{focal_row.gene_id}|{focal_row.variant_key}"
        controls = candidates_by_group.get(group_id)
        if controls is None:
            continue
        controls = controls[
            [
                (str(row.variant_key), str(focal_row.strategy)) not in observed_controls
                for row in controls.itertuples(index=False)
            ]
        ].head(MATCHED_POOL_SIZE)
        if controls.empty:
            continue
        common = {
            "focal_id": focal_row.focal_id,
            "strategy": focal_row.strategy,
            "gene_id": str(focal_row.gene_id),
            "context": focal_row.context,
            "primary_consequence": focal_row.primary_consequence,
        }
        rows.append(
            {
                **common,
                "role": "observed",
                "option": 0,
                "variant_key": focal_row.variant_key,
                "chrom": focal_row.chrom,
                "pos": int(focal_row.pos),
                "target_pos": int(focal_row.target_pos),
                "ref": focal_row.ref,
                "alt": focal_row.alt,
                "vep_consequence_terms": focal_row.vep_consequence_terms,
                "vep_transcript_id": focal_row.vep_transcript_id,
            }
        )
        for option, control in enumerate(controls.itertuples(index=False), start=1):
            rows.append(
                {
                    **common,
                    "role": "control",
                    "option": option,
                    "variant_key": control.variant_key,
                    "chrom": control.chrom,
                    "pos": int(control.pos),
                    "target_pos": int(control.target_pos),
                    "ref": control.ref,
                    "alt": control.alt,
                    "vep_consequence_terms": control.vep_consequence_terms,
                    "vep_transcript_id": control.vep_transcript_id,
                }
            )
    return pd.DataFrame(rows)


def _matching_diagnostics(focal: pd.DataFrame, matched: pd.DataFrame) -> pd.DataFrame:
    eligible = (
        focal.groupby(["strategy", "primary_consequence"], sort=True)["focal_id"]
        .nunique()
        .rename("eligible_focals")
    )
    matched_ids = set(matched.loc[matched["role"] == "observed", "focal_id"].astype(str))
    retained = (
        focal[focal["focal_id"].astype(str).isin(matched_ids)]
        .groupby(["strategy", "primary_consequence"], sort=True)["focal_id"]
        .nunique()
        .rename("matched_focals")
    )
    result = eligible.to_frame().join(retained, how="left").fillna({"matched_focals": 0}).reset_index()
    result["matched_focals"] = result["matched_focals"].astype(int)
    result["match_rate"] = result["matched_focals"] / result["eligible_focals"]
    return result


def _annotate_conservation(
    matched: pd.DataFrame,
    output_path: Path,
) -> tuple[pd.DataFrame, dict]:
    unique = (
        matched[["variant_key", "chrom", "pos", "ref", "alt"]]
        .drop_duplicates("variant_key")
        .sort_values(["chrom", "pos", "variant_key"], kind="mergesort")
    )
    rows = [
        {
            "variant_key": str(row.variant_key),
            "chrom": str(row.chrom),
            "pos": str(int(row.pos)),
            "ref": str(row.ref),
            "alt": str(row.alt),
        }
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
    manifest: dict,
    manifest_path: Path,
    matched_path: Path,
    conservation_path: Path,
    vep_cache_path: Path,
    resamples: int,
    seed: int,
) -> TargetSpaceNullAnalysis:
    matched = matched.copy()
    matched["phyloP100way"] = pd.to_numeric(matched["phyloP100way"], errors="coerce")
    return TargetSpaceNullAnalysis(
        summary=_matched_summary(matched, ["strategy"], resamples, seed),
        consequence_summary=_matched_summary(
            matched,
            ["strategy", "primary_consequence"],
            resamples,
            seed + 1,
        ),
        ecdf=_matched_ecdf(matched),
        manifest=manifest,
        manifest_path=manifest_path,
        matched_path=matched_path,
        conservation_path=conservation_path,
        vep_cache_path=vep_cache_path,
        resamples=resamples,
    )


def _group_key(values: object) -> tuple[object, ...]:
    return values if isinstance(values, tuple) else (values,)


def _paired_values(frame: pd.DataFrame, value_column: str) -> tuple[np.ndarray, list[np.ndarray]]:
    observed = []
    controls = []
    for _focal_id, group in frame.groupby("focal_id", sort=False):
        focal_values = pd.to_numeric(
            group.loc[group["role"] == "observed", value_column], errors="coerce"
        ).dropna()
        control_values = pd.to_numeric(
            group.loc[group["role"] == "control", value_column], errors="coerce"
        ).dropna()
        if focal_values.empty or control_values.empty:
            continue
        observed.append(float(focal_values.iloc[0]))
        controls.append(control_values.to_numpy(dtype=float))
    return np.asarray(observed, dtype=float), controls


def _resampled_statistics(
    observed: np.ndarray,
    controls: list[np.ndarray],
    resamples: int,
    seed: int,
) -> tuple[float, np.ndarray]:
    if observed.size == 0 or not controls:
        return math.nan, np.array([], dtype=float)
    rng = np.random.default_rng(seed)
    null = np.empty(resamples, dtype=float)
    for index in range(resamples):
        draw = np.fromiter(
            (values[int(rng.integers(0, len(values)))] for values in controls),
            dtype=float,
            count=len(controls),
        )
        null[index] = float(np.median(draw))
    return float(np.median(observed)), null


def _matched_summary(
    frame: pd.DataFrame,
    group_columns: list[str],
    resamples: int,
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
    ]
    if frame.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    grouper = group_columns[0] if len(group_columns) == 1 else group_columns
    for raw_key, group in frame.groupby(grouper, sort=True):
        key = _group_key(raw_key)
        observed, controls = _paired_values(group, "phyloP100way")
        observed_median, null = _resampled_statistics(
            observed,
            controls,
            resamples,
            _stable_rank(seed, *key),
        )
        null_median = float(np.median(null)) if null.size else math.nan
        row = dict(zip(group_columns, key))
        row.update(
            {
                "matched_focals": len(controls),
                "observed_median": observed_median,
                "null_median": null_median,
                "null_ci_low": float(np.quantile(null, 0.025)) if null.size else math.nan,
                "null_ci_high": float(np.quantile(null, 0.975)) if null.size else math.nan,
                "median_difference": observed_median - null_median,
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
            ("Matched target-space null", control_fractions),
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
