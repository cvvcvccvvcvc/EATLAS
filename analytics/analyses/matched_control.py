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
import warnings
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from analytics.io.artifacts import (
    directory_metadata,
    path_metadata,
    write_json_atomic,
    write_tsv_atomic,
)
from analytics.io.performance import PerformanceProfile, profile_stage
from .conservation import annotate_track, parse_tracks
from .external_evidence import build_external_evidence
from .observed_variant_store import (
    ObservedVariantStore,
    build_or_load_observed_variant_store,
)
from .target_context import context_at, read_disjoint_contexts
from genomics.variants import (
    changed_target_position,
    normalize_chrom,
    parse_variant_key,
    refseq_accession_to_chrom,
)
from analytics.annotation.vep import annotate_vep_consequences
from analytics.annotation.vep_result_cache import DEFAULT_TILE_SIZE_BP


CONTROL_VERSION = 5
FOCAL_CACHE_VERSION = 1
MATCHED_POOL_SIZE = 5
CANDIDATE_POOL_SIZE = MATCHED_POOL_SIZE * 3
CANDIDATE_FOCAL_CHUNK_SIZE = 2_000
RESAMPLE_BLOCK_SIZE = 16


@dataclass(frozen=True)
class TargetSpaceNullAnalysis:
    summary: pd.DataFrame
    consequence_summary: pd.DataFrame
    ecdf: pd.DataFrame
    gnomad_summary: pd.DataFrame
    clinvar_summary: pd.DataFrame
    clinvar_class_summary: pd.DataFrame
    manifest: dict
    manifest_path: Path
    matched_path: Path
    conservation_path: Path
    vep_cache_path: Path
    external_evidence_path: Path
    external_evidence_manifest_path: Path
    resamples: int
    focal_path: Path | None = None
    focal_manifest_path: Path | None = None


def build_target_space_null(
    *,
    run_dir: Path,
    variant_annotations_tsv: Path,
    target_features_tsv: Path,
    genes_tsv: Path,
    target_sequences_dir: Path,
    clinvar_vcf: Path,
    strategies: list[str],
    sample_size_per_strategy: int = 25_000,
    resamples: int = 1_000,
    seed: int = 20_260_721,
    gnomad_cache_dir: Path | None = None,
    vep_backend: str = "rest",
    vep_release: str | None = None,
    vep_executable: str | Path = "vep",
    vep_cache_dir: Path | None = None,
    vep_forks: int = 1,
    vep_result_cache_dir: Path | None = None,
    vep_result_cache_tile_size_bp: int = DEFAULT_TILE_SIZE_BP,
    performance_profile: PerformanceProfile | None = None,
) -> TargetSpaceNullAnalysis:
    """Build or load the target-space null for one completed run."""

    if sample_size_per_strategy < 1:
        raise ValueError("target-space-null sample size must be >= 1")
    if resamples < 100:
        raise ValueError("target-space-null resamples must be >= 100")

    outdir = run_dir / "analytics" / "negative_control"
    outdir.mkdir(parents=True, exist_ok=True)
    matched_path = outdir / "target_space_null.snv.tsv.gz"
    focal_path = outdir / "target_space_null.focal_snvs.tsv.gz"
    focal_manifest_path = outdir / "target_space_null.focal_snvs.manifest.json"
    conservation_path = outdir / "target_space_null.phyloP100way.tsv.gz"
    vep_cache_path = outdir / "vep_consequences.sqlite"
    manifest_path = outdir / "manifest.json"
    external_evidence_path = outdir / "target_space_null.external_evidence.tsv.gz"
    external_evidence_manifest_path = outdir / "target_space_null.external_evidence.manifest.json"

    expected_inputs = {
        "version": CONTROL_VERSION,
        "variant_annotations": path_metadata(variant_annotations_tsv),
        "target_features": path_metadata(target_features_tsv),
        "genes": path_metadata(genes_tsv),
        "target_sequences": directory_metadata(target_sequences_dir, "*.fa.gz"),
        "strategies": sorted(strategies),
        "sample_size_per_strategy": sample_size_per_strategy,
        "matched_pool_size": MATCHED_POOL_SIZE,
        "candidate_pool_size": CANDIDATE_POOL_SIZE,
        "candidate_focal_chunk_size": CANDIDATE_FOCAL_CHUNK_SIZE,
        "seed": seed,
        "matching": ["gene_id", "target_context", "ref", "alt", "vep_primary_consequence"],
        "vep": {
            "release": str(vep_release) if vep_release is not None else "current",
            "refseq": True,
            "pick_allele_gene": True,
        },
        "conservation_track": "phyloP100way",
    }
    with profile_stage(performance_profile, "Target-null cache lookup") as timing:
        cache_hit = _cache_is_valid(
            manifest_path,
            expected_inputs,
            [matched_path, conservation_path],
        )
        timing["details"] = "cache hit" if cache_hit else "cache miss"
    if cache_hit:
        manifest = json.loads(manifest_path.read_text())
        return _load_analysis(
            matched_path,
            conservation_path,
            vep_cache_path,
            manifest,
            manifest_path,
            resamples,
            seed,
            clinvar_vcf,
            external_evidence_path,
            external_evidence_manifest_path,
            gnomad_cache_dir,
            performance_profile,
        )

    with profile_stage(performance_profile, "Target-null observed store") as timing:
        observed_store = build_or_load_observed_variant_store(
            variant_annotations_tsv=variant_annotations_tsv,
            analytics_dir=run_dir / "analytics",
            strategies=strategies,
        )
        timing["details"] = "cache hit" if observed_store.cache_hit else "cache miss"
        timing["metrics"] = {
            "source_rows": int(observed_store.manifest["source_row_count"]),
            "allele_gene_rows": int(observed_store.manifest["allele_gene_count"]),
            "alleles": int(observed_store.manifest["allele_count"]),
            "store_bytes": int(observed_store.allele_gene_path.stat().st_size)
            + int(observed_store.allele_path.stat().st_size),
        }

    focal_inputs = {
        "version": FOCAL_CACHE_VERSION,
        "variant_annotations": path_metadata(variant_annotations_tsv),
        "target_features": path_metadata(target_features_tsv),
        "genes": path_metadata(genes_tsv),
        "target_sequences": directory_metadata(target_sequences_dir, "*.fa.gz"),
        "strategies": sorted(strategies),
        "sample_size_per_strategy": sample_size_per_strategy,
        "seed": seed,
    }
    with profile_stage(performance_profile, "Target-null focal sampling") as timing:
        genes = _read_genes(genes_tsv)
        contexts = read_disjoint_contexts(
            target_features_tsv,
            {gene_id: int(gene["length"]) for gene_id, gene in genes.items()},
        )
        focal_cache_hit = _cache_is_valid(
            focal_manifest_path,
            focal_inputs,
            [focal_path],
        )
        sequences: dict[str, str] | None = None
        if focal_cache_hit:
            focal_manifest = json.loads(focal_manifest_path.read_text())
            focal = pd.read_csv(
                focal_path,
                sep="\t",
                compression="gzip",
                dtype=str,
                keep_default_na=False,
            )
            sampled_focal_count = int(focal_manifest["sampled_focal_count"])
            reference_valid_focal_count = int(
                focal_manifest["reference_valid_focal_count"]
            )
            reference_mismatch_count = int(
                focal_manifest["reference_mismatch_count"]
            )
        else:
            focal = _sample_focal_snvs(
                observed_store,
                contexts,
                genes,
                strategies,
                sample_size_per_strategy,
                seed,
            )
            if focal.empty:
                raise ValueError(
                    "No normalized GAPH SNVs were available for the target-space null."
                )
            sampled_focal_count = len(focal)
            sequences = _read_target_sequences(
                target_sequences_dir,
                set(focal["gene_id"]),
            )
            focal, reference_mismatch_count = _validate_focal_reference(
                focal,
                genes,
                sequences,
            )
            reference_valid_focal_count = len(focal)
            if focal.empty:
                raise ValueError(
                    "No sampled GAPH SNVs matched the target reference sequence."
                )
            _write_tsv(focal_path, focal)
            write_json_atomic(
                focal_manifest_path,
                {
                    "inputs": focal_inputs,
                    "complete": True,
                    "sampled_focal_count": sampled_focal_count,
                    "reference_valid_focal_count": reference_valid_focal_count,
                    "reference_mismatch_count": reference_mismatch_count,
                    "outputs": {focal_path.name: path_metadata(focal_path)},
                },
            )
        if sequences is None:
            sequences = _read_target_sequences(
                target_sequences_dir,
                set(focal["gene_id"]),
            )
        timing["details"] = "cache hit" if focal_cache_hit else "cache miss"
        timing["metrics"] = {
            "sampled_focals": int(sampled_focal_count),
            "reference_valid_focals": int(reference_valid_focal_count),
        }

    with profile_stage(performance_profile, "Target-null focal VEP") as timing:
        focal_annotations, focal_vep = annotate_vep_consequences(
            focal[["variant_key", "gene_id", "chrom", "pos", "ref", "alt"]],
            vep_cache_path,
            backend=vep_backend,
            release=vep_release,
            vep_executable=vep_executable,
            vep_cache_dir=vep_cache_dir,
            vep_forks=vep_forks,
            vep_result_cache_dir=vep_result_cache_dir,
            vep_result_cache_tile_size_bp=vep_result_cache_tile_size_bp,
        )
        focal = _merge_vep(focal, focal_annotations)
        focal = focal[focal["vep_status"].eq("ok")].reset_index(drop=True)
        if focal.empty:
            raise ValueError("VEP returned no target-gene consequences for sampled GAPH SNVs.")
        timing["metrics"] = {
            "requested": int(focal_vep.get("requested", 0)),
            "queried": int(focal_vep.get("queried", 0)),
            "retained_focals": int(focal["focal_id"].nunique()),
        }

    with profile_stage(performance_profile, "Target-null control VEP") as timing:
        candidates, generated_candidate_count, candidate_vep = _annotate_candidate_controls(
            focal,
            contexts,
            genes,
            sequences,
            vep_cache_path,
            str(focal_vep["release"]),
            seed,
            vep_backend=vep_backend,
            vep_executable=vep_executable,
            vep_cache_dir=vep_cache_dir,
            vep_forks=vep_forks,
            vep_result_cache_dir=vep_result_cache_dir,
            vep_result_cache_tile_size_bp=vep_result_cache_tile_size_bp,
        )
        if candidates.empty:
            raise ValueError(
                "No consequence-matched target-space control candidates were available."
            )
        timing["metrics"] = {
            "generated_candidates": int(generated_candidate_count),
            "consequence_matched_candidates": int(len(candidates)),
            "vep_requested": int(candidate_vep.get("requested", 0)),
            "vep_queried": int(candidate_vep.get("queried", 0)),
        }

    with profile_stage(
        performance_profile,
        "Target-null observed-control exclusion",
    ) as timing:
        observed_controls = _collect_observed_control_keys(
            observed_store,
            candidates,
            strategies,
        )
        timing["metrics"] = {"observed_memberships": int(len(observed_controls))}

    with profile_stage(performance_profile, "Target-null matching") as timing:
        matched = _build_matched_rows(focal, candidates, observed_controls)
        if matched.empty:
            raise ValueError("No consequence-matched target-space controls were available.")
        matching_diagnostics = _matching_diagnostics(focal, matched)
        timing["metrics"] = {
            "matched_focals": int(
                matched.loc[matched["role"] == "observed", "focal_id"].nunique()
            ),
            "matched_controls": int((matched["role"] == "control").sum()),
        }

    with profile_stage(performance_profile, "Target-null phyloP") as timing:
        conservation_rows, conservation_manifest = _annotate_conservation(
            matched,
            conservation_path,
        )
        matched = matched.merge(
            conservation_rows,
            on="variant_key",
            how="left",
            validate="many_to_one",
        )
        timing["metrics"] = {
            "unique_positions": int(conservation_manifest.get("unique_positions", 0)),
            "annotated_positions": int(
                conservation_manifest.get("annotated_positions", 0)
            ),
        }

    with profile_stage(performance_profile, "Target-null artifact write"):
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
        "outputs": {
            matched_path.name: path_metadata(matched_path),
            conservation_path.name: path_metadata(conservation_path),
        },
        "matched_tsv": str(matched_path),
        "conservation_tsv": str(conservation_path),
        "vep_cache": str(vep_cache_path),
    }
    write_json_atomic(manifest_path, manifest)
    return _summarize_analysis(
        matched,
        manifest,
        manifest_path,
        matched_path,
        conservation_path,
        vep_cache_path,
        resamples,
        seed,
        clinvar_vcf,
        external_evidence_path,
        external_evidence_manifest_path,
        gnomad_cache_dir,
        performance_profile,
    )


def _cache_is_valid(manifest_path: Path, expected_inputs: dict, outputs: list[Path]) -> bool:
    if not manifest_path.exists() or not all(path.exists() for path in outputs):
        return False
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    if manifest.get("complete") is not True or manifest.get("inputs") != expected_inputs:
        return False
    expected_outputs = {path.name: path_metadata(path) for path in outputs}
    return manifest.get("outputs") == expected_outputs


def _load_analysis(
    matched_path: Path,
    conservation_path: Path,
    vep_cache_path: Path,
    manifest: dict,
    manifest_path: Path,
    resamples: int,
    seed: int,
    clinvar_vcf: Path,
    external_evidence_path: Path,
    external_evidence_manifest_path: Path,
    gnomad_cache_dir: Path | None,
    performance_profile: PerformanceProfile | None,
) -> TargetSpaceNullAnalysis:
    with profile_stage(performance_profile, "Target-null cached artifact read"):
        matched = pd.read_csv(
            matched_path,
            sep="\t",
            compression="gzip",
            keep_default_na=False,
        )
        matched["phyloP100way"] = pd.to_numeric(
            matched["phyloP100way"],
            errors="coerce",
        )
    return _summarize_analysis(
        matched,
        manifest,
        manifest_path,
        matched_path,
        conservation_path,
        vep_cache_path,
        resamples,
        seed,
        clinvar_vcf,
        external_evidence_path,
        external_evidence_manifest_path,
        gnomad_cache_dir,
        performance_profile,
    )


def _write_tsv(path: Path, frame: pd.DataFrame) -> None:
    write_tsv_atomic(path, frame)


def _read_genes(path: Path) -> dict[str, dict[str, object]]:
    frame = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
    required = {
        "gene_id",
        "genomic_accession",
        "begin",
        "end",
        "sequence_length",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Genes table missing columns: {', '.join(sorted(missing))}")
    genes = {}
    for row in frame.itertuples(index=False):
        chrom = refseq_accession_to_chrom(str(row.genomic_accession))
        if chrom is None:
            raise ValueError(
                "Genes table contains an unsupported genomic accession: "
                f"{row.genomic_accession}"
            )
        genes[str(row.gene_id)] = {
            "chrom": chrom,
            "begin": int(row.begin),
            "end": int(row.end),
            "length": int(row.sequence_length),
        }
    return genes


def _stable_rank(seed: int, *parts: object) -> int:
    payload = "|".join([str(seed), *(str(part) for part in parts)]).encode()
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "big")


def _sample_focal_snvs(
    observed_store: ObservedVariantStore,
    contexts: dict[str, list[tuple[int, int, str]]],
    genes: dict[str, dict[str, object]],
    strategies: list[str],
    limit: int,
    seed: int,
) -> pd.DataFrame:
    strategy_bits = [
        (strategy, 1 << observed_store.strategies.index(strategy))
        for strategy in strategies
    ]
    heaps: dict[str, list[tuple[int, str, dict[str, object]]]] = defaultdict(list)

    for chunk in observed_store.iter_focal_rows(strategies):
        for variant_key, raw_gene_id, _source_ref, _source_alt, raw_mask in chunk:
            gene_id = str(raw_gene_id)
            parsed = parse_variant_key(variant_key)
            gene = genes.get(gene_id)
            if parsed is None or gene is None:
                continue
            chrom, pos, ref, alt = parsed
            target_pos = changed_target_position(parsed, int(gene["begin"]))
            record_base = {
                "gene_id": gene_id,
                "variant_key": str(variant_key),
                "target_pos": target_pos,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "context": context_at(contexts.get(gene_id, []), target_pos),
            }
            strategy_mask = int(raw_mask)
            for strategy, bit in strategy_bits:
                if not strategy_mask & bit:
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
        chrom = normalize_chrom(str(row["chrom"]))
        if (
            chrom != str(genes[gene_id]["chrom"])
            or target_pos < 0
            or target_pos >= len(sequence)
            or sequence[target_pos] != str(row["ref"])
        ):
            mismatches += 1
            continue
        rows.append(row)
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
    *,
    vep_backend: str,
    vep_executable: str | Path,
    vep_cache_dir: Path | None,
    vep_forks: int,
    vep_result_cache_dir: Path | None,
    vep_result_cache_tile_size_bp: int,
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
            backend=vep_backend,
            vep_executable=vep_executable,
            vep_cache_dir=vep_cache_dir,
            vep_forks=vep_forks,
            vep_result_cache_dir=vep_result_cache_dir,
            vep_result_cache_tile_size_bp=vep_result_cache_tile_size_bp,
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
        "backend": first.get("backend", ""),
        "base_url": first.get("base_url", ""),
        "assembly": first.get("assembly", ""),
        "release": release,
        "options": first.get("options", {}),
        "requested": sum(int(item.get("requested", 0)) for item in summaries),
        "cached": sum(int(item.get("cached", 0)) for item in summaries),
        "shared_cached": sum(int(item.get("shared_cached", 0)) for item in summaries),
        "local_cached": sum(int(item.get("local_cached", 0)) for item in summaries),
        "queried": sum(int(item.get("queried", 0)) for item in summaries),
        "batch_count": sum(int(item.get("batch_count", 0)) for item in summaries),
        "status_counts": dict(sorted(status_counts.items())),
        "cache_path": str(cache_path),
        "vep_cache_dir": first.get("vep_cache_dir", ""),
        "vep_executable": first.get("vep_executable", ""),
        "vep_forks": first.get("vep_forks", ""),
    }


def _collect_observed_control_keys(
    observed_store: ObservedVariantStore,
    candidates: pd.DataFrame,
    strategies: list[str],
) -> set[tuple[str, str]]:
    return observed_store.observed_control_keys(candidates["variant_key"], strategies)


def _build_matched_rows(
    focal: pd.DataFrame,
    candidates: pd.DataFrame,
    observed_controls: set[tuple[str, str]],
) -> pd.DataFrame:
    output_columns = [
        "focal_id",
        "strategy",
        "gene_id",
        "context",
        "primary_consequence",
        "role",
        "option",
        "variant_key",
        "chrom",
        "pos",
        "target_pos",
        "ref",
        "alt",
        "vep_consequence_terms",
        "vep_transcript_id",
    ]
    if focal.empty or candidates.empty:
        return pd.DataFrame(columns=output_columns)

    focal_work = focal.copy()
    focal_work["gene_id"] = focal_work["gene_id"].astype(str)
    focal_work["_focal_order"] = np.arange(len(focal_work))
    focal_work["control_group"] = (
        focal_work["gene_id"] + "|" + focal_work["variant_key"].astype(str)
    )
    candidate_columns = [
        "control_group",
        "variant_key",
        "chrom",
        "pos",
        "target_pos",
        "ref",
        "alt",
        "vep_consequence_terms",
        "vep_transcript_id",
    ]
    candidate_work = candidates[candidate_columns].drop_duplicates(
        ["control_group", "variant_key"],
        keep="first",
    )
    candidate_work["_candidate_order"] = candidate_work.groupby(
        "control_group",
        sort=False,
    ).cumcount()
    candidate_work = candidate_work.rename(
        columns={
            column: f"control_{column}"
            for column in candidate_columns
            if column != "control_group"
        }
    )
    common_columns = [
        "focal_id",
        "strategy",
        "gene_id",
        "context",
        "primary_consequence",
        "control_group",
        "_focal_order",
    ]
    controls = focal_work[common_columns].merge(
        candidate_work,
        on="control_group",
        how="inner",
        sort=False,
        validate="many_to_many",
    )
    if observed_controls and not controls.empty:
        observed_index = pd.MultiIndex.from_tuples(
            observed_controls,
            names=["variant_key", "strategy"],
        )
        control_index = pd.MultiIndex.from_arrays(
            [
                controls["control_variant_key"].astype(str),
                controls["strategy"].astype(str),
            ],
            names=observed_index.names,
        )
        controls = controls[~control_index.isin(observed_index)]
    if controls.empty:
        return pd.DataFrame(columns=output_columns)

    controls = controls.sort_values(
        ["_focal_order", "_candidate_order"],
        kind="mergesort",
    )
    controls["option"] = controls.groupby("focal_id", sort=False).cumcount() + 1
    controls = controls[controls["option"] <= MATCHED_POOL_SIZE].copy()
    matched_focal_ids = set(controls["focal_id"].astype(str))

    observed = focal_work[
        focal_work["focal_id"].astype(str).isin(matched_focal_ids)
    ].copy()
    observed["role"] = "observed"
    observed["option"] = 0

    control_output = controls[common_columns[:5] + ["_focal_order", "option"]].copy()
    control_output["role"] = "control"
    for column in [
        "variant_key",
        "chrom",
        "pos",
        "target_pos",
        "ref",
        "alt",
        "vep_consequence_terms",
        "vep_transcript_id",
    ]:
        control_output[column] = controls[f"control_{column}"]

    matched = pd.concat(
        [
            observed[[*output_columns, "_focal_order"]],
            control_output[[*output_columns, "_focal_order"]],
        ],
        ignore_index=True,
    )
    matched = matched.sort_values(
        ["_focal_order", "option"],
        kind="mergesort",
    ).drop(columns="_focal_order")
    matched["pos"] = pd.to_numeric(matched["pos"], errors="raise").astype(int)
    matched["target_pos"] = pd.to_numeric(
        matched["target_pos"],
        errors="raise",
    ).astype(int)
    return matched[output_columns].reset_index(drop=True)


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
    clinvar_vcf: Path,
    external_evidence_path: Path,
    external_evidence_manifest_path: Path,
    gnomad_cache_dir: Path | None,
    performance_profile: PerformanceProfile | None,
) -> TargetSpaceNullAnalysis:
    matched = matched.copy()
    matched["phyloP100way"] = pd.to_numeric(matched["phyloP100way"], errors="coerce")
    with profile_stage(performance_profile, "Target-null external evidence") as timing:
        evidence, evidence_manifest = build_external_evidence(
            matched=matched,
            matched_path=matched_path,
            clinvar_vcf=clinvar_vcf,
            output_path=external_evidence_path,
            manifest_path=external_evidence_manifest_path,
            gnomad_cache_dir=gnomad_cache_dir,
            performance_profile=performance_profile,
        )
        matched = matched.merge(
            evidence,
            on="variant_key",
            how="left",
            validate="many_to_one",
        )
        timing["metrics"] = {
            "unique_alleles": int(matched["variant_key"].nunique()),
            "gnomad_queried_alleles": int(
                evidence_manifest.get("gnomad", {}).get("queried_allele_count", 0)
            ),
        }
        timing["details"] = (
            "cache hit" if evidence_manifest.get("cache_hit") else "cache miss"
        )

    with profile_stage(performance_profile, "Target-null resampling") as timing:
        matched["gnomad_found_value"] = np.where(
            matched["gnomad_status"].eq("ok"),
            matched["gnomad_found"].astype(float),
            np.nan,
        )
        gnomad_af_values = pd.to_numeric(matched["gnomad_af"], errors="coerce")
        matched["gnomad_af_value"] = gnomad_af_values.where(
            matched["gnomad_status"].eq("ok") & gnomad_af_values.gt(0)
        )
        matched["clinvar_found_value"] = matched["clinvar_found"].astype(float)

        gnomad_found = _matched_metric_summary(
            matched,
            ["strategy"],
            "gnomad_found_value",
            "mean",
            resamples,
            seed + 2,
            status_column="gnomad_status",
        )
        gnomad_found.insert(1, "metric", "found_fraction")
        gnomad_af = _matched_metric_summary(
            matched,
            ["strategy"],
            "gnomad_af_value",
            "median",
            resamples,
            seed + 3,
            status_column="gnomad_status",
        )
        gnomad_af.insert(1, "metric", "median_af")
        clinvar_found = _matched_metric_summary(
            matched,
            ["strategy"],
            "clinvar_found_value",
            "mean",
            resamples,
            seed + 4,
        )

        clinvar_classes = []
        for index, category in enumerate(["B/LB", "P/LP", "VUS", "Other"]):
            value_column = f"clinvar_class_{index}"
            matched[value_column] = np.where(
                matched["clinvar_classified"].astype(bool),
                matched["clinvar_class"].eq(category).astype(float),
                np.nan,
            )
            category_summary = _matched_metric_summary(
                matched,
                ["strategy"],
                value_column,
                "mean",
                resamples,
                seed + 10 + index,
            )
            category_summary.insert(1, "clinvar_class", category)
            clinvar_classes.append(category_summary)

        summary = _matched_summary(matched, ["strategy"], resamples, seed)
        consequence_summary = _matched_summary(
            matched,
            ["strategy", "primary_consequence"],
            resamples,
            seed + 1,
        )
        ecdf = _matched_ecdf(matched)
        timing["metrics"] = {
            "resamples": int(resamples),
            "matched_focal_memberships": int(summary["matched_focals"].sum())
            if not summary.empty
            else 0,
        }

    analysis_manifest = {**manifest, "external_evidence": evidence_manifest}
    return TargetSpaceNullAnalysis(
        summary=summary,
        consequence_summary=consequence_summary,
        ecdf=ecdf,
        gnomad_summary=pd.concat([gnomad_found, gnomad_af], ignore_index=True),
        clinvar_summary=clinvar_found,
        clinvar_class_summary=pd.concat(clinvar_classes, ignore_index=True),
        manifest=analysis_manifest,
        manifest_path=manifest_path,
        matched_path=matched_path,
        conservation_path=conservation_path,
        vep_cache_path=vep_cache_path,
        external_evidence_path=external_evidence_path,
        external_evidence_manifest_path=external_evidence_manifest_path,
        resamples=resamples,
        focal_path=manifest_path.parent / "target_space_null.focal_snvs.tsv.gz",
        focal_manifest_path=(
            manifest_path.parent / "target_space_null.focal_snvs.manifest.json"
        ),
    )


def _group_key(values: object) -> tuple[object, ...]:
    return values if isinstance(values, tuple) else (values,)


def _paired_values(frame: pd.DataFrame, value_column: str) -> tuple[np.ndarray, np.ndarray]:
    working = frame[["focal_id", "role", value_column]].copy()
    working["value"] = pd.to_numeric(working.pop(value_column), errors="coerce")

    observed = (
        working.loc[working["role"].eq("observed") & working["value"].notna(), ["focal_id", "value"]]
        .drop_duplicates("focal_id")
        .set_index("focal_id")["value"]
    )
    controls = working.loc[
        working["role"].eq("control") & working["value"].notna(),
        ["focal_id", "value"],
    ].copy()
    if observed.empty or controls.empty:
        return np.array([], dtype=float), np.empty((0, 0), dtype=float)

    controls["control_index"] = controls.groupby("focal_id", sort=False).cumcount()
    control_matrix = controls.pivot(
        index="focal_id",
        columns="control_index",
        values="value",
    )
    paired = observed.rename("observed").to_frame().join(control_matrix, how="inner")
    if paired.empty:
        return np.array([], dtype=float), np.empty((0, 0), dtype=float)
    return (
        paired.pop("observed").to_numpy(dtype=float),
        paired.to_numpy(dtype=float),
    )


def _resampled_statistics(
    observed: np.ndarray,
    controls: np.ndarray,
    resamples: int,
    seed: int,
) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    if observed.size == 0 or controls.size == 0:
        empty = np.array([], dtype=float)
        return math.nan, empty, empty, empty
    control_counts = np.isfinite(controls).sum(axis=1)
    if np.any(control_counts == 0):
        raise ValueError("Every paired focal must have at least one finite control value.")

    rng = np.random.default_rng(seed)
    observed_bootstrap = np.empty(resamples, dtype=float)
    null = np.empty(resamples, dtype=float)
    difference = np.empty(resamples, dtype=float)
    for start in range(0, resamples, RESAMPLE_BLOCK_SIZE):
        stop = min(resamples, start + RESAMPLE_BLOCK_SIZE)
        focal_indices, control_indices = _matched_set_draw_indices(
            rng,
            control_counts,
            stop - start,
        )
        observed_values = np.median(observed[focal_indices], axis=1)
        null_values = np.median(controls[focal_indices, control_indices], axis=1)
        observed_bootstrap[start:stop] = observed_values
        null[start:stop] = null_values
        difference[start:stop] = observed_values - null_values
    return float(np.median(observed)), observed_bootstrap, null, difference


def _matched_set_draw_indices(
    rng: np.random.Generator,
    control_counts: np.ndarray,
    block_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    focal_indices = rng.integers(
        0,
        len(control_counts),
        size=(block_size, len(control_counts)),
    )
    selected_control_counts = control_counts[focal_indices]
    control_indices = rng.integers(0, selected_control_counts)
    return focal_indices, control_indices


def _paired_metric_values(
    frame: pd.DataFrame,
    value_column: str,
    status_column: str | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    columns = ["focal_id", "role", value_column]
    if status_column:
        columns.append(status_column)
    working = frame[columns].copy()
    if status_column:
        working = working[working[status_column].eq("ok")]
    working["value"] = pd.to_numeric(working.pop(value_column), errors="coerce")

    observed = (
        working.loc[working["role"].eq("observed"), ["focal_id", "value"]]
        .drop_duplicates("focal_id")
        .set_index("focal_id")["value"]
    )
    controls = working.loc[working["role"].eq("control"), ["focal_id", "value"]].copy()
    if observed.empty or controls.empty:
        return np.array([], dtype=float), np.empty((0, 0), dtype=float), np.array([], dtype=int)

    controls["control_index"] = controls.groupby("focal_id", sort=False).cumcount()
    control_counts = controls.groupby("focal_id", sort=False).size().rename("control_count")
    control_matrix = controls.pivot(index="focal_id", columns="control_index", values="value")
    paired = (
        observed.rename("observed")
        .to_frame()
        .join(control_counts, how="inner")
        .join(control_matrix, how="inner")
    )
    if paired.empty:
        return np.array([], dtype=float), np.empty((0, 0), dtype=float), np.array([], dtype=int)
    observed_values = paired.pop("observed").to_numpy(dtype=float)
    counts = paired.pop("control_count").to_numpy(dtype=int)
    return observed_values, paired.to_numpy(dtype=float), counts


def _metric_statistic(values: np.ndarray, statistic: str, axis: int | None = None):
    function = np.nanmean if statistic == "mean" else np.nanmedian
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return function(values, axis=axis)


def _resampled_metric_statistics(
    observed: np.ndarray,
    controls: np.ndarray,
    control_counts: np.ndarray,
    statistic: str,
    resamples: int,
    seed: int,
) -> tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if observed.size == 0 or controls.size == 0:
        empty = np.array([], dtype=float)
        return math.nan, empty, empty, empty, empty
    rng = np.random.default_rng(seed)
    observed_bootstrap = np.empty(resamples, dtype=float)
    null = np.empty(resamples, dtype=float)
    difference = np.empty(resamples, dtype=float)
    null_nonmissing = np.empty(resamples, dtype=float)
    for start in range(0, resamples, RESAMPLE_BLOCK_SIZE):
        stop = min(resamples, start + RESAMPLE_BLOCK_SIZE)
        focal_indices, control_indices = _matched_set_draw_indices(
            rng,
            control_counts,
            stop - start,
        )
        observed_draws = observed[focal_indices]
        control_draws = controls[focal_indices, control_indices]
        observed_values = _metric_statistic(observed_draws, statistic, axis=1)
        null_values = _metric_statistic(control_draws, statistic, axis=1)
        observed_bootstrap[start:stop] = observed_values
        null[start:stop] = null_values
        difference[start:stop] = observed_values - null_values
        null_nonmissing[start:stop] = np.isfinite(control_draws).sum(axis=1)
    return (
        float(_metric_statistic(observed, statistic)),
        observed_bootstrap,
        null,
        difference,
        null_nonmissing,
    )


def _percentile_interval(values: np.ndarray) -> tuple[float, float, int]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return math.nan, math.nan, 0
    low, high = np.quantile(finite, [0.025, 0.975])
    return float(low), float(high), int(finite.size)


def _matched_metric_summary(
    frame: pd.DataFrame,
    group_columns: list[str],
    value_column: str,
    statistic: str,
    resamples: int,
    seed: int,
    *,
    status_column: str | None = None,
) -> pd.DataFrame:
    columns = [
        *group_columns,
        "matched_focals",
        "observed_value",
        "observed_ci_low",
        "observed_ci_high",
        "null_value",
        "null_ci_low",
        "null_ci_high",
        "difference",
        "difference_ci_low",
        "difference_ci_high",
        "observed_nonmissing",
        "null_nonmissing_median",
        "valid_resamples",
    ]
    rows = []
    grouper = group_columns[0] if len(group_columns) == 1 else group_columns
    for raw_key, group in frame.groupby(grouper, sort=True):
        key = _group_key(raw_key)
        observed, controls, control_counts = _paired_metric_values(
            group,
            value_column,
            status_column,
        )
        (
            observed_value,
            observed_bootstrap,
            null,
            difference,
            null_nonmissing,
        ) = _resampled_metric_statistics(
            observed,
            controls,
            control_counts,
            statistic,
            resamples,
            _stable_rank(seed, *key),
        )
        finite_null = null[np.isfinite(null)]
        observed_ci_low, observed_ci_high, _ = _percentile_interval(observed_bootstrap)
        null_ci_low, null_ci_high, _ = _percentile_interval(null)
        difference_ci_low, difference_ci_high, valid_resamples = _percentile_interval(difference)
        null_value = float(np.median(finite_null)) if finite_null.size else math.nan
        row = dict(zip(group_columns, key))
        row.update(
            {
                "matched_focals": len(observed),
                "observed_value": observed_value,
                "observed_ci_low": observed_ci_low,
                "observed_ci_high": observed_ci_high,
                "null_value": null_value,
                "null_ci_low": null_ci_low,
                "null_ci_high": null_ci_high,
                "difference": observed_value - null_value,
                "difference_ci_low": difference_ci_low,
                "difference_ci_high": difference_ci_high,
                "observed_nonmissing": int(np.isfinite(observed).sum()),
                "null_nonmissing_median": float(np.median(null_nonmissing)) if null_nonmissing.size else math.nan,
                "valid_resamples": valid_resamples,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


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
        "observed_ci_low",
        "observed_ci_high",
        "null_median",
        "null_ci_low",
        "null_ci_high",
        "median_difference",
        "difference_ci_low",
        "difference_ci_high",
        "valid_resamples",
    ]
    if frame.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    grouper = group_columns[0] if len(group_columns) == 1 else group_columns
    for raw_key, group in frame.groupby(grouper, sort=True):
        key = _group_key(raw_key)
        observed, controls = _paired_values(group, "phyloP100way")
        observed_median, observed_bootstrap, null, difference = _resampled_statistics(
            observed,
            controls,
            resamples,
            _stable_rank(seed, *key),
        )
        finite_null = null[np.isfinite(null)]
        null_median = float(np.median(finite_null)) if finite_null.size else math.nan
        observed_ci_low, observed_ci_high, _ = _percentile_interval(observed_bootstrap)
        null_ci_low, null_ci_high, _ = _percentile_interval(null)
        difference_ci_low, difference_ci_high, valid_resamples = _percentile_interval(difference)
        row = dict(zip(group_columns, key))
        row.update(
            {
                "matched_focals": len(controls),
                "observed_median": observed_median,
                "observed_ci_low": observed_ci_low,
                "observed_ci_high": observed_ci_high,
                "null_median": null_median,
                "null_ci_low": null_ci_low,
                "null_ci_high": null_ci_high,
                "median_difference": observed_median - null_median,
                "difference_ci_low": difference_ci_low,
                "difference_ci_high": difference_ci_high,
                "valid_resamples": valid_resamples,
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
        finite_controls = np.isfinite(controls)
        control_counts = finite_controls.sum(axis=1)
        pooled_controls = controls[finite_controls]
        combined = np.concatenate([observed, pooled_controls])
        low, high = np.quantile(combined, [0.005, 0.995])
        if not np.isfinite(low) or not np.isfinite(high):
            continue
        grid = np.linspace(low, high, 121) if high > low else np.array([low])
        observed_ordered = np.sort(observed)
        observed_fractions = np.searchsorted(observed_ordered, grid, side="right") / len(observed)
        control_weights = np.broadcast_to(
            (1.0 / control_counts)[:, None],
            controls.shape,
        )[finite_controls]
        order = np.argsort(pooled_controls, kind="mergesort")
        ordered_controls = pooled_controls[order]
        cumulative_weights = np.cumsum(control_weights[order])
        positions = np.searchsorted(ordered_controls, grid, side="right")
        control_fractions = np.zeros(len(grid), dtype=float)
        present = positions > 0
        control_fractions[present] = cumulative_weights[positions[present] - 1] / len(observed)
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
