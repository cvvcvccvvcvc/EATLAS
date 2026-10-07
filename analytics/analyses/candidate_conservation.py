"""Candidate-wide phyloP distributions for gnomAD-stratified reporting."""

from __future__ import annotations

import json
import tempfile
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd

from analytics.io.artifacts import path_metadata, write_json_atomic, write_tsv_atomic
from analytics.io.performance import PerformanceProfile, profile_stage
from analytics.io.variant_source import variant_source_blocks
from .statistics import weighted_quantiles as _weighted_quantiles
from .candidate_conservation_aggregation import (
    CandidateAlleleStore,
    build_candidate_allele_store,
    resolve_candidate_aggregation_source,
)
from .conservation import (
    DEFAULT_TRACK_NAMES,
    PositionScores,
    format_chrom,
    parse_tracks,
    read_position_scores,
    score_positions,
    track_identity,
)


CACHE_VERSION = 7
QUANTILES = np.linspace(0.0, 1.0, 101)
MAX_HISTOGRAM_BINS = 80
REQUIRED_COLUMNS = {
    "variant_key",
    "lookup_status",
    "strategies",
    "gnomad_af",
}


@dataclass(frozen=True)
class CandidateConservation:
    distributions_path: Path
    histograms_path: Path
    manifest_path: Path
    distributions: pd.DataFrame
    histograms: pd.DataFrame
    manifest: dict
    position_scores: PositionScores | None = None

    def without_position_scores(self) -> "CandidateConservation":
        return replace(self, position_scores=None)


def build_candidate_conservation(
    *,
    variant_annotations_source: Path | Sequence[Path],
    analytics_dir: Path,
    annotation_failures_tsv: Path | Sequence[Path] | None = None,
    additional_rows: list[dict[str, str]] | None = None,
    track_names: str = DEFAULT_TRACK_NAMES,
    max_block_bp: int = 250_000,
    max_gap_bp: int = 50_000,
    remote_retries: int = 3,
    retry_sleep_seconds: float = 5.0,
    precision: int = 6,
    chunk_size: int = 100_000,
    strategies: list[str] | None = None,
    performance_profile: PerformanceProfile | None = None,
    phylop_bigwig: Path | None = None,
    phylop_identity: dict[str, object] | None = None,
) -> CandidateConservation:
    """Compute compact exact percentile curves with temporary allele-level scores."""
    tracks = parse_tracks(track_names, phylop_bigwig=phylop_bigwig)
    if len(tracks) != 1:
        raise ValueError("Candidate-wide conservation currently requires exactly one track.")
    track = tracks[0]
    variant_source = resolve_candidate_aggregation_source(variant_annotations_source)
    analytics_dir.mkdir(parents=True, exist_ok=True)
    distributions_path = (
        analytics_dir / "candidate_variants.phyloP100way.distributions.tsv.gz"
    )
    histograms_path = analytics_dir / "candidate_variants.phyloP100way.histograms.tsv.gz"
    manifest_path = analytics_dir / "candidate_variants.phyloP100way.manifest.json"
    expected_inputs = {
        "cache_version": CACHE_VERSION,
        "variant_annotations": variant_source.identity,
        "annotation_failures": (
            [path_metadata(path) for path in _paths(annotation_failures_tsv)]
            if annotation_failures_tsv is not None
            else None
        ),
        "track": track_identity(track, local_content=phylop_identity),
        "max_block_bp": max_block_bp,
        "max_gap_bp": max_gap_bp,
        "remote_retries": remote_retries,
        "retry_sleep_seconds": retry_sleep_seconds,
        "precision": precision,
        "strategies": sorted(strategies) if strategies is not None else None,
    }
    cached = _load_cache(distributions_path, histograms_path, manifest_path, expected_inputs)
    if cached is not None:
        return cached

    with tempfile.TemporaryDirectory(
        prefix=".candidate_phylop_duckdb.", dir=analytics_dir
    ) as temporary:
        accumulator = CandidateDistributions()
        scan_summary = Counter()
        score_summary = Counter()
        position_summary = {"status": "complete", "processed_positions": 0, "block_count": 0}
        blocks = variant_source_blocks(variant_source)
        for block_index, block in enumerate(blocks, 1):
            label = f" [{block_index}/{len(blocks)}]" if len(blocks) > 1 else ""
            with profile_stage(performance_profile, "Candidate allele collapse" + label):
                store = build_candidate_allele_store(
                    variant_annotations_source=block, strategies=strategies,
                    annotation_failures_path=annotation_failures_tsv,
                    temp_dir=Path(temporary),
                )
            try:
                with profile_stage(performance_profile, "Candidate position index" + label) as timing:
                    positions_by_chrom, scan, unsupported = _candidate_positions(store, track.chrom_style, chunk_size)
                    store.register_unsupported(unsupported)
                    if len(blocks) == 1:
                        _add_positions(positions_by_chrom, additional_rows or [], track.chrom_style)
                    scan_summary.update(scan)
                    timing["metrics"] = scan
                with profile_stage(performance_profile, "Candidate phyloP position read" + label) as timing:
                    position_scores = read_position_scores(
                        positions_by_chrom=positions_by_chrom, track=track,
                        max_block_bp=max_block_bp, max_gap_bp=max_gap_bp,
                        remote_retries=remote_retries, retry_sleep_seconds=retry_sleep_seconds,
                        precision=precision,
                    )
                    timing["metrics"] = dict(position_scores.summary)
                    position_summary["processed_positions"] += sum(map(len, positions_by_chrom.values()))
                    position_summary["block_count"] += int(position_scores.summary.get("block_count", 0))
                    if position_scores.summary.get("status") != "complete":
                        position_summary["status"] = "partial"
                with profile_stage(performance_profile, "Candidate allele score materialization" + label) as timing:
                    scoring = _materialize_candidate_scores(store, position_scores, chunk_size)
                    score_summary.update(scoring)
                    timing["metrics"] = scoring
                with profile_stage(performance_profile, "Candidate distribution summaries" + label):
                    accumulator.add(store)
            finally:
                store.close()
            if len(blocks) > 1:
                del positions_by_chrom, position_scores
        if len(blocks) > 1:
            clinvar_positions = {}
            _add_positions(clinvar_positions, additional_rows or [], track.chrom_style)
            position_scores = read_position_scores(
                positions_by_chrom=clinvar_positions, track=track,
                max_block_bp=max_block_bp, max_gap_bp=max_gap_bp,
                remote_retries=remote_retries, retry_sleep_seconds=retry_sleep_seconds,
                precision=precision,
            )
        distributions, histograms, groups, membership_summary = accumulator.finish()
    write_tsv_atomic(distributions_path, distributions)
    write_tsv_atomic(histograms_path, histograms)
    manifest = {
        "inputs": expected_inputs,
        "complete": position_summary["status"] == "complete" and position_scores.summary.get("status") == "complete",
        "candidate_scan": dict(scan_summary),
        "position_read": {**position_scores.summary, **position_summary},
        "memberships": membership_summary,
        "score_materialization": dict(score_summary),
        "aggregation": {
            "engine": "duckdb",
            "source_mode": variant_source.mode,
            "source_file_count": len(variant_source.paths),
            "source_identity": variant_source.identity,
        },
        "groups": groups,
        "quantile_count": len(QUANTILES),
        "histogram_rule": "Freedman-Diaconis with an 80-bin display cap",
        "distributions_tsv": str(distributions_path),
        "histograms_tsv": str(histograms_path),
        "outputs": {
            distributions_path.name: path_metadata(distributions_path),
            histograms_path.name: path_metadata(histograms_path),
        },
    }
    write_json_atomic(manifest_path, manifest)
    return CandidateConservation(
        distributions_path,
        histograms_path,
        manifest_path,
        distributions,
        histograms,
        manifest,
        position_scores,
    )


def _paths(value: Path | Sequence[Path]) -> tuple[Path, ...]:
    return (value,) if isinstance(value, Path) else tuple(value)


def _load_cache(
    distributions_path: Path,
    histograms_path: Path,
    manifest_path: Path,
    expected_inputs: dict,
) -> CandidateConservation | None:
    if (
        not distributions_path.exists()
        or not histograms_path.exists()
        or not manifest_path.exists()
    ):
        return None
    try:
        manifest = json.loads(manifest_path.read_text())
        expected_outputs = {
            distributions_path.name: path_metadata(distributions_path),
            histograms_path.name: path_metadata(histograms_path),
        }
        if (
            manifest.get("inputs") != expected_inputs
            or manifest.get("complete") is not True
            or manifest.get("outputs") != expected_outputs
        ):
            return None
        distributions = pd.read_csv(distributions_path, sep="\t", compression="gzip")
        histograms = pd.read_csv(histograms_path, sep="\t", compression="gzip")
        return CandidateConservation(
            distributions_path,
            histograms_path,
            manifest_path,
            distributions,
            histograms,
            manifest,
        )
    except (OSError, json.JSONDecodeError, ValueError):
        return None


def _candidate_positions(
    store: CandidateAlleleStore,
    chrom_style: str,
    chunk_size: int,
) -> tuple[dict[str, set[int]], dict[str, int], list[int]]:
    positions_by_chrom: dict[str, set[int]] = {}
    usable_allele_count = 0
    unsupported_allele_count = 0
    unsupported_allele_ids: list[int] = []
    summary = store.context_summary()
    unsupported_allele_count += summary["position_failed_context_count"]
    for rows in store.iter_position_rows(chunk_size):
        for allele_id, key_valid, chrom, pos, ref, alt, context_count in rows:
            if not key_valid:
                unsupported_allele_count += int(context_count)
                unsupported_allele_ids.append(int(allele_id))
                continue
            positions, _basis = score_positions(int(pos), str(ref), str(alt))
            positions = [position for position in positions if position >= 0]
            if not positions:
                unsupported_allele_count += int(context_count)
                unsupported_allele_ids.append(int(allele_id))
                continue
            usable_allele_count += int(context_count)
            if int(context_count) == 0:
                continue
            formatted_chrom = format_chrom(chrom, chrom_style)
            positions_by_chrom.setdefault(formatted_chrom, set()).update(positions)
    return positions_by_chrom, {
        "variant_context_row_count": summary["variant_context_row_count"],
        "usable_allele_context_count": usable_allele_count,
        "unsupported_allele_context_count": unsupported_allele_count,
        "candidate_unique_position_count": sum(len(values) for values in positions_by_chrom.values()),
    }, unsupported_allele_ids


def _add_positions(
    positions_by_chrom: dict[str, set[int]],
    rows: list[dict[str, str]],
    chrom_style: str,
) -> None:
    for row in rows:
        chrom = format_chrom(row.get("chrom", ""), chrom_style)
        positions, _basis = score_positions(int(row["pos"]), row.get("ref", ""), row.get("alt", ""))
        positions_by_chrom.setdefault(chrom, set()).update(position for position in positions if position >= 0)


def _aggregate_distributions(
    store: CandidateAlleleStore,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]], dict[str, int]]:
    accumulator = CandidateDistributions()
    accumulator.add(store)
    return accumulator.finish()


class CandidateDistributions:
    """Exact score frequencies merge across disjoint allele blocks without row expansion."""

    def __init__(self):
        self.groups = Counter()
        self.scores = defaultdict(Counter)
        self.summary = Counter()

    def add(self, store: CandidateAlleleStore) -> None:
        self.summary.update(store.summary())
        for row in store.group_counts().itertuples(index=False):
            key = (str(row.strategy), str(row.gnomad_status))
            self.groups[key] += int(row.variant_count)
            frame = store.group_score_counts(strategy=key[0], gnomad_status=key[1])
            self.scores[key].update({float(r.score): int(r.weight) for r in frame.itertuples(index=False)})

    def finish(self):
        distribution_rows, histogram_rows, groups = [], [], []
        for strategy in sorted({key[0] for key in self.groups}):
            combined = Counter()
            for key, frequencies in self.scores.items():
                if key[0] == strategy:
                    combined.update(frequencies)
            if combined:
                all_values = np.asarray(sorted(combined), dtype=float)
                all_weights = np.asarray([combined[v] for v in all_values], dtype=np.int64)
                edges = _weighted_histogram_edges(all_values, all_weights)
            for key in sorted(k for k in self.groups if k[0] == strategy):
                frequencies = self.scores[key]
                count = sum(frequencies.values())
                group = {
                    "strategy": strategy, "gnomad_status": key[1], "variant_count": self.groups[key],
                    "scored_count": count, "score_coverage": count / self.groups[key] if self.groups[key] else 0.0,
                }
                if count:
                    values = np.asarray(sorted(frequencies), dtype=float)
                    weights = np.asarray([frequencies[v] for v in values], dtype=np.int64)
                    quantiles = _weighted_quantiles(values, weights, QUANTILES)
                    distribution_rows.extend({
                        "strategy": strategy, "gnomad_status": key[1], "quantile": float(q),
                        "phyloP100way": float(v), "variant_count": self.groups[key], "scored_count": count,
                    } for q, v in zip(QUANTILES, quantiles))
                    q1, median, q3 = _weighted_quantiles(values, weights, np.asarray([.25, .5, .75]))
                    group.update({
                        "q1": float(q1), "median": float(median), "q3": float(q3),
                        "lower_whisker": float(values[np.searchsorted(values, q1 - 1.5 * (q3 - q1))]),
                        "upper_whisker": float(values[np.searchsorted(values, q3 + 1.5 * (q3 - q1), side="right") - 1]),
                    })
                    counts = np.histogram(values, bins=edges, weights=weights)[0]
                    histogram_rows.extend({
                        "strategy": strategy, "gnomad_status": key[1], "bin_left": float(left),
                        "bin_right": float(right), "count": int(n), "fraction": float(n / count),
                    } for left, right, n in zip(edges[:-1], edges[1:], counts))
                groups.append(group)
        keys = ("unique_usable_allele_count", "strategy_variant_membership_count",
                "lookup_failed_allele_context_count", "gnomad_status_conflict_membership_count")
        distributions = pd.DataFrame(distribution_rows, columns=[
            "strategy", "gnomad_status", "quantile", "phyloP100way", "variant_count", "scored_count",
        ])
        histograms = pd.DataFrame(histogram_rows, columns=[
            "strategy", "gnomad_status", "bin_left", "bin_right", "count", "fraction",
        ])
        return distributions, histograms, groups, {k: self.summary[k] for k in keys}


def _weighted_histogram_edges(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    minimum, maximum = float(values[0]), float(values[-1])
    if minimum == maximum:
        padding = max(abs(minimum) * .05, .5)
        return np.asarray([minimum - padding, maximum + padding])
    q1, q3 = _weighted_quantiles(values, weights, np.asarray([.25, .75]))
    width = 2 * (q3 - q1) * int(weights.sum()) ** (-1 / 3)
    bins = min(MAX_HISTOGRAM_BINS, int(np.ceil((maximum - minimum) / width))) if width else 1
    return np.linspace(minimum, maximum, bins + 1)


def _materialize_candidate_scores(
    store: CandidateAlleleStore,
    position_scores: PositionScores,
    chunk_size: int,
) -> dict[str, int]:
    attempted_count = 0
    scored_count = 0
    chrom_names = {}
    for rows in store.iter_scoring_rows(chunk_size):
        scores = []
        for allele_id, chrom, pos, ref, alt in rows:
            attempted_count += 1
            if len(ref) == len(alt) == 1:
                if chrom not in chrom_names:
                    chrom_names[chrom] = format_chrom(chrom, position_scores.track.chrom_style)
                value = position_scores.values.get((chrom_names[chrom], int(pos) - 1))
                if value is not None:
                    scores.append((int(allele_id), float(value)))
                continue
            required = _required_positions(
                chrom,
                pos,
                ref,
                alt,
                position_scores.track.chrom_style,
            )
            values = [position_scores.values.get(position) for position in required]
            if values and all(value is not None for value in values):
                scores.append((int(allele_id), float(np.mean(values))))
        store.append_scores(scores)
        scored_count += len(scores)
    return {
        "attempted_allele_count": attempted_count,
        "scored_allele_count": scored_count,
        "missing_score_allele_count": attempted_count - scored_count,
    }


def _required_positions(
    chrom: object,
    pos: object,
    ref: object,
    alt: object,
    chrom_style: str,
) -> list[tuple[str, int]]:
    positions, _basis = score_positions(int(pos), str(ref), str(alt))
    formatted_chrom = format_chrom(chrom, chrom_style)
    return [(formatted_chrom, position) for position in positions if position >= 0]
