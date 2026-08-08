"""Candidate-wide phyloP distributions for gnomAD-stratified reporting."""

from __future__ import annotations

import json
import tempfile
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd

from analytics.io.artifacts import path_metadata, write_json_atomic, write_tsv_atomic
from analytics.io.performance import PerformanceProfile, profile_stage
from .candidate_conservation_aggregation import (
    CandidateAlleleStore,
    build_candidate_allele_store,
)
from .conservation import (
    DEFAULT_TRACK_NAMES,
    PositionScores,
    format_chrom,
    parse_tracks,
    read_position_scores,
    score_positions,
)
from genomics.variants import parse_variant_key


CACHE_VERSION = 5
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
    variant_annotations_tsv: Path,
    analytics_dir: Path,
    annotation_failures_tsv: Path | None = None,
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
) -> CandidateConservation:
    """Compute compact exact percentile curves without persisting allele-level scores."""
    tracks = parse_tracks(track_names)
    if len(tracks) != 1:
        raise ValueError("Candidate-wide conservation currently requires exactly one track.")
    track = tracks[0]
    analytics_dir.mkdir(parents=True, exist_ok=True)
    distributions_path = (
        analytics_dir / "candidate_variants.phyloP100way.distributions.tsv.gz"
    )
    histograms_path = analytics_dir / "candidate_variants.phyloP100way.histograms.tsv.gz"
    manifest_path = analytics_dir / "candidate_variants.phyloP100way.manifest.json"
    expected_inputs = {
        "cache_version": CACHE_VERSION,
        "variant_annotations": path_metadata(variant_annotations_tsv),
        "annotation_failures": (
            path_metadata(annotation_failures_tsv) if annotation_failures_tsv is not None else None
        ),
        "track": asdict(track),
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
        with profile_stage(performance_profile, "Candidate allele collapse") as timing:
            store = build_candidate_allele_store(
                variant_annotations_tsv=variant_annotations_tsv,
                strategies=strategies,
                annotation_failures_path=annotation_failures_tsv,
                temp_dir=Path(temporary),
            )
            timing["metrics"] = {
                "source_mode": store.source.mode,
                "source_file_count": len(store.source.paths),
                "strategy_count": len(store.strategies),
            }
        try:
            with profile_stage(performance_profile, "Candidate position index") as timing:
                positions_by_chrom, scan_summary, unsupported = _candidate_positions(
                    store,
                    track.chrom_style,
                    chunk_size,
                )
                store.register_unsupported(unsupported)
                _add_positions(positions_by_chrom, additional_rows or [], track.chrom_style)
                timing["metrics"] = {
                    **scan_summary,
                    "position_count_with_clinvar": sum(
                        len(values) for values in positions_by_chrom.values()
                    ),
                }
            with profile_stage(performance_profile, "Candidate phyloP position read") as timing:
                position_scores = read_position_scores(
                    positions_by_chrom=positions_by_chrom,
                    track=track,
                    max_block_bp=max_block_bp,
                    max_gap_bp=max_gap_bp,
                    remote_retries=remote_retries,
                    retry_sleep_seconds=retry_sleep_seconds,
                    precision=precision,
                )
                timing["metrics"] = dict(position_scores.summary)
            with profile_stage(performance_profile, "Candidate distribution summaries") as timing:
                distributions, histograms, groups, membership_summary = _aggregate_distributions(
                    store,
                    position_scores,
                    chunk_size,
                )
                timing["metrics"] = dict(membership_summary)
        finally:
            store.close()
    write_tsv_atomic(distributions_path, distributions)
    write_tsv_atomic(histograms_path, histograms)
    manifest = {
        "inputs": expected_inputs,
        "complete": position_scores.summary.get("status") == "complete",
        "candidate_scan": scan_summary,
        "position_read": position_scores.summary,
        "memberships": membership_summary,
        "aggregation": {
            "engine": "duckdb",
            "source_mode": store.source.mode,
            "source_file_count": len(store.source.paths),
            "source_identity": store.source.identity,
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
) -> tuple[dict[str, set[int]], dict[str, int], list[str]]:
    positions_by_chrom: dict[str, set[int]] = {}
    usable_allele_count = 0
    unsupported_allele_count = 0
    unsupported_keys: list[str] = []
    summary = store.context_summary()
    unsupported_allele_count += summary["position_failed_context_count"]
    for rows in store.iter_position_rows(chunk_size):
        for variant_key, context_count in rows:
            parsed = parse_variant_key(variant_key)
            if parsed is None:
                unsupported_allele_count += int(context_count)
                unsupported_keys.append(str(variant_key))
                continue
            chrom, pos, ref, alt = parsed
            positions, _basis = score_positions(int(pos), str(ref), str(alt))
            positions = [position for position in positions if position >= 0]
            if not positions:
                unsupported_allele_count += int(context_count)
                unsupported_keys.append(str(variant_key))
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
    }, unsupported_keys


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
    position_scores: PositionScores,
    chunk_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]], dict[str, int]]:
    groups_frame = store.group_counts()
    groups_frame["scored_count"] = 0
    distribution_rows = []
    histogram_rows = []
    box_summaries = {}
    for strategy in groups_frame["strategy"].astype(str).unique():
        strategy_groups = groups_frame[groups_frame["strategy"].astype(str).eq(strategy)]
        values_by_status = {}
        for group in strategy_groups.itertuples(index=False):
            status = str(group.gnomad_status)
            values = np.fromiter(
                _iter_group_scores(
                    store=store,
                    strategy=strategy,
                    gnomad_status=status,
                    position_scores=position_scores,
                    chunk_size=chunk_size,
                ),
                dtype=float,
            )
            groups_frame.loc[
                groups_frame["strategy"].eq(strategy)
                & groups_frame["gnomad_status"].eq(status),
                "scored_count",
            ] = int(values.size)
            if values.size == 0:
                continue
            values_by_status[status] = values
            quantile_values = np.quantile(values, QUANTILES)
            distribution_rows.extend(
                {
                    "strategy": strategy,
                    "gnomad_status": status,
                    "quantile": float(quantile),
                    "phyloP100way": float(score),
                    "variant_count": int(group.variant_count),
                    "scored_count": int(values.size),
                }
                for quantile, score in zip(QUANTILES, quantile_values)
            )
            box_summaries[(strategy, status)] = _box_summary(values)

        if not values_by_status:
            continue
        edges = _histogram_edges(np.concatenate(list(values_by_status.values())))
        for status, values in values_by_status.items():
            counts, _ = np.histogram(values, bins=edges)
            histogram_rows.extend(
                {
                    "strategy": strategy,
                    "gnomad_status": status,
                    "bin_left": float(left),
                    "bin_right": float(right),
                    "count": int(count),
                    "fraction": float(count / values.size),
                }
                for left, right, count in zip(edges[:-1], edges[1:], counts)
            )

    groups_frame["scored_count"] = pd.to_numeric(
        groups_frame["scored_count"], errors="coerce"
    ).fillna(0).astype(int)
    summary = store.summary()

    groups = groups_frame.to_dict(orient="records")
    for group in groups:
        variant_count = int(group["variant_count"])
        scored_count = int(group["scored_count"])
        group["variant_count"] = variant_count
        group["scored_count"] = scored_count
        group["score_coverage"] = scored_count / variant_count if variant_count else 0.0
        group.update(box_summaries.get((str(group["strategy"]), str(group["gnomad_status"])), {}))
    return pd.DataFrame(distribution_rows), pd.DataFrame(histogram_rows), groups, {
        "unique_usable_allele_count": summary["unique_usable_allele_count"],
        "strategy_variant_membership_count": summary["strategy_variant_membership_count"],
        "lookup_failed_allele_context_count": summary["lookup_failed_allele_context_count"],
        "gnomad_status_conflict_membership_count": summary[
            "gnomad_status_conflict_membership_count"
        ],
    }


def _iter_group_scores(
    *,
    store: CandidateAlleleStore,
    strategy: str,
    gnomad_status: str,
    position_scores: PositionScores,
    chunk_size: int,
):
    for variant_keys in store.iter_group_variant_keys(
        strategy=strategy,
        gnomad_status=gnomad_status,
        chunk_size=chunk_size,
    ):
        for variant_key in variant_keys:
            parsed = parse_variant_key(variant_key)
            if parsed is None:
                continue
            chrom, pos, ref, alt = parsed
            required = _required_positions(
                chrom,
                pos,
                ref,
                alt,
                position_scores.track.chrom_style,
            )
            values = [position_scores.values.get(position) for position in required]
            if values and all(value is not None for value in values):
                yield float(np.mean(values))


def _histogram_edges(values: np.ndarray) -> np.ndarray:
    if values.size == 0:
        raise ValueError("Cannot calculate histogram bins without values.")
    minimum = float(np.min(values))
    maximum = float(np.max(values))
    if minimum == maximum:
        padding = max(abs(minimum) * 0.05, 0.5)
        return np.asarray([minimum - padding, maximum + padding])
    edges = np.histogram_bin_edges(values, bins="fd")
    if len(edges) - 1 > MAX_HISTOGRAM_BINS:
        edges = np.linspace(minimum, maximum, MAX_HISTOGRAM_BINS + 1)
    return edges


def _box_summary(values: np.ndarray) -> dict[str, float]:
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    iqr = q3 - q1
    lower_candidates = values[values >= q1 - 1.5 * iqr]
    upper_candidates = values[values <= q3 + 1.5 * iqr]
    return {
        "q1": float(q1),
        "median": float(median),
        "q3": float(q3),
        "lower_whisker": float(np.min(lower_candidates)),
        "upper_whisker": float(np.max(upper_candidates)),
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
