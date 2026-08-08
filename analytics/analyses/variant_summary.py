"""Bounded-memory aggregation of variant annotations for the HTML report."""

from __future__ import annotations

import gzip
import json
import math
import tempfile
from collections import defaultdict
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from analytics.io.performance import PerformanceProfile
from genomics.clinvar import CLINVAR_CLASS_ORDER
from .variant_summary_aggregation import (
    VariantGroupedAggregation,
    aggregate_variant_groups,
    resolve_variant_aggregation_source,
)


VARIANT_USECOLS = [
    "variant_key",
    "gene_id",
    "event_type",
    "ref",
    "alt",
    "lookup_status",
    "strategies",
    "support_row_count",
    "support_ortholog_count",
    "clinvar_id",
    "clinvar_allele_id",
    "clinvar_sig",
    "clinvar_revstat",
    "clinvar_review_stars",
    "clinvar_scv_count",
    "clinvar_hgvs",
    "clinvar_disease",
    "clinvar_variant_type",
    "gnomad_af",
    "gnomad_csq",
]
VEP_USECOLS = ["vep_status", "vep_primary_consequence"]
SUMMARY_CACHE_VERSION = 13
SUMMARY_CACHE_NAME = "variant_summary.json.gz"
SPECIAL_FLOAT_KEY = "__gaph_float__"
ORTHOLOG_EVIDENCE_COLUMNS = [
    "strategy",
    "target_context",
    "taxonomic_scope",
    "evidence_unit",
    "quantile_count",
    "depth_bin",
    "alt_bin",
    "depth_label",
    "alt_label",
    "gnomad_found_count",
    "gnomad_eligible_count",
    "gnomad_found_fraction",
]
ORTHOLOG_EVIDENCE_DISTRIBUTION_COLUMNS = [
    "strategy",
    "taxonomic_scope",
    "evidence_unit",
    "metric",
    "value",
    "variant_count",
]


@dataclass(frozen=True)
class StrategyOverlap:
    strategies: list[str]
    intersections: np.ndarray
    unions: np.ndarray
    jaccard: np.ndarray


@dataclass(frozen=True)
class VariantSummary:
    input_row_count: int
    unique_variant_count: int
    all_strategy_variant_count: int
    strategy_record_count: int
    gene_count: int
    clinvar_found: int
    clinvar_classified: int
    gnomad_found: int
    gnomad_lookup_failed: int
    pathogenic_variant_count: int
    consequence_source: str
    strategies: list[str]
    strategy_stats: pd.DataFrame
    unique_contribution: pd.DataFrame
    gene_variant_counts: pd.DataFrame
    event_counts: pd.DataFrame
    target_context_counts: pd.DataFrame
    gnomad_event_counts: pd.DataFrame
    gnomad_context_counts: pd.DataFrame
    overlap: StrategyOverlap | None
    clinvar_counts: pd.DataFrame
    gnomad_af_summary: pd.DataFrame
    pathogenic_star_counts: pd.DataFrame
    consequence_counts: pd.DataFrame
    pathogenic_consequence_counts: pd.DataFrame
    pathogenic_rows: pd.DataFrame
    ortholog_evidence_available: bool
    ortholog_evidence_cells: pd.DataFrame
    ortholog_evidence_distributions: pd.DataFrame
    cache_hit: bool = False


def _input_metadata(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _variant_input_metadata(path: Path) -> dict[str, object]:
    source = resolve_variant_aggregation_source(path)
    return {
        "path": str(path.resolve()),
        "source": source.identity,
    }


def _encode_scalar(value: object) -> object:
    if value is None or value is pd.NA:
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        numeric = float(value)
        if math.isnan(numeric):
            return None
        if math.isinf(numeric):
            return {SPECIAL_FLOAT_KEY: "inf" if numeric > 0 else "-inf"}
        return numeric
    if isinstance(value, str):
        return value
    raise TypeError(f"Unsupported VariantSummary cache value: {type(value).__name__}")


def _decode_scalar(value: object) -> object:
    if isinstance(value, dict) and set(value) == {SPECIAL_FLOAT_KEY}:
        return float(str(value[SPECIAL_FLOAT_KEY]))
    return value


def _frame_payload(frame: pd.DataFrame) -> dict[str, object]:
    return {
        "columns": frame.columns.astype(str).tolist(),
        "data": [
            [_encode_scalar(value) for value in row]
            for row in frame.itertuples(index=False, name=None)
        ],
    }


def _frame_from_payload(payload: dict[str, object]) -> pd.DataFrame:
    columns = [str(value) for value in payload["columns"]]
    data = [
        [_decode_scalar(value) for value in row]
        for row in payload["data"]
    ]
    return pd.DataFrame(data, columns=columns)


def _summary_payload(
    summary: VariantSummary,
    source: Path,
    target_features: Path | None,
    genes: Path | None,
    annotation_failures: Path | None,
    variant_strategy_support: Path | None,
    ortholog_evidence_summary: Path | None,
    strategy_label: Callable[[str], str],
) -> dict[str, object]:
    overlap = None
    if summary.overlap is not None:
        overlap = {
            "strategies": summary.overlap.strategies,
            "intersections": summary.overlap.intersections.tolist(),
            "unions": summary.overlap.unions.tolist(),
            "jaccard": summary.overlap.jaccard.tolist(),
        }
    frame_names = [
        "strategy_stats",
        "unique_contribution",
        "gene_variant_counts",
        "event_counts",
        "target_context_counts",
        "gnomad_event_counts",
        "gnomad_context_counts",
        "clinvar_counts",
        "gnomad_af_summary",
        "pathogenic_star_counts",
        "consequence_counts",
        "pathogenic_consequence_counts",
        "pathogenic_rows",
        "ortholog_evidence_cells",
        "ortholog_evidence_distributions",
    ]
    return {
        "cache_version": SUMMARY_CACHE_VERSION,
        "input": _variant_input_metadata(source),
        "target_features": _input_metadata(target_features) if target_features is not None else None,
        "genes": _input_metadata(genes) if genes is not None else None,
        "annotation_failures": (
            _input_metadata(annotation_failures) if annotation_failures is not None else None
        ),
        "variant_strategy_support": (
            _input_metadata(variant_strategy_support)
            if variant_strategy_support is not None
            else None
        ),
        "ortholog_evidence_summary": (
            _input_metadata(ortholog_evidence_summary)
            if ortholog_evidence_summary is not None
            else None
        ),
        "strategy_labels": {
            strategy: strategy_label(strategy)
            for strategy in summary.strategies
        },
        "summary": {
            "input_row_count": summary.input_row_count,
            "unique_variant_count": summary.unique_variant_count,
            "all_strategy_variant_count": summary.all_strategy_variant_count,
            "strategy_record_count": summary.strategy_record_count,
            "gene_count": summary.gene_count,
            "clinvar_found": summary.clinvar_found,
            "clinvar_classified": summary.clinvar_classified,
            "gnomad_found": summary.gnomad_found,
            "gnomad_lookup_failed": summary.gnomad_lookup_failed,
            "pathogenic_variant_count": summary.pathogenic_variant_count,
            "consequence_source": summary.consequence_source,
            "ortholog_evidence_available": summary.ortholog_evidence_available,
            "strategies": summary.strategies,
            "overlap": overlap,
            "frames": {
                name: _frame_payload(getattr(summary, name))
                for name in frame_names
            },
        },
    }


def _summary_from_payload(payload: dict[str, object]) -> VariantSummary:
    summary = payload["summary"]
    overlap_payload = summary["overlap"]
    overlap = None
    if overlap_payload is not None:
        overlap = StrategyOverlap(
            strategies=[str(value) for value in overlap_payload["strategies"]],
            intersections=np.asarray(overlap_payload["intersections"], dtype=np.int64),
            unions=np.asarray(overlap_payload["unions"], dtype=np.int64),
            jaccard=np.asarray(overlap_payload["jaccard"], dtype=float),
        )
    frames = summary["frames"]
    return VariantSummary(
        input_row_count=int(summary["input_row_count"]),
        unique_variant_count=int(summary["unique_variant_count"]),
        all_strategy_variant_count=int(summary["all_strategy_variant_count"]),
        strategy_record_count=int(summary["strategy_record_count"]),
        gene_count=int(summary["gene_count"]),
        clinvar_found=int(summary["clinvar_found"]),
        clinvar_classified=int(summary["clinvar_classified"]),
        gnomad_found=int(summary["gnomad_found"]),
        gnomad_lookup_failed=int(summary["gnomad_lookup_failed"]),
        pathogenic_variant_count=int(summary["pathogenic_variant_count"]),
        consequence_source=str(summary["consequence_source"]),
        ortholog_evidence_available=bool(summary["ortholog_evidence_available"]),
        strategies=[str(value) for value in summary["strategies"]],
        strategy_stats=_frame_from_payload(frames["strategy_stats"]),
        unique_contribution=_frame_from_payload(frames["unique_contribution"]),
        gene_variant_counts=_frame_from_payload(frames["gene_variant_counts"]),
        event_counts=_frame_from_payload(frames["event_counts"]),
        target_context_counts=_frame_from_payload(frames["target_context_counts"]),
        gnomad_event_counts=_frame_from_payload(frames["gnomad_event_counts"]),
        gnomad_context_counts=_frame_from_payload(frames["gnomad_context_counts"]),
        overlap=overlap,
        clinvar_counts=_frame_from_payload(frames["clinvar_counts"]),
        gnomad_af_summary=_frame_from_payload(frames["gnomad_af_summary"]),
        pathogenic_star_counts=_frame_from_payload(frames["pathogenic_star_counts"]),
        consequence_counts=_frame_from_payload(frames["consequence_counts"]),
        pathogenic_consequence_counts=_frame_from_payload(frames["pathogenic_consequence_counts"]),
        pathogenic_rows=_frame_from_payload(frames["pathogenic_rows"]),
        ortholog_evidence_cells=_frame_from_payload(frames["ortholog_evidence_cells"]),
        ortholog_evidence_distributions=_frame_from_payload(
            frames["ortholog_evidence_distributions"]
        ),
        cache_hit=True,
    )


def _load_summary_cache(
    cache_path: Path,
    source: Path,
    target_features: Path | None,
    genes: Path | None,
    annotation_failures: Path | None,
    variant_strategy_support: Path | None,
    ortholog_evidence_summary: Path | None,
    strategy_label: Callable[[str], str],
) -> VariantSummary | None:
    if not cache_path.exists():
        return None
    try:
        with gzip.open(cache_path, "rt") as handle:
            payload = json.load(handle)
        if payload.get("cache_version") != SUMMARY_CACHE_VERSION:
            return None
        if payload.get("input") != _variant_input_metadata(source):
            return None
        expected_features = _input_metadata(target_features) if target_features is not None else None
        if payload.get("target_features") != expected_features:
            return None
        expected_genes = _input_metadata(genes) if genes is not None else None
        if payload.get("genes") != expected_genes:
            return None
        expected_failures = (
            _input_metadata(annotation_failures) if annotation_failures is not None else None
        )
        if payload.get("annotation_failures") != expected_failures:
            return None
        expected_support = (
            _input_metadata(variant_strategy_support)
            if variant_strategy_support is not None
            else None
        )
        if payload.get("variant_strategy_support") != expected_support:
            return None
        expected_ortholog_evidence = (
            _input_metadata(ortholog_evidence_summary)
            if ortholog_evidence_summary is not None
            else None
        )
        if payload.get("ortholog_evidence_summary") != expected_ortholog_evidence:
            return None
        labels = payload.get("strategy_labels", {})
        if any(strategy_label(strategy) != label for strategy, label in labels.items()):
            return None
        return _summary_from_payload(payload)
    except (OSError, EOFError, json.JSONDecodeError, KeyError, TypeError, ValueError):
        return None


def _write_summary_cache(
    cache_path: Path,
    summary: VariantSummary,
    source: Path,
    target_features: Path | None,
    genes: Path | None,
    annotation_failures: Path | None,
    variant_strategy_support: Path | None,
    ortholog_evidence_summary: Path | None,
    strategy_label: Callable[[str], str],
) -> None:
    payload = _summary_payload(
        summary,
        source,
        target_features,
        genes,
        annotation_failures,
        variant_strategy_support,
        ortholog_evidence_summary,
        strategy_label,
    )
    with tempfile.NamedTemporaryFile(
        prefix=f".{cache_path.name}.",
        suffix=".tmp",
        dir=cache_path.parent,
        delete=False,
    ) as handle:
        temporary_path = Path(handle.name)
    try:
        with gzip.open(temporary_path, "wt") as handle:
            json.dump(payload, handle, allow_nan=False, separators=(",", ":"))
        temporary_path.chmod(0o644)
        temporary_path.replace(cache_path)
    finally:
        temporary_path.unlink(missing_ok=True)


def _categorize_clinvar(values: pd.Series, record_presence: pd.Series) -> pd.Categorical:
    text = values.fillna("").astype(str).str.lower()
    category = pd.Series("Other", index=values.index, dtype="object")
    if pd.api.types.is_bool_dtype(record_presence.dtype):
        found = record_presence.fillna(False)
    else:
        found = record_presence.fillna("").astype(str).ne("")
    category[~found] = "Not in ClinVar"
    category[found & text.eq("")] = "Unclassified"
    conflicting = text.str.contains("conflicting", na=False)
    uncertain = text.str.contains("uncertain|vus", regex=True, na=False)
    benign = text.str.contains("benign", na=False)
    pathogenic = text.str.contains("pathogenic", na=False)
    category[found & conflicting] = "Other"
    category[found & uncertain & ~conflicting] = "VUS"
    category[found & pathogenic & ~benign & ~uncertain & ~conflicting] = "P/LP"
    category[found & benign & ~pathogenic & ~uncertain & ~conflicting] = "B/LB"
    return pd.Categorical(category, categories=CLINVAR_CLASS_ORDER, ordered=True)


def _add_pathogenic_strategy_support(path: Path | None, variants: pd.DataFrame) -> pd.DataFrame:
    variants = variants.copy()
    for column in ["support_ortholog_mean", "support_ortholog_min", "support_ortholog_max"]:
        variants[column] = np.nan
    if path is None or variants.empty:
        return variants

    keys = {str(value).encode() for value in variants["variant_id"]}
    support: dict[str, list[int]] = defaultdict(list)
    with gzip.open(path, "rb") as handle:
        header = handle.readline().rstrip(b"\r\n").split(b"\t")
        required = [b"variant_key", b"alt_support_ortholog_count"]
        if any(column not in header for column in required):
            raise ValueError(
                "Variant strategy support table needs variant_key and alt_support_ortholog_count."
            )
        key_index = header.index(b"variant_key")
        count_index = header.index(b"alt_support_ortholog_count")
        for line in handle:
            fields = line.rstrip(b"\r\n").split(b"\t")
            if len(fields) <= max(key_index, count_index) or fields[key_index] not in keys:
                continue
            try:
                count = int(fields[count_index])
            except ValueError:
                continue
            support[fields[key_index].decode()].append(count)

    for index, variant_id in variants["variant_id"].items():
        values = support.get(str(variant_id), [])
        if not values:
            continue
        variants.at[index, "support_ortholog_mean"] = float(np.mean(values))
        variants.at[index, "support_ortholog_min"] = min(values)
        variants.at[index, "support_ortholog_max"] = max(values)
    return variants


def _weighted_quantile_bins(
    values: pd.Series,
    weights: pd.Series,
    quantile_count: int,
) -> np.ndarray:
    grouped = (
        pd.DataFrame({"value": values, "weight": weights})
        .groupby("value", as_index=False, sort=True)["weight"]
        .sum()
    )
    cumulative = grouped["weight"].cumsum().to_numpy()
    total = int(grouped["weight"].sum())
    boundaries = []
    for index in range(1, quantile_count):
        position = int(np.searchsorted(cumulative, total * index / quantile_count, side="left"))
        boundaries.append(float(grouped.iloc[position]["value"]))
    return np.searchsorted(boundaries, values.to_numpy(dtype=float), side="left")


def _bin_labels(values: pd.Series, bins: np.ndarray, count: int, *, percent: bool) -> dict[int, str]:
    labels = {}
    numeric = values.to_numpy(dtype=float)
    for index in range(count):
        observed = numeric[bins == index]
        if observed.size == 0:
            labels[index] = f"Q{index + 1} (empty)"
            continue
        low = float(observed.min())
        high = float(observed.max())
        if percent:
            low_text = f"{low:.0%}"
            high_text = f"{high:.0%}"
        else:
            low_text = f"{int(low):,}"
            high_text = f"{int(high):,}"
        labels[index] = low_text if low == high else f"{low_text}-{high_text}"
    return labels


def _ortholog_evidence_distributions(
    frame: pd.DataFrame,
    *,
    count_column: str,
    site_column: str,
    alt_column: str,
) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=ORTHOLOG_EVIDENCE_DISTRIBUTION_COLUMNS)
    group_columns = ["strategy", "taxonomic_scope", "evidence_unit"]
    distributions = []
    for metric, value_column in (
        ("site_aligned", site_column),
        ("exact_alt", alt_column),
    ):
        grouped = (
            frame.groupby([*group_columns, value_column], as_index=False, sort=True)[
                count_column
            ]
            .sum()
            .rename(columns={value_column: "value", count_column: "variant_count"})
        )
        grouped["metric"] = metric
        distributions.append(grouped)
    return pd.concat(distributions, ignore_index=True)[
        ORTHOLOG_EVIDENCE_DISTRIBUTION_COLUMNS
    ]



def read_taxonomic_ortholog_evidence(
    path: Path,
) -> tuple[bool, pd.DataFrame, pd.DataFrame]:
    """Bin a compact pipeline histogram for interactive report heatmaps."""
    empty = pd.DataFrame(columns=ORTHOLOG_EVIDENCE_COLUMNS)
    empty_distributions = pd.DataFrame(columns=ORTHOLOG_EVIDENCE_DISTRIBUTION_COLUMNS)
    required = {
        "strategy",
        "target_context",
        "taxonomic_scope",
        "evidence_unit",
        "site_aligned_count",
        "alt_support_count",
        "gnomad_found_count",
        "gnomad_not_found_count",
        "gnomad_lookup_failed_count",
    }
    frame = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(
            f"Ortholog evidence summary {path} missing columns: {', '.join(sorted(missing))}"
        )
    if frame.empty:
        return True, empty, empty_distributions
    for column in (
        "site_aligned_count",
        "alt_support_count",
        "gnomad_found_count",
        "gnomad_not_found_count",
        "gnomad_lookup_failed_count",
    ):
        frame[column] = pd.to_numeric(frame[column], errors="raise").astype("int64")
    invalid = (
        frame["site_aligned_count"].le(0)
        | frame["alt_support_count"].lt(0)
        | frame["alt_support_count"].gt(frame["site_aligned_count"])
    )
    if invalid.any():
        row = frame.loc[invalid].iloc[0]
        raise ValueError(
            "Invalid taxonomic ortholog evidence: "
            f"ALT={row['alt_support_count']}, site={row['site_aligned_count']}"
        )
    frame["variant_count"] = (
        frame["gnomad_found_count"]
        + frame["gnomad_not_found_count"]
        + frame["gnomad_lookup_failed_count"]
    )
    frame = frame[frame["target_context"].isin(("cds", "utr", "intron"))]
    distributions = _ortholog_evidence_distributions(
        frame,
        count_column="variant_count",
        site_column="site_aligned_count",
        alt_column="alt_support_count",
    )
    frame["gnomad_eligible_count"] = (
        frame["gnomad_found_count"] + frame["gnomad_not_found_count"]
    )
    frame = frame[frame["gnomad_eligible_count"].gt(0)]
    if frame.empty:
        return True, empty, distributions

    cells = []
    group_columns = ["strategy", "target_context", "taxonomic_scope", "evidence_unit"]
    for group, subset in frame.groupby(group_columns, sort=True):
        strategy, context, scope, unit = group
        subset = subset.copy()
        for quantile_count in (2, 4, 10):
            depth_bins = _weighted_quantile_bins(
                subset["site_aligned_count"], subset["gnomad_eligible_count"], quantile_count
            )
            alt_bins = _weighted_quantile_bins(
                subset["alt_support_count"], subset["gnomad_eligible_count"], quantile_count
            )
            depth_labels = _bin_labels(
                subset["site_aligned_count"], depth_bins, quantile_count, percent=False
            )
            alt_labels = _bin_labels(
                subset["alt_support_count"], alt_bins, quantile_count, percent=False
            )
            binned = subset.assign(depth_bin=depth_bins, alt_bin=alt_bins)
            aggregated = binned.groupby(["depth_bin", "alt_bin"], as_index=False)[
                ["gnomad_found_count", "gnomad_eligible_count"]
            ].sum()
            for row in aggregated.itertuples(index=False):
                eligible = int(row.gnomad_eligible_count)
                found = int(row.gnomad_found_count)
                cells.append(
                    {
                        "strategy": str(strategy),
                        "target_context": str(context),
                        "taxonomic_scope": str(scope),
                        "evidence_unit": str(unit),
                        "quantile_count": quantile_count,
                        "depth_bin": int(row.depth_bin),
                        "alt_bin": int(row.alt_bin),
                        "depth_label": depth_labels[int(row.depth_bin)],
                        "alt_label": alt_labels[int(row.alt_bin)],
                        "gnomad_found_count": found,
                        "gnomad_eligible_count": eligible,
                        "gnomad_found_fraction": found / eligible,
                    }
                )
    return True, pd.DataFrame(cells, columns=ORTHOLOG_EVIDENCE_COLUMNS), distributions


def _expand_strategy_masks(
    frame: pd.DataFrame,
    strategies: tuple[str, ...],
) -> pd.DataFrame:
    if frame.empty:
        return frame.assign(strategy=pd.Series(dtype="object"))
    rows = []
    for row in frame.to_dict(orient="records"):
        mask = int(row.pop("strategy_mask"))
        for index, strategy in enumerate(strategies):
            if mask & (1 << index):
                rows.append({**row, "strategy": strategy})
    return pd.DataFrame(rows)


def _ortholog_evidence_from_grouped(
    grouped: pd.DataFrame,
    distribution_source: pd.DataFrame,
) -> tuple[bool, pd.DataFrame, pd.DataFrame]:
    empty = pd.DataFrame(columns=ORTHOLOG_EVIDENCE_COLUMNS)
    distributions = _ortholog_evidence_distributions(
        distribution_source.assign(taxonomic_scope="all", evidence_unit="ortholog"),
        count_column="variant_count",
        site_column="site_depth",
        alt_column="alt_support",
    )
    if grouped.empty:
        return not distribution_source.empty, empty, distributions
    cells = []
    for (strategy, context), subset in grouped.groupby(
        ["strategy", "target_context"], sort=True
    ):
        for quantile_count in (2, 4, 10):
            depth_bins = _weighted_quantile_bins(
                subset["site_depth"], subset["gnomad_eligible_count"], quantile_count
            )
            alt_bins = _weighted_quantile_bins(
                subset["alt_support"], subset["gnomad_eligible_count"], quantile_count
            )
            depth_labels = _bin_labels(
                subset["site_depth"], depth_bins, quantile_count, percent=False
            )
            alt_labels = _bin_labels(
                subset["alt_support"], alt_bins, quantile_count, percent=False
            )
            binned = subset.assign(depth_bin=depth_bins, alt_bin=alt_bins)
            aggregated = binned.groupby(["depth_bin", "alt_bin"], as_index=False)[
                ["gnomad_found_count", "gnomad_eligible_count"]
            ].sum()
            for row in aggregated.itertuples(index=False):
                eligible = int(row.gnomad_eligible_count)
                found = int(row.gnomad_found_count)
                cells.append(
                    {
                        "strategy": str(strategy),
                        "target_context": str(context),
                        "taxonomic_scope": "all",
                        "evidence_unit": "ortholog",
                        "quantile_count": quantile_count,
                        "depth_bin": int(row.depth_bin),
                        "alt_bin": int(row.alt_bin),
                        "depth_label": depth_labels[int(row.depth_bin)],
                        "alt_label": alt_labels[int(row.alt_bin)],
                        "gnomad_found_count": found,
                        "gnomad_eligible_count": eligible,
                        "gnomad_found_fraction": found / eligible,
                    }
                )
    return True, pd.DataFrame(cells, columns=ORTHOLOG_EVIDENCE_COLUMNS), distributions


def _summary_from_grouped_aggregation(
    grouped: VariantGroupedAggregation,
    strategy_label: Callable[[str], str],
    variant_strategy_support_path: Path | None,
    ortholog_evidence_summary_path: Path | None,
) -> VariantSummary:
    strategies = list(grouped.masks.strategies)
    global_rows = _expand_strategy_masks(grouped.global_groups, grouped.masks.strategies)
    allele_gene_rows = _expand_strategy_masks(
        grouped.allele_gene_groups,
        grouped.masks.strategies,
    )
    strategy_counts = grouped.masks.strategy_counts()

    stats_rows = []
    af_by_strategy = (
        grouped.gnomad_af_summary.set_index("strategy").to_dict(orient="index")
        if not grouped.gnomad_af_summary.empty
        else {}
    )
    for strategy in strategies:
        rows = global_rows[global_rows["strategy"].eq(strategy)]
        gene_rows = allele_gene_rows[allele_gene_rows["strategy"].eq(strategy)]
        total = int(strategy_counts[strategy])
        ti = int(rows.loc[rows["titv_kind"].eq("ti"), "variant_count"].sum())
        tv = int(rows.loc[rows["titv_kind"].eq("tv"), "variant_count"].sum())
        gnomad_found = int(
            rows.loc[rows["gnomad_status"].eq("found"), "variant_count"].sum()
        )
        gnomad_eligible = int(
            rows.loc[rows["gnomad_status"].isin(["found", "not_found"]), "variant_count"].sum()
        )
        clinvar_found = int(
            rows.loc[rows["clinvar_found"].astype(bool), "variant_count"].sum()
        )
        clinvar_classified = int(
            rows.loc[rows["clinvar_classified"].astype(bool), "variant_count"].sum()
        )
        af = af_by_strategy.get(strategy, {})
        stats_rows.append(
            {
                "Strategy": strategy_label(strategy),
                "Unique Variants": total,
                "Genes": int(gene_rows["gene_id"].nunique()),
                "Ti/Tv": np.nan if ti == 0 and tv == 0 else (float("inf") if tv == 0 else round(ti / tv, 3)),
                "Found in ClinVar": clinvar_found,
                "ClinVar found %": clinvar_found / total if total else np.nan,
                "ClinVar classified": clinvar_classified,
                "ClinVar classified %": clinvar_classified / total if total else np.nan,
                "gnomAD Found": gnomad_found,
                "gnomAD Eligible": gnomad_eligible,
                "gnomAD lookup failed": int(
                    rows.loc[rows["gnomad_status"].eq("lookup_failed"), "variant_count"].sum()
                ),
                "gnomAD found %": gnomad_found / gnomad_eligible if gnomad_eligible else np.nan,
                "P/LP": int(rows.loc[rows["clinvar_category"].eq("P/LP"), "variant_count"].sum()),
                "B/LB": int(rows.loc[rows["clinvar_category"].eq("B/LB"), "variant_count"].sum()),
                "VUS": int(rows.loc[rows["clinvar_category"].eq("VUS"), "variant_count"].sum()),
                "Other ClinVar": int(rows.loc[rows["clinvar_category"].eq("Other"), "variant_count"].sum()),
                "Median gnomAD AF": af.get("Median gnomAD AF", np.nan),
            }
        )
    strategy_stats = pd.DataFrame(stats_rows).sort_values("Strategy").reset_index(drop=True)

    totals = grouped.masks.strategy_counts()
    ordered = sorted(strategies, key=lambda strategy: (-totals[strategy], strategy))
    if len(ordered) < 2:
        overlap = None
    else:
        lexical_intersections = np.asarray(grouped.masks.intersections(), dtype=np.int64)
        indices = [strategies.index(strategy) for strategy in ordered]
        intersections = lexical_intersections[np.ix_(indices, indices)]
        ordered_totals = np.asarray([totals[strategy] for strategy in ordered], dtype=np.int64)
        unions = ordered_totals[:, None] + ordered_totals[None, :] - intersections
        jaccard = np.divide(
            intersections,
            unions,
            out=np.zeros_like(intersections, dtype=float),
            where=unions != 0,
        )
        overlap = StrategyOverlap(ordered, intersections, unions, jaccard)

    def grouped_counts(
        frame: pd.DataFrame,
        value_columns: list[str],
        *,
        where: pd.Series | None = None,
    ) -> pd.DataFrame:
        selected = frame if where is None else frame.loc[where]
        if selected.empty:
            return pd.DataFrame(columns=["strategy", *value_columns, "Variant_Count"])
        result = (
            selected.groupby(["strategy", *value_columns], as_index=False, sort=True)[
                "variant_count"
            ]
            .sum()
            .rename(columns={"variant_count": "Variant_Count"})
        )
        result["strategy"] = result["strategy"].map(strategy_label)
        return result

    event_counts = grouped_counts(global_rows, ["event_type"])
    target_context_counts = grouped_counts(allele_gene_rows, ["target_context"])
    gene_variant_counts = grouped_counts(allele_gene_rows, ["gene_id"]).sort_values(
        ["strategy", "Variant_Count", "gene_id"],
        ascending=[True, False, True],
    )
    gnomad_event_counts = grouped_counts(
        global_rows,
        ["gnomad_status", "event_type"],
        where=global_rows["gnomad_status"].isin(["found", "not_found"]),
    )
    gnomad_context_counts = grouped_counts(
        allele_gene_rows,
        ["gnomad_status", "target_context"],
        where=allele_gene_rows["gnomad_status"].isin(["found", "not_found"]),
    )
    clinvar_counts = grouped_counts(global_rows, ["clinvar_category"])
    star_counts = grouped_counts(
        global_rows,
        ["review_stars"],
        where=global_rows["clinvar_category"].eq("P/LP"),
    ).rename(columns={"strategy": "Strategy", "review_stars": "Review stars"})
    consequence_counts = grouped_counts(
        allele_gene_rows,
        ["consequence"],
        where=allele_gene_rows["consequence"].ne(""),
    ).rename(columns={"consequence": "value"})
    pathogenic_consequence_counts = grouped_counts(
        allele_gene_rows,
        ["consequence"],
        where=(
            allele_gene_rows["consequence"].ne("")
            & allele_gene_rows["clinvar_category"].eq("P/LP")
        ),
    ).rename(columns={"consequence": "value"})

    af_summary = grouped.gnomad_af_summary.drop(columns=["Median gnomAD AF"], errors="ignore").copy()
    if not af_summary.empty:
        af_summary["Strategy"] = af_summary.pop("strategy").map(strategy_label)
        af_summary = af_summary[["Strategy", "Count", "Q05", "Q25", "Median", "Q75", "Q95"]]
    unique_contribution = pd.DataFrame(
        {
            "Strategy": [strategy_label(strategy) for strategy in strategies],
            "Unique To Strategy": [
                grouped.masks.unique_strategy_counts()[strategy] for strategy in strategies
            ],
        }
    ).sort_values("Strategy")

    if ortholog_evidence_summary_path is not None:
        ortholog_available, ortholog_cells, ortholog_distributions = read_taxonomic_ortholog_evidence(
            ortholog_evidence_summary_path
        )
    else:
        ortholog_available, ortholog_cells, ortholog_distributions = (
            _ortholog_evidence_from_grouped(
                grouped.ortholog_evidence_grouped,
                grouped.ortholog_distribution_source,
            )
        )

    pathogenic_rows = _add_pathogenic_strategy_support(
        variant_strategy_support_path,
        grouped.pathogenic_rows,
    )
    return VariantSummary(
        input_row_count=grouped.masks.input_row_count,
        unique_variant_count=grouped.masks.unique_variant_count,
        all_strategy_variant_count=grouped.masks.all_strategy_variant_count,
        strategy_record_count=grouped.masks.strategy_record_count,
        gene_count=grouped.gene_count,
        clinvar_found=int(grouped.global_groups.loc[grouped.global_groups["clinvar_found"].astype(bool), "variant_count"].sum()),
        clinvar_classified=int(grouped.global_groups.loc[grouped.global_groups["clinvar_classified"].astype(bool), "variant_count"].sum()),
        gnomad_found=int(grouped.global_groups.loc[grouped.global_groups["gnomad_status"].eq("found"), "variant_count"].sum()),
        gnomad_lookup_failed=int(grouped.global_groups.loc[grouped.global_groups["gnomad_status"].eq("lookup_failed"), "variant_count"].sum()),
        pathogenic_variant_count=int(grouped.global_groups.loc[grouped.global_groups["clinvar_category"].eq("P/LP"), "variant_count"].sum()),
        consequence_source=grouped.consequence_source,
        strategies=strategies,
        strategy_stats=strategy_stats,
        unique_contribution=unique_contribution,
        gene_variant_counts=gene_variant_counts,
        event_counts=event_counts,
        target_context_counts=target_context_counts,
        gnomad_event_counts=gnomad_event_counts,
        gnomad_context_counts=gnomad_context_counts,
        overlap=overlap,
        clinvar_counts=clinvar_counts,
        gnomad_af_summary=af_summary,
        pathogenic_star_counts=star_counts,
        consequence_counts=consequence_counts,
        pathogenic_consequence_counts=pathogenic_consequence_counts,
        pathogenic_rows=pathogenic_rows,
        ortholog_evidence_available=ortholog_available,
        ortholog_evidence_cells=ortholog_cells,
        ortholog_evidence_distributions=ortholog_distributions,
    )


def build_variant_summary(
    path: Path,
    work_dir: Path,
    strategy_label: Callable[[str], str],
    target_features_path: Path | None = None,
    genes_path: Path | None = None,
    annotation_failures_path: Path | None = None,
    variant_strategy_support_path: Path | None = None,
    ortholog_evidence_summary_path: Path | None = None,
    chunk_size: int = 100_000,
    performance_profile: PerformanceProfile | None = None,
) -> VariantSummary:
    """Aggregate a variant annotation table without retaining row-level data in memory."""
    work_dir.mkdir(parents=True, exist_ok=True)
    cache_path = work_dir / SUMMARY_CACHE_NAME
    cache_stage = (
        performance_profile.stage("Variant summary cache lookup")
        if performance_profile is not None
        else nullcontext()
    )
    with cache_stage:
        cached = _load_summary_cache(
            cache_path,
            path,
            target_features_path,
            genes_path,
            annotation_failures_path,
            variant_strategy_support_path,
            ortholog_evidence_summary_path,
            strategy_label,
        )
    if cached is not None:
        return cached

    aggregation_stage = (
        performance_profile.stage("Variant summary aggregation")
        if performance_profile is not None
        else nullcontext()
    )
    with aggregation_stage:
        summary = _compute_variant_summary(
            path,
            work_dir,
            target_features_path,
            genes_path,
            annotation_failures_path,
            variant_strategy_support_path,
            ortholog_evidence_summary_path,
            strategy_label,
            chunk_size,
            performance_profile,
        )
        if performance_profile is not None:
            performance_profile.add_metric("input_rows", summary.input_row_count)
            performance_profile.add_metric(
                "strategy_membership_rows", summary.strategy_record_count
            )

    write_stage = (
        performance_profile.stage("Variant summary cache write")
        if performance_profile is not None
        else nullcontext()
    )
    with write_stage:
        _write_summary_cache(
            cache_path,
            summary,
            path,
            target_features_path,
            genes_path,
            annotation_failures_path,
            variant_strategy_support_path,
            ortholog_evidence_summary_path,
            strategy_label,
        )
    return summary


def _compute_variant_summary(
    path: Path,
    work_dir: Path,
    target_features_path: Path | None,
    genes_path: Path | None,
    annotation_failures_path: Path | None,
    variant_strategy_support_path: Path | None,
    ortholog_evidence_summary_path: Path | None,
    strategy_label: Callable[[str], str],
    chunk_size: int,
    performance_profile: PerformanceProfile | None,
) -> VariantSummary:
    del chunk_size
    source = resolve_variant_aggregation_source(path)
    with tempfile.TemporaryDirectory(
        prefix=".variant_summary_duckdb.",
        dir=work_dir,
    ) as temporary_dir:
        grouped = aggregate_variant_groups(
            source,
            genes_path=genes_path,
            target_features_path=target_features_path,
            annotation_failures_path=annotation_failures_path,
            variant_strategy_support_path=(
                None if ortholog_evidence_summary_path is not None else variant_strategy_support_path
            ),
            temp_dir=Path(temporary_dir),
        )
    if performance_profile is not None:
        for name, seconds in grouped.timings.items():
            performance_profile.add_metric(f"duckdb_{name}_seconds", seconds)
    return _summary_from_grouped_aggregation(
        grouped,
        strategy_label,
        variant_strategy_support_path,
        ortholog_evidence_summary_path,
    )
