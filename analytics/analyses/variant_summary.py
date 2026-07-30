"""Bounded-memory aggregation of variant annotations for the HTML report."""

from __future__ import annotations

import gzip
import json
import math
import sqlite3
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from genomics.clinvar import CLINVAR_CLASS_ORDER
from .target_context import context_at, read_disjoint_contexts
from genomics.variants import (
    RegionIndex,
    changed_target_position,
    gnomad_lookup_status,
    parse_variant_key,
    read_failed_regions,
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
VARIANT_REQUIRED = {"variant_key", "gene_id", "event_type", "strategies"}
SUMMARY_CACHE_VERSION = 10
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
        "input": _input_metadata(source),
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
    strategy_label: Callable[[str], str],
) -> VariantSummary | None:
    if not cache_path.exists():
        return None
    try:
        with gzip.open(cache_path, "rt") as handle:
            payload = json.load(handle)
        if payload.get("cache_version") != SUMMARY_CACHE_VERSION:
            return None
        if payload.get("input") != _input_metadata(source):
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
    strategy_label: Callable[[str], str],
) -> None:
    payload = _summary_payload(
        summary,
        source,
        target_features,
        genes,
        annotation_failures,
        variant_strategy_support,
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


def _header(path: Path) -> list[str]:
    compression = "gzip" if path.suffix == ".gz" else None
    return pd.read_csv(path, sep="\t", compression=compression, nrows=0).columns.tolist()


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


def _normalize_chunk(
    df: pd.DataFrame,
    gene_begins: dict[str, int] | None = None,
    gnomad_failed_regions: RegionIndex | None = None,
) -> pd.DataFrame:
    has_vep_consequences = {"vep_status", "vep_primary_consequence"}.issubset(df.columns)
    for column in VARIANT_USECOLS:
        if column not in df.columns:
            df[column] = ""

    df["variant_id"] = df["variant_key"].astype(str)
    missing_key = df["variant_id"].eq("")
    if missing_key.any():
        df.loc[missing_key, "variant_id"] = (
            df.loc[missing_key, "gene_id"].astype(str)
            + ":"
            + df.loc[missing_key, "event_type"].astype(str)
            + ":"
            + df.loc[missing_key, "ref"].astype(str)
            + ">"
            + df.loc[missing_key, "alt"].astype(str)
        )

    df["gnomad_af"] = pd.to_numeric(df["gnomad_af"], errors="coerce")
    parsed_keys = df["variant_key"].map(parse_variant_key)
    failed_regions = gnomad_failed_regions or {}
    df["gnomad_status"] = [
        gnomad_lookup_status(
            key=key,
            lookup_status=status,
            found=not pd.isna(af),
            failed_regions=failed_regions,
        )
        for key, status, af in zip(parsed_keys, df["lookup_status"], df["gnomad_af"])
    ]
    if gene_begins is not None:
        df["target_start0"] = [
            changed_target_position(key, gene_begins[str(gene_id)])
            if key is not None and str(gene_id) in gene_begins
            else pd.NA
            for key, gene_id in zip(parsed_keys, df["gene_id"])
        ]
    elif "target_start0" in df.columns:
        df["target_start0"] = pd.to_numeric(df["target_start0"], errors="coerce")
    else:
        df["target_start0"] = pd.NA
    for column in ["support_row_count", "support_ortholog_count", "clinvar_scv_count"]:
        df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0).astype("int64")
    clinvar_evidence_columns = [
        "clinvar_id",
        "clinvar_allele_id",
        "clinvar_sig",
        "clinvar_revstat",
        "clinvar_hgvs",
        "clinvar_disease",
        "clinvar_variant_type",
    ]
    df["clinvar_found"] = df[clinvar_evidence_columns].fillna("").astype(str).ne("").any(axis=1)
    df["clinvar_classified"] = df["clinvar_sig"].astype(str) != ""
    df["clinvar_category"] = _categorize_clinvar(df["clinvar_sig"], df["clinvar_found"])
    if has_vep_consequences:
        df["consequence"] = df["vep_primary_consequence"].where(
            df["vep_status"].astype(str).eq("ok"),
            "",
        )
    else:
        df["consequence"] = df["gnomad_csq"].astype(str)

    ref = df["ref"].astype(str).str.upper()
    alt = df["alt"].astype(str).str.upper()
    valid_snv = df["event_type"].astype(str).eq("snv") & ref.str.len().eq(1) & alt.str.len().eq(1)
    transition = ref.str.cat(alt, sep=">").isin(["A>G", "G>A", "C>T", "T>C"])
    df["titv_kind"] = ""
    df.loc[valid_snv & transition, "titv_kind"] = "ti"
    df.loc[valid_snv & ~transition, "titv_kind"] = "tv"
    df["review_stars"] = df["clinvar_review_stars"].astype(str).where(
        df["clinvar_review_stars"].astype(str).isin(["0", "1", "2", "3", "4"]),
        "Unmapped",
    )
    return df


def _create_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.executescript(
        """
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous = OFF;
        PRAGMA temp_store = FILE;
        CREATE TABLE memberships (
            variant_id TEXT NOT NULL,
            strategy TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            event_type TEXT NOT NULL,
            target_context TEXT NOT NULL,
            clinvar_found INTEGER NOT NULL,
            clinvar_classified INTEGER NOT NULL,
            clinvar_category TEXT NOT NULL,
            gnomad_af REAL,
            gnomad_status TEXT NOT NULL,
            titv_kind TEXT NOT NULL,
            review_stars TEXT NOT NULL,
            consequence TEXT NOT NULL,
            PRIMARY KEY (variant_id, strategy)
        ) WITHOUT ROWID;
        """
    )
    return connection


def _keep_top_pathogenic(current: pd.DataFrame, chunk: pd.DataFrame, limit: int = 100) -> pd.DataFrame:
    pathogenic = chunk[chunk["clinvar_category"] == "P/LP"].copy()
    if pathogenic.empty:
        return current
    pathogenic["_star_rank"] = pd.to_numeric(pathogenic["review_stars"], errors="coerce").fillna(-1)
    keep = pd.concat([current, pathogenic], ignore_index=True) if not current.empty else pathogenic
    keep = keep.sort_values(
        ["_star_rank", "clinvar_scv_count", "variant_id"],
        ascending=[False, False, True],
        kind="mergesort",
    ).drop_duplicates("variant_id", keep="first")
    return keep.head(limit)


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


def _strategy_stats(connection: sqlite3.Connection, strategy_label: Callable[[str], str]) -> pd.DataFrame:
    stats = pd.read_sql_query(
        """
        SELECT strategy,
               COUNT(*) AS "Unique Variants",
               COUNT(DISTINCT gene_id) AS "Genes",
               SUM(titv_kind = 'ti') AS ti,
               SUM(titv_kind = 'tv') AS tv,
               SUM(clinvar_found) AS "Found in ClinVar",
               SUM(clinvar_classified) AS "ClinVar classified",
               SUM(gnomad_status = 'found') AS "gnomAD Found",
               SUM(gnomad_status IN ('found', 'not_found')) AS "gnomAD Eligible",
               SUM(gnomad_status = 'lookup_failed') AS "gnomAD lookup failed",
               SUM(clinvar_category = 'P/LP') AS "P/LP",
               SUM(clinvar_category = 'B/LB') AS "B/LB",
               SUM(clinvar_category = 'VUS') AS "VUS",
               SUM(clinvar_category = 'Other') AS "Other ClinVar"
        FROM memberships
        GROUP BY strategy
        ORDER BY strategy
        """,
        connection,
    )
    if stats.empty:
        return pd.DataFrame()
    medians = pd.read_sql_query(
        """
        WITH ranked AS (
            SELECT strategy, gnomad_af,
                   ROW_NUMBER() OVER (PARTITION BY strategy ORDER BY gnomad_af) AS row_number,
                   COUNT(*) OVER (PARTITION BY strategy) AS row_count
            FROM memberships
            WHERE gnomad_af IS NOT NULL
        )
        SELECT strategy, AVG(gnomad_af) AS "Median gnomAD AF"
        FROM ranked
        WHERE row_number IN ((row_count + 1) / 2, (row_count + 2) / 2)
        GROUP BY strategy
        """,
        connection,
    )
    stats = stats.merge(medians, on="strategy", how="left")
    stats["Ti/Tv"] = [
        np.nan if tv == 0 and ti == 0 else (float("inf") if tv == 0 else round(ti / tv, 3))
        for ti, tv in zip(stats.pop("ti"), stats.pop("tv"))
    ]
    variant_denominator = stats["Unique Variants"].replace(0, np.nan)
    gnomad_denominator = stats["gnomAD Eligible"].replace(0, np.nan)
    stats["ClinVar found %"] = stats["Found in ClinVar"] / variant_denominator
    stats["ClinVar classified %"] = stats["ClinVar classified"] / variant_denominator
    stats["gnomAD found %"] = stats["gnomAD Found"] / gnomad_denominator
    stats["Strategy"] = stats.pop("strategy").map(strategy_label)
    columns = [
        "Strategy",
        "Unique Variants",
        "Genes",
        "Ti/Tv",
        "Found in ClinVar",
        "ClinVar found %",
        "ClinVar classified",
        "ClinVar classified %",
        "gnomAD Found",
        "gnomAD Eligible",
        "gnomAD lookup failed",
        "gnomAD found %",
        "P/LP",
        "B/LB",
        "VUS",
        "Other ClinVar",
        "Median gnomAD AF",
    ]
    return stats[columns].sort_values("Strategy").reset_index(drop=True)


def _grouped_counts(
    connection: sqlite3.Connection,
    value_column: str,
    strategy_label: Callable[[str], str],
    where: str = "",
) -> pd.DataFrame:
    condition = f"WHERE {where}" if where else ""
    counts = pd.read_sql_query(
        f"""
        SELECT strategy, {value_column} AS value, COUNT(*) AS Variant_Count
        FROM memberships
        {condition}
        GROUP BY strategy, {value_column}
        ORDER BY strategy, {value_column}
        """,
        connection,
    )
    if not counts.empty:
        counts["strategy"] = counts["strategy"].map(strategy_label)
    return counts


def _overlap_summary(
    connection: sqlite3.Connection,
    strategies: list[str],
    totals: dict[str, int],
) -> tuple[int, dict[str, int], int, StrategyOverlap | None]:
    intersections: Counter[tuple[str, str]] = Counter()
    unique_counts: Counter[str] = Counter()
    unique_variant_count = 0
    all_strategy_variant_count = 0
    current_variant = None
    current_strategies: list[str] = []

    def consume_group(group: list[str]) -> None:
        nonlocal all_strategy_variant_count, unique_variant_count
        if not group:
            return
        unique_variant_count += 1
        if len(group) == len(strategies):
            all_strategy_variant_count += 1
        if len(group) == 1:
            unique_counts[group[0]] += 1
        intersections.update(combinations_with_replacement(group, 2))

    cursor = connection.execute("SELECT variant_id, strategy FROM memberships ORDER BY variant_id, strategy")
    for variant_id, strategy in cursor:
        if current_variant is not None and variant_id != current_variant:
            consume_group(current_strategies)
            current_strategies = []
        current_variant = variant_id
        current_strategies.append(strategy)
    consume_group(current_strategies)

    ordered = sorted(strategies, key=lambda strategy: (-totals.get(strategy, 0), strategy))
    if len(ordered) < 2:
        return unique_variant_count, dict(unique_counts), all_strategy_variant_count, None
    size = len(ordered)
    shared = np.zeros((size, size), dtype=np.int64)
    unions = np.zeros((size, size), dtype=np.int64)
    jaccard = np.zeros((size, size), dtype=float)
    for row_index, row_strategy in enumerate(ordered):
        for col_index, col_strategy in enumerate(ordered):
            key = tuple(sorted((row_strategy, col_strategy)))
            count = intersections[key]
            union = totals[row_strategy] + totals[col_strategy] - count
            shared[row_index, col_index] = count
            unions[row_index, col_index] = union
            jaccard[row_index, col_index] = count / union if union else 0.0
    return (
        unique_variant_count,
        dict(unique_counts),
        all_strategy_variant_count,
        StrategyOverlap(ordered, shared, unions, jaccard),
    )


def _gnomad_af_summary(
    connection: sqlite3.Connection,
    strategy_label: Callable[[str], str],
) -> pd.DataFrame:
    rows = []
    strategies = [row[0] for row in connection.execute(
        "SELECT DISTINCT strategy FROM memberships WHERE gnomad_af > 0 ORDER BY strategy"
    )]
    for strategy in strategies:
        values = np.asarray(
            [row[0] for row in connection.execute(
                "SELECT gnomad_af FROM memberships WHERE strategy = ? AND gnomad_af > 0",
                (strategy,),
            )],
            dtype=float,
        )
        if values.size == 0:
            continue
        quantiles = np.quantile(np.log10(values), [0.05, 0.25, 0.5, 0.75, 0.95])
        rows.append(
            {
                "Strategy": strategy_label(strategy),
                "Count": int(values.size),
                "Q05": float(quantiles[0]),
                "Q25": float(quantiles[1]),
                "Median": float(quantiles[2]),
                "Q75": float(quantiles[3]),
                "Q95": float(quantiles[4]),
            }
        )
    return pd.DataFrame(rows)


def _gene_lengths_from_features(path: Path) -> dict[str, int]:
    frame = pd.read_csv(
        path,
        sep="\t",
        compression="gzip" if path.suffix == ".gz" else None,
        keep_default_na=False,
        usecols=["gene_id", "feature_type", "target_start0", "target_end0"],
    )
    genes = frame[frame["feature_type"].astype(str).str.lower().eq("gene")]
    if genes.empty:
        raise ValueError("Target features contain no gene rows for target-context assignment.")
    return {
        str(row.gene_id): int(row.target_end0) - int(row.target_start0)
        for row in genes.itertuples(index=False)
    }


def _gene_begins(path: Path) -> dict[str, int]:
    frame = pd.read_csv(
        path,
        sep="\t",
        compression="gzip" if path.suffix == ".gz" else None,
        keep_default_na=False,
        usecols=["gene_id", "begin"],
    )
    return {
        str(row.gene_id): int(row.begin)
        for row in frame.itertuples(index=False)
    }


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


def _ortholog_evidence_summary(
    connection: sqlite3.Connection,
    support_path: Path | None,
    chunk_size: int,
) -> tuple[bool, pd.DataFrame, pd.DataFrame]:
    empty = pd.DataFrame(columns=ORTHOLOG_EVIDENCE_COLUMNS)
    empty_distributions = pd.DataFrame(columns=ORTHOLOG_EVIDENCE_DISTRIBUTION_COLUMNS)
    if support_path is None or "site_aligned_ortholog_count" not in _header(support_path):
        return False, empty, empty_distributions

    connection.executescript(
        """
        CREATE TABLE ortholog_support (
            variant_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            strategy TEXT NOT NULL,
            alt_support INTEGER NOT NULL,
            site_depth INTEGER NOT NULL,
            PRIMARY KEY (variant_id, gene_id, strategy)
        ) WITHOUT ROWID;
        """
    )
    columns = [
        "variant_key",
        "gene_id",
        "strategy",
        "alt_support_ortholog_count",
        "site_aligned_ortholog_count",
    ]
    insert_sql = "INSERT INTO ortholog_support VALUES (?, ?, ?, ?, ?)"
    for chunk in pd.read_csv(
        support_path,
        sep="\t",
        compression="gzip" if support_path.suffix == ".gz" else None,
        usecols=columns,
        keep_default_na=False,
        chunksize=chunk_size,
    ):
        alt = pd.to_numeric(chunk["alt_support_ortholog_count"], errors="coerce")
        depth = pd.to_numeric(chunk["site_aligned_ortholog_count"], errors="coerce")
        keep = depth.notna() & depth.gt(0) & alt.notna()
        chunk = chunk.loc[keep].copy()
        chunk["alt_support_ortholog_count"] = alt.loc[keep].astype("int64")
        chunk["site_aligned_ortholog_count"] = depth.loc[keep].astype("int64")
        invalid = (
            chunk["alt_support_ortholog_count"].lt(0)
            | chunk["alt_support_ortholog_count"].gt(chunk["site_aligned_ortholog_count"])
        )
        if invalid.any():
            row = chunk.loc[invalid].iloc[0]
            raise ValueError(
                "Invalid ortholog evidence for "
                f"{row['variant_key']} / {row['strategy']}: "
                f"ALT={row['alt_support_ortholog_count']}, depth={row['site_aligned_ortholog_count']}"
            )
        connection.executemany(
            insert_sql,
            chunk[columns].itertuples(index=False, name=None),
        )
    connection.commit()

    grouped = pd.read_sql_query(
        """
        SELECT m.strategy,
               m.target_context,
               s.site_depth,
               s.alt_support,
               SUM(m.gnomad_status = 'found') AS gnomad_found_count,
               COUNT(*) AS gnomad_eligible_count
        FROM memberships AS m
        JOIN ortholog_support AS s
          ON s.variant_id = m.variant_id
         AND s.gene_id = m.gene_id
         AND s.strategy = m.strategy
        WHERE m.event_type = 'snv'
          AND m.target_context IN ('cds', 'utr', 'intron')
          AND m.gnomad_status IN ('found', 'not_found')
        GROUP BY m.strategy, m.target_context, s.site_depth, s.alt_support
        """,
        connection,
    )
    distribution_source = pd.read_sql_query(
        """
        SELECT m.strategy,
               'all' AS taxonomic_scope,
               'ortholog' AS evidence_unit,
               s.site_depth,
               s.alt_support,
               COUNT(*) AS variant_count
        FROM memberships AS m
        JOIN ortholog_support AS s
          ON s.variant_id = m.variant_id
         AND s.gene_id = m.gene_id
         AND s.strategy = m.strategy
        WHERE m.event_type = 'snv'
          AND m.target_context IN ('cds', 'utr', 'intron')
        GROUP BY m.strategy, s.site_depth, s.alt_support
        """,
        connection,
    )
    distributions = _ortholog_evidence_distributions(
        distribution_source,
        count_column="variant_count",
        site_column="site_depth",
        alt_column="alt_support",
    )
    if grouped.empty:
        return True, empty, distributions
    cells = []
    for (strategy, context), subset in grouped.groupby(
        ["strategy", "target_context"], sort=True
    ):
        subset = subset.copy()
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
            binned = subset.assign(
                depth_bin=depth_bins,
                alt_bin=alt_bins,
            )
            aggregated = binned.groupby(
                ["depth_bin", "alt_bin"], as_index=False
            )[["gnomad_found_count", "gnomad_eligible_count"]].sum()
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


def build_variant_summary(
    path: Path,
    work_dir: Path,
    strategy_label: Callable[[str], str],
    target_features_path: Path | None = None,
    genes_path: Path | None = None,
    annotation_failures_path: Path | None = None,
    variant_strategy_support_path: Path | None = None,
    chunk_size: int = 100_000,
) -> VariantSummary:
    """Aggregate a variant annotation table without retaining row-level data in memory."""
    work_dir.mkdir(parents=True, exist_ok=True)
    cache_path = work_dir / SUMMARY_CACHE_NAME
    cached = _load_summary_cache(
        cache_path,
        path,
        target_features_path,
        genes_path,
        annotation_failures_path,
        variant_strategy_support_path,
        strategy_label,
    )
    if cached is not None:
        return cached

    summary = _compute_variant_summary(
        path,
        work_dir,
        target_features_path,
        genes_path,
        annotation_failures_path,
        variant_strategy_support_path,
        strategy_label,
        chunk_size,
    )
    _write_summary_cache(
        cache_path,
        summary,
        path,
        target_features_path,
        genes_path,
        annotation_failures_path,
        variant_strategy_support_path,
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
    strategy_label: Callable[[str], str],
    chunk_size: int,
) -> VariantSummary:
    header = _header(path)
    missing = VARIANT_REQUIRED - set(header)
    if missing:
        raise ValueError(f"Variant annotations missing required columns: {', '.join(sorted(missing))}")
    usecols = [column for column in VARIANT_USECOLS if column in header]
    has_vep_consequences = {"vep_status", "vep_primary_consequence"}.issubset(header)
    if has_vep_consequences:
        usecols.extend(VEP_USECOLS)
    consequence_source = "Ensembl VEP" if has_vep_consequences else "gnomAD CSQ (legacy)"
    if "target_start0" in header:
        usecols.append("target_start0")
    contexts: dict[str, list[tuple[int, int, str]]] = {}
    context_starts: dict[str, list[int]] = {}
    gene_begins = _gene_begins(genes_path) if genes_path is not None else None
    gnomad_failed_regions = read_failed_regions(annotation_failures_path, "gnomad")
    if target_features_path is not None:
        if gene_begins is None and "target_start0" not in header:
            raise ValueError(
                "Target-context reporting needs genes.tsv.gz or a legacy target_start0 column."
            )
        contexts = read_disjoint_contexts(
            target_features_path,
            _gene_lengths_from_features(target_features_path),
        )
        context_starts = {
            gene_id: [interval[0] for interval in intervals]
            for gene_id, intervals in contexts.items()
        }
    with tempfile.NamedTemporaryFile(
        prefix=".strategy_report_variants.",
        suffix=".sqlite3",
        dir=work_dir,
        delete=False,
    ) as temporary_database:
        database_path = Path(temporary_database.name)

    input_row_count = 0
    clinvar_found = 0
    clinvar_classified = 0
    gnomad_found = 0
    gnomad_lookup_failed = 0
    genes: set[str] = set()
    pathogenic_rows = pd.DataFrame(
        columns=[*VARIANT_USECOLS, *VEP_USECOLS, "variant_id", "clinvar_category", "review_stars"]
    )

    connection = _create_database(database_path)
    try:
        reader = pd.read_csv(
            path,
            sep="\t",
            compression="gzip" if path.suffix == ".gz" else None,
            usecols=usecols,
            keep_default_na=False,
            low_memory=False,
            chunksize=chunk_size,
        )
        insert_sql = """
            INSERT OR IGNORE INTO memberships (
                variant_id, strategy, gene_id, event_type, target_context, clinvar_found,
                clinvar_classified, clinvar_category, gnomad_af, gnomad_status, titv_kind,
                review_stars, consequence
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        for chunk in reader:
            chunk = _normalize_chunk(chunk, gene_begins, gnomad_failed_regions)
            chunk["target_context"] = [
                context_at(
                    contexts.get(str(gene_id), []),
                    int(position),
                    context_starts.get(str(gene_id), []),
                )
                if not pd.isna(position) else "unknown"
                for gene_id, position in zip(chunk["gene_id"], chunk["target_start0"])
            ]
            input_row_count += len(chunk)
            clinvar_found += int(chunk["clinvar_found"].sum())
            clinvar_classified += int(chunk["clinvar_classified"].sum())
            gnomad_found += int(chunk["gnomad_af"].notna().sum())
            gnomad_lookup_failed += int(chunk["gnomad_status"].eq("lookup_failed").sum())
            genes.update(chunk["gene_id"].astype(str))
            pathogenic_rows = _keep_top_pathogenic(pathogenic_rows, chunk)

            records = []
            record_columns = [
                "variant_id",
                "strategies",
                "gene_id",
                "event_type",
                "target_context",
                "clinvar_found",
                "clinvar_classified",
                "clinvar_category",
                "gnomad_af",
                "gnomad_status",
                "titv_kind",
                "review_stars",
                "consequence",
            ]
            for (
                variant_id,
                strategy_text,
                gene_id,
                event_type,
                target_context,
                clinvar_found_value,
                clinvar_classified_value,
                clinvar_category,
                gnomad_af,
                gnomad_status,
                titv_kind,
                review_stars,
                consequence,
            ) in chunk[record_columns].itertuples(index=False, name=None):
                strategies = [item.strip() for item in str(strategy_text).split(",") if item.strip()]
                for strategy in strategies:
                    records.append(
                        (
                            str(variant_id),
                            strategy,
                            str(gene_id),
                            str(event_type),
                            str(target_context),
                            int(clinvar_found_value),
                            int(clinvar_classified_value),
                            str(clinvar_category),
                            None if pd.isna(gnomad_af) else float(gnomad_af),
                            str(gnomad_status),
                            str(titv_kind),
                            str(review_stars),
                            str(consequence),
                        )
                    )
            connection.executemany(insert_sql, records)
            connection.commit()

        raw_stats = pd.read_sql_query(
            "SELECT strategy, COUNT(*) AS count FROM memberships GROUP BY strategy ORDER BY strategy",
            connection,
        )
        strategies = raw_stats["strategy"].astype(str).tolist()
        totals = dict(zip(raw_stats["strategy"], raw_stats["count"].astype(int)))
        unique_variant_count, unique_counts, all_strategy_variant_count, overlap = _overlap_summary(
            connection,
            strategies,
            totals,
        )
        strategy_stats = _strategy_stats(connection, strategy_label)
        unique_contribution = pd.DataFrame(
            {
                "Strategy": [strategy_label(strategy) for strategy in strategies],
                "Unique To Strategy": [unique_counts.get(strategy, 0) for strategy in strategies],
            }
        ).sort_values("Strategy")

        event_counts = _grouped_counts(connection, "event_type", strategy_label).rename(columns={"value": "event_type"})
        target_context_counts = _grouped_counts(connection, "target_context", strategy_label).rename(
            columns={"value": "target_context"}
        )
        gnomad_event_counts = pd.read_sql_query(
            """
            SELECT strategy, gnomad_status,
                   event_type, COUNT(*) AS Variant_Count
            FROM memberships
            WHERE gnomad_status IN ('found', 'not_found')
            GROUP BY strategy, gnomad_status, event_type
            """,
            connection,
        )
        gnomad_context_counts = pd.read_sql_query(
            """
            SELECT strategy, gnomad_status,
                   target_context, COUNT(*) AS Variant_Count
            FROM memberships
            WHERE gnomad_status IN ('found', 'not_found')
            GROUP BY strategy, gnomad_status, target_context
            """,
            connection,
        )
        for frame in (gnomad_event_counts, gnomad_context_counts):
            if not frame.empty:
                frame["strategy"] = frame["strategy"].map(strategy_label)
        clinvar_counts = _grouped_counts(connection, "clinvar_category", strategy_label).rename(
            columns={"value": "clinvar_category"}
        )
        star_counts = _grouped_counts(
            connection,
            "review_stars",
            strategy_label,
            "clinvar_category = 'P/LP'",
        ).rename(columns={"strategy": "Strategy", "value": "Review stars"})
        consequence_counts = _grouped_counts(
            connection,
            "consequence",
            strategy_label,
            "consequence != ''",
        )
        pathogenic_consequence_counts = _grouped_counts(
            connection,
            "consequence",
            strategy_label,
            "consequence != '' AND clinvar_category = 'P/LP'",
        )
        strategy_record_count = int(sum(totals.values()))
        pathogenic_variant_count = int(
            connection.execute(
                "SELECT COUNT(DISTINCT variant_id) FROM memberships "
                "WHERE clinvar_category = 'P/LP'"
            ).fetchone()[0]
        )
        gnomad_af_summary = _gnomad_af_summary(connection, strategy_label)
        (
            ortholog_evidence_available,
            ortholog_evidence_cells,
            ortholog_evidence_distributions,
        ) = _ortholog_evidence_summary(
            connection,
            variant_strategy_support_path,
            chunk_size,
        )
    finally:
        connection.close()
        database_path.unlink(missing_ok=True)

    pathogenic_rows = _add_pathogenic_strategy_support(
        variant_strategy_support_path,
        pathogenic_rows,
    )

    return VariantSummary(
        input_row_count=input_row_count,
        unique_variant_count=unique_variant_count,
        all_strategy_variant_count=all_strategy_variant_count,
        strategy_record_count=strategy_record_count,
        gene_count=len(genes),
        clinvar_found=clinvar_found,
        clinvar_classified=clinvar_classified,
        gnomad_found=gnomad_found,
        gnomad_lookup_failed=gnomad_lookup_failed,
        pathogenic_variant_count=pathogenic_variant_count,
        consequence_source=consequence_source,
        strategies=strategies,
        strategy_stats=strategy_stats,
        unique_contribution=unique_contribution,
        event_counts=event_counts,
        target_context_counts=target_context_counts,
        gnomad_event_counts=gnomad_event_counts,
        gnomad_context_counts=gnomad_context_counts,
        overlap=overlap,
        clinvar_counts=clinvar_counts,
        gnomad_af_summary=gnomad_af_summary,
        pathogenic_star_counts=star_counts,
        consequence_counts=consequence_counts,
        pathogenic_consequence_counts=pathogenic_consequence_counts,
        pathogenic_rows=pathogenic_rows.drop(columns=["_star_rank"], errors="ignore"),
        ortholog_evidence_available=ortholog_evidence_available,
        ortholog_evidence_cells=ortholog_evidence_cells,
        ortholog_evidence_distributions=ortholog_evidence_distributions,
    )
