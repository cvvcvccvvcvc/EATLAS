"""Bounded-memory aggregation of variant annotations for the HTML report."""

from __future__ import annotations

import gzip
import json
import math
import sqlite3
import tempfile
from collections import Counter
from dataclasses import dataclass
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .target_context import context_at, read_disjoint_contexts


VARIANT_USECOLS = [
    "variant_key",
    "gene_id",
    "target_start0",
    "event_type",
    "ref",
    "alt",
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
VARIANT_REQUIRED = {"variant_key", "gene_id", "event_type", "strategies"}
SUMMARY_CACHE_VERSION = 3
SUMMARY_CACHE_NAME = "variant_summary.json.gz"
SPECIAL_FLOAT_KEY = "__gaph_float__"


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
    ]
    return {
        "cache_version": SUMMARY_CACHE_VERSION,
        "input": _input_metadata(source),
        "target_features": _input_metadata(target_features) if target_features is not None else None,
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
        cache_hit=True,
    )


def _load_summary_cache(
    cache_path: Path,
    source: Path,
    target_features: Path | None,
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
    strategy_label: Callable[[str], str],
) -> None:
    payload = _summary_payload(summary, source, target_features, strategy_label)
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


def _categorize_clinvar(values: pd.Series) -> pd.Series:
    text = values.fillna("").astype(str).str.lower()
    category = pd.Series("Other", index=values.index, dtype="object")
    category[text.eq("")] = "Not Found"
    conflicting = text.str.contains("conflicting", na=False)
    uncertain = text.str.contains("uncertain|vus", regex=True, na=False)
    benign = text.str.contains("benign", na=False)
    pathogenic = text.str.contains("pathogenic", na=False)
    category[conflicting] = "Other"
    category[uncertain & ~conflicting] = "VUS"
    category[pathogenic & ~benign & ~uncertain & ~conflicting] = "P/LP"
    category[benign & ~pathogenic & ~uncertain & ~conflicting] = "B/LB"
    return category


def _normalize_chunk(df: pd.DataFrame) -> pd.DataFrame:
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
    df["target_start0"] = pd.to_numeric(df["target_start0"], errors="coerce")
    for column in ["support_row_count", "support_ortholog_count", "clinvar_scv_count"]:
        df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0).astype("int64")
    df["clinvar_found"] = df["clinvar_id"].astype(str) != ""
    df["clinvar_classified"] = df["clinvar_sig"].astype(str) != ""
    df["clinvar_category"] = _categorize_clinvar(df["clinvar_sig"])

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
            titv_kind TEXT NOT NULL,
            review_stars TEXT NOT NULL,
            gnomad_csq TEXT NOT NULL,
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
    return keep.sort_values(
        ["_star_rank", "support_ortholog_count", "support_row_count"],
        ascending=False,
        kind="mergesort",
    ).head(limit)


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
               SUM(gnomad_af IS NOT NULL) AS "gnomAD Found",
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
    denominator = stats["Unique Variants"].replace(0, np.nan)
    stats["ClinVar found %"] = stats["Found in ClinVar"] / denominator
    stats["ClinVar classified %"] = stats["ClinVar classified"] / denominator
    stats["gnomAD found %"] = stats["gnomAD Found"] / denominator
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


def build_variant_summary(
    path: Path,
    work_dir: Path,
    strategy_label: Callable[[str], str],
    target_features_path: Path | None = None,
    chunk_size: int = 100_000,
) -> VariantSummary:
    """Aggregate a variant annotation table without retaining row-level data in memory."""
    work_dir.mkdir(parents=True, exist_ok=True)
    cache_path = work_dir / SUMMARY_CACHE_NAME
    cached = _load_summary_cache(cache_path, path, target_features_path, strategy_label)
    if cached is not None:
        return cached

    summary = _compute_variant_summary(path, work_dir, target_features_path, strategy_label, chunk_size)
    _write_summary_cache(cache_path, summary, path, target_features_path, strategy_label)
    return summary


def _compute_variant_summary(
    path: Path,
    work_dir: Path,
    target_features_path: Path | None,
    strategy_label: Callable[[str], str],
    chunk_size: int,
) -> VariantSummary:
    header = _header(path)
    missing = VARIANT_REQUIRED - set(header)
    if missing:
        raise ValueError(f"Variant annotations missing required columns: {', '.join(sorted(missing))}")
    usecols = [column for column in VARIANT_USECOLS if column in header]
    contexts: dict[str, list[tuple[int, int, str]]] = {}
    context_starts: dict[str, list[int]] = {}
    if target_features_path is not None:
        if "target_start0" not in header:
            raise ValueError("Variant annotations need target_start0 for target-context reporting.")
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
    genes: set[str] = set()
    pathogenic_rows = pd.DataFrame(columns=[*VARIANT_USECOLS, "variant_id", "clinvar_category", "review_stars"])

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
                clinvar_classified, clinvar_category, gnomad_af, titv_kind,
                review_stars, gnomad_csq
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        for chunk in reader:
            chunk = _normalize_chunk(chunk)
            if target_features_path is not None and chunk["target_start0"].isna().any():
                raise ValueError("Variant annotations contain missing target_start0 values.")
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
                "titv_kind",
                "review_stars",
                "gnomad_csq",
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
                titv_kind,
                review_stars,
                gnomad_csq,
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
                            str(titv_kind),
                            str(review_stars),
                            str(gnomad_csq),
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
            SELECT strategy, CASE WHEN gnomad_af IS NULL THEN 'not_found' ELSE 'found' END AS gnomad_status,
                   event_type, COUNT(*) AS Variant_Count
            FROM memberships GROUP BY strategy, gnomad_status, event_type
            """,
            connection,
        )
        gnomad_context_counts = pd.read_sql_query(
            """
            SELECT strategy, CASE WHEN gnomad_af IS NULL THEN 'not_found' ELSE 'found' END AS gnomad_status,
                   target_context, COUNT(*) AS Variant_Count
            FROM memberships GROUP BY strategy, gnomad_status, target_context
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
            "gnomad_csq",
            strategy_label,
            "gnomad_af IS NOT NULL",
        )
        pathogenic_consequence_counts = _grouped_counts(
            connection,
            "gnomad_csq",
            strategy_label,
            "gnomad_af IS NOT NULL AND clinvar_category = 'P/LP'",
        )
        strategy_record_count = int(sum(totals.values()))
        gnomad_af_summary = _gnomad_af_summary(connection, strategy_label)
    finally:
        connection.close()
        database_path.unlink(missing_ok=True)

    return VariantSummary(
        input_row_count=input_row_count,
        unique_variant_count=unique_variant_count,
        all_strategy_variant_count=all_strategy_variant_count,
        strategy_record_count=strategy_record_count,
        gene_count=len(genes),
        clinvar_found=clinvar_found,
        clinvar_classified=clinvar_classified,
        gnomad_found=gnomad_found,
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
    )
