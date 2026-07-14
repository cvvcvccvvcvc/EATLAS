"""Bounded-memory aggregation of variant annotations for the HTML report."""

from __future__ import annotations

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


VARIANT_USECOLS = [
    "variant_key",
    "gene_id",
    "event_type",
    "ref",
    "alt",
    "lookup_ref",
    "lookup_alt",
    "strategies",
    "support_row_count",
    "support_ortholog_count",
    "support_strategy_count",
    "clinvar_id",
    "clinvar_allele_id",
    "clinvar_sig",
    "clinvar_revstat",
    "clinvar_review_stars",
    "clinvar_review_stars_status",
    "clinvar_sig_conflict",
    "clinvar_scv_count",
    "clinvar_hgvs",
    "clinvar_geneinfo",
    "clinvar_disease",
    "clinvar_variant_type",
    "clinvar_origin",
    "clinvar_rs",
    "gnomad_af",
    "gnomad_af_source",
    "gnomad_csq",
    "gnomad_hgvsc",
    "gnomad_hgvsp",
]
VARIANT_REQUIRED = {"variant_key", "gene_id", "event_type", "strategies"}


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
    strategy_record_count: int
    gene_count: int
    clinvar_found: int
    clinvar_classified: int
    gnomad_found: int
    strategies: list[str]
    strategy_stats: pd.DataFrame
    unique_contribution: pd.DataFrame
    event_counts: pd.DataFrame
    overlap: StrategyOverlap | None
    clinvar_counts: pd.DataFrame
    gnomad_bins: pd.DataFrame
    pathogenic_star_counts: pd.DataFrame
    consequence_counts: pd.DataFrame
    pathogenic_consequence_counts: pd.DataFrame
    pathogenic_rows: pd.DataFrame


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
    for column in ["support_row_count", "support_ortholog_count", "support_strategy_count", "clinvar_scv_count"]:
        df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0).astype("int64")
    df["clinvar_found"] = df["clinvar_id"].astype(str) != ""
    df["clinvar_classified"] = df["clinvar_sig"].astype(str) != ""
    df["clinvar_category"] = _categorize_clinvar(df["clinvar_sig"])

    ref = df["lookup_ref"].where(df["lookup_ref"].astype(str) != "", df["ref"]).astype(str).str.upper()
    alt = df["lookup_alt"].where(df["lookup_alt"].astype(str) != "", df["alt"]).astype(str).str.upper()
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
) -> tuple[int, dict[str, int], StrategyOverlap | None]:
    intersections: Counter[tuple[str, str]] = Counter()
    unique_counts: Counter[str] = Counter()
    unique_variant_count = 0
    current_variant = None
    current_strategies: list[str] = []

    def consume_group(group: list[str]) -> None:
        nonlocal unique_variant_count
        if not group:
            return
        unique_variant_count += 1
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
        return unique_variant_count, dict(unique_counts), None
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
    return unique_variant_count, dict(unique_counts), StrategyOverlap(ordered, shared, unions, jaccard)


def _gnomad_bins(
    connection: sqlite3.Connection,
    strategy_label: Callable[[str], str],
    bin_count: int = 10,
) -> pd.DataFrame:
    limits = connection.execute(
        "SELECT MIN(gnomad_af), MAX(gnomad_af) FROM memberships WHERE gnomad_af > 0"
    ).fetchone()
    if not limits or limits[0] is None:
        return pd.DataFrame(columns=["strategy", "bin_mid", "Variant_Count", "Density"])
    min_value, max_value = (math.log10(float(value)) for value in limits)
    if min_value == max_value:
        min_value -= 0.5
        max_value += 0.5
    edges = np.linspace(min_value, max_value, bin_count + 1)
    counts: Counter[tuple[str, int]] = Counter()
    totals: Counter[str] = Counter()
    for strategy, af in connection.execute("SELECT strategy, gnomad_af FROM memberships WHERE gnomad_af > 0"):
        value = math.log10(af)
        index = min(bin_count - 1, max(0, int((value - min_value) / (max_value - min_value) * bin_count)))
        counts[(strategy, index)] += 1
        totals[strategy] += 1
    rows = [
        {
            "strategy": strategy_label(strategy),
            "bin_mid": float((edges[index] + edges[index + 1]) / 2),
            "Variant_Count": count,
            "Density": count / totals[strategy],
        }
        for (strategy, index), count in sorted(counts.items())
    ]
    return pd.DataFrame(rows)


def build_variant_summary(
    path: Path,
    work_dir: Path,
    strategy_label: Callable[[str], str],
    chunk_size: int = 100_000,
) -> VariantSummary:
    """Aggregate a variant annotation table without retaining row-level data in memory."""
    header = _header(path)
    missing = VARIANT_REQUIRED - set(header)
    if missing:
        raise ValueError(f"Variant annotations missing required columns: {', '.join(sorted(missing))}")
    usecols = [column for column in VARIANT_USECOLS if column in header]
    work_dir.mkdir(parents=True, exist_ok=True)
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
                variant_id, strategy, gene_id, event_type, clinvar_found,
                clinvar_classified, clinvar_category, gnomad_af, titv_kind,
                review_stars, gnomad_csq
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        for chunk in reader:
            chunk = _normalize_chunk(chunk)
            input_row_count += len(chunk)
            clinvar_found += int(chunk["clinvar_found"].sum())
            clinvar_classified += int(chunk["clinvar_classified"].sum())
            gnomad_found += int(chunk["gnomad_af"].notna().sum())
            genes.update(chunk["gene_id"].astype(str))
            pathogenic_rows = _keep_top_pathogenic(pathogenic_rows, chunk)

            records = []
            for row in chunk.itertuples(index=False):
                strategies = [item.strip() for item in str(row.strategies).split(",") if item.strip()]
                for strategy in strategies:
                    records.append(
                        (
                            str(row.variant_id),
                            strategy,
                            str(row.gene_id),
                            str(row.event_type),
                            int(row.clinvar_found),
                            int(row.clinvar_classified),
                            str(row.clinvar_category),
                            None if pd.isna(row.gnomad_af) else float(row.gnomad_af),
                            str(row.titv_kind),
                            str(row.review_stars),
                            str(row.gnomad_csq),
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
        unique_variant_count, unique_counts, overlap = _overlap_summary(connection, strategies, totals)
        strategy_stats = _strategy_stats(connection, strategy_label)
        unique_contribution = pd.DataFrame(
            {
                "Strategy": [strategy_label(strategy) for strategy in strategies],
                "Unique To Strategy": [unique_counts.get(strategy, 0) for strategy in strategies],
            }
        ).sort_values("Strategy")

        event_counts = _grouped_counts(connection, "event_type", strategy_label).rename(columns={"value": "event_type"})
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
        gnomad_bins = _gnomad_bins(connection, strategy_label)
    finally:
        connection.close()
        database_path.unlink(missing_ok=True)

    return VariantSummary(
        input_row_count=input_row_count,
        unique_variant_count=unique_variant_count,
        strategy_record_count=strategy_record_count,
        gene_count=len(genes),
        clinvar_found=clinvar_found,
        clinvar_classified=clinvar_classified,
        gnomad_found=gnomad_found,
        strategies=strategies,
        strategy_stats=strategy_stats,
        unique_contribution=unique_contribution,
        event_counts=event_counts,
        overlap=overlap,
        clinvar_counts=clinvar_counts,
        gnomad_bins=gnomad_bins,
        pathogenic_star_counts=star_counts,
        consequence_counts=consequence_counts,
        pathogenic_consequence_counts=pathogenic_consequence_counts,
        pathogenic_rows=pathogenic_rows.drop(columns=["_star_rank"], errors="ignore"),
    )
