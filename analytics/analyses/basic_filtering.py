"""Simple support-threshold analyses for candidate variant filtering."""

from __future__ import annotations

import csv
import gzip
import json
import math
import tempfile
from dataclasses import dataclass
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from analytics.vep.consequences import (
    VALIDATION_CONSEQUENCE_OPTIONS as CONSEQUENCE_OPTIONS,
)
from analytics.io.artifacts import path_metadata, write_json_atomic
from analytics.io.variant_source import (
    resolve_pre_vep_variant_source,
    sql_string,
    variant_source_sql,
)
from genomics.variants import read_failed_regions
from .conservation_validation import (
    PHYLOP_BANDS,
    SCORE_COLUMN,
    TARGET_CONTEXT_OPTIONS,
    VARIANT_TYPE_OPTIONS,
    ConservationCohort,
    add_grouped_bh,
    assign_phylop_band,
    filter_consequence,
    filter_target_context,
    filter_variant_type,
)
from .statistics import EnrichmentResult, enrichment_result, mantel_haenszel_adjusted
from .variant_summary_aggregation import (
    _configure_duckdb_memory,
    available_cpu_count,
)


FILTER_SCORE_SCHEMA_VERSION = 1
MAX_CLINVAR_THRESHOLDS = 20
FILTER_SCORE_COLUMNS = [
    "variant_key",
    "strategy",
    "variant_type",
    "gnomad_status",
    "ortholog_support",
    "strategy_support",
    "genus_support",
]
FILTER_OPTIONS = [
    ("ortholog", "Exact-ALT ortholog support", "ortholog_support"),
    ("strategy", "Strategy support", "strategy_support"),
    ("genus", "Independent genus support", "genus_support"),
]


@dataclass(frozen=True)
class BasicFilteringAnalysis:
    score_path: Path
    score_manifest_path: Path
    candidate_curves: pd.DataFrame
    clinvar_curves: pd.DataFrame
    cache_hit: bool = False


def build_basic_filtering_analysis(
    *,
    variant_annotations_tsv: Path,
    variant_strategy_support_tsv: Path,
    annotation_failures_tsv: Path,
    analytics_dir: Path,
    cohort: ConservationCohort,
    strategies: list[str],
    eligible_gene_ids_by_strategy: dict[str, set[str]],
) -> BasicFilteringAnalysis:
    """Build candidate-retention, gnomAD, and ClinVar threshold curves."""

    score_path, manifest_path, cache_hit = build_or_load_filter_score_store(
        variant_annotations_tsv=variant_annotations_tsv,
        variant_strategy_support_tsv=variant_strategy_support_tsv,
        annotation_failures_tsv=annotation_failures_tsv,
        analytics_dir=analytics_dir,
        strategies=strategies,
    )
    histograms = read_filter_score_histograms(score_path)
    candidate_curves = candidate_curves_from_histograms(histograms)
    clinvar_scores = read_clinvar_filter_scores(score_path, cohort.variants)
    clinvar_curves = compute_clinvar_filter_curves(
        cohort=cohort,
        clinvar_scores=clinvar_scores,
        candidate_curves=candidate_curves,
        strategies=strategies,
        eligible_gene_ids_by_strategy=eligible_gene_ids_by_strategy,
    )
    return BasicFilteringAnalysis(
        score_path,
        manifest_path,
        candidate_curves,
        clinvar_curves,
        cache_hit,
    )


def build_or_load_filter_score_store(
    *,
    variant_annotations_tsv: Path,
    variant_strategy_support_tsv: Path,
    annotation_failures_tsv: Path,
    analytics_dir: Path,
    strategies: list[str],
) -> tuple[Path, Path, bool]:
    """Materialize one bounded-width allele/strategy table for all filter views."""

    source = resolve_pre_vep_variant_source(
        variant_annotations_tsv,
        required_columns={"variant_key", "event_type", "lookup_status", "gnomad_af"},
    )
    support_columns = _read_header(variant_strategy_support_tsv)
    required_support = {
        "variant_key",
        "gene_id",
        "strategy",
        "alt_support_ortholog_count",
        "alt_support_genus_count",
    }
    missing = required_support - set(support_columns)
    if missing:
        raise ValueError(
            f"Variant strategy support {variant_strategy_support_tsv} missing columns: "
            + ", ".join(sorted(missing))
        )
    outdir = analytics_dir / "basic_filtering"
    score_path = outdir / "filter_scores.parquet"
    manifest_path = outdir / "manifest.json"
    expected_inputs = {
        "schema_version": FILTER_SCORE_SCHEMA_VERSION,
        "variant_source": source.identity,
        "strategy_support": path_metadata(variant_strategy_support_tsv),
        "annotation_failures": path_metadata(annotation_failures_tsv),
        "strategies": sorted(str(value) for value in strategies),
    }
    if score_path.exists() and manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
            if (
                manifest.get("status") == "complete"
                and manifest.get("inputs") == expected_inputs
                and manifest.get("output") == path_metadata(score_path)
            ):
                _validate_score_store(score_path, int(manifest["row_count"]))
                return score_path, manifest_path, True
        except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
            pass

    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".filter_scores.", dir=outdir) as temporary:
        row_count = _build_filter_score_store(
            source=source,
            support_path=variant_strategy_support_tsv,
            support_columns=support_columns,
            annotation_failures_tsv=annotation_failures_tsv,
            score_path=score_path,
            temp_dir=Path(temporary),
        )
    manifest = {
        "status": "complete",
        "inputs": expected_inputs,
        "row_count": row_count,
        "output": path_metadata(score_path),
    }
    write_json_atomic(manifest_path, manifest)
    return score_path, manifest_path, False


def _build_filter_score_store(
    *,
    source,
    support_path: Path,
    support_columns: list[str],
    annotation_failures_tsv: Path,
    score_path: Path,
    temp_dir: Path,
) -> int:
    thread_count = available_cpu_count()
    with tempfile.NamedTemporaryFile(
        dir=score_path.parent,
        prefix=f".{score_path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary_path = Path(handle.name)
    temporary_path.unlink(missing_ok=True)
    connection = duckdb.connect()
    try:
        connection.execute(f"SET threads={thread_count}")
        connection.execute("SET preserve_insertion_order=false")
        connection.execute(f"SET temp_directory={sql_string(temp_dir)}")
        _configure_duckdb_memory(connection, thread_count)
        _register_gnomad_failures(connection, annotation_failures_tsv)
        connection.execute(
            f"CREATE VIEW variant_source_rows AS SELECT * FROM {variant_source_sql(source)}"
        )
        support_sql = _support_source_sql(support_path, support_columns)
        connection.execute(f"CREATE VIEW strategy_support_rows AS SELECT * FROM {support_sql}")

        failed = (
            "EXISTS (SELECT 1 FROM gnomad_failures f WHERE f.chrom = key_chrom "
            "AND key_pos BETWEEN f.start1 AND f.end1)"
        )
        connection.execute(
            "CREATE TEMP TABLE global_annotations AS WITH parsed AS ("
            "SELECT variant_key, lower(event_type) AS event_type, lookup_status, gnomad_af, "
            "regexp_replace(split_part(variant_key, ':', 1), '^chr', '') AS key_chrom, "
            "try_cast(split_part(variant_key, ':', 2) AS BIGINT) AS key_pos "
            "FROM variant_source_rows WHERE variant_key <> ''"
            "), classified AS (SELECT *, "
            "CASE WHEN try_cast(nullif(gnomad_af, '') AS DOUBLE) IS NOT NULL THEN 'found' "
            "WHEN coalesce(lookup_status, '') NOT IN ('', 'ok') OR key_pos IS NULL "
            "THEN 'lookup_failed' "
            f"WHEN {failed} THEN 'lookup_failed' ELSE 'not_found' END AS row_gnomad_status "
            "FROM parsed) SELECT variant_key, "
            "CASE WHEN bool_or(event_type = 'snv') THEN 'snv' ELSE 'indel' END AS variant_type, "
            "CASE WHEN bool_or(row_gnomad_status = 'found') THEN 'found' "
            "WHEN bool_or(row_gnomad_status = 'lookup_failed') THEN 'lookup_failed' "
            "ELSE 'not_found' END AS gnomad_status "
            "FROM classified GROUP BY variant_key"
        )
        connection.execute(
            "CREATE TEMP TABLE compact_support AS SELECT variant_key, trim(strategy) AS strategy, "
            "max(try_cast(alt_support_ortholog_count AS BIGINT)) AS ortholog_support, "
            "max(try_cast(nullif(alt_support_genus_count, '') AS BIGINT)) AS genus_support "
            "FROM strategy_support_rows WHERE variant_key <> '' AND trim(strategy) <> '' "
            "GROUP BY variant_key, trim(strategy)"
        )
        invalid = int(
            connection.execute(
                "SELECT count(*) FROM compact_support WHERE ortholog_support IS NULL "
                "OR ortholog_support < 1 OR genus_support < 0 "
                "OR genus_support > ortholog_support"
            ).fetchone()[0]
        )
        if invalid:
            raise ValueError(f"Variant strategy support contains {invalid} invalid filter scores")
        missing_annotations = int(
            connection.execute(
                "SELECT count(*) FROM compact_support s "
                "LEFT JOIN global_annotations a USING (variant_key) "
                "WHERE a.variant_key IS NULL"
            ).fetchone()[0]
        )
        if missing_annotations:
            raise ValueError(
                "Variant strategy support contains "
                f"{missing_annotations} alleles absent from annotations"
            )
        connection.execute(
            "CREATE TEMP TABLE filter_scores AS SELECT s.variant_key, s.strategy, "
            "a.variant_type, a.gnomad_status, s.ortholog_support, "
            "count(*) OVER (PARTITION BY s.variant_key) AS strategy_support, "
            "s.genus_support FROM compact_support s JOIN global_annotations a USING (variant_key)"
        )
        row_count = int(connection.execute("SELECT count(*) FROM filter_scores").fetchone()[0])
        connection.execute(
            f"COPY filter_scores TO {sql_string(temporary_path)} "
            "(FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)"
        )
        _validate_score_store(temporary_path, row_count)
        temporary_path.chmod(0o644)
        temporary_path.replace(score_path)
        return row_count
    finally:
        connection.close()
        temporary_path.unlink(missing_ok=True)


def read_filter_score_histograms(score_path: Path) -> pd.DataFrame:
    unions = []
    for filter_key, _label, column in FILTER_OPTIONS:
        eligibility = " AND genus_support IS NOT NULL" if filter_key == "genus" else ""
        unions.append(
            "SELECT strategy, variant_type, "
            f"{sql_string(filter_key)} AS filter_key, {column} AS score, gnomad_status, "
            "count(*) AS variant_count FROM read_parquet(?) WHERE "
            f"{column} IS NOT NULL{eligibility} GROUP BY ALL"
        )
    with duckdb.connect() as connection:
        return connection.execute(
            " UNION ALL ".join(unions),
            [str(score_path)] * len(unions),
        ).fetchdf()


def candidate_curves_from_histograms(histograms: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "strategy",
        "variant_type",
        "filter_key",
        "threshold",
        "retained_variant_count",
        "total_variant_count",
        "retained_fraction",
        "gnomad_found_count",
        "gnomad_eligible_count",
        "gnomad_lookup_failed_count",
        "gnomad_found_fraction",
    ]
    if histograms.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, object]] = []
    group_columns = ["strategy", "variant_type", "filter_key"]
    for keys, subset in histograms.groupby(group_columns, sort=True):
        strategy, variant_type, filter_key = keys
        subset = subset.copy()
        subset["score"] = pd.to_numeric(subset["score"], errors="raise").astype(int)
        subset["variant_count"] = pd.to_numeric(
            subset["variant_count"], errors="raise"
        ).astype(int)
        max_score = int(subset["score"].max())
        totals = {
            status: np.bincount(
                subset.loc[subset["gnomad_status"].eq(status), "score"],
                weights=subset.loc[subset["gnomad_status"].eq(status), "variant_count"],
                minlength=max_score + 1,
            )
            for status in ("found", "not_found", "lookup_failed")
        }
        cumulative = {
            status: np.cumsum(values[::-1])[::-1]
            for status, values in totals.items()
        }
        total = int(subset["variant_count"].sum())
        for threshold in range(1, max_score + 1):
            found = int(cumulative["found"][threshold])
            not_found = int(cumulative["not_found"][threshold])
            failed = int(cumulative["lookup_failed"][threshold])
            retained = found + not_found + failed
            eligible = found + not_found
            rows.append(
                {
                    "strategy": str(strategy),
                    "variant_type": str(variant_type),
                    "filter_key": str(filter_key),
                    "threshold": threshold,
                    "retained_variant_count": retained,
                    "total_variant_count": total,
                    "retained_fraction": retained / total if total else float("nan"),
                    "gnomad_found_count": found,
                    "gnomad_eligible_count": eligible,
                    "gnomad_lookup_failed_count": failed,
                    "gnomad_found_fraction": found / eligible if eligible else float("nan"),
                }
            )
    return pd.DataFrame(rows, columns=columns)


def read_clinvar_filter_scores(score_path: Path, cohort: pd.DataFrame) -> pd.DataFrame:
    keys = pd.DataFrame(
        {"variant_key": sorted(set(cohort["variant_key"].astype(str)))}
    )
    with duckdb.connect() as connection:
        connection.register("clinvar_keys", keys)
        return connection.execute(
            "SELECT s.variant_key, s.strategy, s.variant_type, s.ortholog_support, "
            "s.strategy_support, s.genus_support FROM read_parquet(?) s "
            "JOIN clinvar_keys c USING (variant_key)",
            [str(score_path)],
        ).fetchdf()


def compute_clinvar_filter_curves(
    *,
    cohort: ConservationCohort,
    clinvar_scores: pd.DataFrame,
    candidate_curves: pd.DataFrame,
    strategies: list[str],
    eligible_gene_ids_by_strategy: dict[str, set[str]],
) -> pd.DataFrame:
    """Compute unadjusted and fixed-band OR curves without repeated cohort scans."""

    score_maps = _clinvar_score_maps(clinvar_scores)
    threshold_map = {
        (str(strategy), str(variant_type), str(filter_key)): _select_clinvar_thresholds(
            group
        )
        for (strategy, variant_type, filter_key), group in candidate_curves.groupby(
            ["strategy", "variant_type", "filter_key"], sort=False
        )
    }
    rows: list[dict[str, object]] = []
    for strategy in strategies:
        eligible = eligible_gene_ids_by_strategy.get(strategy, set())
        strategy_frame = cohort.variants[
            cohort.variants["gene_ids"].map(
                lambda value: bool(eligible.intersection(str(value).split("|")))
            )
        ].copy()
        for filter_key, _filter_label, score_column in FILTER_OPTIONS:
            score_map = score_maps.get((strategy, score_column), {})
            strategy_frame["filter_score"] = (
                strategy_frame["variant_key"].astype(str).map(score_map).fillna(0).astype(int)
            )
            for variant_type, _variant_label in VARIANT_TYPE_OPTIONS:
                thresholds = threshold_map.get((strategy, variant_type, filter_key), [])
                if not thresholds:
                    continue
                type_frame = filter_variant_type(strategy_frame, variant_type)
                for target_context, _context_label in TARGET_CONTEXT_OPTIONS:
                    context_frame = filter_target_context(type_frame, target_context)
                    for consequence, _consequence_label in CONSEQUENCE_OPTIONS:
                        working = filter_consequence(context_frame, consequence)
                        rows.extend(
                            _unadjusted_curve_rows(
                                working,
                                thresholds,
                                strategy=strategy,
                                variant_type=variant_type,
                                target_context=target_context,
                                consequence=consequence,
                                filter_key=filter_key,
                            )
                        )
                        rows.extend(
                            _fixed_curve_rows(
                                working,
                                thresholds,
                                strategy=strategy,
                                variant_type=variant_type,
                                target_context=target_context,
                                consequence=consequence,
                                filter_key=filter_key,
                            )
                        )
    results = pd.DataFrame(rows)
    if results.empty:
        return results
    results = results[results["status"].astype(str).eq("estimated")].reset_index(drop=True)
    if results.empty:
        return results
    add_grouped_bh(
        results,
        "result_p",
        "result_q",
        [
            "mode",
            "filter_key",
            "threshold",
            "variant_type",
            "target_context",
            "consequence",
        ],
    )
    return results


def _select_clinvar_thresholds(group: pd.DataFrame) -> list[int]:
    """Keep at most 20 informative thresholds across the retention curve."""

    ordered = group.sort_values("threshold")
    if "retained_variant_count" in ordered:
        ordered = ordered.drop_duplicates("retained_variant_count", keep="first")
    thresholds = ordered["threshold"].astype(int).drop_duplicates().tolist()
    if len(thresholds) <= MAX_CLINVAR_THRESHOLDS:
        return thresholds
    indices = np.linspace(0, len(thresholds) - 1, MAX_CLINVAR_THRESHOLDS)
    return [thresholds[index] for index in sorted(set(np.rint(indices).astype(int)))]


def _unadjusted_curve_rows(
    working: pd.DataFrame,
    thresholds: list[int],
    **dimensions: object,
) -> list[dict[str, object]]:
    count_rows = _count_curves(working, thresholds)
    return [
        _association_row_from_counts(
            len(working),
            threshold,
            counts,
            mode="unadjusted",
            **dimensions,
        )
        for threshold, counts in zip(thresholds, count_rows)
    ]


def _fixed_curve_rows(
    working: pd.DataFrame,
    thresholds: list[int],
    **dimensions: object,
) -> list[dict[str, object]]:
    scored = working[np.isfinite(pd.to_numeric(working[SCORE_COLUMN], errors="coerce"))].copy()
    scored["band"] = assign_phylop_band(scored[SCORE_COLUMN])
    band_count_curves = {
        band_key: _count_curves(scored[scored["band"] == band_key], thresholds)
        for band_key, _band_label, _range_text in PHYLOP_BANDS
    }
    rows = []
    for threshold_index, threshold in enumerate(thresholds):
        strata = []
        for band_key, _band_label, _range_text in PHYLOP_BANDS:
            counts = band_count_curves[band_key][threshold_index]
            strata.append(
                EnrichmentResult(
                    str(band_key),
                    *counts,
                    float("nan"),
                    float("nan"),
                    float("nan"),
                    float("nan"),
                )
            )
        benign_total = int((scored["label_class"].astype(str) == "benign").sum())
        pathogenic_total = int((scored["label_class"].astype(str) == "pathogenic").sum())
        observed_total = int((scored["filter_score"].astype(int) >= threshold).sum())
        reason = _estimability_reason(
            len(scored), benign_total, pathogenic_total, observed_total
        )
        if scored["band"].nunique() < 2:
            reason = reason or "At least two populated phyloP bands are required."
        adjusted = mantel_haenszel_adjusted(strata) if not reason else None
        if adjusted is None or not math.isfinite(adjusted.cmh_p):
            reason = reason or "The fixed-band tables do not provide an estimable common OR."
        elif not (
            math.isfinite(adjusted.odds_ratio_mh)
            and adjusted.odds_ratio_mh > 0
            and math.isfinite(adjusted.ci_low)
            and math.isfinite(adjusted.ci_high)
        ):
            reason = "The fixed-band common OR or confidence interval is not finite."
        rows.append(
            {
                **dimensions,
                "mode": "fixed",
                "threshold": threshold,
                "usable_rows": int(len(scored)),
                "benign_observed": sum(item.benign_observed for item in strata),
                "pathogenic_observed": sum(item.pathogenic_observed for item in strata),
                "benign_not_observed": sum(item.benign_not_observed for item in strata),
                "pathogenic_not_observed": sum(item.pathogenic_not_observed for item in strata),
                "result_or": adjusted.odds_ratio_mh if adjusted and not reason else float("nan"),
                "ci_low": adjusted.ci_low if adjusted and not reason else float("nan"),
                "ci_high": adjusted.ci_high if adjusted and not reason else float("nan"),
                "result_p": adjusted.cmh_p if adjusted and not reason else float("nan"),
                "status": "estimated" if adjusted and not reason else "not_estimable",
                "reason": reason,
            }
        )
    return rows


def _association_row_from_counts(
    row_count: int,
    threshold: int,
    counts: tuple[int, int, int, int],
    *,
    mode: str,
    **dimensions: object,
) -> dict[str, object]:
    observed_total = counts[0] + counts[1]
    reason = _estimability_reason(
        row_count, counts[0] + counts[2], counts[1] + counts[3], observed_total
    )
    result = enrichment_result(str(threshold), *counts) if not reason else None
    return {
        **dimensions,
        "mode": mode,
        "threshold": threshold,
        "usable_rows": int(row_count),
        "benign_observed": counts[0],
        "pathogenic_observed": counts[1],
        "benign_not_observed": counts[2],
        "pathogenic_not_observed": counts[3],
        "result_or": result.odds_ratio if result else float("nan"),
        "ci_low": result.ci_low if result else float("nan"),
        "ci_high": result.ci_high if result else float("nan"),
        "result_p": result.fisher_p if result else float("nan"),
        "status": "not_estimable" if reason else "estimated",
        "reason": reason,
    }


def _count_curves(
    working: pd.DataFrame,
    thresholds: list[int],
) -> list[tuple[int, int, int, int]]:
    """Return 2x2 counts for all thresholds from two cumulative histograms."""

    if not thresholds:
        return []
    max_score = max(max(thresholds), int(working["filter_score"].max()) if len(working) else 0)
    labels = working["label_class"].astype(str)
    scores = working["filter_score"].astype(int)
    observed_by_label = {}
    totals = {}
    for label in ("benign", "pathogenic"):
        values = scores[labels.eq(label)].to_numpy(dtype=int)
        totals[label] = int(len(values))
        histogram = np.bincount(values, minlength=max_score + 1)
        observed_by_label[label] = np.cumsum(histogram[::-1])[::-1]
    rows = []
    for threshold in thresholds:
        benign_observed = int(observed_by_label["benign"][threshold])
        pathogenic_observed = int(observed_by_label["pathogenic"][threshold])
        rows.append(
            (
                benign_observed,
                pathogenic_observed,
                totals["benign"] - benign_observed,
                totals["pathogenic"] - pathogenic_observed,
            )
        )
    return rows


def _estimability_reason(
    row_count: int,
    benign_total: int,
    pathogenic_total: int,
    observed_total: int,
) -> str:
    if row_count == 0:
        return "No ClinVar alleles in this selector."
    if benign_total == 0 or pathogenic_total == 0:
        return "Both B/LB and P/LP alleles are required."
    if observed_total == 0 or observed_total == row_count:
        return "Both ALT-observed and ALT-not-observed alleles are required."
    return ""


def _clinvar_score_maps(
    scores: pd.DataFrame,
) -> dict[tuple[str, str], dict[str, int]]:
    maps: dict[tuple[str, str], dict[str, int]] = {}
    for strategy, group in scores.groupby("strategy", sort=False):
        for _filter_key, _label, column in FILTER_OPTIONS:
            valid = group.dropna(subset=[column])
            maps[(str(strategy), column)] = dict(
                zip(valid["variant_key"].astype(str), valid[column].astype(int))
            )
    return maps


def _register_gnomad_failures(connection, path: Path) -> None:
    rows = []
    for chrom, (_starts, intervals) in read_failed_regions(path, "gnomad").items():
        rows.extend(
            {"chrom": str(chrom), "start1": int(start), "end1": int(end)}
            for start, end in intervals
        )
    frame = pd.DataFrame(rows, columns=["chrom", "start1", "end1"])
    connection.register("gnomad_failures_input", frame)
    connection.execute(
        "CREATE VIEW gnomad_failures AS SELECT cast(chrom AS VARCHAR) AS chrom, "
        "cast(start1 AS BIGINT) AS start1, cast(end1 AS BIGINT) AS end1 "
        "FROM gnomad_failures_input"
    )


def _support_source_sql(path: Path, columns: list[str]) -> str:
    schema = "{" + ",".join(f"{sql_string(column)}:'VARCHAR'" for column in columns) + "}"
    return (
        f"read_csv({sql_string(path)}, delim='\\t', header=true, columns={schema}, "
        "auto_detect=false, compression='auto', parallel=true, "
        "nullstr='__GAPH_NULL_SENTINEL__')"
    )


def _read_header(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def _validate_score_store(path: Path, expected_rows: int) -> None:
    with duckdb.connect() as connection:
        relation = connection.read_parquet(str(path))
        if relation.columns != FILTER_SCORE_COLUMNS:
            raise ValueError(
                f"Basic-filter score columns changed in {path}: "
                + ", ".join(relation.columns)
            )
        observed = int(relation.count("*").fetchone()[0])
        if observed != expected_rows:
            raise ValueError(
                f"Basic-filter score row count changed: {observed} != {expected_rows}"
            )
