"""DuckDB-backed allele collapse for candidate-wide conservation summaries."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterator

import pandas as pd

from analytics.io.variant_source import (
    VariantTableSource,
    resolve_pre_vep_variant_source,
    sql_string,
    variant_source_sql,
)
from genomics.variants import read_failed_regions


REQUIRED_COLUMNS = {"variant_key", "lookup_status", "strategies", "gnomad_af"}
MAX_STRATEGIES = 63


CandidateAggregationSource = VariantTableSource


class CandidateAlleleStore:
    """Temporary one-row-per-allele relation used throughout phyloP annotation."""

    def __init__(
        self,
        *,
        source: CandidateAggregationSource,
        strategies: tuple[str, ...],
        annotation_failures_path: Path | None,
        temp_dir: Path,
    ) -> None:
        duckdb = _import_duckdb()
        self.source = source
        self.strategies = strategies
        self.connection = duckdb.connect()
        try:
            self.connection.execute(f"SET threads={available_cpu_count()}")
            self.connection.execute("SET preserve_insertion_order=false")
            self.connection.execute("SET enable_progress_bar=false")
            temp_dir.mkdir(parents=True, exist_ok=True)
            self.connection.execute(f"SET temp_directory={sql_string(temp_dir)}")
            self._register_failed_regions(annotation_failures_path)
            self._build_alleles()
            observed_rows = int(
                self.connection.execute(
                    "SELECT coalesce(sum(context_count), 0) FROM candidate_alleles"
                ).fetchone()[0]
            )
            if source.row_count is not None and observed_rows != source.row_count:
                raise ValueError(
                    "Candidate source row count changed: "
                    f"observed {observed_rows}, expected {source.row_count}"
                )
            self.connection.execute(
                "CREATE TEMP TABLE unsupported_variants (variant_key VARCHAR PRIMARY KEY)"
            )
        except BaseException:
            self.connection.close()
            raise

    def close(self) -> None:
        self.connection.close()

    def iter_position_rows(self, chunk_size: int) -> Iterator[list[tuple[object, ...]]]:
        cursor = self.connection.execute(
            "SELECT variant_key, eligible_position_context_count "
            "FROM candidate_alleles WHERE nonfailed_context_count > 0"
        )
        while rows := cursor.fetchmany(chunk_size):
            yield rows

    def register_unsupported(self, variant_keys: list[str]) -> None:
        if not variant_keys:
            return
        frame = pd.DataFrame({"variant_key": sorted(set(variant_keys))})
        self.connection.register("unsupported_variants_input", frame)
        try:
            self.connection.execute(
                "INSERT OR IGNORE INTO unsupported_variants "
                "SELECT cast(variant_key AS VARCHAR) FROM unsupported_variants_input"
            )
        finally:
            self.connection.unregister("unsupported_variants_input")

    def group_counts(self) -> pd.DataFrame:
        rows = []
        for index, strategy in enumerate(self.strategies):
            bit = 1 << index
            frame = self.connection.execute(
                "SELECT status.gnomad_status, count(*) AS variant_count "
                "FROM candidate_alleles a, "
                "LATERAL (VALUES "
                "  ('found', a.found_strategy_mask), "
                "  ('not_found', a.not_found_strategy_mask)"
                ") status(gnomad_status, strategy_mask) "
                "WHERE status.strategy_mask & ? != 0 "
                "AND NOT EXISTS (SELECT 1 FROM unsupported_variants u "
                "                WHERE u.variant_key = a.variant_key) "
                "GROUP BY status.gnomad_status ORDER BY status.gnomad_status",
                [bit],
            ).fetchdf()
            if frame.empty:
                continue
            frame.insert(0, "strategy", strategy)
            rows.append(frame)
        if not rows:
            return pd.DataFrame(columns=["strategy", "gnomad_status", "variant_count"])
        return pd.concat(rows, ignore_index=True)

    def iter_group_variant_keys(
        self,
        *,
        strategy: str,
        gnomad_status: str,
        chunk_size: int,
    ) -> Iterator[list[str]]:
        try:
            bit = 1 << self.strategies.index(strategy)
        except ValueError as exc:
            raise ValueError(f"Unknown candidate strategy {strategy!r}") from exc
        cursor = self.connection.execute(
            "SELECT a.variant_key FROM candidate_alleles a "
            "WHERE CASE WHEN ? = 'found' THEN a.found_strategy_mask "
            "           WHEN ? = 'not_found' THEN a.not_found_strategy_mask "
            "           ELSE 0 END & ? != 0 "
            "AND NOT EXISTS (SELECT 1 FROM unsupported_variants u "
            "                WHERE u.variant_key = a.variant_key)",
            [gnomad_status, gnomad_status, bit],
        )
        while rows := cursor.fetchmany(chunk_size):
            yield [str(row[0]) for row in rows]

    def summary(self) -> dict[str, int]:
        row = self.connection.execute(
            "SELECT "
            "sum(context_count), "
            "sum(position_failed_context_count), "
            "sum(lookup_failed_context_count), "
            "sum(bit_count(found_strategy_mask & raw_not_found_strategy_mask)), "
            "count(*) FILTER (WHERE found_strategy_mask | not_found_strategy_mask != 0 "
            "  AND NOT EXISTS (SELECT 1 FROM unsupported_variants u "
            "                  WHERE u.variant_key = candidate_alleles.variant_key)), "
            "sum(bit_count(found_strategy_mask | not_found_strategy_mask)) "
            "  FILTER (WHERE found_strategy_mask | not_found_strategy_mask != 0 "
            "  AND NOT EXISTS (SELECT 1 FROM unsupported_variants u "
            "                  WHERE u.variant_key = candidate_alleles.variant_key)) "
            "FROM candidate_alleles"
        ).fetchone()
        return {
            "variant_context_row_count": int(row[0] or 0),
            "position_failed_context_count": int(row[1] or 0),
            "lookup_failed_allele_context_count": int(row[2] or 0),
            "gnomad_status_conflict_membership_count": int(row[3] or 0),
            "unique_usable_allele_count": int(row[4] or 0),
            "strategy_variant_membership_count": int(row[5] or 0),
        }

    def context_summary(self) -> dict[str, int]:
        row = self.connection.execute(
            "SELECT sum(context_count), sum(position_failed_context_count), "
            "sum(lookup_failed_context_count) FROM candidate_alleles"
        ).fetchone()
        return {
            "variant_context_row_count": int(row[0] or 0),
            "position_failed_context_count": int(row[1] or 0),
            "lookup_failed_allele_context_count": int(row[2] or 0),
        }

    def _register_failed_regions(self, path: Path | None) -> None:
        rows = []
        for chrom, (_starts, intervals) in read_failed_regions(path, "gnomad").items():
            rows.extend(
                {"chrom": chrom, "start1": start, "end1": end}
                for start, end in intervals
            )
        failures = pd.DataFrame(rows, columns=["chrom", "start1", "end1"])
        self.connection.register("gnomad_failures_input", failures)
        self.connection.execute(
            "CREATE TEMP VIEW gnomad_failures AS "
            "SELECT cast(chrom AS VARCHAR) AS chrom, cast(start1 AS BIGINT) AS start1, "
            "cast(end1 AS BIGINT) AS end1 FROM gnomad_failures_input"
        )

    def _build_alleles(self) -> None:
        mask_sql = " + ".join(
            "CASE WHEN list_contains(list_transform(string_split(strategies, ','), "
            f"item -> trim(item)), {sql_string(strategy)}) THEN {1 << index} ELSE 0 END"
            for index, strategy in enumerate(self.strategies)
        ) or "0"
        raw_chrom = "trim(split_part(variant_key, ':', 1))"
        normalized_chrom = (
            f"CASE WHEN {raw_chrom} IN ('M', 'chrM') THEN 'MT' "
            f"WHEN starts_with({raw_chrom}, 'chr') THEN substr({raw_chrom}, 4) "
            f"ELSE {raw_chrom} END"
        )
        source_sql = variant_source_sql(self.source)
        self.connection.execute(
            "CREATE TEMP TABLE candidate_alleles AS WITH parsed AS ("
            "SELECT variant_key, lookup_status, strategies, "
            "try_cast(nullif(gnomad_af, '') AS DOUBLE) AS gnomad_af_value, "
            f"{normalized_chrom} AS key_chrom, "
            "split_part(variant_key, ':', 2) AS key_pos_text, "
            "try_cast(split_part(variant_key, ':', 2) AS BIGINT) AS key_pos, "
            "upper(split_part(split_part(variant_key, ':', 3), '>', 1)) AS key_ref, "
            "upper(split_part(split_part(variant_key, ':', 3), '>', 2)) AS key_alt "
            f"FROM {source_sql}"
            "), normalized AS ("
            "SELECT *, "
            f"({mask_sql})::UBIGINT AS strategy_mask, "
            "CASE WHEN gnomad_af_value IS NOT NULL AND NOT isnan(gnomad_af_value) "
            "THEN 'found' "
            "WHEN coalesce(lookup_status, '') NOT IN ('', 'ok') "
            "  OR key_pos IS NULL OR key_pos <= 0 OR key_chrom = '' "
            "  OR NOT regexp_full_match(key_pos_text, '[0-9]+') "
            "  OR NOT regexp_full_match(key_ref, '[ACGT]+') "
            "  OR NOT regexp_full_match(key_alt, '[ACGT]+') "
            "THEN 'lookup_failed' "
            "WHEN EXISTS (SELECT 1 FROM gnomad_failures f "
            "             WHERE f.chrom = key_chrom AND key_pos BETWEEN f.start1 AND f.end1) "
            "THEN 'lookup_failed' ELSE 'not_found' END AS row_gnomad_status "
            "FROM parsed"
            "), collapsed AS (SELECT variant_key, "
            "bit_or(CASE WHEN row_gnomad_status = 'found' THEN strategy_mask ELSE 0 END) "
            "  AS found_strategy_mask, "
            "bit_or(CASE WHEN row_gnomad_status = 'not_found' THEN strategy_mask ELSE 0 END) "
            "  AS raw_not_found_strategy_mask, "
            "count(*) AS context_count, "
            "count_if(row_gnomad_status <> 'lookup_failed') AS nonfailed_context_count, "
            "count_if(lookup_status = 'ok' AND row_gnomad_status <> 'lookup_failed') "
            "  AS eligible_position_context_count, "
            "count_if(lookup_status = 'ok' AND row_gnomad_status = 'lookup_failed') "
            "  AS position_failed_context_count, "
            "count_if(row_gnomad_status = 'lookup_failed') AS lookup_failed_context_count "
            "FROM normalized GROUP BY variant_key"
            ") SELECT *, raw_not_found_strategy_mask & ~found_strategy_mask "
            "  AS not_found_strategy_mask FROM collapsed"
        )


def build_candidate_allele_store(
    *,
    variant_annotations_tsv: Path,
    strategies: list[str] | tuple[str, ...] | None,
    annotation_failures_path: Path | None,
    temp_dir: Path,
) -> CandidateAlleleStore:
    source = resolve_candidate_aggregation_source(variant_annotations_tsv)
    selected = tuple(sorted({str(value).strip() for value in (strategies or []) if str(value).strip()}))
    if not selected:
        selected = _read_strategies(source)
    if len(selected) > MAX_STRATEGIES:
        raise ValueError(
            f"Candidate conservation supports at most {MAX_STRATEGIES} strategies, found {len(selected)}"
        )
    return CandidateAlleleStore(
        source=source,
        strategies=selected,
        annotation_failures_path=annotation_failures_path,
        temp_dir=temp_dir,
    )


def resolve_candidate_aggregation_source(path: Path) -> CandidateAggregationSource:
    return resolve_pre_vep_variant_source(
        path,
        required_columns=REQUIRED_COLUMNS,
    )


def available_cpu_count() -> int:
    allocated = os.environ.get("SLURM_CPUS_PER_TASK")
    if allocated and allocated.isdigit() and int(allocated) > 0:
        return int(allocated)
    return os.cpu_count() or 1


def _read_strategies(source: CandidateAggregationSource) -> tuple[str, ...]:
    duckdb = _import_duckdb()
    connection = duckdb.connect()
    try:
        return tuple(
            str(row[0])
            for row in connection.execute(
                "SELECT DISTINCT trim(strategy) FROM ("
                "SELECT unnest(string_split(strategies, ',')) AS strategy "
                f"FROM {variant_source_sql(source)}"
                ") WHERE trim(strategy) <> '' ORDER BY 1"
            ).fetchall()
        )
    finally:
        connection.close()


def _import_duckdb():
    try:
        import duckdb
    except ImportError as exc:  # pragma: no cover - analytics environment contract
        raise RuntimeError(
            "Candidate conservation aggregation requires the python-duckdb package"
        ) from exc
    return duckdb
