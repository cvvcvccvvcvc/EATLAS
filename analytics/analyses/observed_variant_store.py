"""Compact observed-variant memberships shared by report analyses."""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import pandas as pd

from analytics.io.artifacts import path_metadata, write_json_atomic
from analytics.io.variant_source import (
    VariantTableSource,
    resolve_pre_vep_variant_source,
    sql_string,
    variant_source_sql,
)


STORE_SCHEMA_VERSION = 1
MAX_STRATEGIES = 63
FOCAL_RANK_METHOD = "duckdb_md5_topk_v1"
REQUIRED_COLUMNS = {
    "variant_key",
    "gene_id",
    "event_type",
    "ref",
    "alt",
    "lookup_status",
    "strategies",
}
ALLELE_GENE_COLUMNS = [
    "variant_key",
    "gene_id",
    "event_type",
    "ref",
    "alt",
    "lookup_status",
    "strategy_mask",
]
ALLELE_COLUMNS = ["variant_key", "strategy_mask"]


@dataclass(frozen=True)
class ObservedVariantStore:
    allele_gene_path: Path
    allele_path: Path
    manifest_path: Path
    manifest: dict[str, object]
    strategies: tuple[str, ...]
    cache_hit: bool = False

    def iter_sampled_focal_rows(
        self,
        selected_strategies: list[str] | tuple[str, ...],
        eligible_gene_ids: set[str],
        *,
        limit: int,
        seed: int,
        chunk_size: int = 200_000,
    ) -> Iterator[list[tuple[object, ...]]]:
        """Yield deterministic per-strategy top-K focal SNV memberships."""

        if limit < 1:
            raise ValueError("Observed-variant focal sample limit must be >= 1")
        if chunk_size < 1:
            raise ValueError("Observed-variant focal chunk size must be >= 1")
        selected = tuple(dict.fromkeys(str(value) for value in selected_strategies))
        self.strategy_mask(selected)
        genes = pd.DataFrame(
            {"gene_id": sorted({str(value) for value in eligible_gene_ids})}
        )
        if not selected or genes.empty:
            return
        duckdb = _import_duckdb()
        with tempfile.TemporaryDirectory(
            prefix=".focal_sampling.",
            dir=self.allele_gene_path.parent,
        ) as temporary:
            with duckdb.connect() as connection:
                connection.execute(f"SET threads={available_cpu_count()}")
                connection.execute("SET preserve_insertion_order=false")
                connection.execute("SET enable_progress_bar=false")
                connection.execute(f"SET temp_directory={sql_string(Path(temporary))}")
                connection.register("eligible_focal_genes", genes)
                connection.execute(
                    "CREATE TEMP TABLE eligible_focal_rows AS WITH parsed AS ("
                    "SELECT a.variant_key, a.gene_id, a.strategy_mask, "
                    "trim(split_part(a.variant_key, ':', 1)) AS key_chrom, "
                    "split_part(a.variant_key, ':', 2) AS key_pos_text, "
                    "try_cast(split_part(a.variant_key, ':', 2) AS BIGINT) AS key_pos, "
                    "upper(split_part(split_part(a.variant_key, ':', 3), '>', 1)) "
                    "  AS key_ref, "
                    "upper(split_part(split_part(a.variant_key, ':', 3), '>', 2)) "
                    "  AS key_alt "
                    f"FROM read_parquet({sql_string(self.allele_gene_path)}) a "
                    "JOIN eligible_focal_genes g USING (gene_id) "
                    "WHERE a.event_type = 'snv' AND a.lookup_status = 'ok' "
                    "AND length(a.ref) = 1 AND length(a.alt) = 1 "
                    "AND upper(a.ref) IN ('A','C','G','T') "
                    "AND upper(a.alt) IN ('A','C','G','T')"
                    ") SELECT variant_key, gene_id, strategy_mask FROM parsed "
                    "WHERE CASE WHEN starts_with(key_chrom, 'chr') "
                    "           THEN substr(key_chrom, 4) ELSE key_chrom END <> '' "
                    "AND key_pos > 0 "
                    "AND length(variant_key) - length(replace(variant_key, ':', '')) = 2 "
                    "AND length(variant_key) - length(replace(variant_key, '>', '')) = 1 "
                    "AND regexp_full_match(key_pos_text, '[0-9]+') "
                    "AND regexp_full_match(key_ref, '[ACGT]+') "
                    "AND regexp_full_match(key_alt, '[ACGT]+')"
                )
                for strategy in selected:
                    bit = 1 << self.strategies.index(strategy)
                    cursor = connection.execute(
                        "SELECT variant_key, gene_id, ? AS strategy "
                        "FROM eligible_focal_rows WHERE strategy_mask & ? != 0 "
                        "ORDER BY md5(? || '|' || ? || '|' || gene_id || ':' || variant_key), "
                        "gene_id || ':' || variant_key LIMIT ?",
                        [strategy, bit, str(seed), strategy, limit],
                    )
                    while rows := cursor.fetchmany(chunk_size):
                        yield rows

    def observed_strategy_keys(
        self,
        variant_keys: pd.Series,
        selected_strategies: list[str] | tuple[str, ...],
    ) -> set[tuple[str, str]]:
        """Return observed variant/strategy pairs for the requested alleles."""

        selected = tuple(dict.fromkeys(str(value) for value in selected_strategies))
        selected_mask = self.strategy_mask(selected)
        keys = pd.DataFrame(
            {"variant_key": sorted(set(variant_keys.astype(str)))}
        )
        if keys.empty or selected_mask == 0:
            return set()
        duckdb = _import_duckdb()
        with duckdb.connect() as connection:
            connection.register("requested_variants", keys)
            rows = connection.execute(
                "SELECT a.variant_key, a.strategy_mask "
                f"FROM read_parquet({sql_string(self.allele_path)}) a "
                "JOIN requested_variants r USING (variant_key) "
                f"WHERE a.strategy_mask & {selected_mask} != 0"
            ).fetchall()
        strategy_bits = [
            (strategy, 1 << self.strategies.index(strategy))
            for strategy in selected
            if strategy in self.strategies
        ]
        return {
            (str(variant_key), strategy)
            for variant_key, mask in rows
            for strategy, bit in strategy_bits
            if int(mask) & bit
        }

    def strategy_mask(self, selected_strategies: list[str] | tuple[str, ...]) -> int:
        unknown = sorted(set(selected_strategies) - set(self.strategies))
        if unknown:
            raise ValueError(
                "Observed-variant store does not contain strategies: "
                + ", ".join(unknown)
            )
        return sum(1 << self.strategies.index(strategy) for strategy in set(selected_strategies))


def build_or_load_observed_variant_store(
    *,
    variant_annotations_tsv: Path,
    analytics_dir: Path,
    strategies: list[str] | tuple[str, ...],
) -> ObservedVariantStore:
    selected_strategies = tuple(
        sorted({str(value).strip() for value in strategies if str(value).strip()})
    )
    if not selected_strategies:
        raise ValueError("Observed-variant store requires at least one strategy")
    if len(selected_strategies) > MAX_STRATEGIES:
        raise ValueError(
            f"Observed-variant store supports at most {MAX_STRATEGIES} strategies, "
            f"found {len(selected_strategies)}"
        )
    source = resolve_pre_vep_variant_source(
        variant_annotations_tsv,
        required_columns=REQUIRED_COLUMNS,
    )
    outdir = analytics_dir / "observed_variants"
    allele_gene_path = outdir / "allele_gene_memberships.parquet"
    allele_path = outdir / "allele_memberships.parquet"
    manifest_path = outdir / "manifest.json"
    expected_inputs = {
        "schema_version": STORE_SCHEMA_VERSION,
        "source": source.identity,
        "strategies": list(selected_strategies),
    }
    cached = _load_cache(
        allele_gene_path=allele_gene_path,
        allele_path=allele_path,
        manifest_path=manifest_path,
        expected_inputs=expected_inputs,
    )
    if cached is not None:
        return cached

    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".observed_variants.", dir=outdir) as temporary:
        manifest = _build_store(
            source=source,
            allele_gene_path=allele_gene_path,
            allele_path=allele_path,
            temp_dir=Path(temporary),
            expected_inputs=expected_inputs,
            strategies=selected_strategies,
        )
    write_json_atomic(manifest_path, manifest)
    return ObservedVariantStore(
        allele_gene_path=allele_gene_path,
        allele_path=allele_path,
        manifest_path=manifest_path,
        manifest=manifest,
        strategies=tuple(str(value) for value in manifest["strategies"]),
        cache_hit=False,
    )


def _build_store(
    *,
    source: VariantTableSource,
    allele_gene_path: Path,
    allele_path: Path,
    temp_dir: Path,
    expected_inputs: dict[str, object],
    strategies: tuple[str, ...],
) -> dict[str, object]:
    duckdb = _import_duckdb()
    with duckdb.connect() as connection:
        connection.execute(f"SET threads={available_cpu_count()}")
        connection.execute("SET preserve_insertion_order=false")
        connection.execute("SET enable_progress_bar=false")
        connection.execute(f"SET temp_directory={sql_string(temp_dir)}")
        mask_sql = " + ".join(
            "CASE WHEN list_contains(list_transform(string_split(strategies, ','), "
            f"item -> trim(item)), {sql_string(strategy)}) THEN {1 << index} ELSE 0 END"
            for index, strategy in enumerate(strategies)
        )
        connection.execute(
            "CREATE TEMP TABLE allele_gene_source AS SELECT "
            "cast(variant_key AS VARCHAR) AS variant_key, "
            "cast(gene_id AS VARCHAR) AS gene_id, "
            "cast(event_type AS VARCHAR) AS event_type, "
            "cast(ref AS VARCHAR) AS ref, cast(alt AS VARCHAR) AS alt, "
            "cast(lookup_status AS VARCHAR) AS lookup_status, "
            "cast(strategies AS VARCHAR) AS strategies, "
            f"({mask_sql})::UBIGINT AS strategy_mask "
            f"FROM {variant_source_sql(source)}"
        )
        source_row_count, missing_strategies = connection.execute(
            "SELECT count(*), count_if(trim(strategies) = '') FROM allele_gene_source"
        ).fetchone()
        if source.row_count is not None and int(source_row_count) != source.row_count:
            raise ValueError(
                "Observed-variant source row count changed: "
                f"observed {source_row_count}, expected {source.row_count}"
            )
        if int(missing_strategies) > 0:
            raise ValueError(
                f"Variant annotations contain {missing_strategies} rows without a strategy"
            )
        observed_strategies = tuple(
            str(row[0])
            for row in connection.execute(
                "SELECT DISTINCT trim(strategy) FROM ("
                "SELECT unnest(string_split(strategies, ',')) AS strategy "
                "FROM allele_gene_source"
                ") WHERE trim(strategy) <> '' ORDER BY 1"
            ).fetchall()
        )
        if observed_strategies != strategies:
            raise ValueError(
                "Observed-variant source strategies differ from the report contract: "
                f"observed {', '.join(observed_strategies)}, "
                f"expected {', '.join(strategies)}"
            )
        connection.execute(
            "CREATE TEMP VIEW allele_gene_memberships AS SELECT "
            "variant_key, gene_id, event_type, ref, alt, lookup_status, "
            "strategy_mask FROM allele_gene_source"
        )
        connection.execute(
            "CREATE TEMP TABLE allele_memberships AS SELECT variant_key, "
            "bit_or(strategy_mask) AS strategy_mask "
            "FROM allele_gene_source WHERE variant_key <> '' GROUP BY variant_key"
        )
        allele_gene_count = int(
            connection.execute("SELECT count(*) FROM allele_gene_memberships").fetchone()[0]
        )
        allele_count = int(
            connection.execute("SELECT count(*) FROM allele_memberships").fetchone()[0]
        )
        _write_parquet_atomic(
            connection,
            relation="allele_gene_memberships",
            destination=allele_gene_path,
            expected_columns=ALLELE_GENE_COLUMNS,
            expected_rows=allele_gene_count,
        )
        _write_parquet_atomic(
            connection,
            relation="allele_memberships",
            destination=allele_path,
            expected_columns=ALLELE_COLUMNS,
            expected_rows=allele_count,
        )
    return {
        "complete": True,
        "inputs": expected_inputs,
        "source_mode": source.mode,
        "source_row_count": int(source_row_count),
        "allele_gene_count": allele_gene_count,
        "allele_count": allele_count,
        "strategies": list(strategies),
        "outputs": {
            allele_gene_path.name: path_metadata(allele_gene_path),
            allele_path.name: path_metadata(allele_path),
        },
    }


def _load_cache(
    *,
    allele_gene_path: Path,
    allele_path: Path,
    manifest_path: Path,
    expected_inputs: dict[str, object],
) -> ObservedVariantStore | None:
    if not allele_gene_path.exists() or not allele_path.exists() or not manifest_path.exists():
        return None
    try:
        manifest = json.loads(manifest_path.read_text())
        expected_outputs = {
            allele_gene_path.name: path_metadata(allele_gene_path),
            allele_path.name: path_metadata(allele_path),
        }
        if (
            manifest.get("complete") is not True
            or manifest.get("inputs") != expected_inputs
            or manifest.get("outputs") != expected_outputs
        ):
            return None
        strategies = tuple(str(value) for value in manifest.get("strategies", []))
        if not strategies or len(strategies) > MAX_STRATEGIES:
            return None
        _validate_parquet(
            allele_gene_path,
            ALLELE_GENE_COLUMNS,
            int(manifest["allele_gene_count"]),
        )
        _validate_parquet(
            allele_path,
            ALLELE_COLUMNS,
            int(manifest["allele_count"]),
        )
        return ObservedVariantStore(
            allele_gene_path=allele_gene_path,
            allele_path=allele_path,
            manifest_path=manifest_path,
            manifest=manifest,
            strategies=strategies,
            cache_hit=True,
        )
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return None


def _write_parquet_atomic(
    connection,
    *,
    relation: str,
    destination: Path,
    expected_columns: list[str],
    expected_rows: int,
) -> None:
    temporary = _temporary_path(destination)
    try:
        connection.execute(
            f"COPY {relation} TO {sql_string(temporary)} "
            "(FORMAT PARQUET, COMPRESSION ZSTD)"
        )
        _validate_parquet(temporary, expected_columns, expected_rows)
        temporary.chmod(0o644)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def _validate_parquet(path: Path, expected_columns: list[str], expected_rows: int) -> None:
    duckdb = _import_duckdb()
    with duckdb.connect() as connection:
        relation = connection.read_parquet(str(path))
        if relation.columns != expected_columns:
            raise ValueError(
                f"Observed-variant store columns changed in {path}: "
                + ", ".join(relation.columns)
            )
        observed_rows = int(relation.count("*").fetchone()[0])
        if observed_rows != expected_rows:
            raise ValueError(
                f"Observed-variant store row count changed in {path}: "
                f"{observed_rows} != {expected_rows}"
            )


def _temporary_path(destination: Path) -> Path:
    with tempfile.NamedTemporaryFile(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".parquet.tmp",
        delete=False,
    ) as handle:
        path = Path(handle.name)
    path.unlink()
    return path


def available_cpu_count() -> int:
    allocated = os.environ.get("SLURM_CPUS_PER_TASK")
    if allocated and allocated.isdigit() and int(allocated) > 0:
        return int(allocated)
    return os.cpu_count() or 1


def _import_duckdb():
    try:
        import duckdb
    except ImportError as exc:  # pragma: no cover - analytics environment contract
        raise RuntimeError(
            "Observed-variant store requires the python-duckdb package"
        ) from exc
    return duckdb
