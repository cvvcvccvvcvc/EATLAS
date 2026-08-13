"""DuckDB aggregation primitives for large variant-annotation tables."""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from analytics.analyses.target_context import read_disjoint_contexts
from analytics.annotation.consequences import UNANNOTATED_CONSEQUENCE
from analytics.io.artifacts import file_identity, path_metadata
from analytics.io.variant_source import cohort_variant_paths
from genomics.variants import read_failed_regions


REQUIRED_COLUMNS = {
    "variant_key",
    "gene_id",
    "event_type",
    "strategies",
    "vep_status",
    "vep_primary_consequence",
}
MAX_STRATEGIES = 63
DUCKDB_MEMORY_LIMIT_ENV = "GAPH_DUCKDB_MEMORY_LIMIT"
DUCKDB_MEMORY_FRACTION = 0.5

_MEMORY_SETTING = re.compile(r"^([0-9]+(?:\.[0-9]+)?)\s*([KMGT]i?B)$", re.IGNORECASE)
_MEMORY_UNITS = {
    "KB": 1000,
    "MB": 1000**2,
    "GB": 1000**3,
    "TB": 1000**4,
    "KIB": 1024,
    "MIB": 1024**2,
    "GIB": 1024**3,
    "TIB": 1024**4,
}


@dataclass(frozen=True)
class VariantAggregationSource:
    paths: tuple[Path, ...]
    columns: tuple[str, ...]
    row_count: int | None
    partitioned: bool
    identity: dict[str, object]


@dataclass(frozen=True)
class StrategyMaskAggregation:
    input_row_count: int
    missing_variant_key_count: int
    missing_strategy_count: int
    strategies: tuple[str, ...]
    allele_gene_mask_counts: dict[int, int]
    allele_mask_counts: dict[int, int]

    @property
    def unique_variant_count(self) -> int:
        return sum(self.allele_mask_counts.values())

    @property
    def strategy_record_count(self) -> int:
        return sum(mask.bit_count() * count for mask, count in self.allele_mask_counts.items())

    @property
    def all_strategy_variant_count(self) -> int:
        all_mask = (1 << len(self.strategies)) - 1
        return self.allele_mask_counts.get(all_mask, 0)

    def strategy_counts(self) -> dict[str, int]:
        return {
            strategy: sum(
                count
                for mask, count in self.allele_mask_counts.items()
                if mask & (1 << index)
            )
            for index, strategy in enumerate(self.strategies)
        }

    def unique_strategy_counts(self) -> dict[str, int]:
        return {
            strategy: self.allele_mask_counts.get(1 << index, 0)
            for index, strategy in enumerate(self.strategies)
        }

    def intersections(self) -> list[list[int]]:
        return [
            [
                sum(
                    count
                    for mask, count in self.allele_mask_counts.items()
                    if mask & (1 << left_index) and mask & (1 << right_index)
                )
                for right_index in range(len(self.strategies))
            ]
            for left_index in range(len(self.strategies))
        ]


@dataclass(frozen=True)
class VariantGroupedAggregation:
    masks: StrategyMaskAggregation
    consequence_source: str
    gene_count: int
    global_groups: pd.DataFrame
    allele_gene_groups: pd.DataFrame
    gnomad_af_summary: pd.DataFrame
    pathogenic_rows: pd.DataFrame
    ortholog_evidence_grouped: pd.DataFrame
    ortholog_distribution_source: pd.DataFrame
    timings: dict[str, float]
    diagnostics: dict[str, object]


def available_cpu_count() -> int:
    """Return CPUs available to this process, respecting Slurm allocation."""

    allocated = os.environ.get("SLURM_CPUS_PER_TASK")
    if allocated and allocated.isdigit() and int(allocated) > 0:
        return int(allocated)
    return os.cpu_count() or 1


def _configure_duckdb_memory(connection, thread_count: int) -> dict[str, object]:
    override = os.environ.get(DUCKDB_MEMORY_LIMIT_ENV, "").strip()
    if override:
        requested = override
        source = DUCKDB_MEMORY_LIMIT_ENV
    else:
        slurm_bytes, source = _slurm_memory_bytes(thread_count)
        if slurm_bytes is not None:
            requested = _memory_limit_setting(slurm_bytes * DUCKDB_MEMORY_FRACTION)
        else:
            current = str(
                connection.execute("SELECT current_setting('memory_limit')").fetchone()[0]
            )
            requested = _memory_limit_setting(
                _parse_memory_setting(current) * DUCKDB_MEMORY_FRACTION
            )
            source = "duckdb_default_fraction"

    connection.execute(f"SET memory_limit={_sql_string(requested)}")
    return {
        "memory_limit": str(
            connection.execute("SELECT current_setting('memory_limit')").fetchone()[0]
        ),
        "memory_limit_source": source,
    }


def _slurm_memory_bytes(thread_count: int) -> tuple[int | None, str]:
    per_node_mb = _positive_int_environment("SLURM_MEM_PER_NODE")
    if per_node_mb is not None:
        return per_node_mb * 1024**2, "SLURM_MEM_PER_NODE"

    per_cpu_mb = _positive_int_environment("SLURM_MEM_PER_CPU")
    if per_cpu_mb is None:
        return None, "duckdb_default_fraction"
    allocated_cpus = (
        _positive_int_environment("SLURM_CPUS_PER_TASK") or thread_count
    )
    return per_cpu_mb * allocated_cpus * 1024**2, "SLURM_MEM_PER_CPU"


def _positive_int_environment(name: str) -> int | None:
    value = os.environ.get(name, "").strip()
    if not value.isdigit() or int(value) < 1:
        return None
    return int(value)


def _parse_memory_setting(value: str) -> int:
    match = _MEMORY_SETTING.fullmatch(value.strip())
    if match is None:
        raise ValueError(f"Could not parse DuckDB memory limit: {value!r}")
    amount, unit = match.groups()
    return int(float(amount) * _MEMORY_UNITS[unit.upper()])


def _memory_limit_setting(value: float) -> str:
    mebibytes = max(128, int(value) // (1024**2))
    return f"{mebibytes}MiB"


def resolve_variant_aggregation_source(path: Path) -> VariantAggregationSource:
    """Resolve a validated finalized VEP artifact."""

    path = path.resolve()
    cohort_paths = cohort_variant_paths(path)
    if cohort_paths is not None:
        sources = [resolve_variant_aggregation_source(member) for member in cohort_paths]
        columns = sources[0].columns
        if any(source.columns != columns for source in sources[1:]):
            raise ValueError("Cohort VEP outputs have different table columns")
        partitioned = sources[0].partitioned
        if any(source.partitioned != partitioned for source in sources[1:]):
            raise ValueError("Cohort VEP outputs mix partitioned and merged source modes")
        row_count = (
            sum(int(source.row_count) for source in sources)
            if all(source.row_count is not None for source in sources)
            else None
        )
        return VariantAggregationSource(
            paths=tuple(item for source in sources for item in source.paths),
            columns=columns,
            row_count=row_count,
            partitioned=partitioned,
            identity={
                "cohort_descriptor": path_metadata(path),
                "members": [source.identity for source in sources],
            },
        )
    artifact_dir = path.parent
    plan_path = artifact_dir / "plan.json"
    manifest_path = artifact_dir / "manifest.json"
    partitions_dir = artifact_dir / "partitions"
    if path.name == "variant_annotations.vep.tsv.gz" and plan_path.exists():
        if not manifest_path.exists():
            raise ValueError(f"Incomplete partitioned VEP artifact under {artifact_dir}")
        plan = _read_json(plan_path)
        manifest = _read_json(manifest_path)
        if plan.get("status") != "complete" or manifest.get("status") != "complete":
            raise ValueError(f"Incomplete partitioned VEP artifact under {artifact_dir}")
        columns = tuple(str(column) for column in plan.get("output_columns", []))
        _require_columns(columns, path)
        if list(columns) != list(manifest.get("columns", [])):
            raise ValueError("VEP plan and final manifest columns differ")
        row_count = int(plan.get("row_count", -1))
        if row_count < 0 or row_count != int(manifest.get("row_count", -2)):
            raise ValueError("VEP plan and final manifest row counts differ")

        paths = []
        partition_identities = []
        observed_rows = 0
        for entry in plan.get("partitions", []):
            partition_id = str(entry.get("partition_id", ""))
            partition_manifest_path = partitions_dir / f"{partition_id}.json"
            partition_path = partitions_dir / f"{partition_id}.tsv.gz"
            partition_manifest = _read_json(partition_manifest_path)
            if partition_manifest.get("status") != "complete":
                raise ValueError(f"Incomplete VEP partition: {partition_id}")
            if partition_manifest.get("input") != entry:
                raise ValueError(f"VEP partition input contract changed: {partition_id}")
            if list(partition_manifest.get("output_columns", [])) != list(columns):
                raise ValueError(f"VEP partition columns changed: {partition_id}")
            expected_output = dict(partition_manifest.get("output", {}))
            if not partition_path.exists() or file_identity(partition_path) != expected_output:
                raise ValueError(f"VEP partition output changed: {partition_id}")
            partition_rows = int(partition_manifest.get("row_count", -1))
            if partition_rows != int(entry.get("row_count", -2)):
                raise ValueError(f"VEP partition row count changed: {partition_id}")
            observed_rows += partition_rows
            paths.append(partition_path)
            partition_identities.append(path_metadata(partition_manifest_path))
        if not paths or observed_rows != row_count:
            raise ValueError("VEP partition set does not match the declared row count")
        return VariantAggregationSource(
            paths=tuple(paths),
            columns=columns,
            row_count=row_count,
            partitioned=True,
            identity={
                "plan": path_metadata(plan_path),
                "manifest": path_metadata(manifest_path),
                "partition_manifests": partition_identities,
            },
        )

    if not path.exists():
        raise FileNotFoundError(path)
    columns = tuple(_read_header(path))
    _require_columns(columns, path)
    return VariantAggregationSource(
        paths=(path,),
        columns=columns,
        row_count=None,
        partitioned=False,
        identity={"input": path_metadata(path)},
    )


def aggregate_strategy_masks(
    source: VariantAggregationSource,
    *,
    threads: int | None = None,
    temp_dir: Path | None = None,
) -> StrategyMaskAggregation:
    """Aggregate global alleles and allele-gene rows without strategy expansion."""

    duckdb = _import_duckdb()
    thread_count = available_cpu_count() if threads is None else threads
    if thread_count < 1:
        raise ValueError("DuckDB thread count must be >= 1")
    connection = duckdb.connect()
    try:
        connection.execute(f"SET threads={int(thread_count)}")
        connection.execute("SET preserve_insertion_order=false")
        if temp_dir is not None:
            temp_dir.mkdir(parents=True, exist_ok=True)
            connection.execute(f"SET temp_directory={_sql_string(temp_dir)}")
        connection.execute(f"CREATE VIEW source_rows AS SELECT * FROM {_source_sql(source)}")
        input_row_count, missing_keys, missing_strategies = connection.execute(
            "SELECT count(*), count_if(variant_key = ''), count_if(strategies = '') "
            "FROM source_rows"
        ).fetchone()
        if source.row_count is not None and int(input_row_count) != source.row_count:
            raise ValueError(
                "VEP source row count changed: "
                f"observed {input_row_count}, expected {source.row_count}"
            )
        strategies = tuple(
            str(row[0])
            for row in connection.execute(
                "SELECT DISTINCT trim(strategy) FROM ("
                "SELECT unnest(string_split(strategies, ',')) AS strategy FROM source_rows"
                ") WHERE trim(strategy) <> '' ORDER BY 1"
            ).fetchall()
        )
        if not strategies:
            raise ValueError("Variant annotations contain no strategy memberships")
        if len(strategies) > MAX_STRATEGIES:
            raise ValueError(
                f"Variant summary supports at most {MAX_STRATEGIES} strategies, "
                f"found {len(strategies)}"
            )
        if int(missing_strategies) > 0:
            raise ValueError(
                f"Variant annotations contain {missing_strategies} rows without a strategy"
            )
        mask_sql = " + ".join(
            "CASE WHEN list_contains(list_transform(string_split(strategies, ','), "
            "item -> trim(item)), "
            f"{_sql_string(strategy)}) THEN {1 << index} ELSE 0 END"
            for index, strategy in enumerate(strategies)
        )
        ref = "ref" if "ref" in source.columns else "''"
        alt = "alt" if "alt" in source.columns else "''"
        variant_id = (
            "CASE WHEN variant_key <> '' THEN variant_key ELSE "
            f"gene_id || ':' || event_type || ':' || {ref} || '>' || {alt} END"
        )
        connection.execute(
            "CREATE VIEW masked_rows AS SELECT "
            f"{variant_id} AS variant_id, gene_id, ({mask_sql})::UBIGINT AS strategy_mask "
            "FROM source_rows"
        )
        allele_gene_rows = connection.execute(
            "SELECT strategy_mask, count(*) FROM ("
            "SELECT variant_id, gene_id, bit_or(strategy_mask) AS strategy_mask "
            "FROM masked_rows GROUP BY variant_id, gene_id"
            ") GROUP BY strategy_mask ORDER BY strategy_mask"
        ).fetchall()
        allele_rows = connection.execute(
            "SELECT strategy_mask, count(*) FROM ("
            "SELECT variant_id, bit_or(strategy_mask) AS strategy_mask "
            "FROM masked_rows GROUP BY variant_id"
            ") GROUP BY strategy_mask ORDER BY strategy_mask"
        ).fetchall()
        return StrategyMaskAggregation(
            input_row_count=int(input_row_count),
            missing_variant_key_count=int(missing_keys),
            missing_strategy_count=int(missing_strategies),
            strategies=strategies,
            allele_gene_mask_counts={int(mask): int(count) for mask, count in allele_gene_rows},
            allele_mask_counts={int(mask): int(count) for mask, count in allele_rows},
        )
    finally:
        connection.close()


def aggregate_variant_groups(
    source: VariantAggregationSource,
    *,
    genes_path: Path | None = None,
    target_features_path: Path | None = None,
    annotation_failures_path: Path | None = None,
    variant_strategy_support_path: Path | None = None,
    threads: int | None = None,
    temp_dir: Path | None = None,
) -> VariantGroupedAggregation:
    """Build compact global-allele and allele-gene grouped relations."""

    duckdb = _import_duckdb()
    thread_count = available_cpu_count() if threads is None else threads
    if thread_count < 1:
        raise ValueError("DuckDB thread count must be >= 1")
    timings: dict[str, float] = {}
    diagnostics: dict[str, object] = {}
    connection = duckdb.connect()
    try:
        connection.execute(f"SET threads={int(thread_count)}")
        connection.execute("SET preserve_insertion_order=false")
        if temp_dir is not None:
            temp_dir.mkdir(parents=True, exist_ok=True)
            connection.execute(f"SET temp_directory={_sql_string(temp_dir)}")
        diagnostics.update(_configure_duckdb_memory(connection, thread_count))
        diagnostics["max_temp_directory_size"] = str(
            connection.execute(
                "SELECT current_setting('max_temp_directory_size')"
            ).fetchone()[0]
        )
        connection.execute(f"CREATE VIEW source_rows AS SELECT * FROM {_source_sql(source)}")

        started = time.perf_counter()
        input_row_count, missing_keys, missing_strategies, gene_count = connection.execute(
            "SELECT count(*), count_if(variant_key = ''), count_if(strategies = ''), "
            "count(DISTINCT gene_id) FROM source_rows"
        ).fetchone()
        strategies = _read_strategies(connection)
        timings["source_scan"] = time.perf_counter() - started
        if source.row_count is not None and int(input_row_count) != source.row_count:
            raise ValueError(
                "VEP source row count changed: "
                f"observed {input_row_count}, expected {source.row_count}"
            )
        if int(missing_strategies) > 0:
            raise ValueError(
                f"Variant annotations contain {missing_strategies} rows without a strategy"
            )

        _register_reference_tables(
            connection,
            genes_path=genes_path,
            target_features_path=target_features_path,
            annotation_failures_path=annotation_failures_path,
        )
        started = time.perf_counter()
        _create_normalized_views(connection, source, strategies, genes_path is not None)
        timings["normalization_setup"] = time.perf_counter() - started

        started = time.perf_counter()
        allele_gene_row_count, global_allele_row_count = _materialize_compacted_relations(
            connection
        )
        timings["compacted_relations"] = time.perf_counter() - started
        diagnostics.update(
            {
                "allele_gene_row_count": allele_gene_row_count,
                "global_allele_row_count": global_allele_row_count,
                "temp_storage_bytes_after_materialization": _directory_size(temp_dir),
            }
        )

        started = time.perf_counter()
        allele_gene_mask_rows = connection.execute(
            "SELECT strategy_mask, count(*) FROM allele_gene_rows "
            "GROUP BY strategy_mask ORDER BY strategy_mask"
        ).fetchall()
        global_groups = connection.execute(
            "SELECT strategy_mask, event_type, clinvar_found, clinvar_classified, "
            "clinvar_category, gnomad_status, titv_kind, review_stars, count(*) AS variant_count "
            "FROM global_alleles GROUP BY ALL ORDER BY ALL"
        ).fetchdf()
        allele_mask_counts = {
            int(row.strategy_mask): int(row.variant_count)
            for row in global_groups.groupby("strategy_mask", as_index=False)["variant_count"]
            .sum()
            .itertuples(index=False)
        }
        timings["global_aggregates"] = time.perf_counter() - started

        started = time.perf_counter()
        allele_gene_groups = connection.execute(
            "SELECT strategy_mask, gene_id, event_type, target_context, gnomad_status, "
            "consequence, clinvar_category, count(*) AS variant_count "
            "FROM allele_gene_rows GROUP BY ALL ORDER BY ALL"
        ).fetchdf()
        timings["allele_gene_aggregates"] = time.perf_counter() - started

        started = time.perf_counter()
        gnomad_af_summary = _query_gnomad_af_summary(connection, strategies)
        timings["gnomad_af_quantiles"] = time.perf_counter() - started

        started = time.perf_counter()
        pathogenic_rows = _query_pathogenic_rows(connection, strategies)
        timings["pathogenic_rows"] = time.perf_counter() - started

        started = time.perf_counter()
        ortholog_grouped, ortholog_distributions = _query_ortholog_evidence(
            connection,
            variant_strategy_support_path,
            strategies,
        )
        timings["ortholog_evidence"] = time.perf_counter() - started
        diagnostics["temp_storage_bytes_final"] = _directory_size(temp_dir)

        masks = StrategyMaskAggregation(
            input_row_count=int(input_row_count),
            missing_variant_key_count=int(missing_keys),
            missing_strategy_count=int(missing_strategies),
            strategies=strategies,
            allele_gene_mask_counts={
                int(mask): int(count) for mask, count in allele_gene_mask_rows
            },
            allele_mask_counts=allele_mask_counts,
        )
        return VariantGroupedAggregation(
            masks=masks,
            consequence_source="Ensembl VEP",
            gene_count=int(gene_count),
            global_groups=global_groups,
            allele_gene_groups=allele_gene_groups,
            gnomad_af_summary=gnomad_af_summary,
            pathogenic_rows=pathogenic_rows,
            ortholog_evidence_grouped=ortholog_grouped,
            ortholog_distribution_source=ortholog_distributions,
            timings=timings,
            diagnostics=diagnostics,
        )
    finally:
        connection.close()


def _read_strategies(connection) -> tuple[str, ...]:
    strategies = tuple(
        str(row[0])
        for row in connection.execute(
            "SELECT DISTINCT trim(strategy) FROM ("
            "SELECT unnest(string_split(strategies, ',')) AS strategy FROM source_rows"
            ") WHERE trim(strategy) <> '' ORDER BY 1"
        ).fetchall()
    )
    if not strategies:
        raise ValueError("Variant annotations contain no strategy memberships")
    if len(strategies) > MAX_STRATEGIES:
        raise ValueError(
            f"Variant summary supports at most {MAX_STRATEGIES} strategies, "
            f"found {len(strategies)}"
        )
    return strategies


def _register_reference_tables(
    connection,
    *,
    genes_path: Path | None,
    target_features_path: Path | None,
    annotation_failures_path: Path | None,
) -> None:
    if genes_path is not None:
        genes = pd.read_csv(
            genes_path,
            sep="\t",
            compression="gzip" if genes_path.suffix == ".gz" else None,
            usecols=["gene_id", "begin"],
            dtype={"gene_id": str, "begin": "int64"},
        ).rename(columns={"begin": "gene_begin"})
    else:
        genes = pd.DataFrame(columns=["gene_id", "gene_begin"])
    connection.register("gene_begins_input", genes)
    connection.execute(
        "CREATE VIEW gene_begins AS SELECT cast(gene_id AS VARCHAR) AS gene_id, "
        "cast(gene_begin AS BIGINT) AS gene_begin FROM gene_begins_input"
    )

    context_rows: list[dict[str, object]] = []
    if target_features_path is not None:
        lengths = _gene_lengths_from_features(target_features_path)
        for gene_id, intervals in read_disjoint_contexts(target_features_path, lengths).items():
            context_rows.extend(
                {
                    "gene_id": gene_id,
                    "start0": start,
                    "end0": end,
                    "target_context": context,
                }
                for start, end, context in intervals
            )
    contexts = pd.DataFrame(
        context_rows,
        columns=["gene_id", "start0", "end0", "target_context"],
    )
    connection.register("target_contexts_input", contexts)
    connection.execute(
        "CREATE VIEW target_contexts AS SELECT cast(gene_id AS VARCHAR) AS gene_id, "
        "cast(start0 AS BIGINT) AS start0, cast(end0 AS BIGINT) AS end0, "
        "cast(target_context AS VARCHAR) AS target_context FROM target_contexts_input"
    )

    failure_rows = []
    for chrom, (_starts, intervals) in read_failed_regions(
        annotation_failures_path, "gnomad"
    ).items():
        failure_rows.extend(
            {"chrom": chrom, "start1": start, "end1": end}
            for start, end in intervals
        )
    failures = pd.DataFrame(failure_rows, columns=["chrom", "start1", "end1"])
    connection.register("gnomad_failures_input", failures)
    connection.execute(
        "CREATE VIEW gnomad_failures AS SELECT cast(chrom AS VARCHAR) AS chrom, "
        "cast(start1 AS BIGINT) AS start1, cast(end1 AS BIGINT) AS end1 "
        "FROM gnomad_failures_input"
    )


def _create_normalized_views(
    connection,
    source: VariantAggregationSource,
    strategies: tuple[str, ...],
    has_gene_begins: bool,
) -> None:
    column = lambda name, default="''": f"s.{name}" if name in source.columns else default
    mask_sql = " + ".join(
        "CASE WHEN list_contains(list_transform(string_split(s.strategies, ','), "
        f"item -> trim(item)), {_sql_string(strategy)}) THEN {1 << index} ELSE 0 END"
        for index, strategy in enumerate(strategies)
    )
    ref = column("ref")
    alt = column("alt")
    variant_id = (
        "CASE WHEN s.variant_key <> '' THEN s.variant_key ELSE "
        f"s.gene_id || ':' || s.event_type || ':' || {ref} || '>' || {alt} END"
    )
    selected = ", ".join(
        f"{column(name)} AS {name}"
        for name in _pathogenic_source_columns()
    )
    connection.execute(
        "CREATE VIEW keyed_rows AS SELECT "
        f"{selected}, {variant_id} AS variant_id, ({mask_sql})::UBIGINT AS strategy_mask, "
        "split_part(s.variant_key, ':', 1) AS key_chrom, "
        "try_cast(split_part(s.variant_key, ':', 2) AS BIGINT) AS key_pos, "
        "split_part(split_part(s.variant_key, ':', 3), '>', 1) AS key_ref, "
        "split_part(split_part(s.variant_key, ':', 3), '>', 2) AS key_alt "
        "FROM source_rows s"
    )
    prefix = _common_prefix_sql("k.key_ref", "k.key_alt")
    if has_gene_begins:
        affected_start = f"k.key_pos - g.gene_begin + ({prefix})"
        gene_join = "LEFT JOIN gene_begins g USING (gene_id)"
    elif "target_start0" in source.columns:
        affected_start = "try_cast(k.target_start0 AS BIGINT)"
        gene_join = ""
    else:
        affected_start = "NULL::BIGINT"
        gene_join = ""
    connection.execute(
        "CREATE VIEW positioned_rows AS SELECT k.*, "
        f"{affected_start} AS affected_start0 FROM keyed_rows k {gene_join}"
    )

    evidence_columns = [
        name
        for name in (
            "clinvar_id",
            "clinvar_allele_id",
            "clinvar_sig",
            "clinvar_revstat",
            "clinvar_hgvs",
            "clinvar_disease",
            "clinvar_variant_type",
        )
        if name in source.columns
    ]
    clinvar_found = (
        " OR ".join(f"coalesce(p.{name}, '') <> ''" for name in evidence_columns)
        if evidence_columns
        else "false"
    )
    sig = "lower(coalesce(p.clinvar_sig, ''))"
    consequence = (
        "CASE WHEN p.vep_status = 'ok' "
        "AND coalesce(p.vep_primary_consequence, '') <> '' "
        "THEN p.vep_primary_consequence ELSE "
        f"{_sql_string(UNANNOTATED_CONSEQUENCE)} END"
    )
    failed = (
        "EXISTS (SELECT 1 FROM gnomad_failures f WHERE f.chrom = p.key_chrom "
        "AND p.key_pos BETWEEN f.start1 AND f.end1)"
    )
    context_join = (
        "LEFT JOIN target_contexts c ON c.gene_id = p.gene_id "
        "AND p.affected_start0 >= c.start0 AND p.affected_start0 < c.end0"
    )
    connection.execute(
        "CREATE VIEW normalized_rows AS SELECT p.*, "
        f"({clinvar_found}) AS clinvar_found, "
        "coalesce(p.clinvar_sig, '') <> '' AS clinvar_classified, "
        "CASE "
        f"WHEN NOT ({clinvar_found}) THEN 'Not in ClinVar' "
        f"WHEN {sig} = '' THEN 'Unclassified' "
        f"WHEN contains({sig}, 'conflicting') THEN 'Other' "
        f"WHEN regexp_matches({sig}, 'uncertain|vus') THEN 'VUS' "
        f"WHEN contains({sig}, 'pathogenic') AND NOT contains({sig}, 'benign') THEN 'P/LP' "
        f"WHEN contains({sig}, 'benign') AND NOT contains({sig}, 'pathogenic') THEN 'B/LB' "
        "ELSE 'Other' END AS clinvar_category, "
        "try_cast(nullif(p.gnomad_af, '') AS DOUBLE) AS gnomad_af_value, "
        "CASE WHEN try_cast(nullif(p.gnomad_af, '') AS DOUBLE) IS NOT NULL THEN 'found' "
        "WHEN coalesce(p.lookup_status, '') NOT IN ('', 'ok') OR p.key_pos IS NULL THEN 'lookup_failed' "
        f"WHEN {failed} THEN 'lookup_failed' ELSE 'not_found' END AS gnomad_status, "
        "CASE WHEN p.event_type = 'snv' AND length(p.ref) = 1 AND length(p.alt) = 1 "
        "THEN CASE WHEN p.ref || '>' || p.alt IN ('A>G','G>A','C>T','T>C') "
        "THEN 'ti' ELSE 'tv' END ELSE '' END AS titv_kind, "
        "CASE WHEN p.clinvar_review_stars IN ('0','1','2','3','4') "
        "THEN p.clinvar_review_stars ELSE 'Unmapped' END AS review_stars, "
        f"{consequence} AS consequence, CASE WHEN p.affected_start0 IS NULL THEN 'unknown' "
        "ELSE coalesce(c.target_context, 'other') END AS target_context "
        f"FROM positioned_rows p {context_join}"
    )


def _materialize_compacted_relations(connection) -> tuple[int, int]:
    connection.execute(
        "CREATE TEMP TABLE allele_gene_rows AS SELECT variant_id, gene_id, "
        "bit_or(strategy_mask) AS strategy_mask, first(event_type) AS event_type, "
        "first(target_context) AS target_context, first(gnomad_status) AS gnomad_status, "
        "first(consequence) AS consequence, first(clinvar_category) AS clinvar_category, "
        "bool_or(clinvar_found) AS global_clinvar_found, "
        "bool_or(clinvar_classified) AS global_clinvar_classified, "
        "max(gnomad_af_value) AS global_gnomad_af, "
        "bool_or(gnomad_status = 'found') AS global_gnomad_found, "
        "bool_or(gnomad_status = 'lookup_failed') AS global_gnomad_failed, "
        "first(titv_kind) AS global_titv_kind, "
        "first(review_stars) AS global_review_stars "
        "FROM normalized_rows GROUP BY variant_id, gene_id"
    )
    connection.execute(
        "CREATE TEMP TABLE global_alleles AS SELECT variant_id, "
        "bit_or(strategy_mask) AS strategy_mask, first(event_type) AS event_type, "
        "bool_or(global_clinvar_found) AS clinvar_found, "
        "bool_or(global_clinvar_classified) AS clinvar_classified, "
        "first(clinvar_category) AS clinvar_category, "
        "max(global_gnomad_af) AS gnomad_af, "
        "CASE WHEN bool_or(global_gnomad_found) THEN 'found' "
        "WHEN bool_or(global_gnomad_failed) THEN 'lookup_failed' ELSE 'not_found' END "
        "AS gnomad_status, first(global_titv_kind) AS titv_kind, "
        "first(global_review_stars) AS review_stars "
        "FROM allele_gene_rows GROUP BY variant_id"
    )
    allele_gene_count, global_allele_count = connection.execute(
        "SELECT (SELECT count(*) FROM allele_gene_rows), "
        "(SELECT count(*) FROM global_alleles)"
    ).fetchone()
    return int(allele_gene_count), int(global_allele_count)


def _query_gnomad_af_summary(connection, strategies: tuple[str, ...]) -> pd.DataFrame:
    expressions = []
    for index, strategy in enumerate(strategies):
        bit = 1 << index
        membership = f"strategy_mask & {bit} != 0"
        positive = f"{membership} AND gnomad_af > 0"
        expressions.extend(
            [
                f"count(*) FILTER (WHERE {positive}) AS {_sql_identifier(f'{strategy}_count')}",
                "median(gnomad_af) FILTER (WHERE "
                f"{membership} AND gnomad_af IS NOT NULL) AS "
                f"{_sql_identifier(f'{strategy}_median_af')}",
                "quantile_cont(log10(nullif(gnomad_af, 0)), [0.05,0.25,0.5,0.75,0.95]) "
                f"FILTER (WHERE {positive}) AS {_sql_identifier(f'{strategy}_quantiles')}",
            ]
        )
    row = connection.execute("SELECT " + ",".join(expressions) + " FROM global_alleles").fetchone()
    output = []
    offset = 0
    for strategy in strategies:
        count, median_af, quantiles = row[offset : offset + 3]
        offset += 3
        if not count:
            continue
        output.append(
            {
                "strategy": strategy,
                "Count": int(count),
                "Median gnomAD AF": float(median_af),
                "Q05": float(quantiles[0]),
                "Q25": float(quantiles[1]),
                "Median": float(quantiles[2]),
                "Q75": float(quantiles[3]),
                "Q95": float(quantiles[4]),
            }
        )
    return pd.DataFrame(output)


def _query_pathogenic_rows(
    connection,
    strategies: tuple[str, ...],
) -> pd.DataFrame:
    columns = _pathogenic_source_columns()
    aggregate_columns = []
    for name in columns:
        if name == "gene_id":
            aggregate_columns.append(
                "string_agg(DISTINCT gene_id, ', ' ORDER BY gene_id) AS gene_id"
            )
        elif name == "strategies":
            continue
        elif name in {"variant_key", "clinvar_scv_count", "clinvar_review_stars"}:
            aggregate_columns.append(f"max({name}) AS {name}")
        else:
            aggregate_columns.append(f"first({name}) AS {name}")
    rows = connection.execute(
        "SELECT " + ",".join(aggregate_columns) + ", variant_id, "
        "bit_or(strategy_mask) AS strategy_mask, first(clinvar_category) AS clinvar_category, "
        "first(review_stars) AS review_stars "
        "FROM normalized_rows WHERE clinvar_category = 'P/LP' GROUP BY variant_id "
        "ORDER BY try_cast(first(review_stars) AS INTEGER) DESC NULLS LAST, "
        "try_cast(max(clinvar_scv_count) AS BIGINT) DESC NULLS LAST, variant_id LIMIT 100"
    ).fetchdf()
    if rows.empty:
        return pd.DataFrame(columns=[*columns, "variant_id", "clinvar_category", "review_stars"])
    rows["strategies"] = [
        ",".join(
            strategy
            for index, strategy in enumerate(strategies)
            if int(mask) & (1 << index)
        )
        for mask in rows.pop("strategy_mask")
    ]
    for name in columns:
        if name not in rows:
            rows[name] = ""
    return rows[[*columns, "variant_id", "clinvar_category", "review_stars"]]


def _query_ortholog_evidence(
    connection,
    path: Path | None,
    strategies: tuple[str, ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    grouped_columns = [
        "strategy",
        "target_context",
        "site_depth",
        "alt_support",
        "gnomad_found_count",
        "gnomad_eligible_count",
    ]
    distribution_columns = ["strategy", "site_depth", "alt_support", "variant_count"]
    if path is None:
        return pd.DataFrame(columns=grouped_columns), pd.DataFrame(columns=distribution_columns)
    header = _read_header(path)
    required = {
        "variant_key",
        "gene_id",
        "strategy",
        "alt_support_ortholog_count",
        "site_aligned_ortholog_count",
    }
    if not required.issubset(header):
        return pd.DataFrame(columns=grouped_columns), pd.DataFrame(columns=distribution_columns)
    support_source = VariantAggregationSource(
        paths=(path.resolve(),),
        columns=tuple(header),
        row_count=None,
        partitioned=False,
        identity={"input": path_metadata(path)},
    )
    bits = "CASE " + " ".join(
        f"WHEN strategy = {_sql_string(strategy)} THEN {1 << index}"
        for index, strategy in enumerate(strategies)
    ) + " ELSE 0 END"
    connection.execute(
        "CREATE VIEW ortholog_support_rows AS SELECT variant_key AS variant_id, gene_id, strategy, "
        "try_cast(alt_support_ortholog_count AS BIGINT) AS alt_support, "
        "try_cast(site_aligned_ortholog_count AS BIGINT) AS site_depth, "
        f"({bits})::UBIGINT AS strategy_bit FROM {_source_sql(support_source)}"
    )
    invalid = connection.execute(
        "SELECT count(*) FROM ortholog_support_rows WHERE site_depth IS NULL OR site_depth <= 0 "
        "OR alt_support IS NULL OR alt_support < 0 OR alt_support > site_depth"
    ).fetchone()[0]
    if int(invalid) > 0:
        raise ValueError(f"Variant strategy support contains {invalid} invalid ortholog counts")
    grouped = connection.execute(
        "SELECT s.strategy, a.target_context, s.site_depth, s.alt_support, "
        "count_if(a.gnomad_status = 'found') AS gnomad_found_count, "
        "count(*) AS gnomad_eligible_count "
        "FROM ortholog_support_rows s JOIN allele_gene_rows a "
        "ON a.variant_id = s.variant_id AND a.gene_id = s.gene_id "
        "AND a.strategy_mask & s.strategy_bit != 0 "
        "WHERE a.event_type = 'snv' AND a.target_context IN ('cds','utr','intron') "
        "AND a.gnomad_status IN ('found','not_found') "
        "GROUP BY ALL ORDER BY ALL"
    ).fetchdf()
    distributions = connection.execute(
        "SELECT s.strategy, s.site_depth, s.alt_support, count(*) AS variant_count "
        "FROM ortholog_support_rows s JOIN allele_gene_rows a "
        "ON a.variant_id = s.variant_id AND a.gene_id = s.gene_id "
        "AND a.strategy_mask & s.strategy_bit != 0 "
        "WHERE a.event_type = 'snv' AND a.target_context IN ('cds','utr','intron') "
        "GROUP BY ALL ORDER BY ALL"
    ).fetchdf()
    return grouped, distributions


def _pathogenic_source_columns() -> list[str]:
    return [
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
        "vep_status",
        "vep_primary_consequence",
    ]


def _common_prefix_sql(ref: str, alt: str) -> str:
    mismatches = (
        "list_transform(range(1, least(length("
        f"{ref}), length({alt})) + 1), i -> substr({ref}, i, 1) <> substr({alt}, i, 1))"
    )
    return (
        f"CASE WHEN list_position({mismatches}, true) IS NULL "
        f"THEN least(length({ref}), length({alt})) "
        f"ELSE list_position({mismatches}, true) - 1 END"
    )


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
        raise ValueError("Target features contain no gene rows for target-context assignment")
    return {
        str(row.gene_id): int(row.target_end0) - int(row.target_start0)
        for row in genes.itertuples(index=False)
    }


def _sql_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def _source_sql(source: VariantAggregationSource) -> str:
    columns = "{" + ",".join(
        f"{_sql_string(column)}: 'VARCHAR'" for column in source.columns
    ) + "}"
    paths = "[" + ",".join(_sql_string(path) for path in source.paths) + "]"
    header = "false" if source.partitioned else "true"
    return (
        f"read_csv({paths}, delim='\\t', header={header}, columns={columns}, "
        "auto_detect=false, compression='auto', parallel=true, "
        "nullstr='__GAPH_NULL_SENTINEL__')"
    )


def _read_header(path: Path) -> list[str]:
    handle = gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")
    with handle:
        return next(csv.reader(handle, delimiter="\t"))


def _require_columns(columns: tuple[str, ...], path: Path) -> None:
    missing = REQUIRED_COLUMNS - set(columns)
    if missing:
        raise ValueError(
            f"Variant annotations {path} missing columns: {', '.join(sorted(missing))}"
        )


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _directory_size(path: Path | None) -> int:
    if path is None or not path.exists():
        return 0
    total = 0
    for item in path.rglob("*"):
        try:
            if item.is_file():
                total += item.stat().st_size
        except FileNotFoundError:
            continue
    return total


def _sql_string(value: object) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def _import_duckdb():
    try:
        import duckdb
    except ImportError as exc:  # pragma: no cover - environment contract
        raise RuntimeError(
            "Variant summary aggregation requires the python-duckdb package"
        ) from exc
    return duckdb
