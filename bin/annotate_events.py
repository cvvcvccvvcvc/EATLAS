#!/usr/bin/env python3
"""Annotate TSV variants with ClinVar and gnomAD."""

import argparse
import csv
import json
import gzip
import logging
import os
import shutil
import sys
import time
import concurrent.futures
from collections.abc import Iterable, Iterator
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pysam

if __package__ in {None, ""}:
    runtime_root = Path.cwd()
    if not (runtime_root / "genomics").is_dir():
        runtime_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(runtime_root))

from feature_coverage import load_snv_site_depth
from genomics.clinvar import review_stars as clinvar_review_stars
from genomics.gnomad import GNOMAD_API_URL, fetch_region_variants_recursive, select_af_metrics
from genomics.gnomad_cache import GnomadRegionCache
from genomics.variants import (
    add_context_normalized_record,
    build_context_index,
    event_vcf_key,
    load_target_contexts,
    normalize_chrom,
    refseq_accession_to_chrom,
    variant_key_text,
)
from ortholog_evidence_summary import write_ortholog_evidence_summary


csv.field_size_limit(sys.maxsize)


logger = logging.getLogger(__name__)


def start_phase(name: str) -> float:
    logger.info("Starting phase %s", name)
    return time.perf_counter()


def finish_phase(timings: dict[str, float], name: str, started_at: float) -> None:
    elapsed = round(time.perf_counter() - started_at, 3)
    timings[name] = elapsed
    logger.info("Timing %s: %.3f seconds", name, elapsed)

CLINVAR_COLUMNS = [
    "clinvar_sig",
    "clinvar_revstat",
    "clinvar_review_stars",
    "clinvar_review_stars_status",
    "clinvar_id",
    "clinvar_allele_id",
    "clinvar_scv_count",
    "clinvar_hgvs",
    "clinvar_disease",
    "clinvar_variant_type",
]

GNOMAD_COLUMNS = [
    "gnomad_af",
    "gnomad_af_source",
    "gnomad_csq",
]

ANNOTATION_COLUMNS = CLINVAR_COLUMNS + GNOMAD_COLUMNS

VARIANT_ANNOTATION_FIELDS = [
    "variant_key",
    "gene_id",
    "event_type",
    "ref",
    "alt",
    "lookup_status",
    "support_row_count",
    "support_ortholog_count",
    "strategies",
    *ANNOTATION_COLUMNS,
]

VARIANT_STRATEGY_SUPPORT_FIELDS = [
    "variant_key",
    "gene_id",
    "strategy",
    "alt_support_row_count",
    "alt_support_ortholog_count",
    "alt_support_genus_count",
    "site_aligned_ortholog_count",
]

ALT_SUPPORT_GENUS_COLUMN = "all__genus"

VARIANT_ORTHOLOG_SUPPORT_FIELDS = [
    "variant_key",
    "gene_id",
    "strategy",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "mapq",
    "native_alignment_type",
    "support_row_count",
]

PARTITION_TSV_SHARD_FORMAT = "headerless_gzip_member_v1"
PARTITION_TSV_SHARD_FIELDS = {
    "variant_annotations.tsv.gz": VARIANT_ANNOTATION_FIELDS,
    "variant_strategy_support.tsv.gz": VARIANT_STRATEGY_SUPPORT_FIELDS,
}

FAILURE_FIELDS = ["source", "scope", "chrom", "start", "end", "failure_type", "message"]
GNOMAD_DATASET = "gnomad_r4"

@dataclass(slots=True)
class StrategySupport:
    row_count: int = 0
    orthologs: set[str] = field(default_factory=set)
    ortholog_count_hint: int = 0


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-manifest", required=True, type=Path)
    parser.add_argument("--events-tsv", required=True, type=Path)
    parser.add_argument("--event-ortholog-support-tsv", required=True, type=Path)
    parser.add_argument("--snv-site-depth-tsv", required=True, type=Path)
    parser.add_argument("--snv-taxonomic-depth-tsv", required=True, type=Path)
    parser.add_argument("--snv-alt-taxonomic-support-tsv", required=True, type=Path)
    parser.add_argument("--genes-tsv", required=True, type=Path)
    parser.add_argument("--target-sequences-dir", required=True, type=Path)
    parser.add_argument("--target-features", required=True, type=Path)
    parser.add_argument("--clinvar-vcf", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--partition-id", default="")
    parser.add_argument(
        "--gnomad-cache-dir",
        type=Path,
        default=os.environ.get("GAPH_GNOMAD_CACHE_DIR") or None,
        help="Optional shared directory for resumable gnomAD regional responses.",
    )
    return parser.parse_args()

def cluster_positions(positions: list[int], max_gap: int = 100000) -> list[tuple[int, int]]:
    if not positions: return []
    positions = sorted(positions)
    clusters = []
    start = positions[0]
    last = positions[0]
    for p in positions[1:]:
        if p - last > max_gap:
            clusters.append((start, last))
            start = p
        last = p
    clusters.append((start, last))
    return clusters


def resolve_target_feature_paths(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    if path.is_dir():
        paths = sorted(path.glob("*.tsv.gz"))
        if paths:
            return paths
        raise ValueError(f"No target feature tables found in {path}")
    raise FileNotFoundError(f"Target features input not found: {path}")


def manifest_count(manifest: dict, field: str) -> int:
    value = manifest.get(field)
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"Alignment manifest has invalid {field}: {value!r}")
    return value


def load_alignment_manifest(path: Path, partition_id: str) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Alignment manifest not found: {path}")
    manifest = json.loads(path.read_text())
    if manifest.get("stage") != "alignment":
        raise ValueError(f"Alignment manifest has invalid stage: {manifest.get('stage')!r}")
    if manifest.get("alignment_event_mode") != "compact_support":
        raise ValueError(
            "Annotation requires compact_support alignment events, observed "
            f"{manifest.get('alignment_event_mode')!r}"
        )
    if manifest.get("event_ortholog_support_format") != "event_group_id_v1":
        raise ValueError(
            "Alignment manifest has unsupported event ortholog support format: "
            f"{manifest.get('event_ortholog_support_format')!r}"
        )
    observed_partition_id = str(manifest.get("partition_id") or "")
    if observed_partition_id != partition_id:
        raise ValueError(
            "Alignment manifest partition mismatch: "
            f"expected {partition_id!r}, observed {observed_partition_id!r}"
        )
    expected_profile = "annotation-input" if partition_id else "full"
    if manifest.get("output_profile") != expected_profile:
        raise ValueError(
            "Alignment manifest output profile mismatch: "
            f"expected {expected_profile!r}, observed {manifest.get('output_profile')!r}"
        )
    for field in (
        "alignment_event_count",
        "event_ortholog_support_count",
        "snv_site_depth_count",
        "snv_taxonomic_depth_count",
        "snv_alt_taxonomic_support_count",
    ):
        manifest_count(manifest, field)
    return manifest


def count_tsv_gz_rows(path: Path) -> int:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) is None:
            raise ValueError(f"TSV has no header: {path}")
        return sum(1 for _row in reader)


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def write_tsv_gz(
    path: Path,
    fields: list[str],
    rows: list[dict],
    *,
    include_header: bool = True,
) -> int:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        if include_header:
            writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return len(rows)


class EventOrthologSupportStream:
    """Read the compact merge handoff one event group at a time."""

    REQUIRED_FIELDS = {
        "event_group_id",
        "ortholog_gene_id",
        "tax_id",
        "taxname",
        "mapq",
        "native_alignment_type",
        "support_row_count",
    }

    def __init__(self, path: Path):
        self.path = path
        self.handle = None
        self.reader: Iterator[dict[str, str]] | None = None
        self.current: dict[str, str] | None = None
        self.row_count = 0

    def __enter__(self) -> "EventOrthologSupportStream":
        self.handle = open_text(self.path)
        reader = csv.DictReader(self.handle, delimiter="\t")
        missing = self.REQUIRED_FIELDS - set(reader.fieldnames or [])
        if missing:
            self.handle.close()
            raise ValueError(
                "Event ortholog support table missing required columns: "
                + ", ".join(sorted(missing))
            )
        self.reader = iter(reader)
        self.current = next(self.reader, None)
        return self

    def __exit__(self, *_args) -> None:
        if self.handle is not None:
            self.handle.close()

    @staticmethod
    def group_id(row: dict[str, str]) -> int:
        try:
            group_id = int(row["event_group_id"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("event_group_id must be a positive integer") from exc
        if group_id < 1:
            raise ValueError("event_group_id must be a positive integer")
        return group_id

    def take(self, expected_group_id: int) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        if self.current is None:
            return rows
        observed_group_id = self.group_id(self.current)
        if observed_group_id < expected_group_id:
            raise ValueError(
                "Event ortholog support is out of order or has no matching event: "
                f"event_group_id={observed_group_id}"
            )
        if observed_group_id > expected_group_id:
            return rows
        while self.current is not None and self.group_id(self.current) == expected_group_id:
            rows.append(self.current)
            self.row_count += 1
            if self.reader is None:
                raise RuntimeError("Event ortholog support stream is not open")
            self.current = next(self.reader, None)
        return rows

    def finish(self) -> None:
        if self.current is not None:
            raise ValueError(
                "Event ortholog support has no matching compact event: "
                f"event_group_id={self.group_id(self.current)}"
            )


class ExactSupportSpool:
    """Encode exact supporters as narrow local integer IDs for DuckDB."""

    def __init__(self, outdir: Path):
        self.path = outdir / ".variant_ortholog_support_rows.tsv"
        self.handle = self.path.open("w", buffering=1024 * 1024)
        self.strategy_ids: dict[str, int] = {}
        self.strategy_names: list[str] = []
        self.ortholog_ids: dict[tuple[str, str], int] = {}
        self.ortholog_rows: list[dict[str, str | int]] = []
        self.used_variant_ids: set[int] = set()
        self.input_edge_count = 0
        self.missing_key_count = 0

    def close(self) -> None:
        if not self.handle.closed:
            self.handle.close()

    def strategy_id(self, strategy: str) -> int:
        strategy_id = self.strategy_ids.get(strategy)
        if strategy_id is None:
            strategy_id = len(self.strategy_names) + 1
            self.strategy_ids[strategy] = strategy_id
            self.strategy_names.append(strategy)
        return strategy_id

    def ortholog_id(self, gene_id: str, row: dict[str, str]) -> int:
        ortholog_gene_id = str(row.get("ortholog_gene_id") or "")
        if not ortholog_gene_id:
            raise ValueError("Ortholog support row requires ortholog_gene_id")
        key = (gene_id, ortholog_gene_id)
        ortholog_id = self.ortholog_ids.get(key)
        if ortholog_id is None:
            ortholog_id = len(self.ortholog_rows) + 1
            self.ortholog_ids[key] = ortholog_id
            self.ortholog_rows.append(
                {
                    "ortholog_id": ortholog_id,
                    "gene_id": gene_id,
                    "ortholog_gene_id": ortholog_gene_id,
                    "tax_id": str(row.get("tax_id") or ""),
                    "taxname": str(row.get("taxname") or ""),
                }
            )
            return ortholog_id

        current = self.ortholog_rows[ortholog_id - 1]
        for field in ("tax_id", "taxname"):
            observed = str(row.get(field) or "")
            previous = str(current.get(field) or "")
            if previous and observed and previous != observed:
                raise ValueError(
                    f"Conflicting {field} for gene_id={gene_id}, "
                    f"ortholog_gene_id={ortholog_gene_id}: {previous!r} != {observed!r}"
                )
            if not previous and observed:
                current[field] = observed
        return ortholog_id

    def add_group(
        self,
        aggregate: dict,
        event_row: dict[str, str],
        support_rows: list[dict[str, str]],
    ) -> None:
        if not support_rows:
            return
        if not aggregate.get("variant_key"):
            self.missing_key_count += len(support_rows)
            return
        strategy = str(event_row.get("strategy") or "")
        if not strategy:
            raise ValueError("Exact ortholog support requires one event strategy")
        variant_context_id = int(aggregate["_variant_context_id"])
        strategy_id = self.strategy_id(strategy)
        gene_id = str(aggregate.get("gene_id") or "")
        self.used_variant_ids.add(variant_context_id)
        for row in support_rows:
            support_row_count = int_or_default(row.get("support_row_count"), 0)
            if support_row_count < 1:
                raise ValueError("Ortholog support_row_count must be positive")
            ortholog_id = self.ortholog_id(gene_id, row)
            self.input_edge_count += 1
            self.handle.write(
                f"{variant_context_id}\t{strategy_id}\t{ortholog_id}\t{support_row_count}"
                f"\t{row.get('mapq', '')}\t{row.get('native_alignment_type', '')}"
                f"\t{self.input_edge_count}\n"
            )


def int_or_default(value, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def add_strategy_support(aggregate: dict, row: dict[str, str]) -> None:
    strategy = str(row.get("strategy") or "")
    if not strategy:
        raise ValueError("Event row requires one alignment strategy")
    support_row_count = int_or_default(row.get("support_row_count"), 1)
    ortholog_gene_id = row.get("ortholog_gene_id", "")
    ortholog_count_hint = int_or_default(row.get("support_ortholog_count"), 0)
    support_by_strategy = aggregate["_support_by_strategy"]
    support = support_by_strategy.get(strategy)
    if support is None:
        support = StrategySupport()
        support_by_strategy[strategy] = support
    support.row_count += support_row_count
    if ortholog_gene_id:
        support.orthologs.add(ortholog_gene_id)
    else:
        support.ortholog_count_hint += ortholog_count_hint


def sql_string(value: Path | str) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def write_dimension_tsv(path: Path, fields: list[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_empty_exact_support_parquet(connection, output: Path) -> None:
    connection.execute(
        f"""
        COPY (
            SELECT
                CAST(NULL AS VARCHAR) AS variant_key,
                CAST(NULL AS VARCHAR) AS gene_id,
                CAST(NULL AS VARCHAR) AS strategy,
                CAST(NULL AS VARCHAR) AS ortholog_gene_id,
                CAST(NULL AS VARCHAR) AS tax_id,
                CAST(NULL AS VARCHAR) AS taxname,
                CAST(NULL AS USMALLINT) AS mapq,
                CAST(NULL AS VARCHAR) AS native_alignment_type,
                CAST(NULL AS UBIGINT) AS support_row_count
            WHERE FALSE
        ) TO {sql_string(output)} (
            FORMAT PARQUET,
            COMPRESSION ZSTD
        )
        """
    )


def aggregate_exact_support(
    spool: ExactSupportSpool,
    aggregates_by_id: list[dict | None],
    output_dir: Path,
) -> int:
    """Aggregate local integer edges and write one partition Parquet file."""

    try:
        import duckdb
    except ImportError as exc:
        raise RuntimeError(
            "DuckDB is required to aggregate exact ortholog support; "
            "run annotation in the declared alignment environment"
        ) from exc

    spool.close()
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / "part-00000.parquet"
    output.unlink(missing_ok=True)
    temp_dir = output_dir.parent / ".exact_support_duckdb"
    temp_dir.mkdir(parents=True, exist_ok=True)
    variant_dim = output_dir.parent / ".variant_context_dimension.tsv"
    strategy_dim = output_dir.parent / ".strategy_dimension.tsv"
    ortholog_dim = output_dir.parent / ".ortholog_dimension.tsv"
    write_dimension_tsv(
        variant_dim,
        ["variant_context_id", "variant_key", "gene_id"],
        (
            {
                "variant_context_id": variant_context_id,
                "variant_key": aggregates_by_id[variant_context_id]["variant_key"],
                "gene_id": aggregates_by_id[variant_context_id]["gene_id"],
            }
            for variant_context_id in sorted(spool.used_variant_ids)
        ),
    )
    write_dimension_tsv(
        strategy_dim,
        ["strategy_id", "strategy"],
        (
            {"strategy_id": index, "strategy": strategy}
            for index, strategy in enumerate(spool.strategy_names, start=1)
        ),
    )
    write_dimension_tsv(
        ortholog_dim,
        ["ortholog_id", "gene_id", "ortholog_gene_id", "tax_id", "taxname"],
        spool.ortholog_rows,
    )

    memory_limit = os.environ.get("GAPH_ANNOTATION_DUCKDB_MEMORY_LIMIT", "4GB")
    threads = max(1, int_or_default(os.environ.get("GAPH_ANNOTATION_DUCKDB_THREADS"), 1))
    connection = duckdb.connect(
        config={
            "memory_limit": memory_limit,
            "threads": str(threads),
            "temp_directory": str(temp_dir),
            "preserve_insertion_order": "false",
        }
    )
    try:
        if spool.input_edge_count == 0:
            write_empty_exact_support_parquet(connection, output)
            return 0

        connection.execute(
            f"""
            CREATE TEMP TABLE exact_support AS
            SELECT
                variant_context_id,
                strategy_id,
                ortholog_id,
                arg_min(mapq, edge_order) AS mapq,
                arg_min(native_alignment_type, edge_order) AS native_alignment_type,
                CAST(SUM(support_row_count) AS UBIGINT) AS support_row_count
            FROM read_csv(
                {sql_string(spool.path)},
                delim = '\t',
                header = false,
                columns = {{
                    'variant_context_id': 'UBIGINT',
                    'strategy_id': 'UINTEGER',
                    'ortholog_id': 'UBIGINT',
                    'support_row_count': 'UBIGINT',
                    'mapq': 'VARCHAR',
                    'native_alignment_type': 'VARCHAR',
                    'edge_order': 'UBIGINT'
                }}
            )
            GROUP BY variant_context_id, strategy_id, ortholog_id
            """
        )

        observed: set[tuple[int, str]] = set()
        for variant_context_id, strategy_id, ortholog_count, row_count in connection.execute(
            """
            SELECT
                variant_context_id,
                strategy_id,
                COUNT(*) AS ortholog_count,
                SUM(support_row_count) AS row_count
            FROM exact_support
            GROUP BY variant_context_id, strategy_id
            """
        ).fetchall():
            aggregate = aggregates_by_id[int(variant_context_id)]
            if aggregate is None:
                raise ValueError(f"Unknown variant_context_id: {variant_context_id}")
            strategy = spool.strategy_names[int(strategy_id) - 1]
            support = aggregate["_support_by_strategy"].get(strategy)
            if support is None:
                raise ValueError(
                    "Exact ortholog support identifies a strategy absent from events: "
                    f"variant_context_id={variant_context_id}, strategy={strategy}"
                )
            expected_ortholog_count = max(
                len(support.orthologs),
                support.ortholog_count_hint,
            )
            if (
                int(row_count) != support.row_count
                or int(ortholog_count) < 1
                or int(ortholog_count) > expected_ortholog_count
            ):
                raise ValueError(
                    "Exact ortholog support does not match event totals for "
                    f"variant_context_id={variant_context_id}, strategy={strategy}: "
                    f"orthologs={ortholog_count}/at-most-{expected_ortholog_count}, "
                    f"rows={row_count}/{support.row_count}"
                )
            support.orthologs.clear()
            support.ortholog_count_hint = int(ortholog_count)
            support.row_count = int(row_count)
            observed.add((int(variant_context_id), strategy))

        for aggregate in aggregates_by_id[1:]:
            if aggregate is None or not aggregate.get("variant_key"):
                continue
            variant_context_id = int(aggregate["_variant_context_id"])
            for strategy, support in aggregate["_support_by_strategy"].items():
                if support.row_count > 0 and (variant_context_id, strategy) not in observed:
                    raise ValueError(
                        "Event support has no exact ortholog rows for "
                        f"variant_context_id={variant_context_id}, strategy={strategy}"
                    )

        for variant_context_id, ortholog_count in connection.execute(
            """
            SELECT variant_context_id, COUNT(DISTINCT ortholog_id)
            FROM exact_support
            GROUP BY variant_context_id
            """
        ).fetchall():
            aggregate = aggregates_by_id[int(variant_context_id)]
            if aggregate is not None:
                aggregate["_exact_ortholog_count"] = int(ortholog_count)

        connection.execute(
            f"""
            COPY (
                SELECT
                    v.variant_key,
                    v.gene_id,
                    s.strategy,
                    o.ortholog_gene_id,
                    o.tax_id,
                    o.taxname,
                    CAST(NULLIF(e.mapq, '') AS USMALLINT) AS mapq,
                    NULLIF(e.native_alignment_type, '') AS native_alignment_type,
                    CAST(e.support_row_count AS UBIGINT) AS support_row_count
                FROM exact_support AS e
                JOIN read_csv(
                    {sql_string(variant_dim)},
                    delim = '\t',
                    header = true,
                    columns = {{
                        'variant_context_id': 'UBIGINT',
                        'variant_key': 'VARCHAR',
                        'gene_id': 'VARCHAR'
                    }}
                ) AS v USING (variant_context_id)
                JOIN read_csv(
                    {sql_string(strategy_dim)},
                    delim = '\t',
                    header = true,
                    columns = {{'strategy_id': 'UINTEGER', 'strategy': 'VARCHAR'}}
                ) AS s USING (strategy_id)
                JOIN read_csv(
                    {sql_string(ortholog_dim)},
                    delim = '\t',
                    header = true,
                    columns = {{
                        'ortholog_id': 'UBIGINT',
                        'gene_id': 'VARCHAR',
                        'ortholog_gene_id': 'VARCHAR',
                        'tax_id': 'VARCHAR',
                        'taxname': 'VARCHAR'
                    }}
                ) AS o USING (ortholog_id)
            ) TO {sql_string(output)} (
                FORMAT PARQUET,
                COMPRESSION ZSTD,
                ROW_GROUP_SIZE 100000
            )
            """
        )
        return int(connection.execute("SELECT COUNT(*) FROM exact_support").fetchone()[0])
    finally:
        connection.close()
        for path in (spool.path, variant_dim, strategy_dim, ortholog_dim):
            path.unlink(missing_ok=True)
        shutil.rmtree(temp_dir, ignore_errors=True)


def build_variant_strategy_support(
    aggregates: Iterable[dict],
    site_depths: dict[tuple[str, str, int], int] | None = None,
    genus_supports: dict[tuple[str, str, int, str, str], int] | None = None,
) -> tuple[list[dict[str, object]], int]:
    site_depths = site_depths or {}
    rows: list[dict[str, object]] = []
    missing_key_count = 0
    for aggregate in aggregates:
        variant_key = aggregate.get("variant_key", "")
        support_by_strategy = aggregate["_support_by_strategy"]
        if not variant_key:
            missing_key_count += len(support_by_strategy)
            continue
        for strategy, support in support_by_strategy.items():
            alt_support_count = max(
                len(support.orthologs),
                support.ortholog_count_hint,
            )
            site_depth: int | str = ""
            genus_support: int | str = ""
            if aggregate.get("event_type") == "snv":
                depth_key = (
                    str(aggregate.get("gene_id") or ""),
                    strategy,
                    int(aggregate.get("target_start0") or 0),
                )
                if depth_key not in site_depths:
                    raise ValueError(f"Missing site ortholog depth for SNV {depth_key}")
                site_depth = site_depths[depth_key]
                if alt_support_count > site_depth:
                    raise ValueError(
                        "ALT-support ortholog count exceeds site-aligned ortholog count for "
                        f"{depth_key}: {alt_support_count} > {site_depth}"
                    )
                if genus_supports is not None:
                    genus_key = (
                        str(aggregate.get("gene_id") or ""),
                        strategy,
                        int(aggregate.get("target_start0") or 0),
                        str(aggregate.get("ref") or "").upper(),
                        str(aggregate.get("alt") or "").upper(),
                    )
                    if genus_key not in genus_supports:
                        raise ValueError(
                            f"Missing genus ALT-support count for SNV {genus_key}"
                        )
                    genus_support = genus_supports[genus_key]
                    if genus_support < 0 or genus_support > alt_support_count:
                        raise ValueError(
                            "Genus ALT-support count exceeds ortholog ALT-support count for "
                            f"{genus_key}: {genus_support} > {alt_support_count}"
                        )
            rows.append(
                {
                    "variant_key": variant_key,
                    "gene_id": aggregate.get("gene_id", ""),
                    "strategy": strategy,
                    "alt_support_row_count": support.row_count,
                    "alt_support_ortholog_count": alt_support_count,
                    "alt_support_genus_count": genus_support,
                    "site_aligned_ortholog_count": site_depth,
                }
            )
    rows.sort(
        key=lambda row: (
            int_or_default(row.get("gene_id"), 10**18),
            row["variant_key"],
            row["strategy"],
        )
    )
    return rows, missing_key_count


def load_snv_alt_genus_support(
    path: Path,
) -> dict[tuple[str, str, int, str, str], int]:
    """Load exact-ALT genus counts keyed like the compact SNV support table."""

    required = {
        "gene_id",
        "strategy",
        "target_start0",
        "ref",
        "alt",
        ALT_SUPPORT_GENUS_COLUMN,
    }
    counts: dict[tuple[str, str, int, str, str], int] = {}
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"SNV ALT taxonomic support {path} missing columns: "
                + ", ".join(sorted(missing))
            )
        for row in reader:
            key = (
                str(row["gene_id"]),
                str(row["strategy"]),
                int(row["target_start0"]),
                str(row["ref"]).upper(),
                str(row["alt"]).upper(),
            )
            if key in counts:
                raise ValueError(f"Duplicate SNV ALT taxonomic support row: {key}")
            value = int(row[ALT_SUPPORT_GENUS_COLUMN])
            if value < 0:
                raise ValueError(
                    f"Negative SNV ALT genus support for {key}: {value}"
                )
            counts[key] = value
    return counts


def iter_variant_strategy_snv_sites(
    aggregates: Iterable[dict],
) -> Iterable[dict[str, object]]:
    for aggregate in aggregates:
        if aggregate.get("event_type") != "snv" or not aggregate.get("variant_key"):
            continue
        for strategy in aggregate["_support_by_strategy"]:
            yield {
                "gene_id": aggregate.get("gene_id", ""),
                "strategy": strategy,
                "target_start0": aggregate.get("target_start0", ""),
            }


def variant_aggregate_key(row: dict[str, str], variant_key: str) -> tuple:
    gene_id = row.get("gene_id", "")
    if variant_key:
        return "canonical", gene_id, variant_key
    return (
        "raw",
        gene_id,
        row.get("event_type", ""),
        row.get("target_start0", ""),
        row.get("target_end0", ""),
        row.get("genomic_accession", ""),
        row.get("genomic_start1", ""),
        row.get("genomic_end1", ""),
        row.get("ref", ""),
        row.get("alt", ""),
    )


def path_metadata(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": stat.st_size,
        "mtime": int(stat.st_mtime),
    }


def failure_row(
    source: str,
    scope: str,
    chrom: str | None,
    start: int | str | None,
    end: int | str | None,
    failure_type: str,
    message: str,
) -> dict[str, object]:
    return {
        "source": source,
        "scope": scope,
        "chrom": chrom or "",
        "start": start if start is not None else "",
        "end": end if end is not None else "",
        "failure_type": failure_type,
        "message": message,
    }


def fetch_gnomad_for_cluster(
    region_cache: GnomadRegionCache,
    chrom: str,
    start: int,
    end: int,
) -> list[dict]:
    # pad by 100 bases
    return region_cache.fetch_region(chrom, max(1, start - 100), end + 100)


def format_info_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return "|".join(str(item) for item in value)
    return str(value)


def format_review_status_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ",".join(str(item) for item in value)
    return str(value)


def empty_annotation(columns: list[str]) -> dict[str, str]:
    return {column: "" for column in columns}


def count_pipe_values(value: str) -> str:
    if not value:
        return "0"
    return str(len([item for item in value.split("|") if item]))


def format_float(value) -> str:
    return f"{value:.6g}" if value is not None else ""


def clinvar_annotation_from_record(rec) -> dict[str, str]:
    review_status = format_review_status_value(rec.info.get("CLNREVSTAT"))
    stars, stars_status = clinvar_review_stars(review_status)
    scv_accessions = format_info_value(rec.info.get("CLNSIGSCV"))
    return {
        "clinvar_sig": format_info_value(rec.info.get("CLNSIG")),
        "clinvar_revstat": review_status,
        "clinvar_review_stars": stars,
        "clinvar_review_stars_status": stars_status,
        "clinvar_id": rec.id or "",
        "clinvar_allele_id": format_info_value(rec.info.get("ALLELEID")),
        "clinvar_scv_count": count_pipe_values(scv_accessions),
        "clinvar_hgvs": format_info_value(rec.info.get("CLNHGVS")),
        "clinvar_disease": format_info_value(rec.info.get("CLNDN")),
        "clinvar_variant_type": format_info_value(rec.info.get("CLNVC")),
    }


def gnomad_annotation_from_variant(variant: dict) -> dict[str, str]:
    af, af_source, *_ = select_af_metrics(variant)
    return {
        "gnomad_af": format_float(af),
        "gnomad_af_source": af_source or "",
        "gnomad_csq": str(variant.get("consequence") or ""),
    }


def build_gnomad_statuses(
    aggregates: Iterable[dict],
    gnomad_cache: dict[tuple[str, int, str, str], dict],
    failures: Iterable[dict],
) -> dict[tuple[str, int, str, str], str]:
    failed_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for failure in failures:
        if failure.get("source") != "gnomad" or failure.get("scope") != "region":
            continue
        chrom = normalize_chrom(str(failure.get("chrom") or ""))
        try:
            start = int(failure.get("start") or 0)
            end = int(failure.get("end") or 0)
        except (TypeError, ValueError):
            continue
        if chrom and start > 0 and end >= start:
            failed_by_chrom[chrom].append((start, end))
    for chrom in failed_by_chrom:
        failed_by_chrom[chrom].sort()

    statuses: dict[tuple[str, int, str, str], str] = {}
    for aggregate in aggregates:
        if aggregate.get("event_type") != "snv":
            continue
        target_key = (
            str(aggregate.get("gene_id") or ""),
            int(aggregate.get("target_start0") or 0),
            str(aggregate.get("ref") or "").upper(),
            str(aggregate.get("alt") or "").upper(),
        )
        lookup_key = aggregate.get("_lookup_key")
        found = False
        if lookup_key in gnomad_cache:
            found = bool(gnomad_annotation_from_variant(gnomad_cache[lookup_key])["gnomad_af"])
        if found:
            status = "found"
        elif aggregate.get("lookup_status") != "ok" or lookup_key is None:
            status = "lookup_failed"
        else:
            chrom, position, _ref, _alt = lookup_key
            status = "not_found"
            for start, end in failed_by_chrom.get(normalize_chrom(chrom) or "", []):
                if start <= position <= end:
                    status = "lookup_failed"
                    break
                if start > position:
                    break
        previous = statuses.setdefault(target_key, status)
        if previous != status:
            raise ValueError(f"Conflicting gnomAD status for target SNV {target_key}")
    return statuses


def build_clinvar_cache(
    clinvar,
    accession_positions: dict[str, set[int]],
    contexts: dict[str, dict],
    context_index: dict[str, tuple[list[dict], list[int]]],
    failures: list[dict],
) -> tuple[dict[tuple[str, int, str, str], dict[str, str]], Counter]:
    cache = {}
    status_counts = Counter()
    for acc, positions in accession_positions.items():
        chrom = refseq_accession_to_chrom(acc)
        if not chrom:
            failures.append(
                failure_row(
                    "clinvar",
                    "accession",
                    "",
                    "",
                    "",
                    "unknown_chrom",
                    f"Could not map genomic accession to chromosome: {acc}",
                )
            )
            continue
        for start, end in cluster_positions(list(positions), max_gap=200000):
            try:
                for rec in clinvar.fetch(chrom, max(0, start - 1), end + 1):
                    rec_alts = rec.alts or ()
                    annotation = clinvar_annotation_from_record(rec)
                    for alt in rec_alts:
                        add_context_normalized_record(
                            cache,
                            (chrom, rec.pos, rec.ref, alt),
                            annotation,
                            contexts,
                            context_index,
                            status_counts,
                        )
            except ValueError as exc:
                message = f"ClinVar contig not found for chrom={chrom}; skipping region {start}-{end}."
                logger.warning(message)
                failures.append(
                    failure_row("clinvar", "region", chrom, start, end, "contig_not_found", str(exc))
                )
    if status_counts:
        logger.info(f"ClinVar key normalization status: {dict(status_counts)}")
    return cache, status_counts


def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    alignment_manifest = load_alignment_manifest(args.alignment_manifest, args.partition_id)
    timings_seconds: dict[str, float] = {}
    args.outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = args.outdir / "variant_annotations.tsv.gz"
    support_tsv = args.outdir / "variant_strategy_support.tsv.gz"
    ortholog_support_dir = args.outdir / "variant_ortholog_support"
    ortholog_evidence_tsv = args.outdir / "ortholog_evidence_summary.tsv.gz"
    failures_tsv = args.outdir / "failures.tsv.gz"
    manifest_json = args.outdir / "manifest.json"
    if not args.snv_site_depth_tsv.exists():
        raise FileNotFoundError(f"SNV site-depth TSV not found: {args.snv_site_depth_tsv}")
    if not args.event_ortholog_support_tsv.exists():
        raise FileNotFoundError(
            f"Event ortholog support TSV not found: {args.event_ortholog_support_tsv}"
        )
    required_inputs = [
        args.snv_taxonomic_depth_tsv,
        args.snv_alt_taxonomic_support_tsv,
        args.genes_tsv,
        args.target_sequences_dir,
    ]
    for path in required_inputs:
        if not path.exists():
            raise FileNotFoundError(f"Required annotation input not found: {path}")
    target_feature_paths = resolve_target_feature_paths(args.target_features)
    observed_taxonomic_depth_count = count_tsv_gz_rows(args.snv_taxonomic_depth_tsv)
    expected_taxonomic_depth_count = manifest_count(
        alignment_manifest,
        "snv_taxonomic_depth_count",
    )
    if observed_taxonomic_depth_count != expected_taxonomic_depth_count:
        raise ValueError(
            "SNV taxonomic depth row count does not match alignment manifest: "
            f"rows={observed_taxonomic_depth_count}, manifest={expected_taxonomic_depth_count}"
        )
    observed_taxonomic_alt_support_count = count_tsv_gz_rows(
        args.snv_alt_taxonomic_support_tsv
    )
    expected_taxonomic_alt_support_count = manifest_count(
        alignment_manifest,
        "snv_alt_taxonomic_support_count",
    )
    if observed_taxonomic_alt_support_count != expected_taxonomic_alt_support_count:
        raise ValueError(
            "SNV ALT taxonomic support row count does not match alignment manifest: "
            f"rows={observed_taxonomic_alt_support_count}, "
            f"manifest={expected_taxonomic_alt_support_count}"
        )
    genus_supports = load_snv_alt_genus_support(
        args.snv_alt_taxonomic_support_tsv
    )
    if not args.clinvar_vcf.exists():
        raise FileNotFoundError(f"ClinVar VCF not found: {args.clinvar_vcf}")
    clinvar_tbi = Path(f"{args.clinvar_vcf}.tbi")
    if not clinvar_tbi.exists():
        raise FileNotFoundError(f"ClinVar VCF index not found: {clinvar_tbi}")

    failures: list[dict] = []
    phase_started = start_phase("load_target_context")
    contexts = load_target_contexts(args.genes_tsv, args.target_sequences_dir)
    logger.info("Loaded target context for %s gene(s).", len(contexts))
    context_index = build_context_index(contexts)
    finish_phase(timings_seconds, "load_target_context", phase_started)

    # 1. Read events once, collect lookup regions, and collapse repeated support rows.
    phase_started = start_phase("collapse_events")
    accession_positions = defaultdict(set)
    event_key_status_counts = Counter()
    unique_lookup_status_counts = Counter()
    variant_aggregates: dict[tuple, dict] = {}
    aggregates_by_id: list[dict | None] = [None]
    input_row_count = 0
    support_stream = EventOrthologSupportStream(args.event_ortholog_support_tsv)
    exact_spool = ExactSupportSpool(args.outdir)
    collapse_complete = False
    try:
        support_stream.__enter__()
        with gzip.open(args.events_tsv, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            header = reader.fieldnames
            required = {
                "event_group_id",
                "gene_id",
                "event_type",
                "target_start0",
                "genomic_accession",
                "genomic_start1",
                "ref",
                "alt",
                "strategy",
            }
            missing = required - set(header or [])
            if missing:
                raise ValueError(
                    f"Events table missing required columns: {', '.join(sorted(missing))}"
                )
            for row in reader:
                input_row_count += 1
                event_group_id = int_or_default(row.get("event_group_id"), 0)
                if event_group_id != input_row_count:
                    raise ValueError(
                        "Compact event_group_id values must be consecutive from 1; "
                        f"expected {input_row_count}, observed {event_group_id}"
                    )
                gene_id = str(row.get("gene_id") or "")
                if gene_id not in contexts:
                    raise ValueError(
                        f"Alignment event references gene {gene_id!r} outside the supplied target context"
                    )
                exact_support_rows = support_stream.take(event_group_id)

                acc = row["genomic_accession"]
                lookup_key, status = event_vcf_key(row, contexts)
                event_key_status_counts[status] += 1
                if status == "non_concrete_allele":
                    exact_spool.missing_key_count += len(exact_support_rows)
                    continue
                raw_pos = int_or_default(row.get("genomic_start1"), -1)
                if acc and raw_pos > 0:
                    accession_positions[acc].add(raw_pos)
                if acc and lookup_key:
                    accession_positions[acc].add(int(lookup_key[1]))

                variant_key = variant_key_text(lookup_key)
                aggregate_key = variant_aggregate_key(row, variant_key)
                aggregate = variant_aggregates.get(aggregate_key)
                if aggregate is None:
                    variant_context_id = len(aggregates_by_id)
                    aggregate = {
                        "variant_key": variant_key,
                        "gene_id": row.get("gene_id", ""),
                        "event_type": row.get("event_type", ""),
                        "target_start0": row.get("target_start0", ""),
                        "ref": row.get("ref", ""),
                        "alt": row.get("alt", ""),
                        "lookup_status": status,
                        "support_row_count": 0,
                        "_lookup_key": lookup_key,
                        "_variant_context_id": variant_context_id,
                        "_exact_ortholog_count": 0,
                        "_support_by_strategy": {},
                    }
                    variant_aggregates[aggregate_key] = aggregate
                    aggregates_by_id.append(aggregate)
                    unique_lookup_status_counts[status] += 1

                aggregate["support_row_count"] += int_or_default(
                    row.get("support_row_count"),
                    1,
                )
                add_strategy_support(aggregate, row)
                exact_spool.add_group(aggregate, row, exact_support_rows)
        support_stream.finish()
        expected_event_count = manifest_count(alignment_manifest, "alignment_event_count")
        if input_row_count != expected_event_count:
            raise ValueError(
                "Alignment event row count does not match alignment manifest: "
                f"rows={input_row_count}, manifest={expected_event_count}"
            )
        expected_support_count = manifest_count(
            alignment_manifest,
            "event_ortholog_support_count",
        )
        if support_stream.row_count != expected_support_count:
            raise ValueError(
                "Event ortholog support row count does not match alignment manifest: "
                f"rows={support_stream.row_count}, manifest={expected_support_count}"
            )
        collapse_complete = True
    finally:
        support_stream.__exit__(None, None, None)
        exact_spool.close()
        if not collapse_complete:
            exact_spool.path.unlink(missing_ok=True)
    logger.info(f"Event key normalization status: {dict(event_key_status_counts)}")
    logger.info(f"Collapsed {input_row_count} event row(s) to {len(variant_aggregates)} variant-context row(s).")
    finish_phase(timings_seconds, "collapse_events", phase_started)

    phase_started = start_phase("aggregate_ortholog_support")
    ortholog_support_count = aggregate_exact_support(
        exact_spool,
        aggregates_by_id,
        ortholog_support_dir,
    )
    ortholog_support_missing_key_count = exact_spool.missing_key_count
    finish_phase(timings_seconds, "aggregate_ortholog_support", phase_started)

    phase_started = start_phase("load_site_depth")
    site_depths = load_snv_site_depth(args.snv_site_depth_tsv)
    expected_site_depth_count = manifest_count(alignment_manifest, "snv_site_depth_count")
    if len(site_depths) != expected_site_depth_count:
        raise ValueError(
            "SNV site-depth row count does not match alignment manifest: "
            f"rows={len(site_depths)}, manifest={expected_site_depth_count}"
        )
    logger.info(f"Calculated site-aligned ortholog depth for {len(site_depths)} variant-strategy SNV(s).")
    finish_phase(timings_seconds, "load_site_depth", phase_started)

    # 2. Determine gnomAD clusters
    phase_started = start_phase("gnomad_lookup")
    gnomad_tasks = []
    for acc, positions in accession_positions.items():
        chrom = refseq_accession_to_chrom(acc)
        if not chrom:
            continue
        clusters = cluster_positions(list(positions), max_gap=200000)
        for c_start, c_end in clusters:
            gnomad_tasks.append((chrom, c_start, c_end))

    logger.info(f"Will fetch {len(gnomad_tasks)} region(s) from gnomAD API.")

    # 3. Fetch gnomAD in parallel and cache
    gnomad_region_cache = GnomadRegionCache(
        args.gnomad_cache_dir,
        fetcher=fetch_region_variants_recursive,
    )
    gnomad_cache = {}
    gnomad_key_status_counts = Counter()
    gnomad_region_success_count = 0
    gnomad_raw_variant_count = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        future_to_task = {
            executor.submit(fetch_gnomad_for_cluster, gnomad_region_cache, chrom, start, end): (
                chrom,
                start,
                end,
            )
            for chrom, start, end in gnomad_tasks
        }
        for future in concurrent.futures.as_completed(future_to_task):
            task = future_to_task[future]
            try:
                vars_list = future.result()
                gnomad_region_success_count += 1
                gnomad_raw_variant_count += len(vars_list)
                for v in vars_list:
                    key = (v["chrom"], int(v["pos"]), v["ref"], v["alt"])
                    add_context_normalized_record(
                        gnomad_cache,
                        key,
                        v,
                        contexts,
                        context_index,
                        gnomad_key_status_counts,
                    )
            except Exception as exc:
                logger.error(f"gnomAD fetch failed for {task}: {exc}")
                chrom, start, end = task
                failures.append(
                    failure_row(
                        "gnomad",
                        "region",
                        chrom,
                        start,
                        end,
                        type(exc).__name__,
                        str(exc),
                    )
                )

    logger.info(f"Cached {len(gnomad_cache)} gnomAD variants.")
    if gnomad_key_status_counts:
        logger.info(f"gnomAD key normalization status: {dict(gnomad_key_status_counts)}")
    finish_phase(timings_seconds, "gnomad_lookup", phase_started)

    # 4. Open ClinVar
    phase_started = start_phase("clinvar_lookup")
    clinvar = pysam.VariantFile(str(args.clinvar_vcf))
    clinvar_cache, clinvar_key_status_counts = build_clinvar_cache(
        clinvar,
        accession_positions,
        contexts,
        context_index,
        failures,
    )
    logger.info(f"Cached {len(clinvar_cache)} ClinVar variants.")
    finish_phase(timings_seconds, "clinvar_lookup", phase_started)

    # 5. Annotate unique variant-context rows.
    phase_started = start_phase("write_variant_annotations")
    annotation_value_counts = Counter()
    variant_rows: list[dict] = []
    for aggregate in variant_aggregates.values():
        lookup_key = aggregate["_lookup_key"]
        clinvar_annotation = empty_annotation(CLINVAR_COLUMNS)
        gnomad_annotation = empty_annotation(GNOMAD_COLUMNS)

        if lookup_key:
            clinvar_annotation = clinvar_cache.get(lookup_key, clinvar_annotation)
        if lookup_key and lookup_key in gnomad_cache:
            gnomad_annotation = gnomad_annotation_from_variant(gnomad_cache[lookup_key])

        row = {field: aggregate.get(field, "") for field in VARIANT_ANNOTATION_FIELDS}
        support_by_strategy = aggregate["_support_by_strategy"]
        row["support_ortholog_count"] = aggregate["_exact_ortholog_count"]
        row["strategies"] = ",".join(sorted(support_by_strategy))
        row.update(clinvar_annotation)
        row.update(gnomad_annotation)
        for column in ANNOTATION_COLUMNS:
            if row[column]:
                annotation_value_counts[column] += 1
        variant_rows.append(row)

    variant_rows.sort(
        key=lambda row: (
            int_or_default(row.get("gene_id"), 10**18),
            row.get("variant_key", ""),
            row.get("event_type", ""),
        )
    )
    output_row_count = write_tsv_gz(
        out_tsv,
        VARIANT_ANNOTATION_FIELDS,
        variant_rows,
        include_header=not bool(args.partition_id),
    )
    finish_phase(timings_seconds, "write_variant_annotations", phase_started)

    phase_started = start_phase("write_support_tables")
    strategy_support_rows, strategy_support_missing_key_count = build_variant_strategy_support(
        variant_aggregates.values(),
        site_depths,
        genus_supports,
    )
    strategy_support_count = write_tsv_gz(
        support_tsv,
        VARIANT_STRATEGY_SUPPORT_FIELDS,
        strategy_support_rows,
        include_header=not bool(args.partition_id),
    )
    finish_phase(timings_seconds, "write_support_tables", phase_started)

    phase_started = start_phase("write_ortholog_evidence")
    ortholog_evidence_summary_count = write_ortholog_evidence_summary(
        args.snv_taxonomic_depth_tsv,
        args.snv_alt_taxonomic_support_tsv,
        target_feature_paths,
        build_gnomad_statuses(variant_aggregates.values(), gnomad_cache, failures),
        ortholog_evidence_tsv,
    )
    finish_phase(timings_seconds, "write_ortholog_evidence", phase_started)

    phase_started = start_phase("write_failures")
    failure_count = write_tsv_gz(failures_tsv, FAILURE_FIELDS, failures)
    finish_phase(timings_seconds, "write_failures", phase_started)
    manifest = {
        "output_mode": "unique_variant_context",
        "partition_id": args.partition_id,
        "event_row_count": input_row_count,
        "excluded_non_concrete_event_count": event_key_status_counts["non_concrete_allele"],
        "variant_context_count": len(variant_aggregates),
        "annotated_variant_context_count": output_row_count,
        "variant_strategy_support_count": strategy_support_count,
        "variant_strategy_support_missing_key_count": strategy_support_missing_key_count,
        "variant_ortholog_support_count": ortholog_support_count,
        "variant_ortholog_support_missing_key_count": ortholog_support_missing_key_count,
        "variant_ortholog_support_format": "parquet_dataset",
        "variant_ortholog_support_path": "variant_ortholog_support",
        "variant_ortholog_support_file_count": 1,
        "variant_strategy_site_depth_count": len(site_depths),
        "ortholog_evidence_summary_count": ortholog_evidence_summary_count,
        "target_context_count": len(contexts),
        "clinvar_vcf": path_metadata(args.clinvar_vcf),
        "clinvar_tbi": path_metadata(clinvar_tbi),
        "clinvar_cached_variant_count": len(clinvar_cache),
        "gnomad_api_url": GNOMAD_API_URL,
        "gnomad_dataset": GNOMAD_DATASET,
        "gnomad_region_count": len(gnomad_tasks),
        "gnomad_region_success_count": gnomad_region_success_count,
        "gnomad_region_failure_count": len(gnomad_tasks) - gnomad_region_success_count,
        "gnomad_raw_variant_count": gnomad_raw_variant_count,
        "gnomad_cached_variant_count": len(gnomad_cache),
        "gnomad_shared_cache": gnomad_region_cache.snapshot(),
        "failure_count": failure_count,
        "annotation_nonempty_counts": dict(annotation_value_counts),
        "event_key_status_counts": dict(event_key_status_counts),
        "unique_lookup_status_counts": dict(unique_lookup_status_counts),
        "gnomad_key_status_counts": dict(gnomad_key_status_counts),
        "clinvar_key_status_counts": dict(clinvar_key_status_counts),
        "timings_seconds": timings_seconds,
    }
    if args.partition_id:
        manifest["partition_tsv_shard_format"] = PARTITION_TSV_SHARD_FORMAT
        manifest["partition_tsv_shard_fields"] = PARTITION_TSV_SHARD_FIELDS
    manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    logger.info(f"Saved variant annotations to {out_tsv}")
    logger.info(f"Saved variant-strategy support to {support_tsv}")
    logger.info(f"Saved variant-ortholog support to {ortholog_support_dir}")
    logger.info(f"Saved annotation failures to {failures_tsv}")
    logger.info(f"Saved annotation manifest to {manifest_json}")

if __name__ == "__main__":
    main()
