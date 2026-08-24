"""Derive report support tables from durable event-level evidence."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import tempfile
import time
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path

import duckdb

from analytics.derivations.feature_coverage import (
    iter_snv_event_sites,
    load_snv_site_depth,
    write_snv_site_depth,
    write_snv_taxonomic_depth,
)
from analytics.derivations.ortholog_evidence import (
    ORTHOLOG_EVIDENCE_FIELDS,
    write_ortholog_evidence_summary,
)
from analytics.derivations.support import (
    VARIANT_STRATEGY_SUPPORT_FIELDS,
    EventOrthologSupportStream,
    ExactSupportSpool,
    aggregate_exact_support,
    build_gnomad_statuses,
    build_variant_strategy_support,
    load_snv_alt_genus_support,
    merge_ortholog_evidence,
)
from analytics.derivations.taxonomy import (
    COUNT_KEYS,
    count_member_groups,
    load_taxonomy_profiles,
)
from analytics.io.artifacts import content_identity, file_identity, write_json_atomic
from analytics.io.variant_source import (
    VariantTableSource,
    resolve_variant_table_source,
    variant_source_sql,
)
from genomics.variants import parse_variant_key, variant_aggregate_key


CACHE_SCHEMA_VERSION = 5
CACHE_DIRNAME = "annotation_support"
VARIANT_SUPPORT_FILENAME = "variant_strategy_support.tsv.gz"
ORTHOLOG_EVIDENCE_FILENAME = "ortholog_evidence_summary.tsv.gz"
EVENTS_FILENAME = "alignment_events.tsv.gz"
SEGMENTS_FILENAME = "alignment_segments.tsv.gz"
EVENT_SUPPORT_FILENAME = "event_ortholog_support.tsv.gz"
EVENT_MAP_FILENAME = "event_variant_map.tsv.gz"
EVENT_VARIANT_MAP_FIELDS = [
    "event_group_id",
    "variant_key",
    "normalization_status",
]
COMPACT_EVENT_FIELDS = [
    "event_group_id",
    "gene_id",
    "event_type",
    "target_start0",
    "target_end0",
    "genomic_accession",
    "genomic_start1",
    "genomic_end1",
    "ref",
    "alt",
    "strategy",
    "qc_flags",
]
DNA_BASES = frozenset("ACGT")
SNV_ALT_TAXONOMIC_SUPPORT_FIELDS = [
    "gene_id",
    "strategy",
    "target_start0",
    "ref",
    "alt",
    *COUNT_KEYS,
]


@dataclass(frozen=True)
class AnnotationSupportPaths:
    variant_strategy_support_tsv: Path
    ortholog_evidence_summary_tsv: Path


def resolve_annotation_support_paths(run_dir: Path) -> AnnotationSupportPaths:
    """Require the normalized lineage contract and expose analytics-owned products."""

    annotation_dir = run_dir / "annotation"
    annotation_manifest_path = annotation_dir / "manifest.json"
    if not annotation_manifest_path.is_file():
        raise FileNotFoundError(
            f"Missing annotation manifest required for analytics: {annotation_manifest_path}"
        )
    annotation_manifest = _read_json(annotation_manifest_path)
    if annotation_manifest.get("stage") != "annotation":
        raise ValueError(
            f"Annotation manifest has invalid stage: {annotation_manifest_path}"
        )
    if annotation_manifest.get("schema") != "normalized_annotation_evidence_v3":
        raise ValueError(
            f"Annotation manifest has unsupported schema: {annotation_manifest_path}"
        )
    map_contract = annotation_manifest.get("event_variant_map")
    if not isinstance(map_contract, dict):
        raise ValueError(
            "Annotation manifest does not declare the required partitioned "
            f"event_variant_map contract: {annotation_manifest_path}"
        )
    _validate_map_contract(map_contract, annotation_manifest_path)

    alignment_dir = run_dir / "alignment"
    alignment_manifest_path = alignment_dir / "manifest.json"
    map_root = annotation_dir / "event_variant_map" / "partitions"
    evidence_root = alignment_dir / "evidence" / "partitions"
    variant_annotations_source = annotation_dir / "variant_annotations" / "manifest.json"
    failures = annotation_dir / "failures.tsv.gz"
    taxonomy = run_dir / "fetch" / "taxonomy.tsv.gz"
    target_features = run_dir / "fetch" / "target_features.tsv.gz"
    required_paths = (
        alignment_manifest_path,
        map_root,
        evidence_root,
        variant_annotations_source,
        failures,
        taxonomy,
        target_features,
    )
    missing = [str(path) for path in required_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Incomplete analytics annotation-support contract; missing: "
            + ", ".join(missing)
        )

    alignment_manifest = _read_json(alignment_manifest_path)
    if alignment_manifest.get("stage") != "alignment":
        raise ValueError(f"Alignment manifest has invalid stage: {alignment_manifest_path}")
    if alignment_manifest.get("schema") != "normalized_alignment_evidence_v2":
        raise ValueError(f"Alignment manifest has unsupported schema: {alignment_manifest_path}")
    normalized_evidence = alignment_manifest.get("normalized_evidence")
    if not isinstance(normalized_evidence, dict):
        raise ValueError(
            f"Alignment manifest does not declare normalized_evidence: {alignment_manifest_path}"
        )
    expected_evidence = {
        "layout": "partitioned",
        "format": "tsv_gzip_v1",
        "path": "evidence/partitions",
        "event_group_id_scope": "partition",
        "partition_files": [
            "manifest.json",
            "ortholog_alignment_summary.tsv.gz",
            "alignment_segments.tsv.gz",
            "alignment_events.tsv.gz",
            "event_ortholog_support.tsv.gz",
        ],
    }
    for field, value in expected_evidence.items():
        if normalized_evidence.get(field) != value:
            raise ValueError(
                f"Alignment manifest has invalid normalized_evidence.{field}: "
                f"{alignment_manifest_path}"
            )

    partition_dirs = _resolve_partition_dirs(
        evidence_root,
        map_root,
        annotation_manifest,
        map_contract,
    )
    if normalized_evidence.get("partition_count") != len(partition_dirs):
        raise ValueError(
            "Alignment normalized_evidence.partition_count does not match durable partitions"
        )
    return build_or_load_annotation_support(
        partition_dirs=partition_dirs,
        map_root=map_root,
        taxonomy=taxonomy,
        target_features=target_features,
        variant_annotations_source=variant_annotations_source,
        failures=failures,
        alignment_manifest=alignment_manifest_path,
        annotation_manifest=annotation_manifest_path,
        analytics_dir=run_dir / "analytics",
    )


def build_or_load_annotation_support(
    *,
    partition_dirs: list[Path],
    map_root: Path,
    taxonomy: Path,
    target_features: Path,
    variant_annotations_source: Path,
    failures: Path,
    alignment_manifest: Path,
    annotation_manifest: Path,
    analytics_dir: Path,
) -> AnnotationSupportPaths:
    """Build report schemas without using pipeline-owned report aggregates."""

    if not partition_dirs:
        raise ValueError("Annotation support requires at least one evidence partition")
    partition_dirs = sorted(partition_dirs, key=lambda path: path.name)
    partition_ids = [path.name for path in partition_dirs]
    if len(partition_ids) != len(set(partition_ids)):
        raise ValueError(f"Duplicate alignment evidence partition IDs: {partition_ids}")
    _validate_partition_inputs(partition_dirs, map_root)
    for path in (
        taxonomy,
        target_features,
        variant_annotations_source,
        failures,
        alignment_manifest,
        annotation_manifest,
    ):
        if not path.is_file():
            raise FileNotFoundError(f"Missing annotation support input: {path}")

    cache_dir = analytics_dir / CACHE_DIRNAME
    outputs = AnnotationSupportPaths(
        variant_strategy_support_tsv=cache_dir / VARIANT_SUPPORT_FILENAME,
        ortholog_evidence_summary_tsv=cache_dir / ORTHOLOG_EVIDENCE_FILENAME,
    )
    manifest_path = cache_dir / "manifest.json"
    variant_source = resolve_variant_table_source(
        variant_annotations_source,
        required_columns={"variant_key", "gene_id", "lookup_status", "gnomad_af"},
    )
    inputs = _input_identities(
        partition_dirs=partition_dirs,
        map_root=map_root,
        taxonomy=taxonomy,
        target_features=target_features,
        variant_source=variant_source,
        failures=failures,
        alignment_manifest=alignment_manifest,
        annotation_manifest=annotation_manifest,
    )
    fingerprint = _fingerprint(inputs)
    if _cache_is_valid(manifest_path, outputs, inputs, fingerprint):
        return outputs

    build_started = time.perf_counter()
    timings: dict[str, float] = {}
    support_metrics = {
        "canonical_variant_count": 0,
        "canonical_support_edge_count": 0,
        "exact_support_edge_count": 0,
    }
    cache_dir.mkdir(parents=True, exist_ok=True)
    phase_started = time.perf_counter()
    profiles = load_taxonomy_profiles(taxonomy)
    _record_timing(timings, "load_taxonomy", phase_started)
    phase_started = time.perf_counter()
    failure_rows = list(
        _iter_required_tsv(
            failures,
            {"source", "scope", "chrom", "start", "end"},
        )
    )
    _record_timing(timings, "load_failures", phase_started)
    support_row_count = 0
    ortholog_partition_inputs: list[tuple[Path, dict[str, int]]] = []
    with tempfile.TemporaryDirectory(
        prefix=".annotation_support_",
        dir=cache_dir,
    ) as temporary_name:
        temporary_dir = Path(temporary_name)
        temporary_support = temporary_dir / VARIANT_SUPPORT_FILENAME
        temporary_ortholog = temporary_dir / ORTHOLOG_EVIDENCE_FILENAME
        connection = duckdb.connect(str(temporary_dir / "source_annotations.duckdb"))
        try:
            memory_limit = os.environ.get(
                "GAPH_ANALYTICS_DUCKDB_MEMORY_LIMIT",
                "2GB",
            )
            connection.execute(f"SET memory_limit = {_sql_string(memory_limit)}")
            connection.execute(
                f"SET temp_directory = {_sql_string(temporary_dir / 'duckdb_tmp')}"
            )
            connection.execute("SET preserve_insertion_order = false")
            phase_started = time.perf_counter()
            _load_source_annotations(connection, variant_source)
            _record_timing(timings, "load_source_annotations", phase_started)
            with gzip.open(temporary_support, "wt", newline="") as support_handle:
                support_writer = csv.DictWriter(
                    support_handle,
                    fieldnames=VARIANT_STRATEGY_SUPPORT_FIELDS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                support_writer.writeheader()
                for partition_dir in partition_dirs:
                    partition_id = partition_dir.name
                    partition_work = temporary_dir / partition_id
                    partition_work.mkdir()
                    map_path = map_root / partition_id / EVENT_MAP_FILENAME
                    alt_support_path = partition_work / "snv_alt_taxonomic_support.tsv.gz"
                    exact_support_spool = ExactSupportSpool(
                        partition_work / ".exact_support_rows.tsv"
                    )
                    phase_started = time.perf_counter()
                    with exact_support_spool:
                        aggregates_by_id = _collapse_partition(
                            partition_id,
                            partition_dir / EVENTS_FILENAME,
                            partition_dir / EVENT_SUPPORT_FILENAME,
                            map_path,
                            profiles,
                            alt_support_path,
                            exact_support_spool,
                        )
                    _record_timing(timings, "collapse_events", phase_started)
                    phase_started = time.perf_counter()
                    exact_edge_count = aggregate_exact_support(
                        connection,
                        exact_support_spool,
                        aggregates_by_id,
                    )
                    _record_timing(timings, "aggregate_exact_support", phase_started)
                    aggregates = [
                        aggregate
                        for aggregate in aggregates_by_id[1:]
                        if aggregate is not None
                    ]
                    support_metrics["canonical_variant_count"] += len(aggregates)
                    input_edge_count = exact_support_spool.input_edge_count
                    support_metrics[
                        "canonical_support_edge_count"
                    ] += input_edge_count
                    support_metrics["exact_support_edge_count"] += exact_edge_count
                    site_depth_path = partition_work / "snv_site_depth.tsv.gz"
                    taxonomic_depth_path = partition_work / "snv_taxonomic_depth.tsv.gz"
                    phase_started = time.perf_counter()
                    write_snv_site_depth(
                        [partition_dir / SEGMENTS_FILENAME],
                        iter_snv_event_sites([partition_dir / EVENTS_FILENAME]),
                        site_depth_path,
                        partition_work,
                    )
                    _record_timing(timings, "snv_site_depth", phase_started)
                    phase_started = time.perf_counter()
                    write_snv_taxonomic_depth(
                        [partition_dir / SEGMENTS_FILENAME],
                        iter_snv_event_sites([partition_dir / EVENTS_FILENAME]),
                        taxonomy,
                        taxonomic_depth_path,
                        partition_work,
                    )
                    _record_timing(timings, "snv_taxonomic_depth", phase_started)
                    phase_started = time.perf_counter()
                    support_rows, missing_key_count = build_variant_strategy_support(
                        aggregates,
                        load_snv_site_depth(site_depth_path),
                        load_snv_alt_genus_support(alt_support_path),
                    )
                    if missing_key_count:
                        raise ValueError(
                            f"Canonical event map omitted {missing_key_count} support group(s) "
                            f"in partition {partition_id}"
                        )
                    support_writer.writerows(support_rows)
                    support_row_count += len(support_rows)
                    _record_timing(timings, "write_variant_support", phase_started)

                    phase_started = time.perf_counter()
                    annotation_rows = _source_annotations_for_genes(
                        connection,
                        {str(item.get("gene_id") or "") for item in aggregates},
                    )
                    gnomad_cache: dict[tuple[str, int, str, str], dict[str, object]] = {}
                    for aggregate in aggregates:
                        key = (
                            str(aggregate.get("gene_id") or ""),
                            str(aggregate.get("variant_key") or ""),
                        )
                        annotation = annotation_rows.get(key)
                        if annotation is None:
                            raise ValueError(
                                "Source variant annotations are missing canonical event support "
                                f"for partition {partition_id}: {key}"
                            )
                        if annotation["lookup_status"] != aggregate.get("lookup_status"):
                            raise ValueError(
                                f"Normalization status differs for {key}: "
                                f"map={aggregate.get('lookup_status')!r}, "
                                f"annotation={annotation['lookup_status']!r}"
                            )
                        lookup_key = aggregate.get("_lookup_key")
                        if annotation["gnomad_af"] and lookup_key is not None:
                            gnomad_cache[lookup_key] = {
                                "exome": {"af": annotation["gnomad_af"]}
                            }
                    statuses = build_gnomad_statuses(
                        aggregates,
                        gnomad_cache,
                        failure_rows,
                    )
                    _record_timing(timings, "join_source_annotations", phase_started)
                    ortholog_path = partition_work / ORTHOLOG_EVIDENCE_FILENAME
                    phase_started = time.perf_counter()
                    ortholog_count = write_ortholog_evidence_summary(
                        taxonomic_depth_path,
                        alt_support_path,
                        [target_features],
                        statuses,
                        ortholog_path,
                    )
                    _record_timing(timings, "ortholog_evidence_summary", phase_started)
                    ortholog_partition_inputs.append(
                        (
                            partition_work,
                            {"ortholog_evidence_summary_count": ortholog_count},
                        )
                    )
        finally:
            connection.close()

        phase_started = time.perf_counter()
        ortholog_row_count = merge_ortholog_evidence(
            ortholog_partition_inputs,
            temporary_ortholog,
        )
        _record_timing(timings, "merge_ortholog_evidence", phase_started)
        temporary_support.chmod(0o644)
        temporary_ortholog.chmod(0o644)
        temporary_support.replace(outputs.variant_strategy_support_tsv)
        temporary_ortholog.replace(outputs.ortholog_evidence_summary_tsv)

    timings["total"] = round(time.perf_counter() - build_started, 6)
    _report_timings(timings)
    write_json_atomic(
        manifest_path,
        {
            "schema_version": CACHE_SCHEMA_VERSION,
            "status": "complete",
            "fingerprint": fingerprint,
            "inputs": inputs,
            "timings_seconds": timings,
            "exact_support": support_metrics,
            "variant_strategy_support": {
                "columns": VARIANT_STRATEGY_SUPPORT_FIELDS,
                "row_count": support_row_count,
                "output": file_identity(outputs.variant_strategy_support_tsv),
            },
            "ortholog_evidence_summary": {
                "columns": ORTHOLOG_EVIDENCE_FIELDS,
                "row_count": ortholog_row_count,
                "output": file_identity(outputs.ortholog_evidence_summary_tsv),
            },
        },
    )
    return outputs


def _collapse_partition(
    partition_id: str,
    events_path: Path,
    event_support_path: Path,
    event_map_path: Path,
    profiles: dict,
    alt_support_path: Path,
    exact_support_spool: ExactSupportSpool,
) -> list[dict | None]:
    """Join partition-local lineage streams while holding one exact-support group."""

    aggregates: dict[tuple, dict] = {}
    aggregates_by_id: list[dict | None] = [None]
    with (
        gzip.open(events_path, "rt", newline="") as event_handle,
        gzip.open(event_map_path, "rt", newline="") as map_handle,
        EventOrthologSupportStream(event_support_path) as support_stream,
        gzip.open(alt_support_path, "wt", newline="") as alt_handle,
    ):
        alt_writer = csv.DictWriter(
            alt_handle,
            fieldnames=SNV_ALT_TAXONOMIC_SUPPORT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        alt_writer.writeheader()
        event_reader = csv.DictReader(event_handle, delimiter="\t")
        if event_reader.fieldnames != COMPACT_EVENT_FIELDS:
            raise ValueError(
                f"Alignment events {events_path} have an invalid schema: "
                f"expected {COMPACT_EVENT_FIELDS}, observed {event_reader.fieldnames}"
            )
        map_reader = csv.DictReader(map_handle, delimiter="\t")
        if map_reader.fieldnames != EVENT_VARIANT_MAP_FIELDS:
            raise ValueError(
                f"Unexpected event-variant map fields in {event_map_path}: "
                f"{map_reader.fieldnames}"
            )
        for expected_id, pair in enumerate(
            zip_longest(event_reader, map_reader),
            start=1,
        ):
            event_row, map_row = pair
            if event_row is None or map_row is None:
                raise ValueError(
                    f"Event/map row count mismatch in partition {partition_id}"
                )
            try:
                event_group_id = int(event_row["event_group_id"])
                map_group_id = int(map_row["event_group_id"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid event_group_id in partition {partition_id}"
                ) from exc
            if event_group_id != expected_id or map_group_id != expected_id:
                raise ValueError(
                    "Partition-local event_group_id values must be consecutive and aligned: "
                    f"partition={partition_id}, expected={expected_id}, "
                    f"event={event_group_id}, map={map_group_id}"
                )
            support_rows = support_stream.take(event_group_id)
            _validate_exact_event_support(event_row, support_rows, event_support_path)
            alt_row = _alt_taxonomic_support_row(event_row, support_rows, profiles)
            if alt_row is not None:
                alt_writer.writerow(alt_row)

            variant_key = str(map_row.get("variant_key") or "")
            if not variant_key:
                continue
            lookup_key = parse_variant_key(variant_key)
            if lookup_key is None:
                raise ValueError(
                    f"Invalid canonical variant_key in {event_map_path}: {variant_key!r}"
                )
            aggregate_key = variant_aggregate_key(event_row, variant_key)
            aggregate = aggregates.get(aggregate_key)
            if aggregate is None:
                aggregate = {
                    "variant_key": variant_key,
                    "gene_id": event_row.get("gene_id", ""),
                    "event_type": event_row.get("event_type", ""),
                    "target_start0": event_row.get("target_start0", ""),
                    "ref": event_row.get("ref", ""),
                    "alt": event_row.get("alt", ""),
                    "lookup_status": map_row.get("normalization_status", ""),
                    "_lookup_key": lookup_key,
                    "_support_by_strategy": {},
                    "_variant_context_id": len(aggregates_by_id),
                }
                aggregates[aggregate_key] = aggregate
                aggregates_by_id.append(aggregate)
            exact_support_spool.add_group(
                variant_context_id=int(aggregate["_variant_context_id"]),
                gene_id=str(aggregate.get("gene_id") or ""),
                strategy=str(event_row.get("strategy") or ""),
                support_rows=support_rows,
            )
        support_stream.finish()
    return aggregates_by_id


def _alt_taxonomic_support_row(
    event_row: dict[str, str],
    support_rows: list[dict[str, str]],
    profiles: dict,
) -> dict[str, object] | None:
    ref = str(event_row.get("ref") or "").upper()
    alt = str(event_row.get("alt") or "").upper()
    if (
        event_row.get("event_type") != "snv"
        or len(ref) != 1
        or len(alt) != 1
        or ref not in DNA_BASES
        or alt not in DNA_BASES
    ):
        return None
    counts = count_member_groups(
        (
            (str(row.get("ortholog_gene_id") or ""), str(row.get("tax_id") or ""))
            for row in support_rows
            if row.get("tax_id")
        ),
        profiles,
    )
    if counts["all__ortholog"] < 1:
        return None
    return {
        "gene_id": event_row["gene_id"],
        "strategy": event_row["strategy"],
        "target_start0": event_row["target_start0"],
        "ref": ref,
        "alt": alt,
        **counts,
    }


def _validate_exact_event_support(
    event_row: dict[str, str],
    support_rows: list[dict[str, str]],
    path: Path,
) -> None:
    orthologs: set[str] = set()
    for row in support_rows:
        ortholog_gene_id = str(row.get("ortholog_gene_id") or "")
        if not ortholog_gene_id:
            raise ValueError(f"Exact event support has an empty ortholog_gene_id: {path}")
        if ortholog_gene_id in orthologs:
            raise ValueError(
                "Exact event support contains duplicate ortholog rows for "
                f"event_group_id={event_row.get('event_group_id')}: {ortholog_gene_id}"
            )
        orthologs.add(ortholog_gene_id)
        _positive_int(
            row.get("support_row_count"),
            f"support_row_count in {path}",
        )
    if not orthologs:
        raise ValueError(
            "Compact event has no exact ortholog support for "
            f"event_group_id={event_row.get('event_group_id')}: {path}"
        )


def _load_source_annotations(connection, source: VariantTableSource) -> None:
    connection.execute(
        f"""
        CREATE TABLE source_annotations AS
        SELECT
            CAST(variant_key AS VARCHAR) AS variant_key,
            CAST(gene_id AS VARCHAR) AS gene_id,
            CAST(lookup_status AS VARCHAR) AS lookup_status,
            CAST(gnomad_af AS VARCHAR) AS gnomad_af
        FROM {variant_source_sql(source)}
        """
    )
    connection.execute(
        "CREATE INDEX source_annotation_key_idx ON source_annotations(gene_id, variant_key)"
    )


def _source_annotations_for_genes(
    connection,
    gene_ids: set[str],
) -> dict[tuple[str, str], dict[str, str]]:
    connection.execute("DROP TABLE IF EXISTS requested_genes")
    connection.execute("CREATE TEMP TABLE requested_genes(gene_id VARCHAR PRIMARY KEY)")
    if gene_ids:
        connection.executemany(
            "INSERT INTO requested_genes VALUES (?)",
            [(gene_id,) for gene_id in sorted(gene_ids)],
        )
    rows: dict[tuple[str, str], dict[str, str]] = {}
    for variant_key, gene_id, lookup_status, gnomad_af in connection.execute(
        """
        SELECT a.variant_key, a.gene_id, a.lookup_status, a.gnomad_af
        FROM source_annotations AS a
        JOIN requested_genes AS g USING (gene_id)
        """
    ).fetchall():
        key = (str(gene_id or ""), str(variant_key or ""))
        observed = {
            "lookup_status": str(lookup_status or ""),
            "gnomad_af": str(gnomad_af or ""),
        }
        previous = rows.setdefault(key, observed)
        if previous != observed:
            raise ValueError(
                f"Conflicting source annotation evidence for canonical context {key}"
            )
    return rows


def _resolve_partition_dirs(
    evidence_root: Path,
    map_root: Path,
    annotation_manifest: dict,
    map_contract: dict,
) -> list[Path]:
    if not evidence_root.is_dir() or not map_root.is_dir():
        raise NotADirectoryError(
            f"Evidence and event-map partition roots must be directories: "
            f"{evidence_root}, {map_root}"
        )
    evidence = sorted(path for path in evidence_root.iterdir() if path.is_dir())
    maps = sorted(path for path in map_root.iterdir() if path.is_dir())
    evidence_ids = [path.name for path in evidence]
    map_ids = [path.name for path in maps]
    declared_ids = annotation_manifest.get("partition_ids")
    if not isinstance(declared_ids, list) or not all(
        isinstance(value, str) and value for value in declared_ids
    ):
        raise ValueError("Annotation manifest has invalid partition_ids")
    if evidence_ids != sorted(declared_ids) or map_ids != sorted(declared_ids):
        raise ValueError(
            "Alignment evidence, event maps, and annotation manifest have different "
            f"partition IDs: evidence={evidence_ids}, maps={map_ids}, "
            f"annotation={declared_ids}"
        )
    if map_contract.get("partition_count") != len(evidence):
        raise ValueError("Event-variant map partition_count does not match durable partitions")
    if not evidence:
        raise ValueError(f"Annotation support contains no partitions: {evidence_root}")
    return evidence


def _validate_partition_inputs(partition_dirs: list[Path], map_root: Path) -> None:
    for partition_dir in partition_dirs:
        for filename in (EVENTS_FILENAME, SEGMENTS_FILENAME, EVENT_SUPPORT_FILENAME):
            path = partition_dir / filename
            if not path.is_file():
                raise FileNotFoundError(
                    f"Alignment evidence partition {partition_dir.name} is missing {filename}"
                )
        map_path = map_root / partition_dir.name / EVENT_MAP_FILENAME
        if not map_path.is_file():
            raise FileNotFoundError(
                f"Annotation event-map partition {partition_dir.name} is missing "
                f"{EVENT_MAP_FILENAME}"
            )


def _validate_map_contract(contract: dict, path: Path) -> None:
    expected = {
        "layout": "partitioned",
        "format": "tsv_gzip_v1",
        "path": "event_variant_map/partitions",
        "fields": EVENT_VARIANT_MAP_FIELDS,
        "event_group_id_scope": "partition",
    }
    for field, value in expected.items():
        if contract.get(field) != value:
            raise ValueError(
                f"Annotation manifest has invalid event_variant_map.{field} in {path}: "
                f"{contract.get(field)!r}"
            )


def _input_identities(
    *,
    partition_dirs: list[Path],
    map_root: Path,
    taxonomy: Path,
    target_features: Path,
    variant_source: VariantTableSource,
    failures: Path,
    alignment_manifest: Path,
    annotation_manifest: Path,
) -> dict[str, object]:
    return {
        "alignment_manifest": content_identity(alignment_manifest),
        "annotation_manifest": content_identity(annotation_manifest),
        "taxonomy": content_identity(taxonomy),
        "target_features": content_identity(target_features),
        "variant_annotations": variant_source.identity,
        "failures": file_identity(failures),
        "partitions": [
            {
                "partition_id": partition_dir.name,
                "events": file_identity(partition_dir / EVENTS_FILENAME),
                "segments": file_identity(partition_dir / SEGMENTS_FILENAME),
                "event_ortholog_support": file_identity(
                    partition_dir / EVENT_SUPPORT_FILENAME
                ),
                "event_variant_map": file_identity(
                    map_root / partition_dir.name / EVENT_MAP_FILENAME
                ),
            }
            for partition_dir in partition_dirs
        ],
    }


def _fingerprint(inputs: dict[str, object]) -> str:
    payload = {"schema_version": CACHE_SCHEMA_VERSION, "inputs": inputs}
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _cache_is_valid(
    manifest_path: Path,
    outputs: AnnotationSupportPaths,
    inputs: dict[str, object],
    fingerprint: str,
) -> bool:
    if not manifest_path.is_file() or not all(
        path.is_file()
        for path in (
            outputs.variant_strategy_support_tsv,
            outputs.ortholog_evidence_summary_tsv,
        )
    ):
        return False
    try:
        manifest = json.loads(manifest_path.read_text())
        return (
            manifest.get("schema_version") == CACHE_SCHEMA_VERSION
            and manifest.get("status") == "complete"
            and manifest.get("fingerprint") == fingerprint
            and manifest.get("inputs") == inputs
            and manifest.get("variant_strategy_support", {}).get("columns")
            == VARIANT_STRATEGY_SUPPORT_FIELDS
            and manifest.get("variant_strategy_support", {}).get("output")
            == file_identity(outputs.variant_strategy_support_tsv)
            and manifest.get("ortholog_evidence_summary", {}).get("columns")
            == ORTHOLOG_EVIDENCE_FIELDS
            and manifest.get("ortholog_evidence_summary", {}).get("output")
            == file_identity(outputs.ortholog_evidence_summary_tsv)
        )
    except (OSError, ValueError, json.JSONDecodeError):
        return False


def _iter_required_tsv(path: Path, required: set[str]):
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")
        yield from reader


def _require_tsv_header(path: Path, required: set[str]) -> None:
    with gzip.open(path, "rt", newline="") as handle:
        fields = set(next(csv.reader(handle, delimiter="\t"), []))
    missing = required - fields
    if missing:
        raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")


def _positive_int(value: object, label: str) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {label}: {value!r}") from exc
    if result < 1:
        raise ValueError(f"Invalid {label}: {value!r}")
    return result


def _record_timing(timings: dict[str, float], phase: str, started: float) -> None:
    elapsed = time.perf_counter() - started
    timings[phase] = round(timings.get(phase, 0.0) + elapsed, 6)


def _report_timings(timings: dict[str, float]) -> None:
    details = ", ".join(
        f"{phase}={seconds:.3f}s"
        for phase, seconds in timings.items()
        if phase != "total"
    )
    print(
        f"Annotation support cache built in {timings['total']:.3f}s ({details})",
        flush=True,
    )


def _read_json(path: Path) -> dict:
    try:
        payload = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def _sql_string(value: Path | str) -> str:
    return "'" + str(value).replace("'", "''") + "'"
