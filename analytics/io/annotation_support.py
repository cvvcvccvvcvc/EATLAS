"""Derive report support tables from durable event-level evidence."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import sys
import tempfile
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path

from analytics.io.artifacts import content_identity, file_identity, write_json_atomic
from genomics.variants import parse_variant_key


# These standalone pipeline helpers still use sibling imports. Keep the compatibility
# bridge at this migration boundary so the scientific algorithms remain single-source.
_BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
_ADDED_BIN_PATH = str(_BIN_DIR) not in sys.path
if _ADDED_BIN_PATH:
    sys.path.insert(0, str(_BIN_DIR))
try:
    from annotate_events import (
        EVENT_VARIANT_MAP_FIELDS,
        VARIANT_STRATEGY_SUPPORT_FIELDS,
        EventOrthologSupportStream,
        add_strategy_support,
        build_gnomad_statuses,
        build_variant_strategy_support,
        load_snv_alt_genus_support,
        variant_aggregate_key,
    )
    from feature_coverage import (
        iter_snv_event_sites,
        load_snv_site_depth,
        write_snv_site_depth,
        write_snv_taxonomic_depth,
    )
    from finalize_annotation_partitions import (
        ORTHOLOG_EVIDENCE_FIELDS,
        merge_ortholog_evidence,
    )
    from ortholog_evidence_summary import write_ortholog_evidence_summary
    from taxonomic_evidence import COUNT_KEYS, count_member_groups, load_taxonomy_profiles
finally:
    if _ADDED_BIN_PATH:
        sys.path.remove(str(_BIN_DIR))


CACHE_SCHEMA_VERSION = 2
CACHE_DIRNAME = "annotation_support"
VARIANT_SUPPORT_FILENAME = "variant_strategy_support.tsv.gz"
ORTHOLOG_EVIDENCE_FILENAME = "ortholog_evidence_summary.tsv.gz"
EVENTS_FILENAME = "alignment_events.tsv.gz"
SEGMENTS_FILENAME = "alignment_segments.tsv.gz"
EVENT_SUPPORT_FILENAME = "event_ortholog_support.tsv.gz"
EVENT_MAP_FILENAME = "event_variant_map.tsv.gz"
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
    if annotation_manifest.get("schema") != "normalized_annotation_evidence_v1":
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
    source_annotations = annotation_dir / "variant_annotations.tsv.gz"
    failures = annotation_dir / "failures.tsv.gz"
    taxonomy = run_dir / "fetch" / "taxonomy.tsv.gz"
    target_features = run_dir / "fetch" / "target_features.tsv.gz"
    required_paths = (
        alignment_manifest_path,
        map_root,
        evidence_root,
        source_annotations,
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
    if alignment_manifest.get("schema") != "normalized_alignment_evidence_v1":
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
        source_annotations=source_annotations,
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
    source_annotations: Path,
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
        source_annotations,
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
    inputs = _input_identities(
        partition_dirs=partition_dirs,
        map_root=map_root,
        taxonomy=taxonomy,
        target_features=target_features,
        source_annotations=source_annotations,
        failures=failures,
        alignment_manifest=alignment_manifest,
        annotation_manifest=annotation_manifest,
    )
    fingerprint = _fingerprint(inputs)
    if _cache_is_valid(manifest_path, outputs, inputs, fingerprint):
        return outputs

    try:
        import duckdb
    except ImportError as exc:
        raise RuntimeError(
            "DuckDB is required to derive annotation support; "
            "run analytics in envs/analytics.yml"
        ) from exc

    cache_dir.mkdir(parents=True, exist_ok=True)
    profiles = load_taxonomy_profiles(taxonomy)
    failure_rows = list(_iter_required_tsv(failures, {"source", "scope", "chrom", "start", "end"}))
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
            connection.execute(f"SET temp_directory = {_sql_string(temporary_dir / 'duckdb_tmp')}")
            _load_source_annotations(connection, source_annotations)
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
                    aggregates = _collapse_partition(
                        partition_id,
                        partition_dir / EVENTS_FILENAME,
                        partition_dir / EVENT_SUPPORT_FILENAME,
                        map_path,
                        profiles,
                        alt_support_path,
                    )
                    site_depth_path = partition_work / "snv_site_depth.tsv.gz"
                    taxonomic_depth_path = partition_work / "snv_taxonomic_depth.tsv.gz"
                    write_snv_site_depth(
                        [partition_dir / SEGMENTS_FILENAME],
                        iter_snv_event_sites([partition_dir / EVENTS_FILENAME]),
                        site_depth_path,
                        partition_work,
                    )
                    write_snv_taxonomic_depth(
                        [partition_dir / SEGMENTS_FILENAME],
                        iter_snv_event_sites([partition_dir / EVENTS_FILENAME]),
                        taxonomy,
                        taxonomic_depth_path,
                        partition_work,
                    )
                    support_rows, missing_key_count = build_variant_strategy_support(
                        aggregates.values(),
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

                    annotation_rows = _source_annotations_for_genes(
                        connection,
                        {str(item.get("gene_id") or "") for item in aggregates.values()},
                    )
                    gnomad_cache: dict[tuple[str, int, str, str], dict[str, object]] = {}
                    for aggregate in aggregates.values():
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
                        aggregates.values(),
                        gnomad_cache,
                        failure_rows,
                    )
                    ortholog_path = partition_work / ORTHOLOG_EVIDENCE_FILENAME
                    ortholog_count = write_ortholog_evidence_summary(
                        taxonomic_depth_path,
                        alt_support_path,
                        [target_features],
                        statuses,
                        ortholog_path,
                    )
                    ortholog_partition_inputs.append(
                        (
                            partition_work,
                            {"ortholog_evidence_summary_count": ortholog_count},
                        )
                    )
        finally:
            connection.close()

        ortholog_row_count = merge_ortholog_evidence(
            ortholog_partition_inputs,
            temporary_ortholog,
        )
        temporary_support.chmod(0o644)
        temporary_ortholog.chmod(0o644)
        temporary_support.replace(outputs.variant_strategy_support_tsv)
        temporary_ortholog.replace(outputs.ortholog_evidence_summary_tsv)

    write_json_atomic(
        manifest_path,
        {
            "schema_version": CACHE_SCHEMA_VERSION,
            "status": "complete",
            "fingerprint": fingerprint,
            "inputs": inputs,
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
) -> dict[tuple, dict]:
    """Join partition-local lineage streams while holding one exact-support group."""

    aggregates: dict[tuple, dict] = {}
    event_required = {
        "event_group_id",
        "gene_id",
        "strategy",
        "event_type",
        "target_start0",
        "ref",
        "alt",
        "support_row_count",
        "support_ortholog_count",
    }
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
        event_missing = event_required - set(event_reader.fieldnames or [])
        if event_missing:
            raise ValueError(
                f"Alignment events {events_path} missing columns: "
                + ", ".join(sorted(event_missing))
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
                    "support_row_count": 0,
                    "_lookup_key": lookup_key,
                    "_support_by_strategy": {},
                }
                aggregates[aggregate_key] = aggregate
            aggregate["support_row_count"] += _positive_int(
                event_row.get("support_row_count"),
                f"support_row_count in {events_path}",
            )
            add_strategy_support(aggregate, event_row)
            strategy = str(event_row.get("strategy") or "")
            aggregate["_support_by_strategy"][strategy].orthologs.update(
                str(row["ortholog_gene_id"]) for row in support_rows
            )
        support_stream.finish()
    for aggregate in aggregates.values():
        for support in aggregate["_support_by_strategy"].values():
            # add_strategy_support sees only the compact per-event hint. Exact
            # event support is the canonical source of distinct orthologs after
            # multiple raw events collapse to one normalized variant.
            support.ortholog_count_hint = len(support.orthologs)
    return aggregates


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
    row_count = 0
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
        row_count += _positive_int(
            row.get("support_row_count"),
            f"support_row_count in {path}",
        )
    expected_rows = _positive_int(
        event_row.get("support_row_count"),
        "support_row_count in compact event",
    )
    expected_orthologs = _positive_int(
        event_row.get("support_ortholog_count"),
        "support_ortholog_count in compact event",
    )
    if row_count != expected_rows or len(orthologs) != expected_orthologs:
        raise ValueError(
            "Exact event support does not match compact event totals for "
            f"event_group_id={event_row.get('event_group_id')}: "
            f"rows={row_count}/{expected_rows}, "
            f"orthologs={len(orthologs)}/{expected_orthologs}"
        )


def _load_source_annotations(connection, path: Path) -> None:
    required = {"variant_key", "gene_id", "lookup_status", "gnomad_af"}
    _require_tsv_header(path, required)
    connection.execute(
        f"""
        CREATE TABLE source_annotations AS
        SELECT
            CAST(variant_key AS VARCHAR) AS variant_key,
            CAST(gene_id AS VARCHAR) AS gene_id,
            CAST(lookup_status AS VARCHAR) AS lookup_status,
            CAST(gnomad_af AS VARCHAR) AS gnomad_af
        FROM read_csv(
            {_sql_string(path)},
            delim = '\t',
            header = true,
            all_varchar = true
        )
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
    source_annotations: Path,
    failures: Path,
    alignment_manifest: Path,
    annotation_manifest: Path,
) -> dict[str, object]:
    return {
        "alignment_manifest": content_identity(alignment_manifest),
        "annotation_manifest": content_identity(annotation_manifest),
        "taxonomy": content_identity(taxonomy),
        "target_features": content_identity(target_features),
        "source_annotations": file_identity(source_annotations),
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
