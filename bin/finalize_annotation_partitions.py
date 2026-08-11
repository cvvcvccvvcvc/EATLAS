#!/usr/bin/env python3
"""Stream annotation partitions into the canonical Stage 3 outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


csv.field_size_limit(sys.maxsize)


COUNT_FIELDS = [
    "event_row_count",
    "excluded_non_concrete_event_count",
    "variant_context_count",
    "annotated_variant_context_count",
    "variant_strategy_support_count",
    "variant_strategy_support_missing_key_count",
    "variant_ortholog_support_count",
    "variant_ortholog_support_missing_key_count",
    "variant_strategy_site_depth_count",
    "target_context_count",
    "clinvar_cached_variant_count",
    "gnomad_region_count",
    "gnomad_region_success_count",
    "gnomad_region_failure_count",
    "gnomad_raw_variant_count",
    "gnomad_cached_variant_count",
]
COUNTER_FIELDS = [
    "annotation_nonempty_counts",
    "event_key_status_counts",
    "unique_lookup_status_counts",
    "gnomad_key_status_counts",
    "clinvar_key_status_counts",
]
GNOMAD_SHARED_CACHE_COUNT_FIELDS = [
    "tile_hit_count",
    "tile_miss_count",
    "tile_write_count",
    "corrupt_tile_count",
    "fetch_batch_count",
    "split_count",
]
GNOMAD_SHARED_CACHE_IDENTITY_FIELDS = [
    "enabled",
    "directory",
    "schema_version",
    "dataset",
    "reference_genome",
    "tile_size_bp",
]
ORTHOLOG_EVIDENCE_KEY_FIELDS = [
    "strategy",
    "target_context",
    "taxonomic_scope",
    "evidence_unit",
    "site_aligned_count",
    "alt_support_count",
]
ORTHOLOG_EVIDENCE_COUNT_FIELDS = [
    "gnomad_found_count",
    "gnomad_not_found_count",
    "gnomad_lookup_failed_count",
]
ORTHOLOG_EVIDENCE_FIELDS = [
    *ORTHOLOG_EVIDENCE_KEY_FIELDS,
    *ORTHOLOG_EVIDENCE_COUNT_FIELDS,
]
VARIANT_ORTHOLOG_SUPPORT_FIELDS = [
    "variant_key",
    "gene_id",
    "strategy",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "support_row_count",
]
PARTITION_TSV_SHARD_FORMAT = "headerless_gzip_member_v1"
FINAL_TSV_FORMAT = "concatenated_gzip_members_v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partition-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def required_manifest_value(manifest: dict, field: str, partition: Path):
    if field not in manifest:
        raise ValueError(f"Annotation partition manifest missing {field}: {partition}")
    return manifest[field]


def manifest_string(manifest: dict, field: str, partition: Path) -> str:
    value = required_manifest_value(manifest, field, partition)
    if not isinstance(value, str) or not value:
        raise ValueError(f"Annotation partition has invalid {field}: {partition}")
    return value


def manifest_count(manifest: dict, field: str, partition: Path) -> int:
    value = required_manifest_value(manifest, field, partition)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"Annotation partition has invalid {field}: {partition}")
    if value < 0:
        raise ValueError(f"Annotation partition has negative {field}: {partition}")
    return value


def manifest_counter(manifest: dict, field: str, partition: Path) -> dict[str, int]:
    value = required_manifest_value(manifest, field, partition)
    if not isinstance(value, dict):
        raise ValueError(f"Annotation partition has invalid {field}: {partition}")
    counter: dict[str, int] = {}
    for key, raw_count in value.items():
        if not isinstance(key, str) or not key:
            raise ValueError(f"Annotation partition has invalid {field} key: {partition}")
        if isinstance(raw_count, bool) or not isinstance(raw_count, int):
            raise ValueError(f"Annotation partition has invalid {field}: {partition}")
        if raw_count < 0:
            raise ValueError(f"Annotation partition has negative {field}: {partition}")
        counter[key] = raw_count
    return counter


def load_partitions(root: Path) -> list[tuple[Path, dict]]:
    if not root.is_dir():
        raise NotADirectoryError(f"Annotation partition root is not a directory: {root}")
    partitions = []
    seen_ids = set()
    for path in sorted(item for item in root.iterdir() if item.is_dir()):
        manifest_path = path / "manifest.json"
        if not manifest_path.exists():
            raise FileNotFoundError(f"Annotation partition missing manifest.json: {path}")
        manifest = json.loads(manifest_path.read_text())
        partition_id = manifest_string(manifest, "partition_id", path)
        if partition_id in seen_ids:
            raise ValueError(f"Duplicate annotation partition_id: {partition_id}")
        seen_ids.add(partition_id)
        partitions.append((path, manifest))
    if not partitions:
        raise ValueError(f"No annotation partitions found in {root}")
    return partitions


def validate_partition_manifests(partitions: list[tuple[Path, dict]]) -> None:
    reference_identity = None
    for partition, manifest in partitions:
        if manifest_string(manifest, "output_mode", partition) != "unique_variant_context":
            raise ValueError(f"Annotation partition has unexpected output_mode: {partition}")
        for field in [
            *COUNT_FIELDS,
            "failure_count",
            "ortholog_evidence_summary_count",
            "variant_ortholog_support_file_count",
        ]:
            manifest_count(manifest, field, partition)
        for field in COUNTER_FIELDS:
            manifest_counter(manifest, field, partition)
        if manifest_string(
            manifest,
            "variant_ortholog_support_format",
            partition,
        ) != "parquet_dataset":
            raise ValueError(f"Annotation partition has unexpected support format: {partition}")
        if manifest_string(
            manifest,
            "variant_ortholog_support_path",
            partition,
        ) != "variant_ortholog_support":
            raise ValueError(f"Annotation partition has unexpected support path: {partition}")

        clinvar_vcf = required_manifest_value(manifest, "clinvar_vcf", partition)
        clinvar_tbi = required_manifest_value(manifest, "clinvar_tbi", partition)
        if not isinstance(clinvar_vcf, dict) or not isinstance(clinvar_tbi, dict):
            raise ValueError(f"Annotation partition has invalid ClinVar metadata: {partition}")
        for metadata in (clinvar_vcf, clinvar_tbi):
            manifest_string(metadata, "path", partition)
            manifest_count(metadata, "size_bytes", partition)
            manifest_count(metadata, "mtime", partition)
        identity = {
            "clinvar_vcf": clinvar_vcf,
            "clinvar_tbi": clinvar_tbi,
            "gnomad_api_url": manifest_string(manifest, "gnomad_api_url", partition),
            "gnomad_dataset": manifest_string(manifest, "gnomad_dataset", partition),
        }
        if reference_identity is None:
            reference_identity = identity
        elif identity != reference_identity:
            raise ValueError("Annotation reference metadata differs across partitions")


def merge_tsv_gz(partitions: list[tuple[Path, dict]], filename: str, output: Path) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    expected_header = None
    count = 0
    with gzip.open(output, "wt", newline="") as out_handle:
        writer = None
        for partition, _manifest in partitions:
            path = partition / filename
            if not path.exists():
                raise FileNotFoundError(f"Annotation partition missing {filename}: {partition}")
            with gzip.open(path, "rt", newline="") as in_handle:
                reader = csv.reader(in_handle, delimiter="\t")
                header = next(reader, None)
                if header is None:
                    raise ValueError(f"Annotation partition has no TSV header: {path}")
                if expected_header is None:
                    expected_header = header
                    writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(header)
                elif header != expected_header:
                    raise ValueError(f"Header mismatch in {path}: expected {expected_header}, observed {header}")
                for row in reader:
                    writer.writerow(row)
                    count += 1
    return count


def concatenate_tsv_gz_members(
    partitions: list[tuple[Path, dict]],
    filename: str,
    count_field: str,
    output: Path,
) -> int:
    """Assemble headerless partition gzip members without recompressing rows."""

    expected_fields: list[str] | None = None
    source_paths: list[Path] = []
    row_count = 0
    for partition, manifest in partitions:
        if manifest.get("partition_tsv_shard_format") != PARTITION_TSV_SHARD_FORMAT:
            raise ValueError(
                f"Annotation partition does not declare {PARTITION_TSV_SHARD_FORMAT}: "
                f"{partition}"
            )
        fields_by_file = manifest.get("partition_tsv_shard_fields")
        if not isinstance(fields_by_file, dict):
            raise ValueError(
                f"Annotation partition is missing partition_tsv_shard_fields: {partition}"
            )
        fields = fields_by_file.get(filename)
        if (
            not isinstance(fields, list)
            or not fields
            or not all(isinstance(field, str) and field for field in fields)
            or len(set(fields)) != len(fields)
        ):
            raise ValueError(
                f"Annotation partition has invalid fields for {filename}: {partition}"
            )
        if expected_fields is None:
            expected_fields = fields
        elif fields != expected_fields:
            raise ValueError(
                f"Annotation partition fields differ for {filename}: "
                f"expected {expected_fields}, observed {fields} in {partition}"
            )

        path = partition / filename
        if not path.exists():
            raise FileNotFoundError(f"Annotation partition missing {filename}: {partition}")
        try:
            with gzip.open(path, "rt", newline="") as handle:
                first_row = next(csv.reader(handle, delimiter="\t"), None)
        except (EOFError, OSError, UnicodeError) as exc:
            raise ValueError(f"Invalid gzip TSV shard: {path}") from exc
        if first_row == fields:
            raise ValueError(f"Partition TSV shard unexpectedly contains a header: {path}")
        if first_row is not None and len(first_row) != len(fields):
            raise ValueError(
                f"Partition TSV shard has {len(first_row)} columns, expected "
                f"{len(fields)}: {path}"
            )
        source_paths.append(path)
        try:
            partition_count = int(manifest[count_field])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"Annotation partition has invalid {count_field}: {partition}"
            ) from exc
        if partition_count < 0:
            raise ValueError(
                f"Annotation partition has negative {count_field}: {partition_count}"
            )
        if (partition_count == 0) != (first_row is None):
            raise ValueError(
                f"Annotation partition {count_field} does not match whether {path} is empty"
            )
        row_count += partition_count

    if expected_fields is None:
        raise ValueError(f"No annotation partition fields found for {filename}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerow(expected_fields)
    with output.open("ab") as output_handle:
        for source in source_paths:
            with source.open("rb") as source_handle:
                shutil.copyfileobj(source_handle, output_handle, length=16 * 1024 * 1024)
    return row_count


def sql_string(value: Path | str) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def safe_partition_name(value: str) -> str:
    normalized = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return normalized or "partition"


def merge_ortholog_support_dataset(
    partitions: list[tuple[Path, dict]],
    output: Path,
) -> tuple[int, int]:
    try:
        import duckdb
    except ImportError as exc:
        raise RuntimeError(
            "DuckDB is required to validate exact ortholog support Parquet"
        ) from exc

    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError(f"Variant ortholog support output is not empty: {output}")
    connection = duckdb.connect()
    expected_schema = None
    row_count = 0
    file_count = 0
    try:
        for partition_index, (partition, manifest) in enumerate(partitions, start=1):
            if manifest_string(
                manifest,
                "variant_ortholog_support_format",
                partition,
            ) != "parquet_dataset":
                raise ValueError(
                    f"Annotation partition does not declare Parquet exact support: {partition}"
                )
            source_dir = partition / "variant_ortholog_support"
            source_files = sorted(source_dir.glob("*.parquet"))
            expected_file_count = manifest_count(
                manifest,
                "variant_ortholog_support_file_count",
                partition,
            )
            if len(source_files) != expected_file_count:
                raise ValueError(
                    "Variant ortholog support file count does not match partition manifest: "
                    f"partition={partition}, files={len(source_files)}, "
                    f"manifest={expected_file_count}"
                )
            if not source_files:
                raise ValueError(f"Annotation partition has no exact-support Parquet: {partition}")
            partition_id = safe_partition_name(
                manifest_string(manifest, "partition_id", partition)
            )
            for source_index, source in enumerate(source_files, start=1):
                schema = connection.execute(
                    f"DESCRIBE SELECT * FROM read_parquet({sql_string(source)})"
                ).fetchall()
                field_names = [str(row[0]) for row in schema]
                if field_names != VARIANT_ORTHOLOG_SUPPORT_FIELDS:
                    raise ValueError(
                        f"Unexpected exact-support Parquet fields in {source}: {field_names}"
                    )
                normalized_schema = [(str(row[0]), str(row[1])) for row in schema]
                if expected_schema is None:
                    expected_schema = normalized_schema
                elif normalized_schema != expected_schema:
                    raise ValueError(
                        f"Exact-support Parquet schema differs in {source}: "
                        f"expected {expected_schema}, observed {normalized_schema}"
                    )
                source_count = int(
                    connection.execute(
                        f"SELECT COUNT(*) FROM read_parquet({sql_string(source)})"
                    ).fetchone()[0]
                )
                target = output / (
                    f"part-{partition_index:06d}-{partition_id}-{source_index:02d}.parquet"
                )
                shutil.copy2(source, target)
                row_count += source_count
                file_count += 1
    finally:
        connection.close()
    return row_count, file_count


def merge_gnomad_shared_cache(partitions: list[tuple[Path, dict]]) -> dict[str, object] | None:
    snapshots = [manifest.get("gnomad_shared_cache") for _path, manifest in partitions]
    if not any(snapshot is not None for snapshot in snapshots):
        return None
    if not all(isinstance(snapshot, dict) for snapshot in snapshots):
        raise ValueError("gnomAD shared-cache metadata is missing from some annotation partitions")

    first = snapshots[0]
    identity = {
        field: required_manifest_value(first, field, partitions[0][0])
        for field in GNOMAD_SHARED_CACHE_IDENTITY_FIELDS
    }
    for (partition, _manifest), snapshot in zip(partitions[1:], snapshots[1:]):
        observed = {
            field: required_manifest_value(snapshot, field, partition)
            for field in GNOMAD_SHARED_CACHE_IDENTITY_FIELDS
        }
        if observed != identity:
            raise ValueError("gnomAD shared-cache configuration differs across annotation partitions")
    return {
        **identity,
        **{
            field: sum(
                manifest_count(snapshot, field, partition)
                for (partition, _manifest), snapshot in zip(partitions, snapshots)
            )
            for field in GNOMAD_SHARED_CACHE_COUNT_FIELDS
        },
    }


def merge_partition_timings(
    partitions: list[tuple[Path, dict]],
) -> tuple[dict[str, dict[str, float]], dict[str, float]]:
    by_partition = {}
    totals: Counter = Counter()
    for path, manifest in partitions:
        timings = manifest.get("timings_seconds")
        if not isinstance(timings, dict):
            continue
        partition_id = manifest_string(manifest, "partition_id", path)
        normalized = {str(name): float(value) for name, value in timings.items()}
        if any(value < 0 for value in normalized.values()):
            raise ValueError(f"Negative phase timing in annotation partition {partition_id}")
        by_partition[partition_id] = normalized
        totals.update(normalized)
    return (
        dict(sorted(by_partition.items())),
        {name: round(value, 3) for name, value in sorted(totals.items())},
    )


def merge_ortholog_evidence(
    partitions: list[tuple[Path, dict]],
    output: Path,
) -> int:
    totals: dict[tuple[str, ...], Counter] = {}
    for partition, manifest in partitions:
        path = partition / "ortholog_evidence_summary.tsv.gz"
        if not path.exists():
            raise FileNotFoundError(
                f"Annotation partition missing ortholog_evidence_summary.tsv.gz: {partition}"
            )
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(ORTHOLOG_EVIDENCE_FIELDS) - set(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"Ortholog evidence summary {path} missing columns: "
                    f"{', '.join(sorted(missing))}"
                )
            partition_row_count = 0
            for row in reader:
                partition_row_count += 1
                key = tuple(row[field] for field in ORTHOLOG_EVIDENCE_KEY_FIELDS)
                counter = totals.setdefault(key, Counter())
                counter.update(
                    {
                        field: int(row[field])
                        for field in ORTHOLOG_EVIDENCE_COUNT_FIELDS
                    }
                )
        expected_count = manifest_count(
            manifest,
            "ortholog_evidence_summary_count",
            partition,
        )
        if partition_row_count != expected_count:
            raise ValueError(
                "Ortholog evidence row count does not match partition manifest: "
                f"partition={partition}, rows={partition_row_count}, "
                f"manifest={expected_count}"
            )

    with gzip.open(output, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=ORTHOLOG_EVIDENCE_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for key in sorted(
            totals,
            key=lambda item: (
                item[0],
                item[1],
                item[2],
                item[3],
                int(item[4]),
                int(item[5]),
            ),
        ):
            writer.writerow(
                {
                    **dict(zip(ORTHOLOG_EVIDENCE_KEY_FIELDS, key)),
                    **totals[key],
                }
            )
    return len(totals)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    partitions = load_partitions(args.partition_root)
    validate_partition_manifests(partitions)
    annotation_count = concatenate_tsv_gz_members(
        partitions,
        "variant_annotations.tsv.gz",
        "annotated_variant_context_count",
        args.outdir / "variant_annotations.tsv.gz",
    )
    support_count = concatenate_tsv_gz_members(
        partitions,
        "variant_strategy_support.tsv.gz",
        "variant_strategy_support_count",
        args.outdir / "variant_strategy_support.tsv.gz",
    )
    ortholog_support_count, ortholog_support_file_count = merge_ortholog_support_dataset(
        partitions,
        args.outdir / "variant_ortholog_support",
    )
    ortholog_evidence_count = merge_ortholog_evidence(
        partitions,
        args.outdir / "ortholog_evidence_summary.tsv.gz",
    )
    failure_count = merge_tsv_gz(partitions, "failures.tsv.gz", args.outdir / "failures.tsv.gz")

    counts = {
        field: sum(manifest_count(manifest, field, path) for path, manifest in partitions)
        for field in COUNT_FIELDS
    }
    if counts["annotated_variant_context_count"] != annotation_count:
        raise ValueError(
            "Annotation row count does not match partition manifests: "
            f"rows={annotation_count}, manifests={counts['annotated_variant_context_count']}"
        )
    if counts["variant_strategy_support_count"] != support_count:
        raise ValueError(
            "Variant-strategy support row count does not match partition manifests: "
            f"rows={support_count}, manifests={counts['variant_strategy_support_count']}"
        )
    if counts["variant_ortholog_support_count"] != ortholog_support_count:
        raise ValueError(
            "Variant-ortholog support row count does not match partition manifests: "
            f"rows={ortholog_support_count}, manifests={counts['variant_ortholog_support_count']}"
        )
    manifest_failure_count = sum(
        manifest_count(manifest, "failure_count", path) for path, manifest in partitions
    )
    if manifest_failure_count != failure_count:
        raise ValueError(
            "Failure row count does not match partition manifests: "
            f"rows={failure_count}, manifests={manifest_failure_count}"
        )
    counters = {}
    for field in COUNTER_FIELDS:
        counter = Counter()
        for path, manifest in partitions:
            counter.update(manifest_counter(manifest, field, path))
        counters[field] = dict(counter)

    first_manifest = partitions[0][1]
    gnomad_shared_cache = merge_gnomad_shared_cache(partitions)
    partition_timings, timing_totals = merge_partition_timings(partitions)
    manifest = {
        "created_at": utc_now(),
        "output_mode": "unique_variant_context_partitioned",
        "partition_count": len(partitions),
        "partition_ids": [manifest_string(item, "partition_id", path) for path, item in partitions],
        **counts,
        **counters,
        "annotated_variant_context_count": annotation_count,
        "variant_strategy_support_count": support_count,
        "variant_ortholog_support_count": ortholog_support_count,
        "variant_ortholog_support_format": "parquet_dataset",
        "variant_ortholog_support_path": "variant_ortholog_support",
        "variant_ortholog_support_file_count": ortholog_support_file_count,
        "large_tsv_format": FINAL_TSV_FORMAT,
        "ortholog_evidence_summary_count": ortholog_evidence_count,
        "failure_count": failure_count,
        "clinvar_vcf": first_manifest["clinvar_vcf"],
        "clinvar_tbi": first_manifest["clinvar_tbi"],
        "gnomad_api_url": first_manifest["gnomad_api_url"],
        "gnomad_dataset": first_manifest["gnomad_dataset"],
        "cache_count_semantics": "ClinVar and gnomAD cache counts are sums across partitions.",
    }
    if gnomad_shared_cache is not None:
        manifest["gnomad_shared_cache"] = gnomad_shared_cache
    if partition_timings:
        manifest["partition_timings_seconds"] = partition_timings
        manifest["partition_timing_totals_seconds"] = timing_totals
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
