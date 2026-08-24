#!/usr/bin/env python3
"""Stream annotation partitions into the canonical Stage 3 outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
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
    "event_variant_map_count",
    "variant_context_count",
    "annotated_variant_context_count",
    "clinvar_cached_variant_count",
    "gnomad_region_count",
    "gnomad_region_success_count",
    "gnomad_region_failure_count",
    "gnomad_raw_variant_count",
    "gnomad_cached_variant_count",
]
COUNTER_FIELDS = [
    "event_key_status_counts",
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
EVENT_VARIANT_MAP_FIELDS = [
    "event_group_id",
    "variant_key",
    "normalization_status",
]
EVENT_VARIANT_MAP_FILENAME = "event_variant_map.tsv.gz"
VEP_SHARD_SCHEMA = "normalized_vep_annotation_shard_v1"
VARIANT_DATASET_SCHEMA = "gaph_variant_annotation_dataset_v1"
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
VEP_FIELDS = [
    "vep_status",
    "vep_primary_consequence",
    "vep_consequence_terms",
    "vep_transcript_id",
    "vep_mane_select",
    "vep_canonical",
    "vep_impact",
    "vep_variant_class",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partition-root", required=True, type=Path)
    parser.add_argument("--vep-root", required=True, type=Path)
    parser.add_argument("--clinvar-vcf", required=True, type=Path)
    parser.add_argument("--clinvar-tbi", required=True, type=Path)
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
        if SAFE_ID.fullmatch(partition_id) is None:
            raise ValueError(f"Unsafe annotation partition_id: {partition_id!r}")
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
        if manifest_string(manifest, "stage", partition) != "annotation":
            raise ValueError(f"Annotation partition has unexpected stage: {partition}")
        if (
            manifest_string(manifest, "schema", partition)
            != "normalized_annotation_evidence_partition_v2"
        ):
            raise ValueError(f"Annotation partition has unexpected schema: {partition}")
        for field in [*COUNT_FIELDS, "failure_count"]:
            manifest_count(manifest, field, partition)
        for field in COUNTER_FIELDS:
            manifest_counter(manifest, field, partition)

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


def content_identity_from_partition_input(
    path: Path,
    declared_metadata: dict[str, object],
    label: str,
) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(f"{label} input does not exist: {path}")
    before = path.stat()
    if (
        before.st_size != declared_metadata.get("size_bytes")
        or int(before.st_mtime) != declared_metadata.get("mtime")
    ):
        raise ValueError(
            f"{label} input changed between partition annotation and finalization: {path}"
        )
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"{label} input changed while hashing: {path}")
    return {"size_bytes": before.st_size, "sha256": digest.hexdigest()}


def variant_annotation_shards(
    partition: Path,
    manifest: dict,
) -> list[dict[str, object]]:
    dataset = required_manifest_value(manifest, "variant_annotations", partition)
    if not isinstance(dataset, dict):
        raise ValueError(f"Annotation partition has invalid variant_annotations: {partition}")
    if dataset.get("layout") != "partitioned" or dataset.get("format") != "tsv_gzip_v1":
        raise ValueError(
            f"Annotation partition has unsupported variant_annotations format: {partition}"
        )
    relative_root = str(dataset.get("path") or "")
    fields = dataset.get("fields")
    raw_shards = dataset.get("shards")
    if (
        not relative_root
        or not isinstance(fields, list)
        or not fields
        or not all(isinstance(field, str) and field for field in fields)
        or not isinstance(raw_shards, list)
        or not raw_shards
    ):
        raise ValueError(f"Annotation partition has invalid variant_annotations: {partition}")

    shards = []
    seen_ids = set()
    observed_rows = 0
    for raw in raw_shards:
        if not isinstance(raw, dict):
            raise ValueError(f"Annotation partition has invalid annotation shard: {partition}")
        shard_id = str(raw.get("shard_id") or "")
        relative_path = str(raw.get("path") or "")
        try:
            row_count = int(raw["row_count"])
            size_bytes = int(raw["size_bytes"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"Annotation partition has invalid annotation shard: {partition}"
            ) from exc
        if (
            not shard_id
            or SAFE_ID.fullmatch(shard_id) is None
            or shard_id in seen_ids
            or not relative_path
            or row_count < 0
            or size_bytes < 1
        ):
            raise ValueError(f"Annotation partition has invalid annotation shard: {partition}")
        seen_ids.add(shard_id)
        relative_root_path = Path(relative_root)
        relative_shard_path = Path(relative_path)
        if relative_root_path.is_absolute() or relative_shard_path.is_absolute():
            raise ValueError(f"Annotation shard path must be relative: {partition}")
        partition_root = partition.resolve()
        dataset_root = (partition_root / relative_root_path).resolve()
        source = (dataset_root / relative_shard_path).resolve()
        try:
            dataset_root.relative_to(partition_root)
            source.relative_to(dataset_root)
        except ValueError as exc:
            raise ValueError(f"Annotation shard escapes its dataset: {source}") from exc
        if not source.is_file() or source.stat().st_size != size_bytes:
            raise ValueError(f"Annotation input shard changed: {source}")
        with gzip.open(source, "rt", newline="") as handle:
            header = next(csv.reader(handle, delimiter="\t"), None)
        if header != fields:
            raise ValueError(f"Annotation input shard header changed: {source}")
        observed_rows += row_count
        shards.append(
            {
                "partition_id": manifest_string(manifest, "partition_id", partition),
                "shard_id": shard_id,
                "row_count": row_count,
                "size_bytes": size_bytes,
                "fields": fields,
                "source": source,
            }
        )
    if len(shards) != int(dataset.get("shard_count", -1)):
        raise ValueError(f"Annotation partition shard count changed: {partition}")
    if observed_rows != int(dataset.get("row_count", -1)):
        raise ValueError(f"Annotation partition shard row count changed: {partition}")
    if observed_rows != manifest_count(manifest, "annotated_variant_context_count", partition):
        raise ValueError(
            f"Annotation shard rows do not match annotated_variant_context_count: {partition}"
        )
    return shards


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


def copy_event_variant_map_dataset(
    partitions: list[tuple[Path, dict]],
    output: Path,
) -> int:
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Event-variant map output is not empty: {output}")
    partitions_root = output / "partitions"
    total_count = 0
    for partition, manifest in partitions:
        partition_id = manifest_string(manifest, "partition_id", partition)
        source = partition / EVENT_VARIANT_MAP_FILENAME
        if not source.exists():
            raise FileNotFoundError(
                f"Annotation partition missing {EVENT_VARIANT_MAP_FILENAME}: {partition}"
            )
        with gzip.open(source, "rt", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            fields = next(reader, None)
            if fields != EVENT_VARIANT_MAP_FIELDS:
                raise ValueError(
                    f"Unexpected event-variant map fields in {source}: "
                    f"{fields}"
                )
            row_count = 0
            for row_count, values in enumerate(reader, start=1):
                if len(values) != len(EVENT_VARIANT_MAP_FIELDS):
                    raise ValueError(
                        f"Event-variant map row has {len(values)} columns, expected "
                        f"{len(EVENT_VARIANT_MAP_FIELDS)}: {source}"
                    )
                row = dict(zip(EVENT_VARIANT_MAP_FIELDS, values))
                try:
                    event_group_id = int(row["event_group_id"])
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Invalid event_group_id in event-variant map {source}"
                    ) from exc
                if event_group_id != row_count:
                    raise ValueError(
                        "Event-variant map event_group_id values must be consecutive "
                        f"from 1 in {source}: expected {row_count}, "
                        f"observed {event_group_id}"
                    )
                if not row["normalization_status"]:
                    raise ValueError(
                        f"Event-variant map has empty normalization_status: {source}"
                    )
                if (
                    row["normalization_status"] == "non_concrete_allele"
                    and row["variant_key"]
                ):
                    raise ValueError(
                        f"Non-concrete event has a canonical variant_key in {source}"
                    )
        expected_count = manifest_count(
            manifest,
            "event_variant_map_count",
            partition,
        )
        event_count = manifest_count(manifest, "event_row_count", partition)
        if row_count != expected_count or row_count != event_count:
            raise ValueError(
                "Event-variant map row count does not match partition manifest: "
                f"partition={partition}, rows={row_count}, "
                f"map_manifest={expected_count}, events={event_count}"
            )
        target = partitions_root / partition_id / EVENT_VARIANT_MAP_FILENAME
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        total_count += row_count
    return total_count


def event_variant_map_manifest(row_count: int, partition_count: int) -> dict[str, object]:
    return {
        "layout": "partitioned",
        "format": "tsv_gzip_v1",
        "path": "event_variant_map/partitions",
        "partition_count": partition_count,
        "row_count": row_count,
        "fields": EVENT_VARIANT_MAP_FIELDS,
        "event_group_id_scope": "partition",
    }


def load_vep_results(root: Path) -> dict[tuple[str, str], dict[str, object]]:
    if not root.is_dir():
        raise NotADirectoryError(f"VEP shard root is not a directory: {root}")
    results = {}
    for directory in sorted(item for item in root.iterdir() if item.is_dir()):
        manifest_path = directory / "manifest.json"
        output = directory / "variant_annotations.tsv.gz"
        if not manifest_path.is_file() or not output.is_file():
            raise FileNotFoundError(f"Incomplete VEP shard output: {directory}")
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("stage") != "annotation" or manifest.get("schema") != VEP_SHARD_SCHEMA:
            raise ValueError(f"Unexpected VEP shard contract: {directory}")
        partition_id = str(manifest.get("partition_id") or "")
        shard_id = str(manifest.get("shard_id") or "")
        key = (partition_id, shard_id)
        if not partition_id or not shard_id or key in results:
            raise ValueError(f"Duplicate or invalid VEP shard identity: {directory}")
        try:
            row_count = int(manifest["row_count"])
            output_size = int(manifest["output"]["size_bytes"])
            fields = list(manifest["output"]["fields"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Invalid VEP shard manifest: {directory}") from exc
        if row_count < 0 or output_size < 1 or output.stat().st_size != output_size:
            raise ValueError(f"VEP shard output changed: {output}")
        with gzip.open(output, "rt", newline="") as handle:
            header = next(csv.reader(handle, delimiter="\t"), None)
        if header != fields:
            raise ValueError(f"VEP shard output header changed: {output}")
        config = manifest.get("config")
        status_counts = manifest.get("status_counts")
        if not isinstance(config, dict) or not config:
            raise ValueError(f"VEP shard has no semantic config: {directory}")
        if not isinstance(status_counts, dict):
            raise ValueError(f"VEP shard has invalid status counts: {directory}")
        normalized_status_counts = {
            str(status): int(count) for status, count in status_counts.items()
        }
        if (
            any(count < 0 for count in normalized_status_counts.values())
            or sum(normalized_status_counts.values()) != row_count
        ):
            raise ValueError(f"VEP shard status counts do not match rows: {directory}")
        results[key] = {
            "directory": directory,
            "manifest": manifest,
            "output": output,
            "row_count": row_count,
            "fields": fields,
            "config": config,
            "status_counts": normalized_status_counts,
        }
    if not results:
        raise ValueError(f"No VEP shard outputs found in {root}")
    return results


def copy_variant_annotation_dataset(
    partitions: list[tuple[Path, dict]],
    vep_root: Path,
    output: Path,
) -> dict[str, object]:
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Variant annotation output is not empty: {output}")
    expected = {}
    partition_order = []
    for partition, manifest in partitions:
        partition_id = manifest_string(manifest, "partition_id", partition)
        partition_order.append(partition_id)
        for shard in variant_annotation_shards(partition, manifest):
            key = (partition_id, str(shard["shard_id"]))
            if key in expected:
                raise ValueError(f"Duplicate annotation shard identity: {key}")
            expected[key] = shard

    observed = load_vep_results(vep_root)
    missing = sorted(set(expected) - set(observed))
    unexpected = sorted(set(observed) - set(expected))
    if missing or unexpected:
        raise ValueError(
            f"VEP shard set differs from annotation inputs: missing={missing}, "
            f"unexpected={unexpected}"
        )

    fields: list[str] | None = None
    config: dict[str, object] | None = None
    status_counts: Counter[str] = Counter()
    partitions_manifest = []
    row_count = 0
    shard_count = 0
    for partition_id in partition_order:
        durable_shards = []
        partition_keys = sorted(key for key in expected if key[0] == partition_id)
        for key in partition_keys:
            source = expected[key]
            result = observed[key]
            manifest = result["manifest"]
            declared_input = manifest.get("input")
            if not isinstance(declared_input, dict):
                raise ValueError(f"VEP shard has invalid input contract: {result['directory']}")
            if (
                int(declared_input.get("size_bytes", -1)) != int(source["size_bytes"])
                or list(declared_input.get("fields", [])) != list(source["fields"])
                or int(result["row_count"]) != int(source["row_count"])
            ):
                raise ValueError(f"VEP shard input contract changed: {result['directory']}")
            if list(result["fields"]) != [*source["fields"], *VEP_FIELDS]:
                raise ValueError(f"VEP shard output fields changed: {result['directory']}")
            if fields is None:
                fields = list(result["fields"])
            elif list(result["fields"]) != fields:
                raise ValueError("VEP shard output fields differ")
            if config is None:
                config = dict(result["config"])
            elif result["config"] != config:
                raise ValueError("VEP semantic configuration differs across shards")

            shard_id = key[1]
            relative = Path("partitions") / partition_id / f"{shard_id}.tsv.gz"
            destination = output / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(result["output"], destination)
            durable_shards.append(
                {
                    "shard_id": shard_id,
                    "path": str(relative),
                    "row_count": int(result["row_count"]),
                    "size_bytes": destination.stat().st_size,
                }
            )
            row_count += int(result["row_count"])
            shard_count += 1
            status_counts.update(result["status_counts"])
        partitions_manifest.append(
            {
                "partition_id": partition_id,
                "shard_count": len(durable_shards),
                "row_count": sum(item["row_count"] for item in durable_shards),
                "shards": durable_shards,
            }
        )

    if fields is None or config is None:
        raise ValueError("VEP annotation dataset is empty")
    if sum(status_counts.values()) != row_count:
        raise ValueError(
            "VEP status counts do not match annotation rows: "
            f"statuses={sum(status_counts.values())}, rows={row_count}"
        )
    descriptor = {
        "schema": VARIANT_DATASET_SCHEMA,
        "status": "complete",
        "layout": "partitioned",
        "format": "tsv_gzip_v1",
        "path": "variant_annotations/manifest.json",
        "partition_count": len(partitions_manifest),
        "shard_count": shard_count,
        "row_count": row_count,
        "fields": fields,
        "vep_config": config,
        "vep_status_counts": dict(sorted(status_counts.items())),
        "partitions": partitions_manifest,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "manifest.json").write_text(
        json.dumps(descriptor, indent=2, sort_keys=True) + "\n"
    )
    return descriptor


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
) -> dict[str, dict[str, float]]:
    by_partition = {}
    for path, manifest in partitions:
        timings = manifest.get("timings_seconds")
        if not isinstance(timings, dict):
            continue
        partition_id = manifest_string(manifest, "partition_id", path)
        normalized = {str(name): float(value) for name, value in timings.items()}
        if any(value < 0 for value in normalized.values()):
            raise ValueError(f"Negative phase timing in annotation partition {partition_id}")
        by_partition[partition_id] = normalized
    return dict(sorted(by_partition.items()))


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    partitions = load_partitions(args.partition_root)
    validate_partition_manifests(partitions)
    first_manifest = partitions[0][1]
    clinvar_vcf_identity = content_identity_from_partition_input(
        args.clinvar_vcf,
        first_manifest["clinvar_vcf"],
        "ClinVar VCF",
    )
    clinvar_tbi_identity = content_identity_from_partition_input(
        args.clinvar_tbi,
        first_manifest["clinvar_tbi"],
        "ClinVar index",
    )
    variant_annotations = copy_variant_annotation_dataset(
        partitions,
        args.vep_root,
        args.outdir / "variant_annotations",
    )
    annotation_count = int(variant_annotations["row_count"])
    event_variant_map_count = copy_event_variant_map_dataset(
        partitions,
        args.outdir / "event_variant_map",
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
    if counts["event_variant_map_count"] != event_variant_map_count:
        raise ValueError(
            "Event-variant map row count does not match partition manifests: "
            f"rows={event_variant_map_count}, "
            f"manifests={counts['event_variant_map_count']}"
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

    gnomad_shared_cache = merge_gnomad_shared_cache(partitions)
    partition_timings = merge_partition_timings(partitions)
    manifest = {
        "created_at": utc_now(),
        "stage": "annotation",
        "schema": "normalized_annotation_evidence_v4",
        "partition_count": len(partitions),
        "partition_ids": [manifest_string(item, "partition_id", path) for path, item in partitions],
        **counts,
        **counters,
        "annotated_variant_context_count": annotation_count,
        "variant_annotations": variant_annotations,
        "event_variant_map": event_variant_map_manifest(
            event_variant_map_count,
            len(partitions),
        ),
        "failure_count": failure_count,
        "clinvar_vcf": clinvar_vcf_identity,
        "clinvar_tbi": clinvar_tbi_identity,
        "gnomad_api_url": first_manifest["gnomad_api_url"],
        "gnomad_dataset": first_manifest["gnomad_dataset"],
        "cache_count_semantics": "ClinVar and gnomAD cache counts are sums across partitions.",
    }
    if gnomad_shared_cache is not None:
        manifest["gnomad_shared_cache"] = gnomad_shared_cache
    if partition_timings:
        manifest["partition_timings_seconds"] = partition_timings
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
