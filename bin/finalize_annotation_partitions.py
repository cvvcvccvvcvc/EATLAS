#!/usr/bin/env python3
"""Stream annotation partitions into the canonical Stage 3 outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


COUNT_FIELDS = [
    "event_row_count",
    "excluded_non_concrete_event_count",
    "variant_context_count",
    "annotated_variant_context_count",
    "variant_strategy_support_count",
    "variant_strategy_support_missing_key_count",
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partition-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


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
        partition_id = str(manifest.get("partition_id") or path.name)
        if partition_id in seen_ids:
            raise ValueError(f"Duplicate annotation partition_id: {partition_id}")
        seen_ids.add(partition_id)
        partitions.append((path, manifest))
    if not partitions:
        raise ValueError(f"No annotation partitions found in {root}")
    return partitions


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


def merge_gnomad_shared_cache(partitions: list[tuple[Path, dict]]) -> dict[str, object] | None:
    snapshots = [manifest.get("gnomad_shared_cache") for _path, manifest in partitions]
    if not any(snapshot is not None for snapshot in snapshots):
        return None
    if not all(isinstance(snapshot, dict) for snapshot in snapshots):
        raise ValueError("gnomAD shared-cache metadata is missing from some annotation partitions")

    first = snapshots[0]
    identity = {field: first.get(field) for field in GNOMAD_SHARED_CACHE_IDENTITY_FIELDS}
    for snapshot in snapshots[1:]:
        observed = {field: snapshot.get(field) for field in GNOMAD_SHARED_CACHE_IDENTITY_FIELDS}
        if observed != identity:
            raise ValueError("gnomAD shared-cache configuration differs across annotation partitions")
    return {
        **identity,
        **{
            field: sum(int(snapshot.get(field) or 0) for snapshot in snapshots)
            for field in GNOMAD_SHARED_CACHE_COUNT_FIELDS
        },
    }


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    partitions = load_partitions(args.partition_root)
    annotation_count = merge_tsv_gz(
        partitions,
        "variant_annotations.tsv.gz",
        args.outdir / "variant_annotations.tsv.gz",
    )
    support_count = merge_tsv_gz(
        partitions,
        "variant_strategy_support.tsv.gz",
        args.outdir / "variant_strategy_support.tsv.gz",
    )
    failure_count = merge_tsv_gz(partitions, "failures.tsv.gz", args.outdir / "failures.tsv.gz")

    counts = {
        field: sum(int(manifest.get(field) or 0) for _path, manifest in partitions)
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
    manifest_failure_count = sum(
        int(manifest.get("failure_count") or 0) for _path, manifest in partitions
    )
    if manifest_failure_count != failure_count:
        raise ValueError(
            "Failure row count does not match partition manifests: "
            f"rows={failure_count}, manifests={manifest_failure_count}"
        )
    counters = {}
    for field in COUNTER_FIELDS:
        counter = Counter()
        for _path, manifest in partitions:
            counter.update({key: int(value) for key, value in (manifest.get(field) or {}).items()})
        counters[field] = dict(counter)

    first_manifest = partitions[0][1]
    gnomad_shared_cache = merge_gnomad_shared_cache(partitions)
    manifest = {
        "created_at": utc_now(),
        "output_mode": "unique_variant_context_partitioned",
        "partition_count": len(partitions),
        "partition_ids": [str(item.get("partition_id") or path.name) for path, item in partitions],
        **counts,
        **counters,
        "annotated_variant_context_count": annotation_count,
        "variant_strategy_support_count": support_count,
        "failure_count": failure_count,
        "clinvar_vcf": first_manifest.get("clinvar_vcf", {}),
        "clinvar_tbi": first_manifest.get("clinvar_tbi", {}),
        "gnomad_api_url": first_manifest.get("gnomad_api_url", ""),
        "gnomad_dataset": first_manifest.get("gnomad_dataset", ""),
        "cache_count_semantics": "ClinVar and gnomAD cache counts are sums across partitions.",
    }
    if gnomad_shared_cache is not None:
        manifest["gnomad_shared_cache"] = gnomad_shared_cache
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
