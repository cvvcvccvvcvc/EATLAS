"""Validated access to current partitioned variant-annotation datasets."""

from __future__ import annotations

import csv
import gzip
import json
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from analytics.io.artifacts import path_metadata


VARIANT_DATASET_SCHEMA = "gaph_variant_annotation_dataset_v1"
PREPARED_VARIANT_SCHEMA = "gaph_prepared_variant_dataset_v1"


@dataclass(frozen=True)
class VariantTableSource:
    paths: tuple[Path, ...]
    columns: tuple[str, ...]
    row_count: int | None
    header: bool
    mode: str
    identity: dict[str, object]
    partitions: tuple["VariantTablePartition", ...] = ()
    format: str = "tsv_gzip_v1"


@dataclass(frozen=True)
class VariantTablePartition:
    partition_id: str
    paths: tuple[Path, ...]
    row_count: int
    identity: dict[str, object]


def resolve_variant_table_source(
    source_paths: Path | Sequence[Path] | VariantTableSource,
    *,
    required_columns: set[str],
) -> VariantTableSource:
    """Resolve one or more current datasets without materializing a combined copy."""

    if isinstance(source_paths, VariantTableSource):
        missing = required_columns - set(source_paths.columns)
        if missing:
            raise ValueError(f"Variant source lacks required columns: {sorted(missing)}")
        return source_paths
    if not isinstance(source_paths, Path):
        members = tuple(source_paths)
        if not members:
            raise ValueError("Variant source requires at least one dataset")
        if len(members) == 1:
            return _resolve_one_variant_source(
                members[0],
                required_columns=required_columns,
            )
        sources = tuple(
            _resolve_one_variant_source(path, required_columns=required_columns)
            for path in members
        )
        columns = sources[0].columns
        if any(source.columns != columns for source in sources[1:]):
            raise ValueError("Variant annotation datasets have different columns")
        if any(source.format != sources[0].format for source in sources[1:]):
            raise ValueError("Variant annotation datasets have different storage formats")
        row_count = (
            sum(int(source.row_count) for source in sources)
            if all(source.row_count is not None for source in sources)
            else None
        )
        return VariantTableSource(
            paths=tuple(path for source in sources for path in source.paths),
            columns=columns,
            row_count=row_count,
            header=True,
            mode="multi_run_partitioned" if len(sources) > 1 else sources[0].mode,
            identity={"members": [source.identity for source in sources]},
            format=sources[0].format,
            partitions=tuple(partition for source in sources for partition in source.partitions),
        )
    return _resolve_one_variant_source(
        source_paths,
        required_columns=required_columns,
    )


def _resolve_one_variant_source(
    path: Path,
    *,
    required_columns: set[str],
) -> VariantTableSource:
    """Resolve one pipeline dataset or an explicit TSV used by focused tools/tests."""

    path = path.expanduser().resolve()
    if path.name == "manifest.json" and path.is_file():
        manifest = _read_json(path)
        if manifest.get("schema") == PREPARED_VARIANT_SCHEMA:
            return _resolve_prepared_dataset(path, manifest, required_columns)
        if manifest.get("schema") == VARIANT_DATASET_SCHEMA:
            return _resolve_partitioned_dataset(
                path,
                manifest,
                required_columns=required_columns,
            )

    if not path.is_file():
        raise FileNotFoundError(path)
    columns = tuple(_read_header(path))
    _require_columns(columns, required_columns, path)
    return VariantTableSource(
        paths=(path,),
        columns=columns,
        row_count=None,
        header=True,
        mode="explicit_tsv",
        identity={"input": path_metadata(path)},
    )


def _resolve_partitioned_dataset(
    manifest_path: Path,
    manifest: dict[str, object],
    *,
    required_columns: set[str],
) -> VariantTableSource:
    if (
        manifest.get("status") != "complete"
        or manifest.get("layout") != "partitioned"
        or manifest.get("format") != "tsv_gzip_v1"
    ):
        raise ValueError(f"Incomplete variant annotation dataset: {manifest_path}")
    columns = tuple(str(column) for column in manifest.get("fields", []))
    _require_columns(columns, required_columns, manifest_path)
    raw_partitions = manifest.get("partitions")
    if not isinstance(raw_partitions, list) or not raw_partitions:
        raise ValueError(f"Variant annotation dataset has no partitions: {manifest_path}")

    paths = []
    files = []
    partitions = []
    observed_rows = 0
    observed_shards = 0
    seen_paths = set()
    dataset_root = manifest_path.parent.resolve()
    for raw_partition in raw_partitions:
        if not isinstance(raw_partition, dict):
            raise ValueError(f"Invalid variant annotation partition: {manifest_path}")
        partition_id = str(raw_partition.get("partition_id") or "")
        if not partition_id:
            raise ValueError(f"Variant annotation partition has no ID: {manifest_path}")
        if any(partition.partition_id == partition_id for partition in partitions):
            raise ValueError(
                f"Duplicate variant annotation partition ID {partition_id}: {manifest_path}"
            )
        raw_shards = raw_partition.get("shards")
        if not isinstance(raw_shards, list) or not raw_shards:
            raise ValueError(f"Variant annotation partition has no shards: {manifest_path}")
        partition_rows = 0
        partition_paths = []
        partition_files = []
        for raw_shard in raw_shards:
            if not isinstance(raw_shard, dict):
                raise ValueError(f"Invalid variant annotation shard: {manifest_path}")
            relative = Path(str(raw_shard.get("path") or ""))
            if relative.is_absolute() or not relative.parts:
                raise ValueError(f"Unsafe variant annotation shard path: {manifest_path}")
            shard_path = (dataset_root / relative).resolve()
            try:
                shard_path.relative_to(dataset_root)
            except ValueError as exc:
                raise ValueError(
                    f"Variant annotation shard escapes its dataset: {shard_path}"
                ) from exc
            if shard_path in seen_paths:
                raise ValueError(f"Duplicate variant annotation shard: {shard_path}")
            seen_paths.add(shard_path)
            try:
                row_count = int(raw_shard["row_count"])
                size_bytes = int(raw_shard["size_bytes"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"Invalid variant annotation shard: {shard_path}") from exc
            if row_count < 0 or size_bytes < 1 or not shard_path.is_file():
                raise ValueError(f"Invalid variant annotation shard: {shard_path}")
            if shard_path.stat().st_size != size_bytes:
                raise ValueError(f"Variant annotation shard changed: {shard_path}")
            if tuple(_read_header(shard_path)) != columns:
                raise ValueError(f"Variant annotation shard columns changed: {shard_path}")
            paths.append(shard_path)
            metadata = path_metadata(shard_path)
            files.append(metadata)
            partition_paths.append(shard_path)
            partition_files.append(metadata)
            partition_rows += row_count
            observed_rows += row_count
            observed_shards += 1
        if partition_rows != int(raw_partition.get("row_count", -1)):
            raise ValueError(f"Variant annotation partition row count changed: {manifest_path}")
        if len(raw_shards) != int(raw_partition.get("shard_count", -1)):
            raise ValueError(f"Variant annotation partition shard count changed: {manifest_path}")
        partitions.append(
            VariantTablePartition(
                partition_id=partition_id,
                paths=tuple(partition_paths),
                row_count=partition_rows,
                identity={"partition_id": partition_id, "files": partition_files},
            )
        )

    if observed_rows != int(manifest.get("row_count", -1)):
        raise ValueError(f"Variant annotation dataset row count changed: {manifest_path}")
    if observed_shards != int(manifest.get("shard_count", -1)):
        raise ValueError(f"Variant annotation dataset shard count changed: {manifest_path}")
    if len(raw_partitions) != int(manifest.get("partition_count", -1)):
        raise ValueError(f"Variant annotation dataset partition count changed: {manifest_path}")
    return VariantTableSource(
        paths=tuple(paths),
        columns=columns,
        row_count=observed_rows,
        header=True,
        mode="partitioned",
        identity={
            "manifest": path_metadata(manifest_path),
            "files": files,
        },
        partitions=tuple(partitions),
    )


def variant_sources_by_partition(
    source: VariantTableSource,
    partition_ids: Sequence[str],
) -> dict[str, VariantTableSource]:
    """Return the annotation shards owned by each alignment partition."""

    expected = tuple(partition_ids)
    if len(expected) != len(set(expected)):
        raise ValueError(f"Duplicate requested variant partition IDs: {expected}")
    if source.mode == "explicit_tsv":
        if len(expected) != 1:
            raise ValueError(
                "An explicit variant TSV can support only one evidence partition"
            )
        return {expected[0]: source}
    if source.mode != "partitioned":
        raise ValueError(
            f"Annotation support requires one partitioned variant dataset, got {source.mode}"
        )
    observed = tuple(partition.partition_id for partition in source.partitions)
    if set(observed) != set(expected):
        raise ValueError(
            "Alignment evidence and variant annotations have different partition IDs: "
            f"evidence={sorted(expected)}, annotations={sorted(observed)}"
        )
    return {
        partition.partition_id: VariantTableSource(
            paths=partition.paths,
            columns=source.columns,
            row_count=partition.row_count,
            header=source.header,
            mode="partition",
            identity=partition.identity,
        )
        for partition in source.partitions
    }


def variant_source_sql(source: VariantTableSource) -> str:
    if source.format == "parquet":
        paths = "[" + ",".join(sql_string(path) for path in source.paths) + "]"
        return f"read_parquet({paths}, hive_partitioning=false)"
    columns = "{" + ",".join(
        f"{sql_string(column)}: 'VARCHAR'" for column in source.columns
    ) + "}"
    paths = "[" + ",".join(sql_string(path) for path in source.paths) + "]"
    return (
        f"read_csv({paths}, delim='\\t', header={'true' if source.header else 'false'}, "
        f"columns={columns}, auto_detect=false, compression='auto', parallel=true, "
        "nullstr='__GAPH_NULL_SENTINEL__')"
    )


def variant_source_blocks(source: VariantTableSource) -> tuple[VariantTableSource, ...]:
    """Prepared genomic blocks keep all occurrences of an allele together."""
    if source.format != "parquet":
        return (source,)
    groups: dict[str, list[VariantTablePartition]] = {}
    for partition in source.partitions:
        groups.setdefault(partition.partition_id, []).append(partition)
    return tuple(
        VariantTableSource(
            paths=tuple(path for part in parts for path in part.paths),
            columns=source.columns, row_count=sum(part.row_count for part in parts),
            header=True, mode="prepared_block", format="parquet",
            identity={"parent": source.identity, "block": block},
        )
        for block, parts in sorted(groups.items())
    )


def _resolve_prepared_dataset(
    path: Path, manifest: dict[str, object], required_columns: set[str]
) -> VariantTableSource:
    if manifest.get("status") != "complete" or manifest.get("format") != "parquet":
        raise ValueError(f"Incomplete prepared variant dataset: {path}")
    columns = tuple(manifest["fields"])
    _require_columns(columns, required_columns, path)
    paths = []
    partitions = []
    total_rows = 0
    for item in manifest["partitions"]:
        member = (path.parent / item["path"]).resolve()
        member.relative_to(path.parent.resolve())
        if not member.is_file() or member.stat().st_size != item["size_bytes"]:
            raise ValueError(f"Prepared variant shard changed: {member}")
        if member in paths:
            raise ValueError(f"Duplicate prepared variant shard: {member}")
        paths.append(member)
        total_rows += int(item["row_count"])
        partitions.append(VariantTablePartition(
            str(item["partition_id"]), (member,), int(item["row_count"]),
            {"source": manifest["source_id"], "shard": item},
        ))
    if not paths or total_rows != int(manifest["row_count"]):
        raise ValueError(f"Prepared variant row count mismatch: {path}")
    return VariantTableSource(
        tuple(paths), columns, total_rows, True, "prepared",
        {"source_id": manifest["source_id"], "preparation_version": 1},
        tuple(partitions), "parquet",
    )


def sql_string(value: object) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def _read_header(path: Path) -> list[str]:
    handle = gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")
    with handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if not header:
        raise ValueError(f"Variant annotations have no header: {path}")
    return header


def _require_columns(
    columns: tuple[str, ...],
    required_columns: set[str],
    path: Path,
) -> None:
    missing = required_columns - set(columns)
    if missing:
        raise ValueError(
            f"Variant annotations {path} missing columns: {', '.join(sorted(missing))}"
        )


def _read_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid JSON file: {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value
