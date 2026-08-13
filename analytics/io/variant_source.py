"""Validated access to compact pre-VEP variant annotation rows."""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path

from analytics.io.artifacts import file_identity, path_metadata


COHORT_VARIANT_SOURCE_KIND = "gaph_cohort_variant_source"
COHORT_VARIANT_SOURCE_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class VariantTableSource:
    paths: tuple[Path, ...]
    columns: tuple[str, ...]
    row_count: int | None
    header: bool
    mode: str
    identity: dict[str, object]


def resolve_pre_vep_variant_source(
    path: Path,
    *,
    required_columns: set[str],
) -> VariantTableSource:
    """Prefer validated VEP input partitions over the enriched merged table."""

    path = path.resolve()
    cohort_paths = cohort_variant_paths(path)
    if cohort_paths is not None:
        sources = [
            resolve_pre_vep_variant_source(member, required_columns=required_columns)
            for member in cohort_paths
        ]
        columns = sources[0].columns
        if any(source.columns != columns for source in sources[1:]):
            raise ValueError("Cohort VEP inputs have different table columns")
        mode = sources[0].mode
        if any(source.mode != mode for source in sources[1:]):
            raise ValueError("Cohort VEP inputs mix partitioned and merged source modes")
        row_count = (
            sum(int(source.row_count) for source in sources)
            if all(source.row_count is not None for source in sources)
            else None
        )
        return VariantTableSource(
            paths=tuple(item for source in sources for item in source.paths),
            columns=columns,
            row_count=row_count,
            header=True,
            mode="cohort_" + mode,
            identity={
                "cohort_descriptor": path_metadata(path),
                "members": [source.identity for source in sources],
            },
        )
    plan_path = path.parent / "plan.json"
    if path.name == "variant_annotations.vep.tsv.gz" and plan_path.exists():
        manifest_path = path.parent / "manifest.json"
        if not manifest_path.exists():
            raise ValueError(f"Incomplete partitioned VEP artifact under {path.parent}")
        plan = _read_json(plan_path)
        manifest = _read_json(manifest_path)
        if plan.get("status") != "complete" or manifest.get("status") != "complete":
            raise ValueError(f"Incomplete partitioned VEP artifact under {path.parent}")
        if manifest.get("source") != plan.get("source"):
            raise ValueError("VEP plan and final manifest sources differ")
        if manifest.get("output") != file_identity(path):
            raise ValueError(f"Bulk VEP output metadata changed: {path}")
        columns = tuple(str(column) for column in plan.get("input_columns", []))
        _require_columns(columns, required_columns, path)
        paths = []
        observed_rows = 0
        for entry in plan.get("partitions", []):
            partition = path.parent / str(entry.get("path", ""))
            expected = dict(entry.get("file", {}))
            if not partition.exists() or file_identity(partition) != expected:
                raise ValueError(f"Prepared VEP input changed: {partition}")
            observed_rows += int(entry.get("row_count", -1))
            paths.append(partition.resolve())
        row_count = int(plan.get("row_count", -1))
        if (
            not paths
            or row_count < 0
            or observed_rows != row_count
            or row_count != int(manifest.get("row_count", -2))
        ):
            raise ValueError("Prepared VEP inputs do not match their plan row count")
        return VariantTableSource(
            paths=tuple(paths),
            columns=columns,
            row_count=row_count,
            header=True,
            mode="vep_input_partitions",
            identity={
                "plan": path_metadata(plan_path),
                "manifest": path_metadata(manifest_path),
                "source": plan.get("source", {}),
            },
        )

    if not path.exists():
        raise FileNotFoundError(path)
    columns = tuple(_read_header(path))
    _require_columns(columns, required_columns, path)
    return VariantTableSource(
        paths=(path,),
        columns=columns,
        row_count=None,
        header=True,
        mode="merged_annotation",
        identity={"input": path_metadata(path)},
    )


def variant_source_sql(source: VariantTableSource) -> str:
    columns = "{" + ",".join(
        f"{sql_string(column)}: 'VARCHAR'" for column in source.columns
    ) + "}"
    paths = "[" + ",".join(sql_string(path) for path in source.paths) + "]"
    return (
        f"read_csv({paths}, delim='\\t', header={'true' if source.header else 'false'}, "
        f"columns={columns}, auto_detect=false, compression='auto', parallel=true, "
        "nullstr='__GAPH_NULL_SENTINEL__')"
    )


def sql_string(value: object) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def cohort_variant_paths(path: Path) -> tuple[Path, ...] | None:
    """Resolve the small descriptor used to union finalized run-level VEP sources."""

    if path.suffix != ".json" or not path.is_file():
        return None
    payload = _read_json(path)
    if payload.get("kind") != COHORT_VARIANT_SOURCE_KIND:
        return None
    if payload.get("schema_version") != COHORT_VARIANT_SOURCE_SCHEMA_VERSION:
        raise ValueError(f"Unsupported cohort variant source schema: {path}")
    raw_members = payload.get("members")
    if not isinstance(raw_members, list) or not raw_members:
        raise ValueError(f"Cohort variant source has no members: {path}")
    members = []
    for raw in raw_members:
        if not isinstance(raw, dict) or not str(raw.get("path") or "").strip():
            raise ValueError(f"Invalid cohort variant source member: {path}")
        member = Path(str(raw["path"])).expanduser()
        if not member.is_absolute():
            member = path.parent / member
        members.append(member.resolve())
    if len(set(members)) != len(members):
        raise ValueError(f"Cohort variant source repeats a member: {path}")
    return tuple(members)


def _read_header(path: Path) -> list[str]:
    handle = gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")
    with handle:
        return next(csv.reader(handle, delimiter="\t"))


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
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value
