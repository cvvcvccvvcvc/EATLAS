#!/usr/bin/env python3
"""Build a resumable, partitioned local-VEP annotation artifact."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import tempfile
import time
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from analytics.io.artifacts import path_metadata
from genomics.variants import parse_variant_key
from analytics.annotation.vep import annotate_vep_consequences
from analytics.annotation.vep_result_cache import DEFAULT_TILE_SIZE_BP


SCHEMA_VERSION = 1
DEFAULT_PARTITION_SIZE = 250_000
REQUIRED_COLUMNS = {"variant_key", "gene_id"}
VEP_COLUMNS = [
    "vep_status",
    "vep_primary_consequence",
    "vep_consequence_terms",
    "vep_transcript_id",
    "vep_mane_select",
    "vep_canonical",
    "vep_impact",
    "vep_variant_class",
]
VEP_RENAME = {
    "status": "vep_status",
    "primary_consequence": "vep_primary_consequence",
    "consequence_terms": "vep_consequence_terms",
    "transcript_id": "vep_transcript_id",
    "mane_select": "vep_mane_select",
    "canonical": "vep_canonical",
    "impact": "vep_impact",
    "variant_class": "vep_variant_class",
}


def prepare_partitions(
    *,
    annotation_tsv: Path,
    outdir: Path,
    partition_size: int = DEFAULT_PARTITION_SIZE,
) -> dict[str, object]:
    """Split the annotation table into deterministic, independently runnable inputs."""

    if partition_size < 1:
        raise ValueError("VEP partition size must be >= 1")
    if not annotation_tsv.exists():
        raise FileNotFoundError(annotation_tsv)

    source = path_metadata(annotation_tsv)
    plan_path = outdir / "plan.json"
    existing = _read_json(plan_path)
    if existing is not None:
        _validate_existing_plan(existing, source, partition_size)
        if _prepared_inputs_complete(existing, outdir):
            return {**existing, "cache_hit": True}

    inputs_dir = outdir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    existing_partitions = {
        str(entry["partition_id"]): entry
        for entry in (existing or {}).get("partitions", [])
    }
    partitions = []
    total_rows = 0
    input_columns: list[str] | None = None
    for index, frame in enumerate(
        pd.read_csv(
            annotation_tsv,
            sep="\t",
            compression="gzip",
            dtype=str,
            keep_default_na=False,
            chunksize=partition_size,
        ),
        start=1,
    ):
        if input_columns is None:
            input_columns = list(frame.columns)
            missing = REQUIRED_COLUMNS - set(input_columns)
            if missing:
                raise ValueError(
                    f"Variant annotations missing columns: {', '.join(sorted(missing))}"
                )
            conflicts = set(VEP_COLUMNS) & set(input_columns)
            if conflicts:
                raise ValueError(
                    f"Variant annotations already contain VEP columns: {', '.join(sorted(conflicts))}"
                )
        partition_id = f"partition_{index:06d}"
        relative_path = Path("inputs") / f"{partition_id}.tsv.gz"
        partition_path = outdir / relative_path
        row_count = len(frame)
        previous = existing_partitions.get(partition_id)
        if (
            previous is not None
            and int(previous.get("row_count", -1)) == row_count
            and str(previous.get("path", "")) == str(relative_path)
            and partition_path.exists()
            and _file_identity(partition_path) == previous.get("file")
        ):
            partitions.append(previous)
        else:
            _write_tsv_atomic(frame, partition_path, header=True)
            partitions.append(
                {
                    "partition_id": partition_id,
                    "path": str(relative_path),
                    "row_count": row_count,
                    "file": _file_identity(partition_path),
                }
            )
        total_rows += row_count

    if input_columns is None:
        raise ValueError(f"Variant annotations are empty: {annotation_tsv}")
    _validate_source_row_count(annotation_tsv, total_rows)
    plan = {
        "schema_version": SCHEMA_VERSION,
        "status": "complete",
        "created_at": _timestamp(),
        "source": source,
        "partition_size": partition_size,
        "partition_count": len(partitions),
        "row_count": total_rows,
        "input_columns": input_columns,
        "output_columns": [*input_columns, *VEP_COLUMNS],
        "partitions": partitions,
    }
    _write_json_atomic(plan_path, plan)
    return {**plan, "cache_hit": False}


def annotate_partition(
    *,
    outdir: Path,
    partition_index: int,
    backend: str,
    release: str,
    vep_executable: str | Path = "vep",
    vep_cache_dir: Path | None = None,
    vep_forks: int = 1,
    rest_workers: int = 2,
    vep_result_cache_dir: Path | None = None,
    vep_result_cache_tile_size_bp: int = DEFAULT_TILE_SIZE_BP,
) -> dict[str, object]:
    """Annotate one prepared partition and publish it atomically."""

    if backend not in {"rest", "local"}:
        raise ValueError("VEP backend must be 'rest' or 'local'")
    if not release:
        raise ValueError("VEP release is required")
    if vep_forks < 1 or rest_workers < 1:
        raise ValueError("VEP worker counts must be >= 1")
    if backend == "local" and vep_cache_dir is None:
        raise ValueError("Local VEP requires a cache directory")
    plan = _require_plan(outdir)
    if partition_index < 1 or partition_index > int(plan["partition_count"]):
        raise ValueError(
            f"Partition index must be between 1 and {plan['partition_count']}: {partition_index}"
        )
    entry = plan["partitions"][partition_index - 1]
    partition_id = str(entry["partition_id"])
    input_path = outdir / str(entry["path"])
    _require_file_identity(input_path, dict(entry["file"]), f"input {partition_id}")

    config = {
        "backend": backend,
        "release": str(release),
        "vep_executable": str(vep_executable) if backend == "local" else "",
        "vep_cache_dir": str(vep_cache_dir) if backend == "local" else "",
        "vep_forks": int(vep_forks) if backend == "local" else 0,
    }
    output_dir = outdir / "partitions"
    output_path = output_dir / f"{partition_id}.tsv.gz"
    manifest_path = output_dir / f"{partition_id}.json"
    existing = _read_json(manifest_path)
    if existing is not None:
        if existing.get("config") != config or existing.get("input") != entry:
            raise ValueError(f"Existing VEP partition has a different contract: {partition_id}")
        _require_file_identity(
            output_path,
            dict(existing.get("output", {})),
            f"output {partition_id}",
        )
        return {**existing, "cache_hit": True}

    frame = pd.read_csv(
        input_path,
        sep="\t",
        compression="gzip",
        dtype=str,
        keep_default_na=False,
    )
    if list(frame.columns) != list(plan["input_columns"]):
        raise ValueError(f"Input columns changed for {partition_id}")
    if len(frame) != int(entry["row_count"]):
        raise ValueError(f"Input row count changed for {partition_id}")

    requests, invalid_keys = _vep_requests(frame)
    started = time.perf_counter()
    if requests.empty:
        annotations = _empty_vep_annotations()
        summary = {
            "status": "complete",
            "backend": backend,
            "release": str(release),
            "requested": 0,
            "queried": 0,
            "cached": 0,
            "status_counts": {},
        }
    else:
        # Partition-local VEP input, output, and SQLite are disposable. Keep
        # them on the compute node to avoid NFS cleanup races between VEP forks.
        with tempfile.TemporaryDirectory(prefix="gaph_vep_") as temporary:
            annotations, summary = annotate_vep_consequences(
                requests,
                Path(temporary) / "partition.sqlite",
                backend=backend,
                release=release,
                max_workers=rest_workers,
                vep_executable=vep_executable,
                vep_cache_dir=vep_cache_dir,
                vep_forks=vep_forks,
                vep_result_cache_dir=vep_result_cache_dir,
                vep_result_cache_tile_size_bp=vep_result_cache_tile_size_bp,
            )
    enriched = _merge_annotations(frame, annotations, invalid_keys)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv_atomic(enriched, output_path, header=False)
    status_counts = {
        str(status): int(count)
        for status, count in enriched["vep_status"].value_counts().sort_index().items()
    }
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "status": "complete",
        "created_at": _timestamp(),
        "partition_id": partition_id,
        "config": config,
        "input": entry,
        "row_count": len(enriched),
        "valid_variant_key_count": len(enriched) - invalid_keys,
        "invalid_variant_key_count": invalid_keys,
        "status_counts": status_counts,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "vep": summary,
        "output_columns": list(enriched.columns),
        "output": _file_identity(output_path),
    }
    _write_json_atomic(manifest_path, manifest)
    return {**manifest, "cache_hit": False}


def finalize_annotations(*, outdir: Path) -> dict[str, object]:
    """Validate all partitions and concatenate their gzip members."""

    plan = _require_plan(outdir)
    final_path = outdir / "variant_annotations.vep.tsv.gz"
    manifest_path = outdir / "manifest.json"
    existing = _read_json(manifest_path)

    partition_manifests = []
    partition_manifest_files = []
    status_counts: Counter[str] = Counter()
    config: dict[str, object] | None = None
    total_rows = 0
    output_columns = list(plan["output_columns"])
    for entry in plan["partitions"]:
        partition_id = str(entry["partition_id"])
        partition_manifest_path = outdir / "partitions" / f"{partition_id}.json"
        partition_manifest = _read_json(partition_manifest_path)
        if partition_manifest is None:
            raise FileNotFoundError(f"Missing VEP partition manifest: {partition_manifest_path}")
        if partition_manifest.get("input") != entry:
            raise ValueError(f"VEP partition input contract changed: {partition_id}")
        if partition_manifest.get("output_columns") != output_columns:
            raise ValueError(f"VEP partition columns changed: {partition_id}")
        if config is None:
            config = dict(partition_manifest["config"])
        elif partition_manifest.get("config") != config:
            raise ValueError(f"VEP partition configuration changed: {partition_id}")
        output_path = outdir / "partitions" / f"{partition_id}.tsv.gz"
        _require_file_identity(
            output_path,
            dict(partition_manifest.get("output", {})),
            f"output {partition_id}",
        )
        total_rows += int(partition_manifest["row_count"])
        status_counts.update(
            {
                str(status): int(count)
                for status, count in dict(partition_manifest["status_counts"]).items()
            }
        )
        partition_manifests.append(partition_manifest)
        partition_manifest_files.append(
            {
                "partition_id": partition_id,
                "file": _file_identity(partition_manifest_path),
            }
        )

    if total_rows != int(plan["row_count"]):
        raise ValueError(
            f"Final VEP row count mismatch: observed {total_rows}, expected {plan['row_count']}"
        )
    plan_identity = _file_identity(outdir / "plan.json")
    if (
        existing is not None
        and existing.get("plan") == plan_identity
        and existing.get("partition_manifests") == partition_manifest_files
    ):
        _require_file_identity(final_path, dict(existing.get("output", {})), "final VEP output")
        return {**existing, "cache_hit": True}
    _concatenate_gzip_members(
        destination=final_path,
        columns=output_columns,
        parts=[
            outdir / "partitions" / f"{entry['partition_id']}.tsv.gz"
            for entry in plan["partitions"]
        ],
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "status": "complete",
        "created_at": _timestamp(),
        "source": plan["source"],
        "plan": plan_identity,
        "partition_manifests": partition_manifest_files,
        "config": config or {},
        "partition_count": len(partition_manifests),
        "row_count": total_rows,
        "status_counts": dict(sorted(status_counts.items())),
        "columns": output_columns,
        "output": _file_identity(final_path),
    }
    _write_json_atomic(manifest_path, manifest)
    return {**manifest, "cache_hit": False}


def _vep_requests(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    rows = []
    invalid = 0
    for row in frame[["variant_key", "gene_id"]].itertuples(index=False):
        parsed = parse_variant_key(row.variant_key)
        if parsed is None:
            invalid += 1
            continue
        chrom, pos, ref, alt = parsed
        rows.append(
            {
                "variant_key": str(row.variant_key),
                "gene_id": str(row.gene_id),
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
            }
        )
    return pd.DataFrame(rows, columns=["variant_key", "gene_id", "chrom", "pos", "ref", "alt"]), invalid


def _empty_vep_annotations() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "variant_key",
            "gene_id",
            "status",
            "primary_consequence",
            "consequence_terms",
            "transcript_id",
            "mane_select",
            "canonical",
            "impact",
            "variant_class",
        ]
    )


def _merge_annotations(
    frame: pd.DataFrame,
    annotations: pd.DataFrame,
    invalid_key_count: int,
) -> pd.DataFrame:
    renamed = annotations.rename(columns=VEP_RENAME)
    keep = ["variant_key", "gene_id", *VEP_COLUMNS]
    enriched = frame.assign(_source_order=range(len(frame))).merge(
        renamed[keep],
        on=["variant_key", "gene_id"],
        how="left",
        sort=False,
        validate="many_to_one",
    )
    enriched = enriched.sort_values("_source_order", kind="stable").drop(columns="_source_order")
    missing = enriched["vep_status"].isna()
    if int(missing.sum()) != invalid_key_count:
        raise ValueError("Local VEP output did not cover every valid variant/gene pair")
    enriched.loc[missing, "vep_status"] = "invalid_variant_key"
    for column in VEP_COLUMNS:
        if column == "vep_canonical":
            enriched[column] = enriched[column].map(
                lambda value: bool(value) if pd.notna(value) else False
            )
        else:
            enriched[column] = enriched[column].fillna("").astype(str)
    return enriched


def _validate_existing_plan(
    plan: dict[str, object],
    source: dict[str, object],
    partition_size: int,
) -> None:
    if plan.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Existing VEP plan has a different schema version")
    if plan.get("source") != source or int(plan.get("partition_size", 0)) != partition_size:
        raise ValueError("Existing VEP plan was prepared from a different input contract")


def _prepared_inputs_complete(plan: dict[str, object], outdir: Path) -> bool:
    try:
        for entry in plan.get("partitions", []):
            _require_file_identity(
                outdir / str(entry["path"]),
                dict(entry["file"]),
                f"input {entry['partition_id']}",
            )
    except (FileNotFoundError, ValueError):
        return False
    return True


def _validate_source_row_count(annotation_tsv: Path, observed: int) -> None:
    manifest_path = annotation_tsv.parent / "manifest.json"
    manifest = _read_json(manifest_path)
    if manifest is None or manifest.get("variant_context_count") is None:
        return
    expected = int(manifest["variant_context_count"])
    if observed != expected:
        raise ValueError(
            f"Variant annotation row count does not match its manifest: {observed} != {expected}"
        )


def _require_plan(outdir: Path) -> dict[str, object]:
    plan_path = outdir / "plan.json"
    plan = _read_json(plan_path)
    if plan is None:
        raise FileNotFoundError(f"VEP partition plan not found: {plan_path}")
    if plan.get("schema_version") != SCHEMA_VERSION or plan.get("status") != "complete":
        raise ValueError(f"VEP partition plan is not complete: {plan_path}")
    return plan


def _file_identity(path: Path) -> dict[str, int]:
    stat = path.stat()
    return {"size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _require_file_identity(path: Path, expected: dict[str, object], label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")
    if _file_identity(path) != expected:
        raise ValueError(f"File metadata changed for {label}: {path}")


def _write_tsv_atomic(frame: pd.DataFrame, path: Path, *, header: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(path, ".tsv.gz")
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            header=header,
            lineterminator="\n",
            compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
        )
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _concatenate_gzip_members(
    *,
    destination: Path,
    columns: list[str],
    parts: list[Path],
) -> None:
    temporary = _temporary_path(destination, ".tsv.gz")
    try:
        with temporary.open("wb") as output:
            with gzip.GzipFile(fileobj=output, mode="wb", mtime=0) as header:
                header.write(("\t".join(columns) + "\n").encode())
            for part in parts:
                with part.open("rb") as source:
                    shutil.copyfileobj(source, output, length=16 * 1024 * 1024)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def _write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(path, ".json")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _temporary_path(destination: Path, suffix: str) -> Path:
    handle = tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.",
        suffix=suffix,
        dir=destination.parent,
        delete=False,
    )
    handle.close()
    return Path(handle.name)


def _read_json(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid JSON file: {path}") from exc


def _timestamp() -> str:
    return datetime.now(UTC).isoformat()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare", help="Create deterministic VEP input partitions.")
    _add_paths(prepare)
    prepare.add_argument(
        "--annotation-tsv",
        type=Path,
        help="Source variant_annotations.tsv.gz. Default: <run-dir>/annotation/variant_annotations.tsv.gz.",
    )
    prepare.add_argument("--partition-size", type=int, default=DEFAULT_PARTITION_SIZE)

    annotate = subparsers.add_parser("annotate", help="Annotate one prepared partition.")
    _add_paths(annotate)
    annotate.add_argument("--partition-index", type=int, required=True, help="One-based partition index.")
    annotate.add_argument(
        "--vep-backend",
        choices=("rest", "local"),
        default=os.environ.get("GAPH_VEP_BACKEND", "rest"),
    )
    annotate.add_argument("--vep-release", default=os.environ.get("GAPH_VEP_RELEASE") or None)
    annotate.add_argument(
        "--vep-executable",
        default=os.environ.get("GAPH_VEP_EXECUTABLE", "vep"),
    )
    annotate.add_argument(
        "--vep-cache-dir",
        type=Path,
        default=os.environ.get("GAPH_VEP_CACHE_DIR") or None,
    )
    annotate.add_argument(
        "--vep-forks",
        type=int,
        default=int(os.environ.get("GAPH_VEP_FORKS", "4")),
    )
    annotate.add_argument(
        "--vep-result-cache-dir",
        type=Path,
        default=os.environ.get("GAPH_VEP_RESULT_CACHE_DIR") or None,
    )
    annotate.add_argument(
        "--vep-result-cache-tile-size-bp",
        type=int,
        default=int(
            os.environ.get(
                "GAPH_VEP_RESULT_CACHE_TILE_SIZE_BP",
                str(DEFAULT_TILE_SIZE_BP),
            )
        ),
    )
    annotate.add_argument("--rest-workers", type=int, default=2)

    finalize = subparsers.add_parser("finalize", help="Validate and join completed partitions.")
    _add_paths(finalize)
    return parser.parse_args()


def _add_paths(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument(
        "--outdir",
        type=Path,
        help="Artifact directory. Default: <run-dir>/analytics/vep_consequences.",
    )


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.expanduser().resolve()
    outdir = (
        args.outdir.expanduser().resolve()
        if args.outdir is not None
        else run_dir / "analytics" / "vep_consequences"
    )
    if args.command == "prepare":
        annotation_tsv = (
            args.annotation_tsv.expanduser().resolve()
            if args.annotation_tsv is not None
            else run_dir / "annotation" / "variant_annotations.tsv.gz"
        )
        result = prepare_partitions(
            annotation_tsv=annotation_tsv,
            outdir=outdir,
            partition_size=args.partition_size,
        )
    elif args.command == "annotate":
        if not args.vep_release:
            raise ValueError("--vep-release is required for partitioned VEP annotation")
        result = annotate_partition(
            outdir=outdir,
            partition_index=args.partition_index,
            backend=args.vep_backend,
            release=args.vep_release,
            vep_executable=args.vep_executable,
            vep_cache_dir=args.vep_cache_dir,
            vep_forks=args.vep_forks,
            rest_workers=args.rest_workers,
            vep_result_cache_dir=args.vep_result_cache_dir,
            vep_result_cache_tile_size_bp=args.vep_result_cache_tile_size_bp,
        )
    else:
        result = finalize_annotations(outdir=outdir)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
