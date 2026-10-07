"""Verified Google Drive sources with bounded, exclusively owned restore space."""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
from dataclasses import fields
from pathlib import Path

import duckdb

from analytics.io.alignment_aggregates import resolve_alignment_aggregate_paths
from analytics.io.annotation_support import resolve_annotation_support_paths
from analytics.io.artifacts import content_identity, write_json_atomic
from analytics.io.duckdb import available_cpu_count, configure_duckdb_memory
from analytics.io.performance import PerformanceProfile, profile_stage
from analytics.io.run_inputs import SourceRun, _resolve_source_run
from analytics.io.taxonomy_summary import resolve_taxonomy_summary_path
from analytics.io.variant_source import (
    PREPARED_VARIANT_SCHEMA, resolve_variant_table_source, sql_string, variant_source_sql,
)
from run_archiving.archive import build_snapshot, read_archive_manifest, restore_run
from run_archiving.rclone import RcloneClient, validate_remote_root


PREPARATION_VERSION = 1
BLOCK_BASES = 10_000_000
PACK_SCHEMA = "gaph_prepared_report_source_v1"
PATH_FIELDS = {
    "root_manifest_json", "fetch_manifest_json", "genes_tsv", "target_features_tsv",
    "target_sequences_dir", "taxonomy_tsv", "orthologs_selected_tsv",
    "variant_annotations_source", "annotation_manifest_json", "annotation_failures_tsv",
    "alignment_manifest_json",
}


class ArchiveSources:
    """Never owns user run-dir inputs; only restores below its private workspace."""

    def __init__(
        self, *, analytics_root: Path, remote_root: str, calculation_identity: dict,
        cache_policy: str, disk_budget_bytes: int, workers: int,
        client: RcloneClient | None = None,
        performance_profile: PerformanceProfile | None = None,
    ) -> None:
        self.remote_root = validate_remote_root(remote_root)
        archive_env = os.environ.get("GAPH_ARCHIVE_ENV")
        if not archive_env and os.environ.get("GAPH_ROOT"):
            archive_env = str(Path(os.environ["GAPH_ROOT"]) / "envs/run-archiving")
        self.client = client or RcloneClient(
            executable=str(Path(archive_env) / "bin/rclone") if archive_env else "rclone",
            config_path=Path(os.environ.get("RCLONE_CONFIG", Path.home() / ".config/rclone/rclone.conf")),
            transfers=2, checkers=2,
        )
        self.client.require_google_drive(self.remote_root)
        self.calculation_id = hashlib.sha256(
            json.dumps(calculation_identity, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()[:24]
        self.cache_policy = cache_policy
        analytics_root = analytics_root.resolve()
        self.root = analytics_root / ("cache" if cache_policy == "keep" else "work") / "gdrive"
        self.root = self.root / self.calculation_id
        self.root.resolve().relative_to(analytics_root)
        self.disk_budget_bytes = disk_budget_bytes
        self.workers = workers
        self.owned: list[Path] = []
        self.performance_profile = performance_profile

    def prepare(self, run_ids: list[str]) -> tuple[SourceRun, ...]:
        if not run_ids or len(set(run_ids)) != len(run_ids):
            raise ValueError("GDrive run IDs must be nonempty and unique")
        for run_id in run_ids:
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id):
                raise ValueError(f"Invalid GDrive run ID: {run_id}")
        archives = [read_archive_manifest(self.client, remote_root=self.remote_root, run_id=run_id)
                    for run_id in run_ids]
        sources = []
        for run_id, archive in zip(run_ids, archives):
            owner = self.root / archive["tree_sha256"]
            owner.resolve().relative_to(self.root.resolve())
            pack = owner / "prepared" / run_id
            self.owned.append(owner)
            receipt = pack / "prepared_source.json"
            if receipt.is_file():
                source = load_prepared_source(pack, archive["tree_sha256"])
                self._remove_owned(owner / "restore" / run_id)
                print(f"Prepared source cache hit: {run_id}", flush=True)
                sources.append(source)
                continue
            self._require_capacity(int(archive["total_bytes"]) * 2)
            raw = owner / "restore" / run_id
            print(f"Restoring Google Drive source: {run_id}", flush=True)
            with profile_stage(self.performance_profile, f"Restore and verify [{run_id}]") as timing:
                if raw.exists():
                    if build_snapshot(raw).tree_sha256 != archive["tree_sha256"]:
                        raise ValueError(f"Existing restore differs from its archive: {raw}")
                else:
                    restore_run(self.client, remote_root=self.remote_root, run_id=run_id, destination=raw)
                timing["metrics"] = {"restored_bytes": archive["total_bytes"]}
            source = _resolve_source_run(raw)
            partial = pack.with_name(run_id + ".partial")
            partial.mkdir(parents=True, exist_ok=True)
            with profile_stage(self.performance_profile, f"Prepare report source [{run_id}]"):
                prepare_source_pack(source, partial, workers=self.workers)
            payload = source_payload(source)
            payload.update({
                "schema": PACK_SCHEMA, "preparation_version": PREPARATION_VERSION,
                "archive_tree_sha256": archive["tree_sha256"],
                "archive_uri": f"{self.remote_root}/runs/{run_id}",
                "outputs": {
                    str(p.relative_to(partial)): {**content_identity(p), "mtime_ns": p.stat().st_mtime_ns}
                    for p in sorted(partial.rglob("*"))
                    if p.is_file() and p.name != "prepared_source.json"
                },
            })
            write_json_atomic(partial / "prepared_source.json", payload)
            partial.rename(pack)
            prepared = load_prepared_source(pack, archive["tree_sha256"])
            self._remove_owned(raw)
            self._require_capacity(0)
            sources.append(prepared)
            print(f"Prepared {run_id}; verified restore released", flush=True)
            if self.performance_profile is not None:
                self.performance_profile.checkpoint({
                    "prepared_sources": len(sources), "source_count": len(run_ids),
                    "prepared_bytes": sum(x["size_bytes"] for x in payload["outputs"].values()),
                })
        return tuple(sources)

    def finish(self) -> None:
        if self.cache_policy == "discard":
            for path in self.owned:
                self._remove_owned(path)

    def _require_capacity(self, additional_bytes: int) -> None:
        used = sum(p.stat().st_size for p in self.root.rglob("*") if p.is_file()) if self.root.exists() else 0
        if used + additional_bytes > self.disk_budget_bytes:
            raise ValueError(
                f"GDrive workspace budget exceeded: {used} used + {additional_bytes} reserved "
                f"> {self.disk_budget_bytes} bytes. Complete checkpoints remain resumable."
            )

    def _remove_owned(self, path: Path) -> None:
        resolved = path.resolve()
        resolved.relative_to(self.root.resolve())
        if resolved == self.root.resolve() or path.is_symlink():
            raise ValueError(f"Refusing unsafe staging cleanup: {path}")
        if path.exists():
            shutil.rmtree(path)


def prepare_source_pack(source: SourceRun, destination: Path, *, workers: int) -> None:
    """Retain report inputs and derivations, without raw alignment or event lineage."""
    derived = destination / "derived"
    resolve_alignment_aggregate_paths(source.run_dir, analytics_dir=derived)
    resolve_annotation_support_paths(source.run_dir, analytics_dir=derived, workers=workers)
    resolve_taxonomy_summary_path(source.run_dir, analytics_dir=derived)
    shutil.copytree(source.run_dir / "fetch", destination / "fetch", dirs_exist_ok=True)
    for relative in (
        "run_manifest.json", "evidence_inventory.json", "alignment/manifest.json",
        "annotation/manifest.json", "annotation/failures.tsv.gz",
    ):
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source.run_dir / relative, target)
    prepare_variant_parquet(source, destination / "variant_rows")
    support = resolve_variant_table_source(
        derived / "annotation_support" / "variant_strategy_support.tsv.gz",
        required_columns={"variant_key"},
    )
    prepare_parquet_dataset(support, destination / "filter_support", source.source_id + ":support")


def prepare_variant_parquet(source: SourceRun, destination: Path) -> None:
    variants = resolve_variant_table_source(source.variant_annotations_source, required_columns={"variant_key"})
    prepare_parquet_dataset(variants, destination, source.source_id)


def prepare_parquet_dataset(variants, destination: Path, source_id: str) -> None:
    """Repartition typed-as-text input once, preserving every value and duplicate."""
    manifest_path = destination / "manifest.json"
    if manifest_path.is_file():
        cached = resolve_variant_table_source(manifest_path, required_columns={"variant_key"})
        if cached.identity.get("source_id") == source_id:
            return
        raise ValueError(f"Prepared variant source identity mismatch: {destination}")
    if destination.exists():
        if destination.is_symlink():
            raise ValueError(f"Unsafe variant preparation path: {destination}")
        shutil.rmtree(destination)
    destination.mkdir(parents=True)
    with duckdb.connect() as connection:
        cpus = available_cpu_count()
        connection.execute(f"SET threads={cpus}")
        configure_duckdb_memory(connection, cpus)
        connection.execute("SET preserve_insertion_order=false")
        connection.execute(f"SET temp_directory={sql_string(destination / '.tmp')}")
        connection.execute(
            f"COPY (SELECT *, CASE WHEN regexp_full_match(split_part(variant_key, ':', 1), '[A-Za-z0-9]+') "
            "AND try_cast(split_part(variant_key, ':', 2) AS BIGINT)>0 THEN "
            "split_part(variant_key, ':', 1)||'_'||cast(floor((cast(split_part(variant_key, ':', 2) AS BIGINT)-1)/"
            f"{BLOCK_BASES}) AS BIGINT) ELSE 'invalid' END AS report_block FROM {variant_source_sql(variants)}) "
            f"TO {sql_string(destination)} (FORMAT PARQUET, COMPRESSION ZSTD, PARTITION_BY (report_block), "
            "ROW_GROUP_SIZE 100000)"
        )
        partitions = []
        for path in sorted(destination.glob("report_block=*/*.parquet")):
            rows = connection.execute(f"SELECT count(*) FROM read_parquet({sql_string(path)}, hive_partitioning=false)").fetchone()[0]
            partitions.append({
                "partition_id": path.parent.name.removeprefix("report_block="),
                "path": str(path.relative_to(destination)), "row_count": rows,
                "size_bytes": path.stat().st_size,
            })
        if not partitions:
            empty_path = destination / "empty.parquet"
            connection.execute(
                f"COPY (SELECT * FROM {variant_source_sql(variants)} WHERE false) "
                f"TO {sql_string(empty_path)} (FORMAT PARQUET, COMPRESSION ZSTD)"
            )
            partitions.append({"partition_id": "invalid", "path": empty_path.name,
                               "row_count": 0, "size_bytes": empty_path.stat().st_size})
    row_count = sum(x["row_count"] for x in partitions)
    if variants.row_count is not None and row_count != variants.row_count:
        raise ValueError("Prepared variant row count differs from source evidence")
    write_json_atomic(destination / "manifest.json", {
        "schema": PREPARED_VARIANT_SCHEMA, "status": "complete", "format": "parquet",
        "source_id": source_id, "fields": list(variants.columns),
        "row_count": row_count, "partitions": partitions,
    })


def source_payload(source: SourceRun) -> dict:
    values = {}
    for field in fields(SourceRun):
        value = getattr(source, field.name)
        if field.name in PATH_FIELDS:
            value = str(value.relative_to(source.run_dir))
        elif isinstance(value, frozenset):
            value = sorted(value)
        elif field.name in {"run_dir", "prepared_cache_dir", "archive_uri"}:
            continue
        values[field.name] = value
    values["variant_annotations_source"] = "variant_rows/manifest.json"
    return {"source": values}


def load_prepared_source(directory: Path, archive_tree_sha256: str) -> SourceRun:
    payload = json.loads((directory / "prepared_source.json").read_text())
    if (
        payload.get("schema") != PACK_SCHEMA
        or payload.get("preparation_version") != PREPARATION_VERSION
        or payload.get("archive_tree_sha256") != archive_tree_sha256
    ):
        raise ValueError(f"Prepared source identity mismatch: {directory}")
    observed = {str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file()}
    if directory.is_symlink() or any(p.is_symlink() for p in directory.rglob("*")):
        raise ValueError(f"Prepared source must not contain symlinks: {directory}")
    if observed != set(payload["outputs"]) | {"prepared_source.json"}:
        raise ValueError(f"Prepared source membership changed: {directory}")
    for relative, expected in payload["outputs"].items():
        path = directory / relative
        path.resolve().relative_to(directory.resolve())
        stat = path.stat()
        if stat.st_size != expected["size_bytes"]:
            raise ValueError(f"Prepared source file changed: {path}")
        if stat.st_mtime_ns != expected["mtime_ns"] and content_identity(path) != {
            k: expected[k] for k in ("size_bytes", "sha256")
        }:
            raise ValueError(f"Prepared source checksum mismatch: {path}")
    values = dict(payload["source"])
    for key in PATH_FIELDS:
        values[key] = directory / values[key]
        values[key].resolve().relative_to(directory.resolve())
    for key in ("requested_gene_ids", "target_gene_ids"):
        values[key] = frozenset(values[key])
    return SourceRun(
        run_dir=directory, prepared_cache_dir=directory / "derived",
        archive_uri=payload["archive_uri"], **values,
    )
