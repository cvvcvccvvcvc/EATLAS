from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any

from run_archiving import __version__
from run_archiving.rclone import RcloneClient, remote_join


SCHEMA_VERSION = 1
READ_SIZE = 8 * 1024 * 1024
RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
REQUIRED_MANIFESTS = (
    "fetch/manifest.json",
    "alignment/manifest.json",
    "annotation/manifest.json",
)


class ArchiveError(RuntimeError):
    """Raised when a run cannot be archived or verified safely."""


@dataclass(frozen=True)
class FileRecord:
    path: str
    size: int
    mtime_ns: int
    md5: str
    sha256: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "size": self.size,
            "mtime_ns": self.mtime_ns,
            "md5": self.md5,
            "sha256": self.sha256,
        }


@dataclass(frozen=True)
class RunSnapshot:
    run_dir: Path
    run_id: str
    files: tuple[FileRecord, ...]
    total_bytes: int
    tree_sha256: str
    legacy_run: bool

    @property
    def file_count(self) -> int:
        return len(self.files)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(READ_SIZE):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_run_id(run_id: str) -> str:
    if not RUN_ID_PATTERN.fullmatch(run_id):
        raise ArchiveError(f"Invalid run ID: {run_id!r}")
    return run_id


def _validate_run_dir(
    run_dir: Path, *, allow_legacy_run: bool
) -> tuple[Path, str, bool]:
    if run_dir.is_symlink():
        raise ArchiveError(f"Run directory must not be a symlink: {run_dir}")
    resolved = run_dir.expanduser().resolve(strict=True)
    if not resolved.is_dir():
        raise ArchiveError(f"Run path is not a directory: {resolved}")
    run_id = _validate_run_id(resolved.name)
    run_manifest_path = resolved / "run_manifest.json"
    legacy_run = not run_manifest_path.is_file()
    if legacy_run:
        if not allow_legacy_run:
            raise ArchiveError(
                "run_manifest.json is missing. New runs require a successful root "
                "run manifest; use --allow-legacy-run only for historical runs "
                "created before that contract existed."
            )
    else:
        try:
            run_manifest = json.loads(run_manifest_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise ArchiveError(f"Invalid JSON manifest: {run_manifest_path}") from exc
        if (
            run_manifest.get("schema_version") != 2
            or run_manifest.get("status") != "complete"
            or run_manifest.get("success") is not True
            or run_manifest.get("exit_status") != 0
            or not run_manifest.get("completed_at")
        ):
            raise ArchiveError(
                f"Run is not marked successfully complete: {run_manifest_path}"
            )
    for relative in REQUIRED_MANIFESTS:
        path = resolved / relative
        if not path.is_file():
            raise ArchiveError(f"Required completed-run artifact is missing: {path}")
        try:
            json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise ArchiveError(f"Invalid JSON manifest: {path}") from exc
    for forbidden in ("work", ".nextflow"):
        if (resolved / forbidden).exists():
            raise ArchiveError(
                f"Execution cache must not be inside an archived run: "
                f"{resolved / forbidden}"
            )
    return resolved, run_id, legacy_run


def _hash_file(path: Path) -> tuple[int, int, str, str]:
    before = path.stat(follow_symlinks=False)
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(READ_SIZE):
            md5.update(chunk)
            sha256.update(chunk)
    after = path.stat(follow_symlinks=False)
    if (
        before.st_size != after.st_size
        or before.st_mtime_ns != after.st_mtime_ns
        or before.st_ino != after.st_ino
    ):
        raise ArchiveError(f"File changed while it was being hashed: {path}")
    return after.st_size, after.st_mtime_ns, md5.hexdigest(), sha256.hexdigest()


def build_snapshot(
    run_dir: Path, *, allow_legacy_run: bool = False
) -> RunSnapshot:
    resolved, run_id, legacy_run = _validate_run_dir(
        run_dir, allow_legacy_run=allow_legacy_run
    )
    paths: list[Path] = []
    for root, directory_names, file_names in os.walk(
        resolved, topdown=True, followlinks=False
    ):
        root_path = Path(root)
        directory_names.sort()
        file_names.sort()
        for directory_name in directory_names:
            directory_path = root_path / directory_name
            if directory_path.is_symlink():
                raise ArchiveError(
                    f"Symlinks are not supported in archived runs: {directory_path}"
                )
        for file_name in file_names:
            path = root_path / file_name
            if path.is_symlink():
                raise ArchiveError(
                    f"Symlinks are not supported in archived runs: {path}"
                )
            if not path.is_file():
                raise ArchiveError(
                    f"Only regular files are supported in archived runs: {path}"
                )
            relative = path.relative_to(resolved).as_posix()
            if "\n" in relative or "\r" in relative:
                raise ArchiveError(
                    f"Newlines are not supported in archived file names: {relative!r}"
                )
            paths.append(path)

    records: list[FileRecord] = []
    tree_digest = hashlib.sha256()
    total_bytes = 0
    for path in sorted(paths, key=lambda item: item.relative_to(resolved).as_posix()):
        relative = path.relative_to(resolved).as_posix()
        size, mtime_ns, md5, sha256 = _hash_file(path)
        record = FileRecord(
            path=relative,
            size=size,
            mtime_ns=mtime_ns,
            md5=md5,
            sha256=sha256,
        )
        records.append(record)
        total_bytes += size
        tree_digest.update(relative.encode("utf-8"))
        tree_digest.update(b"\0")
        tree_digest.update(str(size).encode("ascii"))
        tree_digest.update(b"\0")
        tree_digest.update(sha256.encode("ascii"))
        tree_digest.update(b"\n")

    if not records:
        raise ArchiveError(f"Run directory contains no files: {resolved}")
    return RunSnapshot(
        run_dir=resolved,
        run_id=run_id,
        files=tuple(records),
        total_bytes=total_bytes,
        tree_sha256=tree_digest.hexdigest(),
        legacy_run=legacy_run,
    )


def assert_snapshot_unchanged(snapshot: RunSnapshot) -> None:
    expected = {
        record.path: (record.size, record.mtime_ns) for record in snapshot.files
    }
    observed: dict[str, tuple[int, int]] = {}
    for root, directory_names, file_names in os.walk(
        snapshot.run_dir, topdown=True, followlinks=False
    ):
        root_path = Path(root)
        for directory_name in directory_names:
            if (root_path / directory_name).is_symlink():
                raise ArchiveError("Run gained a symlink while it was archived.")
        for file_name in file_names:
            path = root_path / file_name
            if path.is_symlink() or not path.is_file():
                raise ArchiveError("Run gained a non-regular file while it was archived.")
            relative = path.relative_to(snapshot.run_dir).as_posix()
            stat_result = path.stat(follow_symlinks=False)
            observed[relative] = (stat_result.st_size, stat_result.st_mtime_ns)
    if observed != expected:
        raise ArchiveError(
            "Run contents changed after the archive snapshot was created; "
            "the remote archive was not marked complete."
        )


def _manifest_payload(snapshot: RunSnapshot) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "archiver_version": __version__,
        "run_id": snapshot.run_id,
        "file_count": snapshot.file_count,
        "total_bytes": snapshot.total_bytes,
        "tree_sha256": snapshot.tree_sha256,
        "legacy_run": snapshot.legacy_run,
        "files": [record.as_dict() for record in snapshot.files],
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_control_files(snapshot: RunSnapshot, directory: Path) -> dict[str, Path]:
    directory.mkdir(parents=True, exist_ok=True)
    manifest = directory / "manifest.json"
    md5sums = directory / "MD5SUMS"
    sha256sums = directory / "SHA256SUMS"
    _write_json(manifest, _manifest_payload(snapshot))
    md5sums.write_text(
        "".join(f"{record.md5}  {record.path}\n" for record in snapshot.files)
    )
    sha256sums.write_text(
        "".join(f"{record.sha256}  {record.path}\n" for record in snapshot.files)
    )
    return {
        "manifest": manifest,
        "md5sums": md5sums,
        "sha256sums": sha256sums,
    }


def _validate_manifest(payload: dict[str, Any], run_id: str) -> None:
    if payload.get("schema_version") != SCHEMA_VERSION:
        raise ArchiveError(
            f"Unsupported archive schema: {payload.get('schema_version')!r}"
        )
    if payload.get("run_id") != run_id:
        raise ArchiveError(
            f"Remote manifest run ID is {payload.get('run_id')!r}, expected {run_id!r}"
        )
    if not isinstance(payload.get("legacy_run"), bool):
        raise ArchiveError("Remote archive manifest has no valid legacy-run flag.")
    files = payload.get("files")
    if not isinstance(files, list) or not files:
        raise ArchiveError("Remote archive manifest has no file inventory.")
    paths: set[str] = set()
    total_bytes = 0
    tree_digest = hashlib.sha256()
    for item in files:
        if not isinstance(item, dict):
            raise ArchiveError("Remote archive manifest contains an invalid file row.")
        path = item.get("path")
        if (
            not isinstance(path, str)
            or not path
            or path.startswith("/")
            or ".." in Path(path).parts
            or path in paths
        ):
            raise ArchiveError(f"Unsafe or duplicate archived path: {path!r}")
        paths.add(path)
        try:
            size = int(item["size"])
            md5 = str(item["md5"])
            sha256 = str(item["sha256"])
            int(item["mtime_ns"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ArchiveError(f"Invalid inventory values for {path!r}") from exc
        if size < 0:
            raise ArchiveError(f"Negative file size in inventory for {path!r}")
        if not re.fullmatch(r"[0-9a-f]{32}", md5):
            raise ArchiveError(f"Invalid MD5 in inventory for {path!r}")
        if not re.fullmatch(r"[0-9a-f]{64}", sha256):
            raise ArchiveError(f"Invalid SHA-256 in inventory for {path!r}")
        total_bytes += size
        tree_digest.update(path.encode("utf-8"))
        tree_digest.update(b"\0")
        tree_digest.update(str(size).encode("ascii"))
        tree_digest.update(b"\0")
        tree_digest.update(sha256.encode("ascii"))
        tree_digest.update(b"\n")
    if len(files) != int(payload.get("file_count", -1)):
        raise ArchiveError("Remote archive file count does not match its inventory.")
    if total_bytes != int(payload.get("total_bytes", -1)):
        raise ArchiveError("Remote archive byte count does not match its inventory.")
    if tree_digest.hexdigest() != payload.get("tree_sha256"):
        raise ArchiveError("Remote archive tree checksum does not match its inventory.")


def _checksum_file_from_manifest(
    manifest: dict[str, Any], algorithm: str, destination: Path
) -> None:
    destination.write_text(
        "".join(
            f"{item[algorithm]}  {item['path']}\n" for item in manifest["files"]
        )
    )


def _archive_paths(remote_root: str, run_id: str) -> dict[str, str]:
    _validate_run_id(run_id)
    base = remote_join(remote_root, "runs", run_id)
    return {
        "base": base,
        "data": remote_join(base, "data"),
        "control": remote_join(base, "_archive"),
        "manifest": remote_join(base, "_archive", "manifest.json"),
        "complete": remote_join(base, "_archive", "COMPLETE.json"),
    }


def _parse_completion_marker(text: str, run_id: str) -> dict[str, Any]:
    try:
        complete = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ArchiveError("Remote completion marker is invalid JSON.") from exc
    if not isinstance(complete, dict):
        raise ArchiveError("Remote completion marker must be a JSON object.")
    if complete.get("schema_version") != SCHEMA_VERSION:
        raise ArchiveError("Remote completion marker has an unsupported schema.")
    if complete.get("run_id") != run_id:
        raise ArchiveError("Remote completion marker has a different run ID.")
    if (
        not isinstance(complete.get("completed_at"), str)
        or not complete["completed_at"]
    ):
        raise ArchiveError("Remote completion marker has no completion time.")
    for field in ("file_count", "total_bytes"):
        value = complete.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ArchiveError(
                f"Remote completion marker has an invalid {field.replace('_', ' ')}."
            )
    if not re.fullmatch(r"[0-9a-f]{64}", str(complete.get("tree_sha256", ""))):
        raise ArchiveError("Remote completion marker has an invalid tree checksum.")
    return complete


def list_archives(
    client: RcloneClient,
    *,
    remote_root: str,
) -> list[dict[str, Any]]:
    client.preflight(remote_root)
    marker_pattern = "runs/*/_archive/COMPLETE.json"
    marker_paths = client.list_files(remote_root, include=marker_pattern)
    archives: list[dict[str, Any]] = []
    seen_run_ids: set[str] = set()
    for marker_path in marker_paths:
        parts = PurePosixPath(marker_path).parts
        if (
            len(parts) != 4
            or parts[0] != "runs"
            or parts[2] != "_archive"
            or parts[3] != "COMPLETE.json"
        ):
            raise ArchiveError(f"Unexpected archive marker path: {marker_path!r}")
        run_id = _validate_run_id(parts[1])
        if run_id in seen_run_ids:
            raise ArchiveError(f"Duplicate remote archive marker for run {run_id!r}.")
        seen_run_ids.add(run_id)
        marker_text = client.read_text_optional(remote_join(remote_root, marker_path))
        if marker_text is None:
            raise ArchiveError(
                f"Remote archive marker disappeared while listing run {run_id!r}."
            )
        complete = _parse_completion_marker(marker_text, run_id)
        archives.append(
            {
                "run_id": run_id,
                "archived_at": complete["completed_at"],
                "file_count": complete["file_count"],
                "total_bytes": complete["total_bytes"],
                "tree_sha256": complete["tree_sha256"],
            }
        )
    archives.sort(
        key=lambda item: (str(item["archived_at"]), str(item["run_id"])),
        reverse=True,
    )
    return archives


def _load_remote_manifest(
    client: RcloneClient, remote_root: str, run_id: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    paths = _archive_paths(remote_root, run_id)
    complete_text = client.read_text_optional(paths["complete"])
    if complete_text is None:
        raise ArchiveError(f"Remote archive is not marked complete: {paths['base']}")
    complete = _parse_completion_marker(complete_text, run_id)
    with tempfile.TemporaryDirectory(prefix="run-archive-verify-") as temporary:
        manifest_path = Path(temporary) / "manifest.json"
        client.download_file(paths["manifest"], manifest_path)
        manifest_bytes = manifest_path.read_bytes()
        try:
            manifest = json.loads(manifest_bytes)
        except json.JSONDecodeError as exc:
            raise ArchiveError("Remote archive manifest is invalid JSON.") from exc
    _validate_manifest(manifest, run_id)
    if complete.get("file_count") != manifest.get("file_count"):
        raise ArchiveError("Remote completion marker has a different file count.")
    if complete.get("total_bytes") != manifest.get("total_bytes"):
        raise ArchiveError("Remote completion marker has a different byte count.")
    manifest_sha256 = hashlib.sha256(manifest_bytes).hexdigest()
    if complete.get("manifest_sha256") != manifest_sha256:
        raise ArchiveError("Remote archive manifest checksum does not match its marker.")
    if complete.get("tree_sha256") != manifest.get("tree_sha256"):
        raise ArchiveError("Remote archive tree checksum does not match its marker.")
    return manifest, complete


def verify_remote(
    client: RcloneClient,
    *,
    remote_root: str,
    run_id: str,
) -> dict[str, Any]:
    paths = _archive_paths(remote_root, run_id)
    manifest, _ = _load_remote_manifest(client, remote_root, run_id)
    with tempfile.TemporaryDirectory(prefix="run-archive-check-") as temporary:
        md5sums = Path(temporary) / "MD5SUMS"
        _checksum_file_from_manifest(manifest, "md5", md5sums)
        client.verify_checksum("md5", md5sums, paths["data"])
    remote_count, remote_bytes = client.size(paths["data"])
    if remote_count != int(manifest["file_count"]):
        raise ArchiveError(
            f"Remote file count is {remote_count}, expected {manifest['file_count']}."
        )
    if remote_bytes != int(manifest["total_bytes"]):
        raise ArchiveError(
            f"Remote byte count is {remote_bytes}, expected {manifest['total_bytes']}."
        )
    return manifest


def archive_run(
    client: RcloneClient,
    *,
    run_dir: Path,
    remote_root: str,
    dry_run: bool = False,
    allow_legacy_run: bool = False,
) -> dict[str, Any]:
    client.preflight(remote_root)
    snapshot = build_snapshot(run_dir, allow_legacy_run=allow_legacy_run)
    paths = _archive_paths(remote_root, snapshot.run_id)

    complete_text = client.read_text_optional(paths["complete"])
    if complete_text is not None:
        manifest = verify_remote(
            client, remote_root=remote_root, run_id=snapshot.run_id
        )
        if manifest["tree_sha256"] != snapshot.tree_sha256:
            raise ArchiveError(
                "A complete archive with this run ID already exists but contains "
                "different data."
            )
        return {
            "status": "already_archived",
            "run_id": snapshot.run_id,
            "remote": paths["base"],
            "file_count": snapshot.file_count,
            "total_bytes": snapshot.total_bytes,
            "tree_sha256": snapshot.tree_sha256,
            "legacy_run": snapshot.legacy_run,
        }

    with tempfile.TemporaryDirectory(prefix="run-archive-") as temporary:
        control_dir = Path(temporary) / "_archive"
        control_files = write_control_files(snapshot, control_dir)
        client.copy_tree(snapshot.run_dir, paths["data"], dry_run=dry_run)
        if dry_run:
            return {
                "status": "dry_run",
                "run_id": snapshot.run_id,
                "remote": paths["base"],
                "file_count": snapshot.file_count,
                "total_bytes": snapshot.total_bytes,
                "tree_sha256": snapshot.tree_sha256,
                "legacy_run": snapshot.legacy_run,
            }

        client.verify_checksum("md5", control_files["md5sums"], paths["data"])
        remote_count, remote_bytes = client.size(paths["data"])
        if remote_count != snapshot.file_count or remote_bytes != snapshot.total_bytes:
            raise ArchiveError(
                "Remote size does not match the local snapshot: "
                f"{remote_count} files/{remote_bytes} bytes remote, "
                f"{snapshot.file_count} files/{snapshot.total_bytes} bytes local."
            )
        assert_snapshot_unchanged(snapshot)

        for name in ("manifest", "md5sums", "sha256sums"):
            client.copy_file(
                control_files[name],
                remote_join(paths["control"], control_files[name].name),
            )

        manifest_sha256 = sha256_file(control_files["manifest"])
        complete_path = control_dir / "COMPLETE.json"
        _write_json(
            complete_path,
            {
                "schema_version": SCHEMA_VERSION,
                "run_id": snapshot.run_id,
                "completed_at": utc_now(),
                "manifest_sha256": manifest_sha256,
                "tree_sha256": snapshot.tree_sha256,
                "file_count": snapshot.file_count,
                "total_bytes": snapshot.total_bytes,
            },
        )
        client.copy_file(complete_path, paths["complete"])

    verify_remote(client, remote_root=remote_root, run_id=snapshot.run_id)
    return {
        "status": "archived",
        "run_id": snapshot.run_id,
        "remote": paths["base"],
        "file_count": snapshot.file_count,
        "total_bytes": snapshot.total_bytes,
        "tree_sha256": snapshot.tree_sha256,
        "legacy_run": snapshot.legacy_run,
    }


def restore_run(
    client: RcloneClient,
    *,
    remote_root: str,
    run_id: str,
    destination: Path,
) -> dict[str, Any]:
    _validate_run_id(run_id)
    manifest = verify_remote(client, remote_root=remote_root, run_id=run_id)
    resolved_destination = destination.expanduser().resolve()
    if resolved_destination.exists():
        if not resolved_destination.is_dir() or any(resolved_destination.iterdir()):
            raise ArchiveError(
                f"Restore destination must be absent or empty: {resolved_destination}"
            )
        resolved_destination.rmdir()
    partial_destination = resolved_destination.with_name(
        f"{resolved_destination.name}.partial"
    )
    if partial_destination.exists() and not partial_destination.is_dir():
        raise ArchiveError(
            f"Partial restore path is not a directory: {partial_destination}"
        )
    partial_destination.mkdir(parents=True, exist_ok=True)
    client.copy_tree(
        _archive_paths(remote_root, run_id)["data"],
        partial_destination,
    )
    restored = build_snapshot(
        partial_destination,
        allow_legacy_run=bool(manifest["legacy_run"]),
    )
    if restored.tree_sha256 != manifest["tree_sha256"]:
        raise ArchiveError("Restored run does not match the archived tree checksum.")
    if resolved_destination.exists():
        raise ArchiveError(
            f"Restore destination appeared during verification: {resolved_destination}"
        )
    partial_destination.rename(resolved_destination)
    return {
        "status": "restored",
        "run_id": run_id,
        "destination": str(resolved_destination),
        "file_count": restored.file_count,
        "total_bytes": restored.total_bytes,
        "tree_sha256": restored.tree_sha256,
    }


def remove_local_run(
    client: RcloneClient,
    *,
    run_dir: Path,
    results_root: Path,
    remote_root: str,
    confirmation: str,
) -> dict[str, Any]:
    if run_dir.is_symlink():
        raise ArchiveError(f"Run directory must not be a symlink: {run_dir}")
    resolved_run = run_dir.expanduser().resolve(strict=True)
    resolved_root = results_root.expanduser().resolve(strict=True)
    if resolved_run.parent != resolved_root:
        raise ArchiveError(
            f"Refusing to remove a path outside the direct results root: "
            f"{resolved_run} is not a direct child of {resolved_root}"
        )
    if confirmation != resolved_run.name:
        raise ArchiveError(
            f"Confirmation must exactly match the run ID: {resolved_run.name}"
        )

    run_id = resolved_run.name
    remote_manifest = verify_remote(
        client, remote_root=remote_root, run_id=run_id
    )
    quarantine = resolved_root / f"{run_id}.removing"
    if quarantine.exists():
        raise ArchiveError(f"Removal quarantine already exists: {quarantine}")
    resolved_run.rename(quarantine)
    try:
        local_snapshot = build_snapshot(
            quarantine,
            allow_legacy_run=bool(remote_manifest["legacy_run"]),
        )
        if local_snapshot.tree_sha256 != remote_manifest["tree_sha256"]:
            raise ArchiveError(
                "Local run differs from the verified remote archive; "
                "nothing was removed."
            )
        if resolved_run.exists():
            raise ArchiveError(
                "The original run path reappeared during removal; "
                "the quarantined run was not removed."
            )
        assert_snapshot_unchanged(local_snapshot)
        removed_bytes = local_snapshot.total_bytes
        removed_files = local_snapshot.file_count
        shutil.rmtree(quarantine)
    except Exception:
        if quarantine.exists() and not resolved_run.exists():
            quarantine.rename(resolved_run)
        raise
    return {
        "status": "removed_local",
        "run_id": run_id,
        "removed_files": removed_files,
        "removed_bytes": removed_bytes,
    }
