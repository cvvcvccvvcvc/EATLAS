from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pytest

from run_archiving.archive import (
    ArchiveError,
    archive_run,
    build_snapshot,
    remove_local_run,
    restore_run,
    verify_remote,
)


def _digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


class LocalRemote:
    def __init__(self, root: Path) -> None:
        self.root = root

    def _path(self, value: str) -> Path:
        _, separator, relative = value.partition(":")
        if not separator:
            raise ValueError(f"Expected fake remote path: {value}")
        return self.root / relative.lstrip("/")

    def preflight(self, remote_root: str) -> None:
        self._path(remote_root).mkdir(parents=True, exist_ok=True)

    def read_text_optional(self, remote_path: str) -> str | None:
        path = self._path(remote_path)
        return path.read_text() if path.is_file() else None

    @staticmethod
    def _copy_contents(source: Path, destination: Path) -> None:
        destination.mkdir(parents=True, exist_ok=True)
        for path in sorted(source.rglob("*")):
            if path.is_dir():
                continue
            relative = path.relative_to(source)
            target = destination / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists():
                if _digest(path, "sha256") != _digest(target, "sha256"):
                    raise RuntimeError(f"immutable destination differs: {target}")
                continue
            shutil.copy2(path, target)

    def copy_tree(
        self, source: str | Path, destination: str | Path, *, dry_run: bool = False
    ) -> None:
        if dry_run:
            return
        source_path = (
            self._path(source) if isinstance(source, str) and ":" in source else Path(source)
        )
        destination_path = (
            self._path(destination)
            if isinstance(destination, str) and ":" in destination
            else Path(destination)
        )
        self._copy_contents(source_path, destination_path)

    def copy_file(self, source: str | Path, destination: str) -> None:
        source_path = Path(source)
        destination_path = self._path(destination)
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        if destination_path.exists():
            if _digest(source_path, "sha256") != _digest(destination_path, "sha256"):
                raise RuntimeError(f"immutable destination differs: {destination_path}")
            return
        shutil.copy2(source_path, destination_path)

    def download_file(self, source: str, destination: Path) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(self._path(source), destination)

    def verify_checksum(
        self, algorithm: str, checksum_file: Path, destination: str
    ) -> None:
        root = self._path(destination)
        for line in checksum_file.read_text().splitlines():
            expected, relative = line.split("  ", 1)
            assert _digest(root / relative, algorithm) == expected

    def size(self, remote_path: str) -> tuple[int, int]:
        files = [path for path in self._path(remote_path).rglob("*") if path.is_file()]
        return len(files), sum(path.stat().st_size for path in files)


def _make_run(results_root: Path, run_id: str = "run_001") -> Path:
    run_dir = results_root / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "run_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "status": "complete",
                "success": True,
                "exit_status": 0,
                "completed_at": "2026-07-30T00:00:00Z",
            }
        )
        + "\n"
    )
    for stage in ("fetch", "alignment", "annotation"):
        stage_dir = run_dir / stage
        stage_dir.mkdir(parents=True, exist_ok=True)
        (stage_dir / "manifest.json").write_text(
            json.dumps({"stage": stage, "failure_count": 0}) + "\n"
        )
    (run_dir / "fetch" / "genes.tsv.gz").write_bytes(b"genes")
    (run_dir / "alignment" / "strategy_summary.tsv.gz").write_bytes(b"summary")
    (run_dir / "annotation" / "variant_annotations.tsv.gz").write_bytes(
        b"variants"
    )
    return run_dir


def test_archive_verify_and_restore_round_trip(tmp_path: Path) -> None:
    results_root = tmp_path / "results"
    run_dir = _make_run(results_root)
    client = LocalRemote(tmp_path / "remote")

    archived = archive_run(
        client, run_dir=run_dir, remote_root="drive:GAPH", dry_run=False
    )
    assert archived["status"] == "archived"
    assert (
        tmp_path
        / "remote"
        / "GAPH"
        / "runs"
        / "run_001"
        / "_archive"
        / "COMPLETE.json"
    ).is_file()

    verified = verify_remote(client, remote_root="drive:GAPH", run_id="run_001")
    assert verified["tree_sha256"] == build_snapshot(run_dir).tree_sha256

    repeated = archive_run(
        client, run_dir=run_dir, remote_root="drive:GAPH", dry_run=False
    )
    assert repeated["status"] == "already_archived"

    destination = tmp_path / "restored"
    restored = restore_run(
        client,
        remote_root="drive:GAPH",
        run_id="run_001",
        destination=destination,
    )
    assert restored["status"] == "restored"
    assert build_snapshot(destination).tree_sha256 == build_snapshot(run_dir).tree_sha256


def test_archive_requires_successful_root_run_manifest(tmp_path: Path) -> None:
    run_dir = _make_run(tmp_path / "results")
    client = LocalRemote(tmp_path / "remote")
    (run_dir / "run_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "status": "failed",
                "success": False,
                "exit_status": 1,
                "completed_at": "2026-07-30T00:00:00Z",
            }
        )
        + "\n"
    )

    with pytest.raises(ArchiveError, match="not marked successfully complete"):
        archive_run(client, run_dir=run_dir, remote_root="drive:GAPH")


def test_archive_legacy_run_requires_explicit_exception(tmp_path: Path) -> None:
    run_dir = _make_run(tmp_path / "results")
    client = LocalRemote(tmp_path / "remote")
    (run_dir / "run_manifest.json").unlink()

    with pytest.raises(ArchiveError, match="allow-legacy-run"):
        archive_run(client, run_dir=run_dir, remote_root="drive:GAPH")

    result = archive_run(
        client,
        run_dir=run_dir,
        remote_root="drive:GAPH",
        allow_legacy_run=True,
    )
    assert result["status"] == "archived"
    assert result["legacy_run"] is True


def test_restore_resumes_matching_partial_directory(tmp_path: Path) -> None:
    run_dir = _make_run(tmp_path / "results")
    client = LocalRemote(tmp_path / "remote")
    archive_run(client, run_dir=run_dir, remote_root="drive:GAPH")
    destination = tmp_path / "restored"
    partial = tmp_path / "restored.partial"
    (partial / "fetch").mkdir(parents=True)
    shutil.copy2(
        run_dir / "fetch" / "manifest.json",
        partial / "fetch" / "manifest.json",
    )

    result = restore_run(
        client,
        remote_root="drive:GAPH",
        run_id="run_001",
        destination=destination,
    )

    assert result["status"] == "restored"
    assert destination.is_dir()
    assert not partial.exists()


def test_complete_archive_refuses_changed_local_run(tmp_path: Path) -> None:
    run_dir = _make_run(tmp_path / "results")
    client = LocalRemote(tmp_path / "remote")
    archive_run(client, run_dir=run_dir, remote_root="drive:GAPH")
    (run_dir / "annotation" / "variant_annotations.tsv.gz").write_bytes(b"changed")

    with pytest.raises(ArchiveError, match="different data"):
        archive_run(client, run_dir=run_dir, remote_root="drive:GAPH")


def test_remove_local_requires_verified_archive_and_exact_confirmation(
    tmp_path: Path,
) -> None:
    results_root = tmp_path / "results"
    run_dir = _make_run(results_root)
    client = LocalRemote(tmp_path / "remote")
    archive_run(client, run_dir=run_dir, remote_root="drive:GAPH")

    with pytest.raises(ArchiveError, match="Confirmation"):
        remove_local_run(
            client,
            run_dir=run_dir,
            results_root=results_root,
            remote_root="drive:GAPH",
            confirmation="wrong-run",
        )
    assert run_dir.is_dir()

    result = remove_local_run(
        client,
        run_dir=run_dir,
        results_root=results_root,
        remote_root="drive:GAPH",
        confirmation="run_001",
    )
    assert result["status"] == "removed_local"
    assert not run_dir.exists()


def test_remove_local_restores_name_when_local_data_changed(tmp_path: Path) -> None:
    results_root = tmp_path / "results"
    run_dir = _make_run(results_root)
    client = LocalRemote(tmp_path / "remote")
    archive_run(client, run_dir=run_dir, remote_root="drive:GAPH")
    (run_dir / "annotation" / "variant_annotations.tsv.gz").write_bytes(b"changed")

    with pytest.raises(ArchiveError, match="differs"):
        remove_local_run(
            client,
            run_dir=run_dir,
            results_root=results_root,
            remote_root="drive:GAPH",
            confirmation="run_001",
        )

    assert run_dir.is_dir()
    assert not (results_root / "run_001.removing").exists()


def test_remove_local_rejects_path_outside_results_root(tmp_path: Path) -> None:
    expected_root = tmp_path / "results"
    expected_root.mkdir()
    run_dir = _make_run(tmp_path / "other")
    client = LocalRemote(tmp_path / "remote")

    with pytest.raises(ArchiveError, match="outside the direct results root"):
        remove_local_run(
            client,
            run_dir=run_dir,
            results_root=expected_root,
            remote_root="drive:GAPH",
            confirmation="run_001",
        )


def test_snapshot_rejects_symlinks(tmp_path: Path) -> None:
    run_dir = _make_run(tmp_path / "results")
    target = run_dir / "fetch" / "genes.tsv.gz"
    (run_dir / "fetch" / "linked.tsv.gz").symlink_to(target)

    with pytest.raises(ArchiveError, match="Symlinks"):
        build_snapshot(run_dir)


def test_verify_rejects_path_traversal_run_id(tmp_path: Path) -> None:
    client = LocalRemote(tmp_path / "remote")

    with pytest.raises(ArchiveError, match="Invalid run ID"):
        verify_remote(client, remote_root="drive:GAPH", run_id="../run_001")
