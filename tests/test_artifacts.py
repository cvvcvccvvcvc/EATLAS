from __future__ import annotations

import stat
from pathlib import Path

import pandas as pd

from analytics.io.artifacts import (
    cached_content_identity,
    directory_metadata,
    path_metadata,
    write_json_atomic,
    write_tsv_atomic,
)
from analytics.io import artifacts


def test_artifact_helpers_write_and_identify_outputs(tmp_path: Path) -> None:
    table_path = tmp_path / "table.tsv.gz"
    manifest_path = tmp_path / "manifest.json"
    sequence_dir = tmp_path / "sequences"
    sequence_dir.mkdir()
    (sequence_dir / "1.fa.gz").write_bytes(b"one")
    (sequence_dir / "ignored.txt").write_text("ignored")

    write_tsv_atomic(table_path, pd.DataFrame({"value": [1, 2]}))
    write_json_atomic(manifest_path, {"output": path_metadata(table_path)})

    observed = pd.read_csv(table_path, sep="\t", compression="gzip")
    metadata = directory_metadata(sequence_dir, "*.fa.gz")
    assert observed["value"].tolist() == [1, 2]
    assert metadata["file_count"] == 1
    assert metadata["files"][0]["path"] == "1.fa.gz"
    assert path_metadata(table_path)["mtime_ns"] > 0
    assert stat.S_IMODE(table_path.stat().st_mode) == 0o644
    assert stat.S_IMODE(manifest_path.stat().st_mode) == 0o644


def test_cached_content_identity_rehashes_only_after_file_change(
    tmp_path: Path,
    monkeypatch,
) -> None:
    source = tmp_path / "reference.bw"
    source.write_bytes(b"first")
    calls = 0
    original = artifacts.sha256_file

    def count_hashes(path: Path, *, chunk_size: int = 16 * 1024 * 1024) -> str:
        nonlocal calls
        calls += 1
        return original(path, chunk_size=chunk_size)

    monkeypatch.setattr(artifacts, "sha256_file", count_hashes)
    cache_dir = tmp_path / "cache"
    first = cached_content_identity(source, cache_dir=cache_dir)
    assert cached_content_identity(source, cache_dir=cache_dir) == first
    assert calls == 1

    source.write_bytes(b"other")
    second = cached_content_identity(source, cache_dir=cache_dir)
    assert second != first
    assert calls == 2
