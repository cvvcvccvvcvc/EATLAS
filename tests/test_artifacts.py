from __future__ import annotations

from pathlib import Path

import pandas as pd

from analytics.io.artifacts import (
    directory_metadata,
    path_metadata,
    write_json_atomic,
    write_tsv_atomic,
)


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
