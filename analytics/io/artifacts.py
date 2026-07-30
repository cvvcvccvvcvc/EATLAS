"""Atomic artifact writes and stable local file identities."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd


def file_identity(path: Path) -> dict[str, int]:
    stat = path.stat()
    return {"size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def path_metadata(path: Path) -> dict[str, object]:
    if not path.exists():
        raise FileNotFoundError(path)
    return {"path": str(path.resolve()), **file_identity(path)}


def directory_metadata(path: Path, pattern: str = "*") -> dict[str, object]:
    if not path.exists():
        raise FileNotFoundError(path)
    files = sorted(item for item in path.glob(pattern) if item.is_file())
    return {
        "path": str(path.resolve()),
        "file_count": len(files),
        "files": [
            {"path": str(item.relative_to(path)), **file_identity(item)}
            for item in files
        ],
    }


def write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        temporary_path.replace(path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def write_json_atomic(path: Path, payload: Any) -> None:
    write_text_atomic(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_tsv_atomic(
    path: Path,
    frame: pd.DataFrame,
    *,
    header: bool = True,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    compression: str | dict[str, object] | None = None
    if str(path).endswith(".gz"):
        compression = {"method": "gzip", "compresslevel": 6, "mtime": 0}
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
        frame.to_csv(
            temporary_path,
            sep="\t",
            index=False,
            header=header,
            compression=compression,
            lineterminator="\n",
        )
        temporary_path.chmod(0o644)
        temporary_path.replace(path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
