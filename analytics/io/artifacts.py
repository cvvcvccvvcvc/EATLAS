"""Atomic artifact writes and stable local file identities."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd


def sha256_file(path: Path, *, chunk_size: int = 16 * 1024 * 1024) -> str:
    """Return a content identity suitable for reproducible analysis inputs."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def content_identity(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    before = path.stat()
    sha256 = sha256_file(path)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"File changed while hashing: {path}")
    return {
        "size_bytes": before.st_size,
        "sha256": sha256,
    }


def cached_content_identity(path: Path, *, cache_dir: Path) -> dict[str, object]:
    """Return a content hash, reusing it only while the same file is unchanged."""

    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    signature = _stat_signature(resolved)
    marker_name = hashlib.sha256(str(resolved).encode("utf-8")).hexdigest() + ".json"
    marker_path = cache_dir / marker_name
    if marker_path.is_file():
        try:
            marker = json.loads(marker_path.read_text())
        except (OSError, json.JSONDecodeError):
            marker = None
        if (
            isinstance(marker, dict)
            and marker.get("schema_version") == 1
            and marker.get("path") == str(resolved)
            and marker.get("stat") == signature
            and _valid_content_identity(marker.get("content"), signature["size_bytes"])
        ):
            return dict(marker["content"])

    identity = content_identity(resolved)
    after = _stat_signature(resolved)
    if signature != after:
        raise ValueError(f"File changed while hashing: {resolved}")
    write_json_atomic(
        marker_path,
        {
            "schema_version": 1,
            "path": str(resolved),
            "stat": after,
            "content": identity,
        },
    )
    return identity


def _stat_signature(path: Path) -> dict[str, int]:
    observed = path.stat()
    return {
        "device": observed.st_dev,
        "inode": observed.st_ino,
        "size_bytes": observed.st_size,
        "mtime_ns": observed.st_mtime_ns,
        "ctime_ns": observed.st_ctime_ns,
    }


def _valid_content_identity(value: object, size_bytes: int) -> bool:
    digest = value.get("sha256") if isinstance(value, dict) else None
    return (
        isinstance(value, dict)
        and value.get("size_bytes") == size_bytes
        and isinstance(digest, str)
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
    )


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
        temporary_path.chmod(0o644)
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
