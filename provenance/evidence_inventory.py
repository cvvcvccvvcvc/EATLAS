"""Content inventory for immutable pipeline evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path, PurePosixPath


INVENTORY_FILENAME = "evidence_inventory.json"
INVENTORY_SCHEMA_VERSION = 1
EVIDENCE_SCOPES = ("fetch", "alignment", "annotation")
READ_SIZE = 16 * 1024 * 1024
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


class EvidenceInventoryError(ValueError):
    """Raised when evidence does not satisfy its inventory contract."""


def build_evidence_inventory(
    inputs: Mapping[str, Path],
    *,
    scopes: tuple[str, ...] = EVIDENCE_SCOPES,
) -> dict[str, object]:
    """Hash logical files and directories into one deterministic inventory."""

    _validate_scopes(scopes)
    paths: dict[str, Path] = {}
    for logical_prefix, source_path in sorted(inputs.items()):
        prefix = _safe_relative_path(logical_prefix)
        source = source_path.expanduser().absolute()
        if source.is_symlink():
            raise EvidenceInventoryError(f"Evidence input is a symlink: {source}")
        if source.is_file():
            members = [(prefix, source)]
        elif source.is_dir():
            members = [
                (f"{prefix}/{relative}", path)
                for relative, path in _regular_tree_files(source)
            ]
        else:
            raise EvidenceInventoryError(
                f"Evidence input is not a regular file or directory: {source}"
            )

        for logical, path in members:
            if logical in paths:
                raise EvidenceInventoryError(f"Duplicate logical evidence path: {logical}")
            paths[logical] = path
    before = {logical: _file_stat(path) for logical, path in paths.items()}
    rows = [
        {"path": logical, "size_bytes": before[logical][2], "sha256": _sha256_file(path)}
        for logical, path in sorted(paths.items())
    ]
    if before != {logical: _file_stat(path) for logical, path in paths.items()}:
        raise EvidenceInventoryError("Evidence files changed while hashing")
    if sorted(paths) != _logical_input_paths(list(inputs.items())):
        raise EvidenceInventoryError("Evidence file set changed while hashing")
    payload = _inventory_payload(rows, scopes=scopes)
    validate_evidence_inventory(payload, scopes=scopes)
    return payload


def build_run_evidence_inventory(run_dir: Path) -> dict[str, object]:
    """Hash the three durable evidence roots of one run directory."""

    root = run_dir.expanduser().resolve(strict=True)
    return build_evidence_inventory(
        {scope: root / scope for scope in EVIDENCE_SCOPES}
    )


def write_evidence_inventory(
    output: Path,
    inputs: Mapping[str, Path],
    *,
    scopes: tuple[str, ...] = EVIDENCE_SCOPES,
) -> dict[str, object]:
    payload = build_evidence_inventory(inputs, scopes=scopes)
    _write_json_atomic(output, payload)
    return payload


def load_evidence_inventory(
    path: Path,
    *,
    scopes: tuple[str, ...] = EVIDENCE_SCOPES,
) -> dict[str, object]:
    content = _read_stable_regular_file(path, allow_staged_symlink=True)
    return _parse_evidence_inventory(content, path, scopes)


def load_bound_evidence_inventory(
    path: Path,
    descriptor: object,
    *,
    scopes: tuple[str, ...] = EVIDENCE_SCOPES,
) -> tuple[dict[str, object], dict[str, object]]:
    """Load and validate the exact inventory bytes bound by a run manifest."""

    content = _read_stable_regular_file(path)
    observed = _inventory_file_descriptor(content)
    if descriptor != observed:
        raise EvidenceInventoryError(
            f"Run manifest evidence inventory descriptor does not match {path}"
        )
    return _parse_evidence_inventory(content, path, scopes), observed


def _parse_evidence_inventory(
    content: bytes,
    path: Path,
    scopes: tuple[str, ...],
) -> dict[str, object]:
    try:
        payload = json.loads(content)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise EvidenceInventoryError(f"Invalid evidence inventory: {path}") from exc
    validate_evidence_inventory(payload, scopes=scopes)
    return payload


def validate_evidence_inventory(
    payload: object,
    *,
    scopes: tuple[str, ...] = EVIDENCE_SCOPES,
) -> None:
    _validate_scopes(scopes)
    required = {
        "schema_version",
        "scope",
        "file_count",
        "total_bytes",
        "tree_sha256",
        "files",
    }
    if not isinstance(payload, dict) or set(payload) != required:
        raise EvidenceInventoryError("Evidence inventory has an invalid top-level schema")
    if type(payload["schema_version"]) is not int or payload["schema_version"] != INVENTORY_SCHEMA_VERSION:
        raise EvidenceInventoryError(
            f"Unsupported evidence inventory schema: {payload['schema_version']!r}"
        )
    if payload["scope"] != list(scopes):
        raise EvidenceInventoryError("Evidence inventory has an invalid scope")
    raw_files = payload["files"]
    if not isinstance(raw_files, list) or not raw_files:
        raise EvidenceInventoryError("Evidence inventory contains no files")

    normalized = []
    previous_path = ""
    for raw in raw_files:
        if not isinstance(raw, dict) or set(raw) != {
            "path",
            "size_bytes",
            "sha256",
        }:
            raise EvidenceInventoryError("Evidence inventory has an invalid file row")
        path = _safe_relative_path(raw["path"])
        if len(PurePosixPath(path).parts) < 2:
            raise EvidenceInventoryError(f"Evidence file must be below a stage directory: {path}")
        size = raw["size_bytes"]
        sha256 = raw["sha256"]
        if path <= previous_path:
            raise EvidenceInventoryError(
                "Evidence inventory paths must be unique and sorted"
            )
        if isinstance(size, bool) or not isinstance(size, int) or size < 0:
            raise EvidenceInventoryError(f"Invalid evidence file size: {path}")
        if not isinstance(sha256, str) or SHA256_PATTERN.fullmatch(sha256) is None:
            raise EvidenceInventoryError(f"Invalid evidence SHA-256: {path}")
        normalized.append({"path": path, "size_bytes": size, "sha256": sha256})
        previous_path = path

    observed_scopes = {PurePosixPath(row["path"]).parts[0] for row in normalized}
    if observed_scopes != set(scopes):
        raise EvidenceInventoryError("Evidence inventory file scopes differ from its declaration")
    expected = _inventory_payload(normalized, scopes=scopes)
    for field in ("file_count", "total_bytes", "tree_sha256"):
        if type(payload[field]) is not type(expected[field]) or payload[field] != expected[field]:
            raise EvidenceInventoryError(
                f"Evidence inventory {field} does not match its file rows"
            )


def combine_evidence_inventories(
    fragments: Mapping[str, dict[str, object]],
) -> dict[str, object]:
    """Combine exactly one independently validated inventory per evidence stage."""

    if set(fragments) != set(EVIDENCE_SCOPES):
        raise EvidenceInventoryError(
            "Evidence inventory fragments must cover fetch, alignment, and annotation"
        )
    rows = []
    for scope in EVIDENCE_SCOPES:
        fragment = fragments[scope]
        validate_evidence_inventory(fragment, scopes=(scope,))
        rows.extend(fragment["files"])
    payload = _inventory_payload(sorted(rows, key=lambda row: row["path"]))
    validate_evidence_inventory(payload)
    return payload


def inventory_file_descriptor(path: Path) -> dict[str, object]:
    """Describe the small inventory file stored beside the root manifest."""

    return _inventory_file_descriptor(_read_stable_regular_file(path))


def _inventory_file_descriptor(content: bytes) -> dict[str, object]:
    return {
        "path": INVENTORY_FILENAME,
        "schema_version": INVENTORY_SCHEMA_VERSION,
        "size_bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def verify_run_evidence(
    run_dir: Path,
    expected: dict[str, object],
) -> str:
    """Fully hash a run's evidence and return its current stat fingerprint."""

    validate_evidence_inventory(expected)
    before = evidence_stat_fingerprint(run_dir, expected)
    observed = build_run_evidence_inventory(run_dir)
    if observed != expected:
        _raise_inventory_difference(expected, observed)
    after = evidence_stat_fingerprint(run_dir, expected)
    if before != after:
        raise EvidenceInventoryError("Evidence files changed during verification")
    return after


def evidence_stat_fingerprint(
    run_dir: Path,
    inventory: dict[str, object],
) -> str:
    """Fingerprint cheap filesystem metadata for a previously verified tree."""

    validate_evidence_inventory(inventory)
    root = run_dir.expanduser().resolve(strict=True)
    expected_paths = [str(row["path"]) for row in inventory["files"]]
    observed_paths = _run_evidence_paths(root)
    if observed_paths != expected_paths:
        raise EvidenceInventoryError("Run evidence file set differs from its inventory")

    digest = hashlib.sha256()
    for relative in expected_paths:
        path = root / relative
        for value in (relative, *_file_stat(path)):
            digest.update(str(value).encode("utf-8"))
            digest.update(b"\0")
        digest.update(b"\n")
    return digest.hexdigest()


def assert_inventory_matches_records(
    inventory: dict[str, object],
    records: Mapping[str, tuple[int, str]],
) -> None:
    """Compare inventory evidence with hashes already computed by an archiver."""

    validate_evidence_inventory(inventory)
    expected = {
        str(row["path"]): (int(row["size_bytes"]), str(row["sha256"]))
        for row in inventory["files"]
    }
    observed = {
        path: identity
        for path, identity in records.items()
        if PurePosixPath(path).parts[0] in EVIDENCE_SCOPES
    }
    if observed != expected:
        missing = sorted(expected.keys() - observed.keys())
        extra = sorted(observed.keys() - expected.keys())
        changed = sorted(
            path
            for path in expected.keys() & observed.keys()
            if expected[path] != observed[path]
        )
        detail = _difference_detail(missing, extra, changed)
        raise EvidenceInventoryError(f"Archived evidence differs from inventory: {detail}")


def _file_stat(
    path: Path,
    *,
    follow_symlinks: bool = False,
) -> tuple[int, int, int, int, int]:
    metadata = path.stat(follow_symlinks=follow_symlinks)
    if not stat.S_ISREG(metadata.st_mode):
        raise EvidenceInventoryError(f"Evidence path is not a regular file: {path}")
    return (
        metadata.st_dev, metadata.st_ino, metadata.st_size,
        metadata.st_mtime_ns, metadata.st_ctime_ns,
    )


def _regular_tree_files(root: Path) -> list[tuple[str, Path]]:
    files = []
    for directory, directory_names, file_names in os.walk(
        root,
        topdown=True,
        onerror=_raise_walk_error,
        followlinks=False,
    ):
        directory_path = Path(directory)
        directory_names.sort()
        file_names.sort()
        for name in directory_names:
            path = directory_path / name
            if path.is_symlink():
                raise EvidenceInventoryError(f"Evidence tree contains a symlink: {path}")
        for name in file_names:
            path = directory_path / name
            metadata = path.stat(follow_symlinks=False)
            if path.is_symlink() or not stat.S_ISREG(metadata.st_mode):
                raise EvidenceInventoryError(
                    f"Evidence tree contains a non-regular file: {path}"
                )
            files.append((path.relative_to(root).as_posix(), path))
    return files


def _raise_walk_error(error: OSError) -> None:
    path = error.filename or "<unknown>"
    raise EvidenceInventoryError(f"Cannot read evidence directory: {path}") from error


def _run_evidence_paths(root: Path) -> list[str]:
    paths = []
    for scope in EVIDENCE_SCOPES:
        directory = root / scope
        if directory.is_symlink() or not directory.is_dir():
            raise EvidenceInventoryError(f"Missing evidence directory: {directory}")
        paths.extend(
            f"{scope}/{relative}"
            for relative, _path in _regular_tree_files(directory)
        )
    return sorted(paths)


def _logical_input_paths(inputs: list[tuple[str, Path]]) -> list[str]:
    paths = []
    for prefix, source in inputs:
        if source.is_symlink():
            raise EvidenceInventoryError(f"Evidence input is a symlink: {source}")
        if source.is_file():
            paths.append(prefix)
        elif source.is_dir():
            paths.extend(
                f"{prefix}/{relative}"
                for relative, _path in _regular_tree_files(source)
            )
        else:
            raise EvidenceInventoryError(f"Evidence input disappeared: {source}")
    return sorted(paths)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(READ_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_stable_regular_file(
    path: Path,
    *,
    allow_staged_symlink: bool = False,
) -> bytes:
    try:
        before = _file_stat(path, follow_symlinks=allow_staged_symlink)
        content = path.read_bytes()
        after = _file_stat(path, follow_symlinks=allow_staged_symlink)
    except OSError as exc:
        raise EvidenceInventoryError(f"Cannot read evidence inventory: {path}") from exc
    if before != after:
        raise EvidenceInventoryError(f"Evidence inventory changed while reading: {path}")
    return content


def _inventory_payload(
    rows: list[dict[str, object]],
    *,
    scopes: tuple[str, ...] = EVIDENCE_SCOPES,
) -> dict[str, object]:
    tree = hashlib.sha256()
    total_bytes = 0
    for row in rows:
        size = int(row["size_bytes"])
        total_bytes += size
        tree.update(str(row["path"]).encode("utf-8"))
        tree.update(b"\0")
        tree.update(str(size).encode("ascii"))
        tree.update(b"\0")
        tree.update(str(row["sha256"]).encode("ascii"))
        tree.update(b"\n")
    return {
        "schema_version": INVENTORY_SCHEMA_VERSION,
        "scope": list(scopes),
        "file_count": len(rows),
        "total_bytes": total_bytes,
        "tree_sha256": tree.hexdigest(),
        "files": rows,
    }


def _safe_relative_path(raw: object) -> str:
    if not isinstance(raw, str) or not raw or any(char in raw for char in "\\\0\r\n"):
        raise EvidenceInventoryError(f"Unsafe evidence path: {raw!r}")
    path = PurePosixPath(raw)
    if path.is_absolute() or path.as_posix() != raw or ".." in path.parts:
        raise EvidenceInventoryError(f"Unsafe evidence path: {raw!r}")
    if not path.parts or path.parts[0] not in EVIDENCE_SCOPES:
        raise EvidenceInventoryError(f"Evidence path is outside the declared scope: {raw!r}")
    return path.as_posix()


def _validate_scopes(scopes: tuple[str, ...]) -> None:
    canonical = tuple(scope for scope in EVIDENCE_SCOPES if scope in scopes)
    if not scopes or scopes != canonical:
        raise EvidenceInventoryError(f"Invalid evidence inventory scopes: {scopes!r}")


def _raise_inventory_difference(
    expected: dict[str, object],
    observed: dict[str, object],
) -> None:
    expected_rows = {str(row["path"]): row for row in expected["files"]}
    observed_rows = {str(row["path"]): row for row in observed["files"]}
    missing = sorted(expected_rows.keys() - observed_rows.keys())
    extra = sorted(observed_rows.keys() - expected_rows.keys())
    changed = sorted(
        path
        for path in expected_rows.keys() & observed_rows.keys()
        if expected_rows[path] != observed_rows[path]
    )
    detail = _difference_detail(missing, extra, changed)
    raise EvidenceInventoryError(f"Run evidence differs from its inventory: {detail}")


def _difference_detail(
    missing: list[str],
    extra: list[str],
    changed: list[str],
) -> str:
    parts = []
    if missing:
        parts.append(f"missing={missing[:5]}")
    if extra:
        parts.append(f"extra={extra[:5]}")
    if changed:
        parts.append(f"changed={changed[:5]}")
    return "; ".join(parts) or "summary fields changed"


def _write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        temporary.chmod(0o644)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _input_mapping(values: Sequence[str]) -> dict[str, Path]:
    inputs = {}
    for value in values:
        logical, separator, source = value.partition("=")
        if not separator or not source:
            raise EvidenceInventoryError(
                f"Evidence input must be LOGICAL_PATH=SOURCE_PATH: {value!r}"
            )
        logical = _safe_relative_path(logical)
        if logical in inputs:
            raise EvidenceInventoryError(f"Duplicate evidence input: {logical}")
        inputs[logical] = Path(source)
    return inputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("create")
    create.add_argument("--scope", choices=EVIDENCE_SCOPES, required=True)
    create.add_argument("--input", action="append", required=True)
    create.add_argument("--output", type=Path, required=True)
    combine = commands.add_parser("combine")
    for scope in EVIDENCE_SCOPES:
        combine.add_argument(f"--{scope}", type=Path, required=True)
    combine.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "create":
        write_evidence_inventory(
            args.output,
            _input_mapping(args.input),
            scopes=(args.scope,),
        )
    else:
        fragments = {}
        for scope in EVIDENCE_SCOPES:
            fragments[scope] = load_evidence_inventory(
                getattr(args, scope), scopes=(scope,)
            )
        _write_json_atomic(args.output, combine_evidence_inventories(fragments))


if __name__ == "__main__":
    main()
