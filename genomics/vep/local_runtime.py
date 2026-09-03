"""Shared command-line contract for the local Ensembl VEP cache."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


VEP_ASSEMBLY = "GRCh38"
VEP_SPECIES = "homo_sapiens"
LOCAL_VEP_PROBE_TIMEOUT_SECONDS = 120


def local_vep_cache_args(*, release: str, cache_dir: str | Path) -> list[str]:
    """Return the cache-selection flags shared by probes and annotations."""

    normalized_release = str(release).strip()
    if not normalized_release:
        raise ValueError("Local VEP requires an explicit release")
    return [
        "--offline",
        "--cache",
        "--refseq",
        "--use_given_ref",
        "--species",
        VEP_SPECIES,
        "--assembly",
        VEP_ASSEMBLY,
        "--cache_version",
        normalized_release,
        "--dir_cache",
        str(cache_dir),
    ]


def resolve_executable(raw: str | Path) -> str | None:
    """Resolve a command name or explicit executable path."""

    value = str(raw)
    path = Path(value).expanduser()
    if path.parent != Path(".") or value.startswith("."):
        if path.is_file() and os.access(path, os.X_OK):
            return str(path.resolve())
        return None
    return shutil.which(value)


def probe_local_vep(
    *,
    release: str,
    executable: str | Path,
    cache_dir: str | Path | None,
    errors: list[str],
) -> dict[str, object]:
    """Check the configured local VEP command and release-specific cache."""

    initial_error_count = len(errors)
    normalized_release = str(release).strip()
    resolved_executable = resolve_executable(executable)
    cache_value = "" if cache_dir is None else str(cache_dir).strip()
    cache_path = Path(cache_value).expanduser() if cache_value else None
    if not normalized_release:
        errors.append("Local VEP release is empty")
    if resolved_executable is None:
        errors.append(f"Local VEP is not executable or was not found: {executable}")
    if cache_path is None:
        errors.append("Local VEP cache directory is empty")
    elif not cache_path.is_dir():
        errors.append(f"Local VEP cache directory does not exist: {cache_path}")
    if len(errors) > initial_error_count:
        return {
            "backend": "local",
            "cache_dir": cache_value,
            "executable": str(executable),
            "probe": "not_run",
            "release": normalized_release,
        }

    command = [
        resolved_executable,
        *local_vep_cache_args(
            release=normalized_release,
            cache_dir=cache_path,
        ),
        "--show_cache_info",
    ]
    try:
        process = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            timeout=LOCAL_VEP_PROBE_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired:
        errors.append(
            f"Local VEP cache probe exceeded {LOCAL_VEP_PROBE_TIMEOUT_SECONDS} seconds"
        )
    except OSError as exc:
        errors.append(f"Could not start local VEP executable {executable}: {exc}")
    else:
        if process.returncode != 0:
            detail = (process.stderr or process.stdout).strip()[-2_000:]
            errors.append(
                f"Local VEP cache probe failed with exit code {process.returncode}: {detail}"
            )

    return {
        "backend": "local",
        "cache_dir": str(cache_path),
        "executable": str(executable),
        "probe": "passed" if len(errors) == initial_error_count else "failed",
        "release": normalized_release,
        "resolved_executable": resolved_executable,
    }
