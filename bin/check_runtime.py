#!/usr/bin/env python3
"""Fail early when selected alignment and annotation dependencies are unavailable."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from genomics.vep.local_runtime import local_vep_cache_args


LOCAL_VEP_PROBE_TIMEOUT_SECONDS = 120


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-strategies", required=True)
    parser.add_argument("--vep-backend", required=True, choices=("rest", "local"))
    parser.add_argument("--vep-release", required=True)
    parser.add_argument("--vep-executable", required=True)
    parser.add_argument("--vep-cache-dir", required=True)
    parser.add_argument("--out-json", required=True, type=Path)
    return parser.parse_args()


def split_strategies(raw: str) -> set[str]:
    return {item.strip() for item in raw.split(",") if item.strip()}


def resolve_executable(raw: str) -> str | None:
    path = Path(raw).expanduser()
    if path.parent != Path(".") or raw.startswith("."):
        if path.is_file() and os.access(path, os.X_OK):
            return str(path.resolve())
        return None
    return shutil.which(raw)


def require_executable(name: str, raw: str, errors: list[str]) -> None:
    if resolve_executable(raw) is None:
        errors.append(f"{name} is not executable or was not found: {raw}")


def require_python_module(name: str, errors: list[str]) -> None:
    if importlib.util.find_spec(name) is None:
        errors.append(f"Python module not importable in task environment: {name}")


def probe_local_vep(
    *,
    release: str,
    executable: str,
    cache_dir: str,
    errors: list[str],
) -> dict[str, object]:
    """Check the configured local VEP command and cache on the task host."""

    initial_error_count = len(errors)
    resolved_executable = resolve_executable(executable)
    cache_path = Path(cache_dir).expanduser()
    if not release.strip():
        errors.append("Local VEP release is empty")
    if resolved_executable is None:
        errors.append(f"Local VEP is not executable or was not found: {executable}")
    if not cache_path.is_dir():
        errors.append(f"Local VEP cache directory does not exist: {cache_path}")
    if len(errors) > initial_error_count:
        return {
            "backend": "local",
            "cache_dir": str(cache_path),
            "executable": executable,
            "probe": "not_run",
            "release": release,
        }

    command = [
        resolved_executable,
        *local_vep_cache_args(release=release, cache_dir=cache_path),
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
        "executable": executable,
        "probe": "passed" if len(errors) == initial_error_count else "failed",
        "release": release,
        "resolved_executable": resolved_executable,
    }


def main() -> None:
    args = parse_args()
    strategies = split_strategies(args.alignment_strategies)
    errors: list[str] = []

    require_python_module("pysam", errors)
    if any(strategy.startswith("minimap2_") for strategy in strategies):
        require_executable("minimap2", "minimap2", errors)
    if "nucmer" in strategies:
        require_executable("nucmer", "nucmer", errors)
    if any(strategy.startswith("bwa_pseudoreads_") for strategy in strategies):
        require_executable("bwa", "bwa", errors)
        require_executable("samtools", "samtools", errors)

    if args.vep_backend == "local":
        vep = probe_local_vep(
            release=args.vep_release,
            executable=args.vep_executable,
            cache_dir=args.vep_cache_dir,
            errors=errors,
        )
    else:
        vep = {"backend": "rest", "probe": "not_applicable"}

    args.out_json.write_text(
        json.dumps(
            {
                "alignment_strategies": sorted(strategies),
                "python": sys.executable,
                "ok": not errors,
                "errors": errors,
                "vep": vep,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    if errors:
        raise SystemExit("Runtime dependency check failed:\n- " + "\n- ".join(errors))


if __name__ == "__main__":
    main()
