#!/usr/bin/env python3
"""Fail early when selected alignment and annotation dependencies are unavailable."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

from genomics.vep.local_runtime import probe_local_vep, resolve_executable


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


def require_executable(name: str, raw: str, errors: list[str]) -> None:
    if resolve_executable(raw) is None:
        errors.append(f"{name} is not executable or was not found: {raw}")


def require_python_module(name: str, errors: list[str]) -> None:
    if importlib.util.find_spec(name) is None:
        errors.append(f"Python module not importable in task environment: {name}")


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
