#!/usr/bin/env python3
"""Fail early when selected pipeline modes cannot run in the task environment."""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--alignment-strategies", required=True)
    parser.add_argument("--datasets-bin", required=True)
    parser.add_argument("--minimap2-bin", required=True)
    parser.add_argument("--nucmer-bin", required=True)
    parser.add_argument("--bwa-bin", required=True)
    parser.add_argument("--samtools-bin", required=True)
    parser.add_argument("--out-json", required=True, type=Path)
    return parser.parse_args()


def split_strategies(raw: str) -> set[str]:
    return {item.strip() for item in raw.split(",") if item.strip()}


def executable_exists(raw: str) -> bool:
    path = Path(raw)
    if path.parent != Path(".") or raw.startswith("."):
        return path.exists() and path.is_file()
    return shutil.which(raw) is not None


def require_executable(name: str, raw: str, errors: list[str]) -> None:
    if not executable_exists(raw):
        errors.append(f"{name} not found: {raw}")


def require_python_module(name: str, errors: list[str]) -> None:
    if importlib.util.find_spec(name) is None:
        errors.append(f"Python module not importable in task environment: {name}")


def main() -> None:
    args = parse_args()
    strategies = split_strategies(args.alignment_strategies)
    errors: list[str] = []

    if args.stage in {"all", "fetch"}:
        require_executable("datasets", args.datasets_bin, errors)

    if args.stage in {"all", "align"}:
        require_executable("bedtools", "bedtools", errors)
        if strategies & {"minimap2_asm10", "minimap2_asm20", "minimap2_taxonomy_adaptive"}:
            require_executable("minimap2", args.minimap2_bin, errors)
        if "nucmer" in strategies:
            require_python_module("pysam", errors)
            require_executable("nucmer", args.nucmer_bin, errors)
        if "bwa_pseudoreads" in strategies:
            require_python_module("pysam", errors)
            require_python_module("bam_filtering_v1", errors)
            require_executable("bwa", args.bwa_bin, errors)
            require_executable("samtools", args.samtools_bin, errors)

    if args.stage in {"all", "annotate"}:
        require_python_module("pysam", errors)

    args.out_json.write_text(
        json.dumps(
            {
                "stage": args.stage,
                "alignment_strategies": sorted(strategies),
                "python": sys.executable,
                "ok": not errors,
                "errors": errors,
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
