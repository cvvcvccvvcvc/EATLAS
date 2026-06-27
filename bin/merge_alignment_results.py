#!/usr/bin/env python3
"""Merge per-gene alignment evidence outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-tasks", required=True, type=Path)
    parser.add_argument("--taxonomy-presets", required=True, type=Path)
    parser.add_argument("--taxonomy-failures", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--result-dir", action="append", required=True, type=Path)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def copy_or_keep(src: Path, dst: Path) -> None:
    if src.resolve() == dst.resolve():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def count_tsv_gz_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with gzip.open(path, "rt", newline="") as handle:
        next(handle, None)
        return sum(1 for _ in handle)


def merge_tsv_gz(paths: list[Path], output: Path) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    header_written = False
    with gzip.open(output, "wt", newline="") as out_handle:
        writer = None
        for path in paths:
            if not path.exists():
                continue
            with gzip.open(path, "rt", newline="") as in_handle:
                reader = csv.reader(in_handle, delimiter="\t")
                header = next(reader, None)
                if header is None:
                    continue
                if not header_written:
                    writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(header)
                    header_written = True
                for row in reader:
                    writer.writerow(row)
                    count += 1

    if not header_written:
        with gzip.open(output, "wt", newline="") as out_handle:
            out_handle.write("")
    return count


def copy_native(result_dirs: list[Path], outdir: Path) -> int:
    copied = 0
    native_root = outdir / "native"
    for result_dir in result_dirs:
        native_dir = result_dir / "native"
        if not native_dir.exists():
            continue
        strategy_dir = native_root / result_dir.name
        for src in sorted(native_dir.rglob("*")):
            if not src.is_file():
                continue
            dst = strategy_dir / src.relative_to(native_dir)
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            copied += 1
    return copied


def load_manifests(result_dirs: list[Path]) -> list[dict]:
    manifests = []
    for result_dir in result_dirs:
        path = result_dir / "manifest.json"
        if path.exists():
            manifests.append(json.loads(path.read_text()))
    return manifests


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result_dirs = sorted(args.result_dir, key=lambda path: path.name)

    copy_or_keep(args.alignment_tasks, args.outdir / "alignment_tasks.tsv.gz")
    copy_or_keep(args.taxonomy_presets, args.outdir / "taxonomy_presets.tsv.gz")
    copy_or_keep(args.taxonomy_failures, args.outdir / "taxonomy_failures.tsv.gz")

    summary_count = merge_tsv_gz(
        [path / "ortholog_alignment_summary.tsv.gz" for path in result_dirs],
        args.outdir / "ortholog_alignment_summary.tsv.gz",
    )
    segment_count = merge_tsv_gz(
        [path / "alignment_segments.tsv.gz" for path in result_dirs],
        args.outdir / "alignment_segments.tsv.gz",
    )
    event_count = merge_tsv_gz(
        [path / "alignment_events.tsv.gz" for path in result_dirs],
        args.outdir / "alignment_events.tsv.gz",
    )
    failure_count = merge_tsv_gz(
        [path / "failures.tsv.gz" for path in result_dirs],
        args.outdir / "failures.tsv.gz",
    )
    native_file_count = copy_native(result_dirs, args.outdir)
    manifests = load_manifests(result_dirs)
    strategies = sorted({manifest.get("strategy", "") for manifest in manifests if manifest.get("strategy")})
    gene_ids = sorted({str(manifest.get("gene_id", "")) for manifest in manifests if manifest.get("gene_id")})

    manifest = {
        "created_at": utc_now(),
        "stage": "alignment",
        "strategy_count": len(strategies),
        "strategies": strategies,
        "gene_count": len(gene_ids),
        "alignment_task_count": count_tsv_gz_rows(args.alignment_tasks),
        "taxonomy_tax_id_count": count_tsv_gz_rows(args.taxonomy_presets),
        "taxonomy_failure_count": count_tsv_gz_rows(args.taxonomy_failures),
        "ortholog_alignment_summary_count": summary_count,
        "alignment_segment_count": segment_count,
        "alignment_event_count": event_count,
        "failure_count": failure_count,
        "native_file_count": native_file_count,
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
