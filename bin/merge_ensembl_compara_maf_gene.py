#!/usr/bin/env python3
"""Consolidate all Ensembl Compara MAF fragments for one target gene."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

from feature_coverage import summarize_feature_coverage
from run_ensembl_compara_maf_alignment import (
    EVENT_FIELDS,
    FAILURE_FIELDS,
    OUTPUT_GZIP_COMPRESSLEVEL,
    SEGMENT_FIELDS,
    SUMMARY_FIELDS,
    interval_union_length,
    write_tsv_gz,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--task-dir", required=True, type=Path)
    parser.add_argument("--fragment-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def fragment_dirs(root: Path, gene_id: str) -> list[Path]:
    if not root.is_dir():
        raise NotADirectoryError(f"MAF fragment root is not a directory: {root}")
    fragments = sorted(path for path in root.iterdir() if path.is_dir())
    if not fragments:
        raise ValueError(f"No MAF fragments found for gene {gene_id} in {root}")
    for fragment in fragments:
        manifest_path = fragment / "manifest.json"
        if not manifest_path.exists():
            raise FileNotFoundError(f"MAF fragment missing manifest.json: {fragment}")
        manifest = json.loads(manifest_path.read_text())
        if str(manifest.get("gene_id") or "") != gene_id:
            raise ValueError(
                f"MAF fragment gene mismatch in {fragment}: "
                f"expected {gene_id}, observed {manifest.get('gene_id')!r}"
            )
    return fragments


def merge_tsv_gz(paths: list[Path], output: Path, expected_fields: list[str]) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(output, "wt", newline="", compresslevel=OUTPUT_GZIP_COMPRESSLEVEL) as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=expected_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for path in paths:
            with gzip.open(path, "rt", newline="") as in_handle:
                reader = csv.DictReader(in_handle, delimiter="\t")
                if reader.fieldnames != expected_fields:
                    raise ValueError(
                        f"Header mismatch in {path}: expected {expected_fields}, observed {reader.fieldnames}"
                    )
                for row in reader:
                    writer.writerow(row)
                    count += 1
    return count


def split_flags(value: str) -> set[str]:
    return {item for item in re.split(r"[,|]", value or "") if item}


def summary_key(row: dict[str, str]) -> tuple[str, str]:
    return row.get("strategy", ""), row.get("ortholog_gene_id", "")


def consolidate_summaries(
    fragment_paths: list[Path],
    segments_path: Path,
    events_path: Path,
    target_length: int,
) -> list[dict[str, object]]:
    base_rows: dict[tuple[str, str], dict[str, str]] = {}
    flags: dict[tuple[str, str], set[str]] = defaultdict(set)
    for path in fragment_paths:
        with gzip.open(path / "ortholog_alignment_summary.tsv.gz", "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != SUMMARY_FIELDS:
                raise ValueError(f"Unexpected MAF summary header in {path}")
            for row in reader:
                key = summary_key(row)
                base_rows.setdefault(key, row)
                flags[key].update(split_flags(row.get("qc_flags", "")))

    target_intervals: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    query_intervals: dict[tuple[str, str], dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    identities: dict[tuple[str, str], list[float]] = defaultdict(list)
    segment_counts: dict[tuple[str, str], int] = defaultdict(int)
    primary_counts: dict[tuple[str, str], int] = defaultdict(int)
    segment_example: dict[tuple[str, str], dict[str, str]] = {}
    with gzip.open(segments_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = summary_key(row)
            segment_example.setdefault(key, row)
            segment_counts[key] += 1
            primary_counts[key] += int(str(row.get("is_primary") or "").lower() == "true")
            target_intervals[key].append((int(row["target_start0"]), int(row["target_end0"])))
            if row.get("query_start0") and row.get("query_end0"):
                query_intervals[key][row.get("query_id", "")].append(
                    (int(row["query_start0"]), int(row["query_end0"]))
                )
            if row.get("identity"):
                identities[key].append(float(row["identity"]))
            flags[key].update(split_flags(row.get("qc_flags", "")))

    event_counts: dict[tuple[str, str], int] = defaultdict(int)
    with gzip.open(events_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            event_counts[summary_key(row)] += 1

    rows: list[dict[str, object]] = []
    for key in sorted(set(base_rows) | set(segment_example)):
        base = dict(base_rows.get(key) or {})
        example = segment_example.get(key, {})
        strategy, ortholog_gene_id = key
        count = segment_counts[key]
        aligned_target_bp = interval_union_length(target_intervals[key])
        aligned_query_bp = sum(
            interval_union_length(intervals) for intervals in query_intervals[key].values()
        )
        identity_values = identities[key]
        row_flags = flags[key]
        if count:
            row_flags.discard("no_alignment")
        else:
            row_flags.add("no_alignment")
        row_flags.add("maf_query_coverage_not_applicable")
        row = {
            "gene_id": base.get("gene_id") or example.get("gene_id", ""),
            "ortholog_gene_id": ortholog_gene_id,
            "tax_id": base.get("tax_id") or example.get("tax_id", ""),
            "taxname": base.get("taxname") or example.get("taxname", ""),
            "strategy": strategy,
            "tool": base.get("tool") or example.get("tool", ""),
            "preset": base.get("preset") or example.get("preset", ""),
            "status": "aligned" if count else "no_alignment",
            "target_length": target_length,
            "query_length": "",
            "segment_count": count,
            "primary_segment_count": primary_counts[key],
            "secondary_segment_count": count - primary_counts[key],
            "aligned_target_bp": aligned_target_bp,
            "aligned_query_bp": aligned_query_bp,
            "target_coverage": f"{aligned_target_bp / target_length if target_length else 0.0:.6f}",
            "query_coverage": "",
            "best_identity": f"{max(identity_values) if identity_values else 0.0:.6f}",
            "mean_identity": f"{sum(identity_values) / len(identity_values) if identity_values else 0.0:.6f}",
            "event_count": event_counts[key],
            "qc_flags": ",".join(sorted(row_flags)),
        }
        rows.append(row)
    return rows


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    fragments = fragment_dirs(args.fragment_root, args.gene_id)
    task = json.loads((args.task_dir / "task.json").read_text())
    if str(task.get("gene_id") or "") != args.gene_id:
        raise ValueError(f"Task directory gene does not match --gene-id {args.gene_id}")
    target_length = int(task.get("target_length") or 0)

    segments_path = args.outdir / "alignment_segments.tsv.gz"
    events_path = args.outdir / "alignment_events.tsv.gz"
    segment_count = merge_tsv_gz(
        [path / "alignment_segments.tsv.gz" for path in fragments],
        segments_path,
        SEGMENT_FIELDS,
    )
    event_count = merge_tsv_gz(
        [path / "alignment_events.tsv.gz" for path in fragments],
        events_path,
        EVENT_FIELDS,
    )
    failure_count = merge_tsv_gz(
        [path / "failures.tsv.gz" for path in fragments],
        args.outdir / "failures.tsv.gz",
        FAILURE_FIELDS,
    )
    summary_rows = consolidate_summaries(fragments, segments_path, events_path, target_length)
    write_tsv_gz(args.outdir / "ortholog_alignment_summary.tsv.gz", SUMMARY_FIELDS, summary_rows)
    feature_coverage_count = summarize_feature_coverage(
        args.task_dir / str(task.get("target_features", "target_features.tsv.gz")),
        args.outdir / "ortholog_alignment_summary.tsv.gz",
        segments_path,
        args.outdir / "feature_coverage.tsv.gz",
    )

    manifests = [json.loads((path / "manifest.json").read_text()) for path in fragments]
    strategy = manifests[0]["strategy"]
    manifest = {
        "task_type": "maf_gene_consolidated",
        "gene_id": args.gene_id,
        "gene_ids": [args.gene_id],
        "strategy": strategy,
        "strategies": [strategy],
        "tool": manifests[0].get("tool", "ensembl_compara_maf"),
        "fragment_count": len(fragments),
        "source_chunk_ids": sorted(str(item.get("chunk_id") or "") for item in manifests),
        "summary_count": len(summary_rows),
        "segment_count": segment_count,
        "event_count": event_count,
        "feature_coverage_count": feature_coverage_count,
        "failure_count": failure_count,
        "output_gzip_compresslevel": OUTPUT_GZIP_COMPRESSLEVEL,
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
