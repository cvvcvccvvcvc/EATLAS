#!/usr/bin/env python3
"""Merge per-gene alignment evidence outputs."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import shutil
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-tasks", required=True, type=Path)
    parser.add_argument("--taxonomy-presets", required=True, type=Path)
    parser.add_argument("--taxonomy-failures", required=True, type=Path)
    parser.add_argument("--target-features", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--result-dir", action="append", required=True, type=Path)
    parser.add_argument("--compact-events", action="store_true")
    return parser.parse_args()


FEATURE_COVERAGE_FIELDS = [
    "gene_id",
    "strategy",
    "feature_type",
    "feature_id",
    "genomic_accession",
    "genomic_start1",
    "genomic_end1",
    "target_start0",
    "target_end0",
    "length_bp",
    "ortholog_count",
    "orthologs_covered",
    "covered_bases",
    "coverage_breadth",
    "depth_bases",
    "mean_depth",
]

COMPACT_EVENT_FIELDS = [
    "gene_id",
    "event_type",
    "target_start0",
    "target_end0",
    "genomic_accession",
    "genomic_start1",
    "genomic_end1",
    "ref",
    "alt",
    "support_row_count",
    "support_ortholog_count",
    "support_strategy_count",
    "strategies",
    "tools",
    "presets",
    "tax_id_count",
    "taxname_count",
    "qc_flags",
]


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


def read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
            count += 1
    return count


def interval_union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def fmt_fraction(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return "0.000000"
    return f"{numerator / denominator:.6f}"


@dataclass(frozen=True)
class FeatureIntervalIndex:
    features_by_gene: dict[str, list[dict[str, str]]]
    starts_by_gene: dict[str, list[int]]

    @classmethod
    def from_features(cls, features_by_gene: dict[str, list[dict[str, str]]]) -> "FeatureIntervalIndex":
        starts_by_gene: dict[str, list[int]] = {}
        for gene_id, features in features_by_gene.items():
            features.sort(key=lambda row: int(row.get("target_start0") or 0))
            starts_by_gene[gene_id] = [int(row.get("target_start0") or 0) for row in features]
        return cls(features_by_gene, starts_by_gene)

    def overlapping(self, gene_id: str, start0: int, end0: int) -> list[dict[str, str]]:
        features = self.features_by_gene.get(gene_id, [])
        if not features:
            return []
        starts = self.starts_by_gene.get(gene_id, [])
        limit = bisect.bisect_left(starts, end0)
        return [
            feature
            for feature in features[:limit]
            if int(feature.get("target_end0") or 0) > start0
        ]


def summarize_feature_coverage(
    target_features: Path,
    summaries_path: Path,
    segments_path: Path,
    output: Path,
) -> int:
    features_by_gene: dict[str, list[dict[str, str]]] = {}
    for row in read_tsv_gz(target_features):
        gene_id = row.get("gene_id", "")
        if gene_id:
            features_by_gene.setdefault(gene_id, []).append(row)
    feature_index = FeatureIntervalIndex.from_features(features_by_gene)

    orthologs_by_gene_strategy: dict[tuple[str, str], set[str]] = {}
    for row in read_tsv_gz(summaries_path):
        gene_id = row.get("gene_id", "")
        strategy = row.get("strategy", "")
        ortholog_gene_id = row.get("ortholog_gene_id", "")
        if gene_id and strategy and ortholog_gene_id:
            orthologs_by_gene_strategy.setdefault((gene_id, strategy), set()).add(ortholog_gene_id)

    strategies_by_gene: dict[str, set[str]] = {}
    for gene_id, strategy in orthologs_by_gene_strategy:
        strategies_by_gene.setdefault(gene_id, set()).add(strategy)

    overlap_any: dict[tuple[str, str, str], list[tuple[int, int]]] = {}
    overlap_by_ortholog: dict[tuple[str, str, str], dict[str, list[tuple[int, int]]]] = {}

    for segment in read_tsv_gz(segments_path):
        gene_id = segment.get("gene_id", "")
        strategy = segment.get("strategy", "")
        ortholog_gene_id = segment.get("ortholog_gene_id", "")
        if not gene_id or not strategy or not ortholog_gene_id:
            continue
        seg_start = int(segment.get("target_start0") or 0)
        seg_end = int(segment.get("target_end0") or 0)
        if seg_end <= seg_start:
            continue
        strategies_by_gene.setdefault(gene_id, set()).add(strategy)
        for feature in feature_index.overlapping(gene_id, seg_start, seg_end):
            feature_start = int(feature.get("target_start0") or 0)
            feature_end = int(feature.get("target_end0") or 0)
            overlap_start = max(seg_start, feature_start)
            overlap_end = min(seg_end, feature_end)
            if overlap_end <= overlap_start:
                continue
            key = (gene_id, strategy, feature["feature_id"])
            overlap_any.setdefault(key, []).append((overlap_start, overlap_end))
            overlap_by_ortholog.setdefault(key, {}).setdefault(ortholog_gene_id, []).append(
                (overlap_start, overlap_end)
            )

    rows: list[dict[str, object]] = []
    for gene_id in sorted(features_by_gene, key=lambda value: int(value) if value.isdigit() else value):
        features = sorted(features_by_gene[gene_id], key=lambda row: int(row.get("target_start0") or 0))
        for strategy in sorted(strategies_by_gene.get(gene_id, [])):
            ortholog_count = len(orthologs_by_gene_strategy.get((gene_id, strategy), set()))
            for feature in features:
                key = (gene_id, strategy, feature["feature_id"])
                per_ortholog = overlap_by_ortholog.get(key, {})
                per_ortholog_lengths = [
                    interval_union_length(intervals) for intervals in per_ortholog.values()
                ]
                depth_bases = sum(per_ortholog_lengths)
                orthologs_covered = sum(1 for length in per_ortholog_lengths if length > 0)
                covered_bases = interval_union_length(overlap_any.get(key, []))
                length_bp = int(feature.get("length_bp") or 0)
                rows.append(
                    {
                        "gene_id": gene_id,
                        "strategy": strategy,
                        "feature_type": feature.get("feature_type", ""),
                        "feature_id": feature.get("feature_id", ""),
                        "genomic_accession": feature.get("genomic_accession", ""),
                        "genomic_start1": feature.get("genomic_start1", ""),
                        "genomic_end1": feature.get("genomic_end1", ""),
                        "target_start0": feature.get("target_start0", ""),
                        "target_end0": feature.get("target_end0", ""),
                        "length_bp": length_bp,
                        "ortholog_count": ortholog_count,
                        "orthologs_covered": orthologs_covered,
                        "covered_bases": covered_bases,
                        "coverage_breadth": fmt_fraction(covered_bases, length_bp),
                        "depth_bases": depth_bases,
                        "mean_depth": fmt_fraction(depth_bases, length_bp),
                    }
                )

    return write_tsv_gz(output, FEATURE_COVERAGE_FIELDS, rows)


def unique_paths(paths: list[Path]) -> list[Path]:
    seen: set[Path] = set()
    unique: list[Path] = []
    for path in paths:
        key = path.resolve() if path.exists() else path
        if key in seen:
            continue
        seen.add(key)
        unique.append(path)
    return unique


def validate_result_dirs(paths: list[Path]) -> list[Path]:
    if not paths:
        raise ValueError("At least one --result-dir is required")
    result_dirs = unique_paths(paths)
    missing = [str(path) for path in result_dirs if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing alignment result dir(s): " + ", ".join(missing))
    not_dirs = [str(path) for path in result_dirs if not path.is_dir()]
    if not_dirs:
        raise NotADirectoryError("Alignment result path(s) are not directories: " + ", ".join(not_dirs))
    missing_manifest = [str(path) for path in result_dirs if not (path / "manifest.json").exists()]
    if missing_manifest:
        raise FileNotFoundError("Alignment result dir(s) missing manifest.json: " + ", ".join(missing_manifest))
    return sorted(result_dirs, key=lambda path: path.name)


def merge_tsv_gz(paths: list[Path], output: Path) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    header_written = False
    expected_header: list[str] | None = None
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
                if expected_header is None:
                    expected_header = header
                elif header != expected_header:
                    raise ValueError(
                        f"Header mismatch while merging {path}: expected {expected_header}, observed {header}"
                    )
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


def create_compact_event_table(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE events (
            gene_id TEXT,
            event_type TEXT,
            target_start0 TEXT,
            target_end0 TEXT,
            genomic_accession TEXT,
            genomic_start1 TEXT,
            genomic_end1 TEXT,
            ref TEXT,
            alt TEXT,
            ortholog_gene_id TEXT,
            strategy TEXT,
            tool TEXT,
            preset TEXT,
            tax_id TEXT,
            taxname TEXT,
            qc_flags TEXT
        )
        """
    )


def insert_event_rows(conn: sqlite3.Connection, path: Path, batch_size: int = 10000) -> int:
    if not path.exists():
        return 0
    fields = [
        "gene_id",
        "event_type",
        "target_start0",
        "target_end0",
        "genomic_accession",
        "genomic_start1",
        "genomic_end1",
        "ref",
        "alt",
        "ortholog_gene_id",
        "strategy",
        "tool",
        "preset",
        "tax_id",
        "taxname",
        "qc_flags",
    ]
    required = set(fields[:9])
    count = 0
    batch: list[tuple[str, ...]] = []
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Events table {path} missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            batch.append(tuple(row.get(field, "") for field in fields))
            count += 1
            if len(batch) >= batch_size:
                conn.executemany(
                    "INSERT INTO events VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    batch,
                )
                batch.clear()
    if batch:
        conn.executemany(
            "INSERT INTO events VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            batch,
        )
    return count


def compact_event_flags(raw_flags: str) -> str:
    flags = sorted({flag for item in raw_flags.split(",") for flag in item.split("|") if flag})
    return ",".join(flags)


def write_compact_events(paths: list[Path], output: Path) -> tuple[int, int]:
    output.parent.mkdir(parents=True, exist_ok=True)
    db_path = output.parent / "alignment_event_support.sqlite"
    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(db_path)
    try:
        conn.execute("PRAGMA journal_mode = OFF")
        conn.execute("PRAGMA synchronous = OFF")
        conn.execute("PRAGMA temp_store = MEMORY")
        create_compact_event_table(conn)
        raw_count = 0
        for path in paths:
            raw_count += insert_event_rows(conn, path)
        conn.commit()
        conn.execute(
            """
            CREATE INDEX events_key_idx ON events (
                gene_id,
                event_type,
                target_start0,
                target_end0,
                genomic_accession,
                genomic_start1,
                genomic_end1,
                ref,
                alt
            )
            """
        )

        query = """
            SELECT
                gene_id,
                event_type,
                target_start0,
                target_end0,
                genomic_accession,
                genomic_start1,
                genomic_end1,
                ref,
                alt,
                COUNT(*) AS support_row_count,
                COUNT(DISTINCT ortholog_gene_id) AS support_ortholog_count,
                COUNT(DISTINCT strategy) AS support_strategy_count,
                GROUP_CONCAT(DISTINCT strategy) AS strategies,
                GROUP_CONCAT(DISTINCT tool) AS tools,
                GROUP_CONCAT(DISTINCT preset) AS presets,
                COUNT(DISTINCT tax_id) AS tax_id_count,
                COUNT(DISTINCT taxname) AS taxname_count,
                GROUP_CONCAT(DISTINCT qc_flags) AS qc_flags
            FROM events
            GROUP BY
                gene_id,
                event_type,
                target_start0,
                target_end0,
                genomic_accession,
                genomic_start1,
                genomic_end1,
                ref,
                alt
            ORDER BY
                CAST(gene_id AS INTEGER),
                CAST(target_start0 AS INTEGER),
                event_type,
                ref,
                alt
        """
        compact_count = 0
        with gzip.open(output, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=COMPACT_EVENT_FIELDS, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in conn.execute(query):
                record = dict(zip(COMPACT_EVENT_FIELDS, row))
                record["qc_flags"] = compact_event_flags(record.get("qc_flags") or "")
                writer.writerow(record)
                compact_count += 1
        return compact_count, raw_count
    finally:
        conn.close()
        if db_path.exists():
            db_path.unlink()


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


def manifest_strategies(manifests: list[dict]) -> list[str]:
    strategies: set[str] = set()
    for manifest in manifests:
        strategy_list = manifest.get("strategies", []) or []
        if strategy_list:
            for strategy in strategy_list:
                if strategy:
                    strategies.add(str(strategy))
        elif manifest.get("strategy"):
            strategies.add(str(manifest["strategy"]))
    return sorted(strategies)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result_dirs = validate_result_dirs(args.result_dir)

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
    feature_coverage_count = summarize_feature_coverage(
        args.target_features,
        args.outdir / "ortholog_alignment_summary.tsv.gz",
        args.outdir / "alignment_segments.tsv.gz",
        args.outdir / "feature_coverage.tsv.gz",
    )
    event_inputs = [path / "alignment_events.tsv.gz" for path in result_dirs]
    if args.compact_events:
        event_count, raw_event_count = write_compact_events(
            event_inputs,
            args.outdir / "alignment_events.tsv.gz",
        )
    else:
        event_count = merge_tsv_gz(
            event_inputs,
            args.outdir / "alignment_events.tsv.gz",
        )
        raw_event_count = event_count
    failure_count = merge_tsv_gz(
        [path / "failures.tsv.gz" for path in result_dirs],
        args.outdir / "failures.tsv.gz",
    )
    native_file_count = copy_native(result_dirs, args.outdir)
    manifests = load_manifests(result_dirs)
    strategies = manifest_strategies(manifests)
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
        "feature_coverage_count": feature_coverage_count,
        "alignment_event_mode": "compact_support" if args.compact_events else "raw",
        "raw_alignment_event_count": raw_event_count,
        "alignment_event_count": event_count,
        "failure_count": failure_count,
        "native_file_count": native_file_count,
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
