#!/usr/bin/env python3
"""Merge per-gene alignment evidence outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from feature_coverage import summarize_feature_coverage


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-tasks", type=Path)
    parser.add_argument("--taxonomy-presets", type=Path)
    parser.add_argument("--taxonomy-failures", type=Path)
    parser.add_argument("--target-features", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--result-dir", action="append", default=[], type=Path)
    parser.add_argument("--result-root", type=Path)
    parser.add_argument("--partition-id")
    parser.add_argument("--compact-events", action="store_true")
    parser.add_argument("--events-already-compacted", action="store_true")
    return parser.parse_args()


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

STRATEGY_SUMMARY_FIELDS = [
    "strategy",
    "summary_row_count",
    "gene_count",
    "aligned_summary_row_count",
    "event_count",
    "aligned_target_bp",
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


def summarize_alignment_tasks(path: Path) -> tuple[int, int]:
    """Return total task rows and distinct genes ready for alignment."""
    task_count = 0
    ready_gene_ids: set[str] = set()
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", "status"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Alignment tasks {path} missing required columns: "
                + ", ".join(sorted(missing))
            )
        for row in reader:
            task_count += 1
            if row["status"] == "ready" and row["gene_id"]:
                ready_gene_ids.add(row["gene_id"])
    return task_count, len(ready_gene_ids)


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


def resolve_result_dirs(explicit: list[Path], root: Path | None) -> list[Path]:
    paths = list(explicit)
    if root:
        if not root.is_dir():
            raise NotADirectoryError(f"Alignment result root is not a directory: {root}")
        paths.extend(path for path in root.iterdir() if path.is_dir())
    return validate_result_dirs(paths)


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


def write_strategy_summary(
    summaries_path: Path,
    output: Path,
    expected_strategies: list[str],
) -> int:
    """Write small per-strategy aggregates from the canonical summary table."""
    aggregates: dict[str, dict[str, int | set[str]]] = {}
    with gzip.open(summaries_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", "strategy", "status", "event_count", "aligned_target_bp"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Alignment summary {summaries_path} missing required columns: "
                + ", ".join(sorted(missing))
            )
        for row in reader:
            strategy = row["strategy"]
            aggregate = aggregates.setdefault(
                strategy,
                {
                    "summary_row_count": 0,
                    "gene_ids": set(),
                    "aligned_summary_row_count": 0,
                    "event_count": 0,
                    "aligned_target_bp": 0,
                },
            )
            aggregate["summary_row_count"] += 1
            aggregate["gene_ids"].add(row["gene_id"])
            aggregate["aligned_summary_row_count"] += int(row["status"] == "aligned")
            aggregate["event_count"] += int(row["event_count"] or 0)
            aggregate["aligned_target_bp"] += int(row["aligned_target_bp"] or 0)

    for strategy in expected_strategies:
        aggregates.setdefault(
            strategy,
            {
                "summary_row_count": 0,
                "gene_ids": set(),
                "aligned_summary_row_count": 0,
                "event_count": 0,
                "aligned_target_bp": 0,
            },
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STRATEGY_SUMMARY_FIELDS, delimiter="\t")
        writer.writeheader()
        for strategy in sorted(aggregates):
            aggregate = aggregates[strategy]
            writer.writerow(
                {
                    "strategy": strategy,
                    "summary_row_count": aggregate["summary_row_count"],
                    "gene_count": len(aggregate["gene_ids"]),
                    "aligned_summary_row_count": aggregate["aligned_summary_row_count"],
                    "event_count": aggregate["event_count"],
                    "aligned_target_bp": aggregate["aligned_target_bp"],
                }
            )
    return len(aggregates)


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


def manifest_gene_ids(manifests: list[dict]) -> list[str]:
    gene_ids: set[str] = set()
    for manifest in manifests:
        if manifest.get("gene_id"):
            gene_ids.add(str(manifest["gene_id"]))
        for gene_id in manifest.get("gene_ids", []) or []:
            if gene_id:
                gene_ids.add(str(gene_id))
    return sorted(gene_ids, key=lambda value: int(value) if value.isdigit() else value)


def main() -> None:
    args = parse_args()
    if args.events_already_compacted and not args.compact_events:
        raise ValueError("--events-already-compacted requires --compact-events")
    args.outdir.mkdir(parents=True, exist_ok=True)
    result_dirs = resolve_result_dirs(args.result_dir, args.result_root)

    global_inputs = [args.alignment_tasks, args.taxonomy_presets, args.taxonomy_failures, args.target_features]
    if not args.partition_id and any(path is None for path in global_inputs):
        raise ValueError(
            "Final alignment merge requires --alignment-tasks, --taxonomy-presets, "
            "--taxonomy-failures, and --target-features"
        )
    if not args.partition_id:
        copy_or_keep(args.alignment_tasks, args.outdir / "alignment_tasks.tsv.gz")
        copy_or_keep(args.taxonomy_presets, args.outdir / "taxonomy_presets.tsv.gz")
        copy_or_keep(args.taxonomy_failures, args.outdir / "taxonomy_failures.tsv.gz")

    manifests = load_manifests(result_dirs)
    strategies = manifest_strategies(manifests)
    gene_ids = manifest_gene_ids(manifests)
    if args.alignment_tasks:
        alignment_task_count, gene_count = summarize_alignment_tasks(args.alignment_tasks)
    else:
        alignment_task_count = len(gene_ids)
        gene_count = len(gene_ids)

    summary_count = merge_tsv_gz(
        [path / "ortholog_alignment_summary.tsv.gz" for path in result_dirs],
        args.outdir / "ortholog_alignment_summary.tsv.gz",
    )
    strategy_summary_count = write_strategy_summary(
        args.outdir / "ortholog_alignment_summary.tsv.gz",
        args.outdir / "strategy_summary.tsv.gz",
        strategies,
    )
    segment_count = merge_tsv_gz(
        [path / "alignment_segments.tsv.gz" for path in result_dirs],
        args.outdir / "alignment_segments.tsv.gz",
    )
    feature_coverage_inputs = [path / "feature_coverage.tsv.gz" for path in result_dirs]
    missing_feature_coverage = [str(path.parent) for path in feature_coverage_inputs if not path.exists()]
    if missing_feature_coverage:
        if args.partition_id:
            raise FileNotFoundError(
                "Alignment partition inputs missing feature_coverage.tsv.gz: "
                + ", ".join(missing_feature_coverage)
            )
        feature_coverage_mode = "global_fallback"
        feature_coverage_count = summarize_feature_coverage(
            args.target_features,
            args.outdir / "ortholog_alignment_summary.tsv.gz",
            args.outdir / "alignment_segments.tsv.gz",
            args.outdir / "feature_coverage.tsv.gz",
        )
    else:
        feature_coverage_mode = "per_result_merge"
        feature_coverage_count = merge_tsv_gz(
            feature_coverage_inputs,
            args.outdir / "feature_coverage.tsv.gz",
        )
    event_inputs = [path / "alignment_events.tsv.gz" for path in result_dirs]
    if args.compact_events and not args.events_already_compacted:
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
        if args.events_already_compacted:
            raw_event_count = sum(
                int(manifest.get("raw_alignment_event_count") or manifest.get("alignment_event_count") or 0)
                for manifest in manifests
            )
    failure_count = merge_tsv_gz(
        [path / "failures.tsv.gz" for path in result_dirs],
        args.outdir / "failures.tsv.gz",
    )
    native_file_count = copy_native(result_dirs, args.outdir)
    manifest = {
        "created_at": utc_now(),
        "stage": "alignment",
        "partition_id": args.partition_id or "",
        "strategy_count": len(strategies),
        "strategies": strategies,
        "gene_count": gene_count,
        "alignment_task_count": alignment_task_count,
        "taxonomy_tax_id_count": count_tsv_gz_rows(args.taxonomy_presets) if args.taxonomy_presets else 0,
        "taxonomy_failure_count": count_tsv_gz_rows(args.taxonomy_failures) if args.taxonomy_failures else 0,
        "ortholog_alignment_summary_count": summary_count,
        "strategy_summary_count": strategy_summary_count,
        "alignment_segment_count": segment_count,
        "feature_coverage_mode": feature_coverage_mode,
        "feature_coverage_missing_result_count": len(missing_feature_coverage),
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
