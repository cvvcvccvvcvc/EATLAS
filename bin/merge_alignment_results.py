#!/usr/bin/env python3
"""Merge per-gene alignment evidence outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import shutil
import sqlite3
import sys
import time
from datetime import datetime, timezone
from itertools import groupby
from pathlib import Path

from feature_coverage import (
    iter_snv_event_sites,
    write_snv_site_depth,
    write_snv_taxonomic_depth,
)
from taxonomic_evidence import COUNT_KEYS, count_member_groups, load_taxonomy_profiles


csv.field_size_limit(sys.maxsize)


logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-tasks", required=True, type=Path)
    parser.add_argument("--taxonomy", type=Path)
    parser.add_argument("--taxonomy-failures", type=Path)
    parser.add_argument("--taxonomy-summary", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--result-dir", action="append", default=[], type=Path)
    parser.add_argument("--result-root", type=Path)
    parser.add_argument("--partition-id")
    parser.add_argument("--expected-strategies", required=True)
    parser.add_argument("--expected-gene-ids")
    parser.add_argument("--compact-events", action="store_true")
    parser.add_argument(
        "--output-profile",
        choices=["full", "annotation-input", "report-input"],
        default="full",
        help="Select full outputs, partitioned annotation inputs, or final report inputs.",
    )
    return parser.parse_args()


COMPACT_EVENT_FIELDS = [
    "event_group_id",
    "gene_id",
    "event_type",
    "target_start0",
    "target_end0",
    "genomic_accession",
    "genomic_start1",
    "genomic_end1",
    "ref",
    "alt",
    "strategy",
    "support_row_count",
    "support_ortholog_count",
    "tools",
    "presets",
    "tax_id_count",
    "taxname_count",
    "qc_flags",
]
EVENT_ORTHOLOG_SUPPORT_FIELDS = [
    "event_group_id",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "support_row_count",
]
EVENT_KEY_FIELDS = [
    "gene_id",
    "strategy",
    "target_start0",
    "event_type",
    "target_end0",
    "genomic_accession",
    "genomic_start1",
    "genomic_end1",
    "ref",
    "alt",
]
EVENT_STREAM_FIELDS = [
    *EVENT_KEY_FIELDS,
    "ortholog_gene_id",
    "tool",
    "preset",
    "tax_id",
    "taxname",
    "qc_flags",
]
SNV_ALT_TAXONOMIC_SUPPORT_FIELDS = [
    "gene_id",
    "strategy",
    "target_start0",
    "ref",
    "alt",
    *COUNT_KEYS,
]

STRATEGY_SUMMARY_FIELDS = [
    "strategy",
    "summary_row_count",
    "gene_count",
    "aligned_summary_row_count",
    "event_count",
    "aligned_target_bp",
]

ENSEMBL_COMPARA_STRATEGY = "precomputed_ensembl_92_mammals_epo_extended"
DNA_BASES = frozenset("ACGT")
ALIGNMENT_MANIFEST_COUNT_FIELDS = (
    "ortholog_alignment_summary_count",
    "alignment_segment_count",
    "feature_coverage_count",
    "raw_alignment_event_count",
    "alignment_event_count",
    "failure_count",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def start_phase(name: str) -> float:
    logger.info("Starting phase %s", name)
    return time.perf_counter()


def finish_phase(timings: dict[str, float], name: str, started_at: float) -> None:
    elapsed = round(time.perf_counter() - started_at, 3)
    timings[name] = elapsed
    logger.info("Timing %s: %.3f seconds", name, elapsed)


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


def gene_sort_key(value: str) -> tuple[int, int | str]:
    return (0, int(value)) if value.isdigit() else (1, value)


def sorted_gene_ids(values: set[str]) -> list[str]:
    return sorted(values, key=gene_sort_key)


def parse_expected_values(raw: str | None, label: str) -> list[str]:
    values = [value.strip() for value in (raw or "").split(",") if value.strip()]
    if not values:
        raise ValueError(f"{label} must contain at least one value")
    if len(values) != len(set(values)):
        raise ValueError(f"{label} contains duplicate values")
    return values


def read_alignment_capabilities(
    path: Path,
) -> tuple[int, dict[str, tuple[bool, bool]]]:
    """Return target/ortholog readiness for every alignment task gene."""
    task_count = 0
    capabilities: dict[str, tuple[bool, bool]] = {}
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", "target_ready", "ortholog_ready"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Alignment tasks {path} missing required columns: "
                + ", ".join(sorted(missing))
            )
        for row in reader:
            task_count += 1
            gene_id = row.get("gene_id", "")
            if not gene_id:
                continue
            if gene_id in capabilities:
                raise ValueError(f"Alignment tasks {path} contain duplicate gene_id {gene_id}")
            capabilities[gene_id] = (
                str(row["target_ready"]).lower() == "true",
                str(row["ortholog_ready"]).lower() == "true",
            )
    return task_count, capabilities


def expected_gene_strategy_pairs(
    capabilities: dict[str, tuple[bool, bool]],
    gene_ids: list[str],
    strategies: list[str],
) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    missing_genes = set(gene_ids) - set(capabilities)
    if missing_genes:
        raise ValueError(
            "Alignment tasks are missing expected genes: "
            + ", ".join(sorted_gene_ids(missing_genes))
        )
    for gene_id in gene_ids:
        target_ready, ortholog_ready = capabilities[gene_id]
        for strategy in strategies:
            eligible = target_ready if strategy == ENSEMBL_COMPARA_STRATEGY else ortholog_ready
            if eligible:
                pairs.add((gene_id, strategy))
    return pairs


def summarize_alignment_tasks(
    path: Path,
    expected_strategies: list[str],
) -> tuple[int, list[str], dict[str, int]]:
    """Return task count, selected-strategy gene union, and per-strategy eligibility."""
    task_count, capabilities = read_alignment_capabilities(path)
    gene_ids = sorted_gene_ids(set(capabilities))
    pairs = expected_gene_strategy_pairs(capabilities, gene_ids, expected_strategies)
    eligible_gene_ids = sorted_gene_ids({gene_id for gene_id, _strategy in pairs})
    strategy_counts = {
        strategy: sum(pair_strategy == strategy for _gene_id, pair_strategy in pairs)
        for strategy in expected_strategies
    }
    return task_count, eligible_gene_ids, strategy_counts


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
                raise FileNotFoundError(f"Missing required alignment table: {path}")
            with gzip.open(path, "rt", newline="") as in_handle:
                reader = csv.reader(in_handle, delimiter="\t")
                header = next(reader, None)
                if header is None:
                    raise ValueError(f"Alignment table has no TSV header: {path}")
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


def merge_compact_event_handoffs(
    result_dirs: list[Path],
    events_output: Path,
    support_output: Path,
) -> tuple[int, int]:
    """Concatenate compact partitions while rebasing partition-local group IDs."""

    event_count = 0
    support_count = 0
    with (
        gzip.open(events_output, "wt", newline="") as event_handle,
        gzip.open(support_output, "wt", newline="") as support_handle,
    ):
        event_writer = csv.DictWriter(
            event_handle,
            fieldnames=COMPACT_EVENT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        support_writer = csv.DictWriter(
            support_handle,
            fieldnames=EVENT_ORTHOLOG_SUPPORT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        event_writer.writeheader()
        support_writer.writeheader()
        for result_dir in result_dirs:
            events_path = result_dir / "alignment_events.tsv.gz"
            support_path = result_dir / "event_ortholog_support.tsv.gz"
            with (
                gzip.open(events_path, "rt", newline="") as partition_event_handle,
                gzip.open(support_path, "rt", newline="") as partition_support_handle,
            ):
                event_reader = csv.DictReader(partition_event_handle, delimiter="\t")
                support_reader = csv.DictReader(partition_support_handle, delimiter="\t")
                event_missing = set(COMPACT_EVENT_FIELDS) - set(event_reader.fieldnames or [])
                support_missing = set(EVENT_ORTHOLOG_SUPPORT_FIELDS) - set(
                    support_reader.fieldnames or []
                )
                if event_missing:
                    raise ValueError(
                        f"Compact events {events_path} missing columns: "
                        + ", ".join(sorted(event_missing))
                    )
                if support_missing:
                    raise ValueError(
                        f"Event ortholog support {support_path} missing columns: "
                        + ", ".join(sorted(support_missing))
                    )
                current_support = next(support_reader, None)
                local_event_count = 0
                for event_row in event_reader:
                    local_event_count += 1
                    local_group_id = int(event_row["event_group_id"])
                    if local_group_id != local_event_count:
                        raise ValueError(
                            f"Compact event_group_id values in {events_path} must be "
                            f"consecutive from 1; expected {local_event_count}, "
                            f"observed {local_group_id}"
                        )
                    event_count += 1
                    event_row["event_group_id"] = str(event_count)
                    event_writer.writerow(event_row)
                    while current_support is not None:
                        support_group_id = int(current_support["event_group_id"])
                        if support_group_id < local_group_id:
                            raise ValueError(
                                f"Unmatched event ortholog support group {support_group_id} "
                                f"in {support_path}"
                            )
                        if support_group_id > local_group_id:
                            break
                        current_support["event_group_id"] = str(event_count)
                        support_writer.writerow(current_support)
                        support_count += 1
                        current_support = next(support_reader, None)
                if current_support is not None:
                    raise ValueError(
                        f"Unmatched event ortholog support group "
                        f"{current_support['event_group_id']} in {support_path}"
                    )
    return event_count, support_count


def write_strategy_summary(
    summary_paths: list[Path],
    output: Path,
    expected_strategies: list[str],
) -> tuple[int, int]:
    """Write small per-strategy aggregates from the canonical summary table."""
    aggregates: dict[str, dict[str, int | set[str]]] = {}
    summary_row_count = 0
    for summaries_path in summary_paths:
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
                summary_row_count += 1
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
    return summary_row_count, len(aggregates)


def merge_strategy_summaries(
    paths: list[Path],
    output: Path,
    expected_strategies: list[str],
) -> int:
    aggregates: dict[str, dict[str, int]] = {}
    numeric_fields = [field for field in STRATEGY_SUMMARY_FIELDS if field != "strategy"]
    for path in paths:
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(STRATEGY_SUMMARY_FIELDS) - set(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"Strategy summary {path} missing required columns: "
                    + ", ".join(sorted(missing))
                )
            for row in reader:
                strategy = row["strategy"]
                aggregate = aggregates.setdefault(
                    strategy,
                    {field: 0 for field in numeric_fields},
                )
                for field in numeric_fields:
                    aggregate[field] += int(row[field] or 0)

    unexpected = sorted(set(aggregates) - set(expected_strategies))
    if unexpected:
        raise ValueError(f"Unexpected strategies in partition summaries: {unexpected}")
    for strategy in expected_strategies:
        aggregates.setdefault(strategy, {field: 0 for field in numeric_fields})

    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STRATEGY_SUMMARY_FIELDS, delimiter="\t")
        writer.writeheader()
        for strategy in sorted(aggregates):
            writer.writerow({"strategy": strategy, **aggregates[strategy]})
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
        raise FileNotFoundError(f"Missing required alignment table: {path}")
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


def merge_ortholog_metadata(
    orthologs: dict[str, dict[str, object]],
    row: dict[str, str],
) -> None:
    ortholog_gene_id = row["ortholog_gene_id"]
    if not ortholog_gene_id:
        return
    support = orthologs.get(ortholog_gene_id)
    if support is None:
        orthologs[ortholog_gene_id] = {
            "ortholog_gene_id": ortholog_gene_id,
            "tax_id": row["tax_id"],
            "taxname": row["taxname"],
            "support_row_count": 1,
        }
        return

    for field in ("tax_id", "taxname"):
        current = str(support[field] or "")
        observed = row[field]
        if current and observed and current != observed:
            raise ValueError(
                f"Conflicting {field} for ortholog_gene_id={ortholog_gene_id}: "
                f"{current!r} != {observed!r}"
            )
        if not current and observed:
            support[field] = observed
    support["support_row_count"] = int(support["support_row_count"]) + 1


def compact_stream_group(
    event_group_id: int,
    key: tuple[str, ...],
    rows: object,
    taxonomy_profiles: dict | None,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    record: dict[str, object] = {
        "event_group_id": event_group_id,
        **dict(zip(EVENT_KEY_FIELDS, key)),
    }
    orthologs: dict[str, dict[str, object]] = {}
    tools: set[str] = set()
    presets: set[str] = set()
    tax_ids: set[str] = set()
    taxnames: set[str] = set()
    qc_flags: list[str] = []
    support_row_count = 0
    for values in rows:
        row = dict(zip(EVENT_STREAM_FIELDS, values))
        support_row_count += 1
        merge_ortholog_metadata(orthologs, row)
        if row["tool"]:
            tools.add(row["tool"])
        if row["preset"]:
            presets.add(row["preset"])
        if row["tax_id"]:
            tax_ids.add(row["tax_id"])
        if row["taxname"]:
            taxnames.add(row["taxname"])
        if row["qc_flags"]:
            qc_flags.append(row["qc_flags"])

    record.update(
        {
            "support_row_count": support_row_count,
            "support_ortholog_count": len(orthologs),
            "tools": ",".join(sorted(tools)),
            "presets": ",".join(sorted(presets)),
            "tax_id_count": len(tax_ids),
            "taxname_count": len(taxnames),
            "qc_flags": compact_event_flags(",".join(qc_flags)),
        }
    )
    if taxonomy_profiles is not None:
        record.update(
            count_member_groups(
                (
                    (ortholog_gene_id, str(support["tax_id"]))
                    for ortholog_gene_id, support in orthologs.items()
                    if support["tax_id"]
                ),
                taxonomy_profiles,
            )
        )
    support_rows = [orthologs[key] for key in sorted(orthologs)]
    return record, support_rows


def write_compact_events(
    paths: list[Path],
    output: Path,
    ortholog_support_output: Path,
    taxonomy: Path | None = None,
    taxonomic_support_output: Path | None = None,
    timings: dict[str, float] | None = None,
) -> tuple[int, int, int, int]:
    if (taxonomy is None) != (taxonomic_support_output is None):
        raise ValueError("Taxonomic event support requires both taxonomy metadata and an output path")
    output.parent.mkdir(parents=True, exist_ok=True)
    db_path = output.parent / "alignment_event_support.sqlite"
    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(db_path)
    phase_timings = timings if timings is not None else {}
    try:
        phase_started = start_phase("load_events_sqlite")
        conn.execute("PRAGMA journal_mode = OFF")
        conn.execute("PRAGMA synchronous = OFF")
        conn.execute("PRAGMA temp_store = MEMORY")
        create_compact_event_table(conn)
        raw_count = 0
        for path in paths:
            raw_count += insert_event_rows(conn, path)
        conn.commit()
        finish_phase(phase_timings, "load_events_sqlite", phase_started)

        phase_started = start_phase("build_event_index")
        conn.execute(
            """
            CREATE INDEX events_key_idx ON events (
                gene_id,
                strategy,
                CAST(target_start0 AS INTEGER),
                target_start0,
                event_type,
                target_end0,
                genomic_accession,
                genomic_start1,
                genomic_end1,
                ref,
                alt
            )
            """
        )
        finish_phase(phase_timings, "build_event_index", phase_started)

        query = """
            SELECT
                gene_id,
                strategy,
                target_start0,
                event_type,
                target_end0,
                genomic_accession,
                genomic_start1,
                genomic_end1,
                ref,
                alt,
                ortholog_gene_id,
                tool,
                preset,
                tax_id,
                taxname,
                qc_flags
            FROM events INDEXED BY events_key_idx
            ORDER BY
                gene_id,
                strategy,
                CAST(target_start0 AS INTEGER),
                target_start0,
                event_type,
                target_end0,
                genomic_accession,
                genomic_start1,
                genomic_end1,
                ref,
                alt
        """
        taxonomy_profiles = (
            load_taxonomy_profiles(taxonomy)
            if taxonomy is not None
            else None
        )
        compact_count = 0
        taxonomic_support_count = 0
        ortholog_support_count = 0
        phase_started = start_phase("stream_event_groups")
        support_handle = (
            gzip.open(taxonomic_support_output, "wt", newline="")
            if taxonomic_support_output is not None
            else None
        )
        try:
            support_writer = None
            if support_handle is not None:
                support_writer = csv.DictWriter(
                    support_handle,
                    fieldnames=SNV_ALT_TAXONOMIC_SUPPORT_FIELDS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                support_writer.writeheader()
            with (
                gzip.open(output, "wt", newline="") as handle,
                gzip.open(ortholog_support_output, "wt", newline="") as ortholog_handle,
            ):
                writer = csv.DictWriter(
                    handle,
                    fieldnames=COMPACT_EVENT_FIELDS,
                    delimiter="\t",
                    extrasaction="ignore",
                    lineterminator="\n",
                )
                ortholog_writer = csv.DictWriter(
                    ortholog_handle,
                    fieldnames=EVENT_ORTHOLOG_SUPPORT_FIELDS,
                    delimiter="\t",
                    extrasaction="ignore",
                    lineterminator="\n",
                )
                writer.writeheader()
                ortholog_writer.writeheader()
                event_rows = conn.execute(query)
                for key, rows in groupby(
                    event_rows,
                    key=lambda values: values[: len(EVENT_KEY_FIELDS)],
                ):
                    compact_count += 1
                    record, ortholog_rows = compact_stream_group(
                        compact_count,
                        key,
                        rows,
                        taxonomy_profiles,
                    )
                    writer.writerow(record)
                    for ortholog_row in ortholog_rows:
                        ortholog_writer.writerow(
                            {"event_group_id": compact_count, **ortholog_row}
                        )
                        ortholog_support_count += 1
                    if (
                        support_writer is None
                        or record.get("event_type") != "snv"
                        or len(str(record.get("ref") or "")) != 1
                        or len(str(record.get("alt") or "")) != 1
                        or str(record.get("ref") or "").upper() not in DNA_BASES
                        or str(record.get("alt") or "").upper() not in DNA_BASES
                        or int(record.get("all__ortholog") or 0) < 1
                    ):
                        continue
                    support_writer.writerow(
                        {
                            "gene_id": record["gene_id"],
                            "strategy": record["strategy"],
                            "target_start0": record["target_start0"],
                            "ref": str(record["ref"]).upper(),
                            "alt": str(record["alt"]).upper(),
                            **{key: int(record[key]) for key in COUNT_KEYS},
                        }
                    )
                    taxonomic_support_count += 1
        finally:
            if support_handle is not None:
                support_handle.close()
        finish_phase(phase_timings, "stream_event_groups", phase_started)
        return compact_count, raw_count, taxonomic_support_count, ortholog_support_count
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
    return [json.loads((result_dir / "manifest.json").read_text()) for result_dir in result_dirs]


def manifest_values(manifest: dict, field: str) -> list[str]:
    values = manifest.get(field)
    if (
        not isinstance(values, list)
        or not values
        or not all(isinstance(value, str) and value for value in values)
        or len(values) != len(set(values))
    ):
        raise ValueError(f"Alignment manifest has invalid {field}: {values!r}")
    return values


def manifest_strategies(manifests: list[dict]) -> list[str]:
    return sorted(
        {
            strategy
            for manifest in manifests
            for strategy in manifest_values(manifest, "strategies")
        }
    )


def manifest_gene_ids(manifests: list[dict]) -> list[str]:
    return sorted_gene_ids(
        {
            gene_id
            for manifest in manifests
            for gene_id in manifest_values(manifest, "gene_ids")
        }
    )


def merge_strategy_parameters(manifests: list[dict]) -> dict[str, dict]:
    merged: dict[str, dict] = {}
    for manifest in manifests:
        parameters = manifest.get("strategy_parameters")
        if not isinstance(parameters, dict):
            raise ValueError("strategy_parameters must be a JSON object")
        declared_strategies = set(manifest_values(manifest, "strategies"))
        if set(parameters) != declared_strategies:
            raise ValueError(
                "strategy_parameters keys must match manifest strategies: "
                f"strategies={sorted(declared_strategies)}, parameters={sorted(parameters)}"
            )
        for strategy, values in parameters.items():
            if not isinstance(values, dict):
                raise ValueError(f"Strategy parameters for {strategy} must be a JSON object")
            normalized = dict(sorted(values.items()))
            existing = merged.get(strategy)
            if existing is not None and existing != normalized:
                raise ValueError(f"Inconsistent strategy parameters for {strategy}")
            merged[strategy] = normalized

    return dict(sorted(merged.items()))


def require_alignment_tables(
    result_dirs: list[Path],
    output_profile: str,
) -> None:
    if output_profile == "report-input":
        filenames = [
            "strategy_summary.tsv.gz",
            "feature_coverage.tsv.gz",
            "failures.tsv.gz",
        ]
    else:
        filenames = [
            "ortholog_alignment_summary.tsv.gz",
            "alignment_segments.tsv.gz",
            "alignment_events.tsv.gz",
            "feature_coverage.tsv.gz",
            "failures.tsv.gz",
        ]
    missing = [
        str(result_dir / filename)
        for result_dir in result_dirs
        for filename in filenames
        if not (result_dir / filename).exists()
    ]
    if missing:
        raise FileNotFoundError("Missing required alignment table(s): " + ", ".join(missing))


def preview_pairs(pairs: set[tuple[str, str]], limit: int = 10) -> str:
    values = [f"{gene_id}:{strategy}" for gene_id, strategy in sorted(pairs)]
    if len(values) > limit:
        values = values[:limit] + [f"... ({len(pairs) - limit} more)"]
    return ", ".join(values)


def validate_partition_manifests(
    result_dirs: list[Path],
    manifests: list[dict],
    expected_gene_ids: list[str],
    expected_strategies: list[str],
    alignment_tasks: Path,
) -> tuple[list[str], list[str]]:
    _task_count, capabilities = read_alignment_capabilities(alignment_tasks)
    expected_pairs = expected_gene_strategy_pairs(
        capabilities,
        expected_gene_ids,
        expected_strategies,
    )
    genes_without_strategy = set(expected_gene_ids) - {
        gene_id for gene_id, _strategy in expected_pairs
    }
    if genes_without_strategy:
        raise ValueError(
            "Expected genes are not eligible for any selected strategy: "
            + ", ".join(sorted_gene_ids(genes_without_strategy))
        )
    pair_owners: dict[tuple[str, str], Path] = {}
    duplicates: list[str] = []

    for result_dir, manifest in zip(result_dirs, manifests):
        gene_ids = manifest_gene_ids([manifest])
        if len(gene_ids) != 1:
            raise ValueError(
                f"Per-gene alignment result {result_dir} must declare exactly one gene_id; "
                f"observed {gene_ids}"
            )
        strategies = manifest_strategies([manifest])
        if not strategies:
            raise ValueError(f"Alignment result {result_dir} does not declare a strategy")
        gene_id = gene_ids[0]
        for strategy in strategies:
            pair = (gene_id, strategy)
            if pair in pair_owners:
                duplicates.append(f"{gene_id}:{strategy} ({pair_owners[pair].name}, {result_dir.name})")
            else:
                pair_owners[pair] = result_dir

    if duplicates:
        raise ValueError("Duplicate gene-strategy alignment result(s): " + ", ".join(duplicates))

    observed_pairs = set(pair_owners)
    missing_pairs = expected_pairs - observed_pairs
    unexpected_pairs = observed_pairs - expected_pairs
    if missing_pairs or unexpected_pairs:
        details = []
        if missing_pairs:
            details.append("missing=" + preview_pairs(missing_pairs))
        if unexpected_pairs:
            details.append("unexpected=" + preview_pairs(unexpected_pairs))
        raise ValueError("Partition alignment coverage mismatch: " + "; ".join(details))

    observed_strategies = {strategy for _gene_id, strategy in observed_pairs}
    return sorted_gene_ids(set(expected_gene_ids)), sorted(observed_strategies)


def validate_final_manifests(
    result_dirs: list[Path],
    manifests: list[dict],
    expected_gene_ids: list[str],
    expected_strategies: list[str],
    alignment_tasks: Path,
) -> tuple[list[str], list[str]]:
    expected_strategy_set = set(expected_strategies)
    observed_strategy_set: set[str] = set()
    gene_owners: dict[str, str] = {}
    partition_ids: set[str] = set()
    _task_count, capabilities = read_alignment_capabilities(alignment_tasks)

    for result_dir, manifest in zip(result_dirs, manifests):
        partition_id = str(manifest.get("partition_id") or "")
        if not partition_id:
            raise ValueError(f"Alignment partition {result_dir} does not declare partition_id")
        if partition_id in partition_ids:
            raise ValueError(f"Duplicate alignment partition_id: {partition_id}")
        partition_ids.add(partition_id)

        gene_ids = manifest_gene_ids([manifest])
        if not gene_ids:
            raise ValueError(f"Alignment partition {partition_id} does not declare gene_ids")
        declared_gene_count = int(manifest.get("gene_count") or 0)
        if declared_gene_count != len(gene_ids):
            raise ValueError(
                f"Alignment partition {partition_id} gene_count mismatch: "
                f"declared {declared_gene_count}, observed {len(gene_ids)} gene_ids"
            )

        strategies = set(manifest_strategies([manifest]))
        expected_partition_strategies = {
            strategy
            for _gene_id, strategy in expected_gene_strategy_pairs(
                capabilities,
                gene_ids,
                expected_strategies,
            )
        }
        if strategies != expected_partition_strategies:
            missing = sorted(expected_partition_strategies - strategies)
            unexpected = sorted(strategies - expected_partition_strategies)
            raise ValueError(
                f"Alignment partition {partition_id} strategy mismatch: "
                f"missing={missing}; unexpected={unexpected}"
            )
        observed_strategy_set.update(strategies)

        for gene_id in gene_ids:
            if gene_id in gene_owners:
                raise ValueError(
                    f"Gene {gene_id} occurs in multiple alignment partitions: "
                    f"{gene_owners[gene_id]} and {partition_id}"
                )
            gene_owners[gene_id] = partition_id

    expected_gene_set = set(expected_gene_ids)
    observed_gene_set = set(gene_owners)
    if observed_gene_set != expected_gene_set:
        missing = sorted_gene_ids(expected_gene_set - observed_gene_set)
        unexpected = sorted_gene_ids(observed_gene_set - expected_gene_set)
        raise ValueError(
            "Final alignment gene coverage mismatch: "
            f"missing={missing}; unexpected={unexpected}"
        )

    if observed_strategy_set != expected_strategy_set:
        missing = sorted(expected_strategy_set - observed_strategy_set)
        unexpected = sorted(observed_strategy_set - expected_strategy_set)
        raise ValueError(
            "Final alignment strategy coverage mismatch: "
            f"missing={missing}; unexpected={unexpected}"
        )

    return sorted_gene_ids(observed_gene_set), sorted(observed_strategy_set)


def manifest_count(manifest: dict, field: str) -> int:
    if field not in manifest:
        raise ValueError(f"Alignment manifest missing required count: {field}")
    try:
        value = int(manifest[field])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Alignment manifest has invalid {field}: {manifest[field]!r}") from exc
    if value < 0:
        raise ValueError(f"Alignment manifest has negative {field}: {value}")
    return value


def sum_manifest_count(manifests: list[dict], field: str) -> int:
    total = 0
    for manifest in manifests:
        total += manifest_count(manifest, field)
    return total


def validate_alignment_manifests(manifests: list[dict]) -> None:
    for manifest in manifests:
        manifest_values(manifest, "gene_ids")
        manifest_values(manifest, "strategies")
        for field in ALIGNMENT_MANIFEST_COUNT_FIELDS:
            manifest_count(manifest, field)
        mode = manifest.get("alignment_event_mode")
        if mode not in {"raw", "compact_support"}:
            raise ValueError(f"Alignment manifest has invalid alignment_event_mode: {mode!r}")
        raw_count = manifest_count(manifest, "raw_alignment_event_count")
        event_count = manifest_count(manifest, "alignment_event_count")
        if raw_count < event_count:
            raise ValueError(
                "Alignment manifest raw_alignment_event_count cannot be smaller than "
                "alignment_event_count"
            )


def merged_event_mode(manifests: list[dict]) -> str:
    modes = {str(manifest.get("alignment_event_mode") or "") for manifest in manifests}
    if modes - {"raw", "compact_support"} or len(modes) != 1:
        raise ValueError(f"Alignment manifests have inconsistent event modes: {sorted(modes)}")
    return next(iter(modes))


def collect_partition_timings(manifests: list[dict]) -> dict[str, dict[str, float]]:
    collected = {}
    for manifest in manifests:
        partition_id = str(manifest.get("partition_id") or "")
        timings = manifest.get("timings_seconds")
        if partition_id and isinstance(timings, dict):
            collected[partition_id] = {
                str(name): float(value)
                for name, value in timings.items()
            }
    return dict(sorted(collected.items()))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    timings_seconds: dict[str, float] = {}
    if args.partition_id and args.output_profile == "report-input":
        raise ValueError("--output-profile report-input is only valid for the final merge")
    if not args.partition_id and args.output_profile == "annotation-input":
        raise ValueError("--output-profile annotation-input requires --partition-id")
    expected_strategies = sorted(
        parse_expected_values(args.expected_strategies, "--expected-strategies")
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    result_dirs = resolve_result_dirs(args.result_dir, args.result_root)

    global_inputs = [args.taxonomy, args.taxonomy_failures]
    if not args.partition_id and any(path is None for path in global_inputs):
        raise ValueError(
            "Final alignment merge requires --taxonomy and --taxonomy-failures"
        )

    manifests = load_manifests(result_dirs)
    validate_alignment_manifests(manifests)
    strategy_parameters = merge_strategy_parameters(manifests)
    require_alignment_tables(
        result_dirs,
        output_profile=args.output_profile,
    )
    if args.partition_id:
        expected_gene_ids = parse_expected_values(args.expected_gene_ids, "--expected-gene-ids")
        gene_ids, strategies = validate_partition_manifests(
            result_dirs,
            manifests,
            expected_gene_ids,
            expected_strategies,
            args.alignment_tasks,
        )
        alignment_task_count = len(gene_ids)
        strategy_eligible_gene_counts = {
            strategy: len(
                {
                    gene_id
                    for manifest in manifests
                    if strategy in manifest_strategies([manifest])
                    for gene_id in manifest_gene_ids([manifest])
                }
            )
            for strategy in strategies
        }
    else:
        (
            alignment_task_count,
            expected_gene_ids,
            strategy_eligible_gene_counts,
        ) = summarize_alignment_tasks(args.alignment_tasks, expected_strategies)
        if not expected_gene_ids:
            raise ValueError(
                f"Alignment tasks {args.alignment_tasks} contain no genes eligible "
                "for the selected strategies"
            )
        gene_ids, strategies = validate_final_manifests(
            result_dirs,
            manifests,
            expected_gene_ids,
            expected_strategies,
            args.alignment_tasks,
        )
        if args.output_profile == "full":
            copy_or_keep(args.alignment_tasks, args.outdir / "alignment_tasks.tsv.gz")
            copy_or_keep(args.taxonomy, args.outdir / "taxonomy.tsv.gz")
            copy_or_keep(args.taxonomy_failures, args.outdir / "taxonomy_failures.tsv.gz")
        if args.taxonomy_summary is not None:
            copy_or_keep(args.taxonomy_summary, args.outdir / "taxonomy_summary.tsv.gz")
    gene_count = len(gene_ids)

    input_event_mode = merged_event_mode(manifests)
    if args.output_profile == "report-input":
        summary_count = sum_manifest_count(manifests, "ortholog_alignment_summary_count")
        strategy_summary_count = merge_strategy_summaries(
            [path / "strategy_summary.tsv.gz" for path in result_dirs],
            args.outdir / "strategy_summary.tsv.gz",
            strategies,
        )
        segment_count = sum_manifest_count(manifests, "alignment_segment_count")
    else:
        summary_inputs = [
            path / "ortholog_alignment_summary.tsv.gz"
            for path in result_dirs
        ]
        if args.output_profile == "full":
            summary_count = merge_tsv_gz(
                summary_inputs,
                args.outdir / "ortholog_alignment_summary.tsv.gz",
            )
            _, strategy_summary_count = write_strategy_summary(
                [args.outdir / "ortholog_alignment_summary.tsv.gz"],
                args.outdir / "strategy_summary.tsv.gz",
                strategies,
            )
        else:
            summary_count, strategy_summary_count = write_strategy_summary(
                summary_inputs,
                args.outdir / "strategy_summary.tsv.gz",
                strategies,
            )
        if args.output_profile == "full":
            segment_count = merge_tsv_gz(
                [path / "alignment_segments.tsv.gz" for path in result_dirs],
                args.outdir / "alignment_segments.tsv.gz",
            )
        else:
            segment_count = sum_manifest_count(manifests, "alignment_segment_count")

    feature_coverage_inputs = [path / "feature_coverage.tsv.gz" for path in result_dirs]
    feature_coverage_count = merge_tsv_gz(
        feature_coverage_inputs,
        args.outdir / "feature_coverage.tsv.gz",
    )
    if args.output_profile == "report-input":
        event_count = sum_manifest_count(manifests, "alignment_event_count")
        raw_event_count = sum_manifest_count(manifests, "raw_alignment_event_count")
        taxonomic_alt_support_count = sum_manifest_count(
            manifests,
            "snv_alt_taxonomic_support_count",
        )
        event_ortholog_support_count = sum_manifest_count(
            manifests,
            "event_ortholog_support_count",
        )
        alignment_event_mode = input_event_mode
    else:
        event_inputs = [path / "alignment_events.tsv.gz" for path in result_dirs]
        if args.compact_events and input_event_mode == "raw":
            (
                event_count,
                raw_event_count,
                taxonomic_alt_support_count,
                event_ortholog_support_count,
            ) = write_compact_events(
                event_inputs,
                args.outdir / "alignment_events.tsv.gz",
                args.outdir / "event_ortholog_support.tsv.gz",
                args.taxonomy if args.partition_id and args.taxonomy else None,
                args.outdir / "snv_alt_taxonomic_support.tsv.gz"
                if args.partition_id and args.taxonomy
                else None,
                timings=timings_seconds,
            )
        elif args.compact_events and input_event_mode == "compact_support":
            raw_event_count = sum_manifest_count(manifests, "raw_alignment_event_count")
            event_count, event_ortholog_support_count = merge_compact_event_handoffs(
                result_dirs,
                args.outdir / "alignment_events.tsv.gz",
                args.outdir / "event_ortholog_support.tsv.gz",
            )
            taxonomic_alt_support_count = 0
        else:
            if input_event_mode != "raw":
                raise ValueError(
                    "Compact alignment inputs require --compact-events in the final merge"
                )
            event_count = merge_tsv_gz(
                event_inputs,
                args.outdir / "alignment_events.tsv.gz",
            )
            raw_event_count = event_count
            taxonomic_alt_support_count = 0
            event_ortholog_support_count = 0
        alignment_event_mode = "compact_support" if args.compact_events else "raw"

    if args.partition_id:
        phase_started = start_phase("snv_site_depth")
        snv_site_depth_count = write_snv_site_depth(
            [path / "alignment_segments.tsv.gz" for path in result_dirs],
            iter_snv_event_sites([args.outdir / "alignment_events.tsv.gz"]),
            args.outdir / "snv_site_depth.tsv.gz",
            args.outdir,
        )
        finish_phase(timings_seconds, "snv_site_depth", phase_started)
        if args.taxonomy is not None:
            phase_started = start_phase("snv_taxonomic_depth")
            snv_taxonomic_depth_count = write_snv_taxonomic_depth(
                [path / "alignment_segments.tsv.gz" for path in result_dirs],
                iter_snv_event_sites([args.outdir / "alignment_events.tsv.gz"]),
                args.taxonomy,
                args.outdir / "snv_taxonomic_depth.tsv.gz",
                args.outdir,
            )
            finish_phase(timings_seconds, "snv_taxonomic_depth", phase_started)
        else:
            snv_taxonomic_depth_count = 0
    elif args.output_profile == "full":
        snv_site_depth_count = merge_tsv_gz(
            [path / "snv_site_depth.tsv.gz" for path in result_dirs],
            args.outdir / "snv_site_depth.tsv.gz",
        )
        snv_taxonomic_depth_count = sum_manifest_count(
            manifests,
            "snv_taxonomic_depth_count",
        )
    else:
        snv_site_depth_count = sum_manifest_count(manifests, "snv_site_depth_count")
        snv_taxonomic_depth_count = sum_manifest_count(
            manifests,
            "snv_taxonomic_depth_count",
        )

    failure_count = merge_tsv_gz(
        [path / "failures.tsv.gz" for path in result_dirs],
        args.outdir / "failures.tsv.gz",
    )
    native_file_count = (
        copy_native(result_dirs, args.outdir)
        if args.output_profile == "full"
        else 0
    )
    manifest = {
        "created_at": utc_now(),
        "stage": "alignment",
        "partition_id": args.partition_id or "",
        "output_profile": args.output_profile,
        "strategy_count": len(strategies),
        "strategies": strategies,
        "strategy_eligible_gene_counts": strategy_eligible_gene_counts,
        "strategy_parameters": strategy_parameters,
        "gene_count": gene_count,
        "gene_ids": gene_ids,
        "alignment_task_count": alignment_task_count,
        "taxonomy_tax_id_count": count_tsv_gz_rows(args.taxonomy) if args.taxonomy else 0,
        "taxonomy_failure_count": count_tsv_gz_rows(args.taxonomy_failures) if args.taxonomy_failures else 0,
        "ortholog_alignment_summary_count": summary_count,
        "strategy_summary_count": strategy_summary_count,
        "alignment_segment_count": segment_count,
        "feature_coverage_count": feature_coverage_count,
        "alignment_event_mode": alignment_event_mode,
        "event_ortholog_support_format": (
            "event_group_id_v1" if alignment_event_mode == "compact_support" else ""
        ),
        "raw_alignment_event_count": raw_event_count,
        "alignment_event_count": event_count,
        "event_ortholog_support_count": event_ortholog_support_count,
        "snv_site_depth_count": snv_site_depth_count,
        "snv_taxonomic_depth_count": snv_taxonomic_depth_count,
        "snv_alt_taxonomic_support_count": taxonomic_alt_support_count,
        "failure_count": failure_count,
        "native_file_count": native_file_count,
    }
    if timings_seconds:
        manifest["timings_seconds"] = timings_seconds
    partition_timings = collect_partition_timings(manifests)
    if partition_timings:
        manifest["partition_timings_seconds"] = partition_timings
        timing_totals: dict[str, float] = {}
        for timings in partition_timings.values():
            for name, value in timings.items():
                timing_totals[name] = timing_totals.get(name, 0.0) + value
        manifest["partition_timing_totals_seconds"] = {
            name: round(value, 3) for name, value in sorted(timing_totals.items())
        }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
