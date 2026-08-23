#!/usr/bin/env python3
"""Merge per-gene alignment evidence outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import shutil
import sqlite3
import sys
import time
from datetime import datetime, timezone
from itertools import groupby
from pathlib import Path

from bin.alignment_table_schema import ALIGNER_OUTPUT_SCHEMAS


csv.field_size_limit(sys.maxsize)


logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-tasks", required=True, type=Path)
    parser.add_argument("--source-genes", type=Path)
    parser.add_argument("--source-target-features", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--result-dir", action="append", default=[], type=Path)
    parser.add_argument("--result-root", type=Path)
    parser.add_argument("--partition-id")
    parser.add_argument("--expected-strategies", required=True)
    parser.add_argument("--expected-gene-ids")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
    "qc_flags",
]
EVENT_ORTHOLOG_SUPPORT_FIELDS = [
    "event_group_id",
    "ortholog_gene_id",
    "tax_id",
    "mapq",
    "native_alignment_type",
    "support_row_count",
]
PARTITION_EVIDENCE_FILES = (
    "manifest.json",
    "ortholog_alignment_summary.tsv.gz",
    "alignment_segments.tsv.gz",
    "alignment_events.tsv.gz",
    "event_ortholog_support.tsv.gz",
)
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
    "tax_id",
    "mapq",
    "native_alignment_type",
    "qc_flags",
]
ENSEMBL_COMPARA_STRATEGY = "precomputed_ensembl_92_mammals_epo_extended"
ALIGNER_MANIFEST_COUNT_FIELDS = (
    "ortholog_alignment_summary_count",
    "alignment_segment_count",
    "raw_alignment_event_count",
    "alignment_event_count",
    "failure_count",
)
PARTITION_MANIFEST_COUNT_FIELDS = (
    "ortholog_alignment_summary_count",
    "alignment_segment_count",
    "raw_alignment_event_count",
    "alignment_event_count",
    "event_ortholog_support_count",
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


def require_exact_header(
    path: Path,
    observed: list[str] | None,
    expected: list[str],
) -> None:
    observed = observed or []
    if observed != expected:
        raise ValueError(
            f"Alignment table {path} has invalid header: "
            f"expected {expected}, observed {observed}"
        )


def require_exact_tsv_gz_header(path: Path, expected: list[str]) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing required alignment table: {path}")
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        require_exact_header(path, next(reader, None), expected)


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
) -> list[str]:
    """Return the selected-strategy-eligible gene union."""
    _task_count, capabilities = read_alignment_capabilities(path)
    gene_ids = sorted_gene_ids(set(capabilities))
    pairs = expected_gene_strategy_pairs(capabilities, gene_ids, expected_strategies)
    return sorted_gene_ids({gene_id for gene_id, _strategy in pairs})


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
            tax_id TEXT,
            mapq TEXT,
            native_alignment_type TEXT,
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
        "tax_id",
        "mapq",
        "native_alignment_type",
        "qc_flags",
    ]
    count = 0
    batch: list[tuple[str, ...]] = []
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            batch.append(tuple(row[field] for field in fields))
            count += 1
            if len(batch) >= batch_size:
                conn.executemany(
                    "INSERT INTO events VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    batch,
                )
                batch.clear()
    if batch:
        conn.executemany(
            "INSERT INTO events VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
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
            "mapq": row["mapq"],
            "native_alignment_type": row["native_alignment_type"],
            "support_row_count": 1,
        }
        return

    current_tax_id = str(support["tax_id"] or "")
    observed_tax_id = row["tax_id"]
    if current_tax_id and observed_tax_id and current_tax_id != observed_tax_id:
        raise ValueError(
            f"Conflicting tax_id for ortholog_gene_id={ortholog_gene_id}: "
            f"{current_tax_id!r} != {observed_tax_id!r}"
        )
    if not current_tax_id and observed_tax_id:
        support["tax_id"] = observed_tax_id
    for field in ("mapq", "native_alignment_type"):
        if not support[field] and row[field]:
            support[field] = row[field]
    support["support_row_count"] = int(support["support_row_count"]) + 1


def compact_stream_group(
    event_group_id: int,
    key: tuple[str, ...],
    rows: object,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    record: dict[str, object] = {
        "event_group_id": event_group_id,
        **dict(zip(EVENT_KEY_FIELDS, key)),
    }
    orthologs: dict[str, dict[str, object]] = {}
    qc_flags: list[str] = []
    for values in rows:
        row = dict(zip(EVENT_STREAM_FIELDS, values))
        merge_ortholog_metadata(orthologs, row)
        if row["qc_flags"]:
            qc_flags.append(row["qc_flags"])

    record.update(
        {
            "qc_flags": compact_event_flags(",".join(qc_flags)),
        }
    )
    support_rows = [orthologs[key] for key in sorted(orthologs)]
    return record, support_rows


def write_compact_events(
    paths: list[Path],
    output: Path,
    ortholog_support_output: Path,
    timings: dict[str, float] | None = None,
) -> tuple[int, int, int]:
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
                tax_id,
                mapq,
                native_alignment_type,
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
        compact_count = 0
        ortholog_support_count = 0
        phase_started = start_phase("stream_event_groups")
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
                )
                writer.writerow(record)
                for ortholog_row in ortholog_rows:
                    ortholog_writer.writerow(
                        {"event_group_id": compact_count, **ortholog_row}
                    )
                    ortholog_support_count += 1
        finish_phase(phase_timings, "stream_event_groups", phase_started)
        return compact_count, raw_count, ortholog_support_count
    finally:
        conn.close()
        if db_path.exists():
            db_path.unlink()


def copy_partitioned_evidence(result_dirs: list[Path], outdir: Path) -> None:
    partitions_root = outdir / "evidence" / "partitions"
    for result_dir in result_dirs:
        manifest = json.loads((result_dir / "manifest.json").read_text())
        partition_id = str(manifest.get("partition_id") or "")
        if not partition_id:
            raise ValueError(f"Alignment partition {result_dir} does not declare partition_id")
        partition_outdir = partitions_root / partition_id
        for filename in PARTITION_EVIDENCE_FILES:
            copy_or_keep(result_dir / filename, partition_outdir / filename)


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
    input_event_mode: str,
) -> None:
    filenames = [
        "ortholog_alignment_summary.tsv.gz",
        "alignment_segments.tsv.gz",
        "alignment_events.tsv.gz",
        "failures.tsv.gz",
    ]
    if input_event_mode == "compact_support":
        filenames.append("event_ortholog_support.tsv.gz")
    missing = [
        str(result_dir / filename)
        for result_dir in result_dirs
        for filename in filenames
        if not (result_dir / filename).exists()
    ]
    if missing:
        raise FileNotFoundError("Missing required alignment table(s): " + ", ".join(missing))


def validate_alignment_table_schemas(
    result_dirs: list[Path],
    input_event_mode: str,
) -> None:
    filenames = [
        "failures.tsv.gz",
        "ortholog_alignment_summary.tsv.gz",
        "alignment_segments.tsv.gz",
    ]
    if input_event_mode == "raw":
        filenames.append("alignment_events.tsv.gz")

    for result_dir in result_dirs:
        for filename in filenames:
            require_exact_tsv_gz_header(
                result_dir / filename,
                ALIGNER_OUTPUT_SCHEMAS[filename],
            )
        if input_event_mode == "compact_support":
            require_exact_tsv_gz_header(
                result_dir / "alignment_events.tsv.gz",
                COMPACT_EVENT_FIELDS,
            )
            require_exact_tsv_gz_header(
                result_dir / "event_ortholog_support.tsv.gz",
                EVENT_ORTHOLOG_SUPPORT_FIELDS,
            )


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
        mode = manifest.get("alignment_event_mode")
        if mode not in {"raw", "compact_support"}:
            raise ValueError(f"Alignment manifest has invalid alignment_event_mode: {mode!r}")
        count_fields = (
            ALIGNER_MANIFEST_COUNT_FIELDS
            if mode == "raw"
            else PARTITION_MANIFEST_COUNT_FIELDS
        )
        for field in count_fields:
            manifest_count(manifest, field)
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
    expected_strategies = sorted(
        parse_expected_values(args.expected_strategies, "--expected-strategies")
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    result_dirs = resolve_result_dirs(args.result_dir, args.result_root)

    source_context_inputs = [args.source_genes, args.source_target_features]
    if not args.partition_id and any(path is None for path in source_context_inputs):
        raise ValueError(
            "Final alignment merge requires --source-genes and --source-target-features"
        )
    for path in source_context_inputs:
        if path is not None and not path.exists():
            raise FileNotFoundError(f"Missing source target context: {path}")

    manifests = load_manifests(result_dirs)
    validate_alignment_manifests(manifests)
    strategy_parameters = merge_strategy_parameters(manifests)
    input_event_mode = merged_event_mode(manifests)
    require_alignment_tables(
        result_dirs,
        input_event_mode=input_event_mode,
    )
    validate_alignment_table_schemas(
        result_dirs,
        input_event_mode=input_event_mode,
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
    else:
        expected_gene_ids = summarize_alignment_tasks(
            args.alignment_tasks,
            expected_strategies,
        )
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
    gene_count = len(gene_ids)

    expected_event_mode = "raw" if args.partition_id else "compact_support"
    if input_event_mode != expected_event_mode:
        merge_level = "Partition" if args.partition_id else "Final"
        raise ValueError(
            f"{merge_level} alignment merge requires {expected_event_mode} inputs, "
            f"observed {input_event_mode}"
        )
    if not args.partition_id:
        copy_partitioned_evidence(result_dirs, args.outdir)
    if not args.partition_id:
        summary_count = sum_manifest_count(manifests, "ortholog_alignment_summary_count")
        segment_count = sum_manifest_count(manifests, "alignment_segment_count")
    else:
        summary_inputs = [
            path / "ortholog_alignment_summary.tsv.gz"
            for path in result_dirs
        ]
        summary_count = merge_tsv_gz(
            summary_inputs,
            args.outdir / "ortholog_alignment_summary.tsv.gz",
        )
        segment_count = merge_tsv_gz(
            [path / "alignment_segments.tsv.gz" for path in result_dirs],
            args.outdir / "alignment_segments.tsv.gz",
        )

    if not args.partition_id:
        event_count = sum_manifest_count(manifests, "alignment_event_count")
        raw_event_count = sum_manifest_count(manifests, "raw_alignment_event_count")
        event_ortholog_support_count = sum_manifest_count(
            manifests,
            "event_ortholog_support_count",
        )
        alignment_event_mode = input_event_mode
    else:
        event_inputs = [path / "alignment_events.tsv.gz" for path in result_dirs]
        (
            event_count,
            raw_event_count,
            event_ortholog_support_count,
        ) = write_compact_events(
            event_inputs,
            args.outdir / "alignment_events.tsv.gz",
            args.outdir / "event_ortholog_support.tsv.gz",
            timings=timings_seconds,
        )
        alignment_event_mode = "compact_support"
    failure_count = merge_tsv_gz(
        [path / "failures.tsv.gz" for path in result_dirs],
        args.outdir / "failures.tsv.gz",
    )
    manifest = {
        "created_at": utc_now(),
        "stage": "alignment",
        "partition_id": args.partition_id or "",
        "schema": (
            "normalized_alignment_evidence_partition_v2"
            if args.partition_id
            else "normalized_alignment_evidence_v2"
        ),
        "strategy_count": len(strategies),
        "strategies": strategies,
        "strategy_parameters": strategy_parameters,
        "gene_count": gene_count,
        "gene_ids": gene_ids,
        "ortholog_alignment_summary_count": summary_count,
        "alignment_segment_count": segment_count,
        "alignment_event_mode": alignment_event_mode,
        "event_ortholog_support_format": (
            "event_group_id_v2" if alignment_event_mode == "compact_support" else ""
        ),
        "raw_alignment_event_count": raw_event_count,
        "alignment_event_count": event_count,
        "event_ortholog_support_count": event_ortholog_support_count,
        "failure_count": failure_count,
    }
    if not args.partition_id:
        manifest["source_target_context"] = {
            "genes_sha256": sha256_file(args.source_genes),
            "target_features_sha256": sha256_file(args.source_target_features),
        }
    if not args.partition_id:
        manifest["normalized_evidence"] = {
            "layout": "partitioned",
            "format": "tsv_gzip_v1",
            "path": "evidence/partitions",
            "partition_count": len(result_dirs),
            "partition_files": list(PARTITION_EVIDENCE_FILES),
            "event_group_id_scope": "partition",
        }
    if timings_seconds:
        manifest["timings_seconds"] = timings_seconds
    partition_timings = collect_partition_timings(manifests)
    if partition_timings:
        manifest["partition_timings_seconds"] = partition_timings
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
