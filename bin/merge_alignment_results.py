#!/usr/bin/env python3
"""Merge per-gene alignment evidence outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

from feature_coverage import (
    iter_snv_event_sites,
    summarize_feature_coverage,
    write_snv_site_depth,
    write_snv_taxonomic_depth,
)
from taxonomic_evidence import COUNT_KEYS, SCOPE_ORDER, UNIT_ORDER, load_taxonomy_profiles


csv.field_size_limit(sys.maxsize)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-tasks", type=Path)
    parser.add_argument("--taxonomy-presets", type=Path)
    parser.add_argument("--taxonomy-failures", type=Path)
    parser.add_argument("--taxonomy-summary", type=Path)
    parser.add_argument("--target-features", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--result-dir", action="append", default=[], type=Path)
    parser.add_argument("--result-root", type=Path)
    parser.add_argument("--partition-id")
    parser.add_argument("--expected-strategies", required=True)
    parser.add_argument("--expected-gene-ids")
    parser.add_argument("--compact-events", action="store_true")
    parser.add_argument("--events-already-compacted", action="store_true")
    parser.add_argument(
        "--output-profile",
        choices=["full", "annotation-input", "report-input"],
        default="full",
        help="Select full outputs, partitioned annotation inputs, or final report inputs.",
    )
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
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "support_row_count",
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

BWA_STRATEGY = "bwa_pseudoreads"
ENSEMBL_COMPARA_STRATEGY = "precomputed_ensembl_92_mammals_epo_extended"
DNA_BASES = frozenset("ACGT")
BWA_REQUIRED_PARAMETERS = (
    "pseudoread_len",
    "pseudoread_step",
    "pseudoread_phred",
)
BWA_OPTIONAL_PARAMETERS = ("task_cpus", "bwa_threads")


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


def _task_ready(row: dict[str, str], field: str) -> bool:
    if field in row:
        return str(row[field]).lower() == "true"
    return row.get("status") == "ready"


def read_alignment_capabilities(
    path: Path,
) -> tuple[int, dict[str, tuple[bool, bool]]]:
    """Return target/ortholog readiness for every alignment task gene."""
    task_count = 0
    capabilities: dict[str, tuple[bool, bool]] = {}
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
            gene_id = row.get("gene_id", "")
            if not gene_id:
                continue
            if gene_id in capabilities:
                raise ValueError(f"Alignment tasks {path} contain duplicate gene_id {gene_id}")
            capabilities[gene_id] = (
                _task_ready(row, "target_ready"),
                _task_ready(row, "ortholog_ready"),
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


def create_taxonomy_table(conn: sqlite3.Connection, taxonomy_presets: Path) -> None:
    profiles = load_taxonomy_profiles(taxonomy_presets)
    scope_columns = [scope for scope in SCOPE_ORDER if scope != "all"]
    conn.execute(
        "CREATE TABLE taxonomy_profiles ("
        "tax_id TEXT PRIMARY KEY, species_id TEXT, genus_id TEXT, family_id TEXT, order_id TEXT, "
        + ", ".join(f"{scope} INTEGER NOT NULL" for scope in scope_columns)
        + ") WITHOUT ROWID"
    )
    fields = ["tax_id", "species_id", "genus_id", "family_id", "order_id", *scope_columns]
    placeholders = ", ".join("?" for _field in fields)
    rows = []
    for profile in profiles.values():
        scopes = set(profile.scopes())
        rows.append(
            (
                profile.tax_id,
                profile.species_id,
                profile.genus_id,
                profile.family_id,
                profile.order_id,
                *(int(scope in scopes) for scope in scope_columns),
            )
        )
    conn.executemany(
        f"INSERT INTO taxonomy_profiles ({', '.join(fields)}) VALUES ({placeholders})",
        rows,
    )


def taxonomy_count_expressions() -> list[str]:
    expressions = []
    for scope in SCOPE_ORDER:
        scope_condition = "e.tax_id != ''" if scope == "all" else f"p.{scope} = 1"
        for unit in UNIT_ORDER:
            group_value = (
                "e.ortholog_gene_id"
                if unit == "ortholog"
                else f"COALESCE(NULLIF(p.{unit}_id, ''), 'taxon:' || e.tax_id)"
            )
            expressions.append(
                "COUNT(DISTINCT CASE WHEN "
                f"{scope_condition} AND e.ortholog_gene_id != '' THEN {group_value} END) "
                f'AS "{scope}__{unit}"'
            )
    return expressions


def write_compact_events(
    paths: list[Path],
    output: Path,
    ortholog_support_output: Path,
    taxonomy_presets: Path | None = None,
    taxonomic_support_output: Path | None = None,
) -> tuple[int, int, int, int]:
    if (taxonomy_presets is None) != (taxonomic_support_output is None):
        raise ValueError("Taxonomic event support requires both taxonomy presets and an output path")
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
        if taxonomy_presets is not None:
            create_taxonomy_table(conn, taxonomy_presets)
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
                alt,
                strategy
            )
            """
        )

        taxonomic_select = ""
        taxonomic_join = ""
        if taxonomy_presets is not None:
            taxonomic_select = ",\n                " + ",\n                ".join(
                taxonomy_count_expressions()
            )
            taxonomic_join = "LEFT JOIN taxonomy_profiles AS p ON p.tax_id = e.tax_id"
        query = f"""
            SELECT
                e.gene_id,
                e.event_type,
                e.target_start0,
                e.target_end0,
                e.genomic_accession,
                e.genomic_start1,
                e.genomic_end1,
                e.ref,
                e.alt,
                e.strategy,
                COUNT(*) AS support_row_count,
                COUNT(DISTINCT e.ortholog_gene_id) AS support_ortholog_count,
                GROUP_CONCAT(DISTINCT e.tool) AS tools,
                GROUP_CONCAT(DISTINCT e.preset) AS presets,
                COUNT(DISTINCT e.tax_id) AS tax_id_count,
                COUNT(DISTINCT e.taxname) AS taxname_count,
                GROUP_CONCAT(DISTINCT e.qc_flags) AS qc_flags
                {taxonomic_select}
            FROM events AS e
            {taxonomic_join}
            GROUP BY
                e.gene_id,
                e.event_type,
                e.target_start0,
                e.target_end0,
                e.genomic_accession,
                e.genomic_start1,
                e.genomic_end1,
                e.ref,
                e.alt,
                e.strategy
            ORDER BY
                e.gene_id,
                e.strategy,
                CAST(e.target_start0 AS INTEGER),
                e.event_type,
                e.ref,
                e.alt
        """
        compact_count = 0
        taxonomic_support_count = 0
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
            with gzip.open(output, "wt", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=COMPACT_EVENT_FIELDS, delimiter="\t", extrasaction="ignore")
                writer.writeheader()
                for row in conn.execute(query):
                    record = dict(
                        zip(
                            [
                                *COMPACT_EVENT_FIELDS,
                                *(COUNT_KEYS if taxonomy_presets is not None else ()),
                            ],
                            row,
                        )
                    )
                    record["qc_flags"] = compact_event_flags(record.get("qc_flags") or "")
                    writer.writerow(record)
                    compact_count += 1
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
        ortholog_support_query = """
            SELECT
                e.gene_id,
                e.event_type,
                e.target_start0,
                e.target_end0,
                e.genomic_accession,
                e.genomic_start1,
                e.genomic_end1,
                e.ref,
                e.alt,
                e.strategy,
                e.ortholog_gene_id,
                e.tax_id,
                e.taxname,
                COUNT(*) AS support_row_count
            FROM events AS e
            WHERE e.ortholog_gene_id != ''
            GROUP BY
                e.gene_id,
                e.event_type,
                e.target_start0,
                e.target_end0,
                e.genomic_accession,
                e.genomic_start1,
                e.genomic_end1,
                e.ref,
                e.alt,
                e.strategy,
                e.ortholog_gene_id,
                e.tax_id,
                e.taxname
            ORDER BY
                e.gene_id,
                e.strategy,
                CAST(e.target_start0 AS INTEGER),
                e.event_type,
                e.ref,
                e.alt,
                e.ortholog_gene_id
        """
        ortholog_support_count = 0
        with gzip.open(ortholog_support_output, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(EVENT_ORTHOLOG_SUPPORT_FIELDS)
            for row in conn.execute(ortholog_support_query):
                writer.writerow(row)
                ortholog_support_count += 1
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
    return sorted_gene_ids(gene_ids)


def merge_strategy_parameters(manifests: list[dict]) -> dict[str, dict]:
    merged: dict[str, dict] = {}
    for manifest in manifests:
        parameters = manifest.get("strategy_parameters") or {}
        if not isinstance(parameters, dict):
            raise ValueError("strategy_parameters must be a JSON object")
        candidates = dict(parameters)

        if manifest.get("strategy") == BWA_STRATEGY:
            present_required = {
                name: manifest[name]
                for name in BWA_REQUIRED_PARAMETERS
                if manifest.get(name) is not None
            }
            if present_required:
                missing = sorted(set(BWA_REQUIRED_PARAMETERS) - set(present_required))
                if missing:
                    raise ValueError(
                        "BWA manifest has incomplete pseudoread parameters: "
                        + ", ".join(missing)
                    )
                direct = dict(present_required)
                direct.update(
                    {
                        name: manifest[name]
                        for name in BWA_OPTIONAL_PARAMETERS
                        if manifest.get(name) is not None
                    }
                )
                nested = candidates.get(BWA_STRATEGY)
                if nested is not None and nested != direct:
                    raise ValueError(
                        "BWA manifest has conflicting direct and nested strategy parameters"
                    )
                candidates[BWA_STRATEGY] = direct

        for strategy, values in candidates.items():
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
    require_feature_coverage: bool,
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
            "failures.tsv.gz",
        ]
        if require_feature_coverage:
            filenames.append("feature_coverage.tsv.gz")
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
    alignment_tasks: Path | None = None,
) -> tuple[list[str], list[str]]:
    if alignment_tasks is None:
        expected_pairs = {
            (gene_id, strategy)
            for gene_id in expected_gene_ids
            for strategy in expected_strategies
        }
    else:
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


def sum_manifest_count(
    manifests: list[dict],
    field: str,
    *fallback_fields: str,
) -> int:
    total = 0
    for manifest in manifests:
        value = manifest.get(field)
        if value is None:
            value = next(
                (manifest.get(name) for name in fallback_fields if manifest.get(name) is not None),
                0,
            )
        total += int(value or 0)
    return total


def merged_event_mode(manifests: list[dict]) -> str:
    modes = {
        str(manifest.get("alignment_event_mode") or "raw")
        for manifest in manifests
    }
    return next(iter(modes)) if len(modes) == 1 else "mixed"


def main() -> None:
    args = parse_args()
    if args.events_already_compacted and not args.compact_events:
        raise ValueError("--events-already-compacted requires --compact-events")
    if args.partition_id and args.output_profile == "report-input":
        raise ValueError("--output-profile report-input is only valid for the final merge")
    if not args.partition_id and args.output_profile == "annotation-input":
        raise ValueError("--output-profile annotation-input requires --partition-id")
    expected_strategies = sorted(
        parse_expected_values(args.expected_strategies, "--expected-strategies")
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    result_dirs = resolve_result_dirs(args.result_dir, args.result_root)

    global_inputs = [
        args.alignment_tasks,
        args.taxonomy_presets,
        args.taxonomy_failures,
        args.target_features,
    ]
    if not args.partition_id and any(path is None for path in global_inputs):
        raise ValueError(
            "Final alignment merge requires --alignment-tasks, --taxonomy-presets, "
            "--taxonomy-failures, and --target-features"
        )

    manifests = load_manifests(result_dirs)
    require_alignment_tables(
        result_dirs,
        require_feature_coverage=bool(args.partition_id),
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
            copy_or_keep(args.taxonomy_presets, args.outdir / "taxonomy_presets.tsv.gz")
            copy_or_keep(args.taxonomy_failures, args.outdir / "taxonomy_failures.tsv.gz")
        if args.taxonomy_summary is not None:
            copy_or_keep(args.taxonomy_summary, args.outdir / "taxonomy_summary.tsv.gz")
    gene_count = len(gene_ids)

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
            segment_count = sum_manifest_count(
                manifests,
                "alignment_segment_count",
                "segment_count",
            )

    feature_coverage_inputs = [path / "feature_coverage.tsv.gz" for path in result_dirs]
    missing_feature_coverage = [str(path.parent) for path in feature_coverage_inputs if not path.exists()]
    if missing_feature_coverage:
        if args.partition_id or args.output_profile == "report-input":
            raise FileNotFoundError(
                "Alignment inputs missing feature_coverage.tsv.gz: "
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
        alignment_event_mode = merged_event_mode(manifests)
    else:
        event_inputs = [path / "alignment_events.tsv.gz" for path in result_dirs]
        if args.compact_events and not args.events_already_compacted:
            (
                event_count,
                raw_event_count,
                taxonomic_alt_support_count,
                event_ortholog_support_count,
            ) = write_compact_events(
                event_inputs,
                args.outdir / "alignment_events.tsv.gz",
                args.outdir / "event_ortholog_support.tsv.gz",
                args.taxonomy_presets if args.partition_id and args.taxonomy_presets else None,
                args.outdir / "snv_alt_taxonomic_support.tsv.gz"
                if args.partition_id and args.taxonomy_presets
                else None,
            )
        else:
            event_count = merge_tsv_gz(
                event_inputs,
                args.outdir / "alignment_events.tsv.gz",
            )
            raw_event_count = event_count
            taxonomic_alt_support_count = 0
            event_ortholog_support_count = 0
            if args.events_already_compacted:
                raw_event_count = sum(
                    int(manifest.get("raw_alignment_event_count") or manifest.get("alignment_event_count") or 0)
                    for manifest in manifests
                )
                event_ortholog_support_count = merge_tsv_gz(
                    [path / "event_ortholog_support.tsv.gz" for path in result_dirs],
                    args.outdir / "event_ortholog_support.tsv.gz",
                )
        alignment_event_mode = "compact_support" if args.compact_events else "raw"

    if args.partition_id:
        snv_site_depth_count = write_snv_site_depth(
            [path / "alignment_segments.tsv.gz" for path in result_dirs],
            iter_snv_event_sites([args.outdir / "alignment_events.tsv.gz"]),
            args.outdir / "snv_site_depth.tsv.gz",
            args.outdir,
        )
        if args.taxonomy_presets is not None:
            snv_taxonomic_depth_count = write_snv_taxonomic_depth(
                [path / "alignment_segments.tsv.gz" for path in result_dirs],
                iter_snv_event_sites([args.outdir / "alignment_events.tsv.gz"]),
                args.taxonomy_presets,
                args.outdir / "snv_taxonomic_depth.tsv.gz",
                args.outdir,
            )
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
    strategy_parameters = merge_strategy_parameters(manifests)
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
        "taxonomy_tax_id_count": count_tsv_gz_rows(args.taxonomy_presets) if args.taxonomy_presets else 0,
        "taxonomy_failure_count": count_tsv_gz_rows(args.taxonomy_failures) if args.taxonomy_failures else 0,
        "ortholog_alignment_summary_count": summary_count,
        "strategy_summary_count": strategy_summary_count,
        "alignment_segment_count": segment_count,
        "feature_coverage_mode": feature_coverage_mode,
        "feature_coverage_missing_result_count": len(missing_feature_coverage),
        "feature_coverage_count": feature_coverage_count,
        "alignment_event_mode": alignment_event_mode,
        "raw_alignment_event_count": raw_event_count,
        "alignment_event_count": event_count,
        "event_ortholog_support_count": event_ortholog_support_count,
        "snv_site_depth_count": snv_site_depth_count,
        "snv_taxonomic_depth_count": snv_taxonomic_depth_count,
        "snv_alt_taxonomic_support_count": taxonomic_alt_support_count,
        "failure_count": failure_count,
        "native_file_count": native_file_count,
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
