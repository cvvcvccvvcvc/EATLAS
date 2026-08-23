"""Small report aggregates derived from per-ortholog alignment evidence."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path


STRATEGY_SUMMARY_FIELDS = [
    "strategy",
    "summary_row_count",
    "gene_count",
    "aligned_summary_row_count",
    "event_count",
    "aligned_target_bp",
]


def write_strategy_summary(
    summary_paths: list[Path],
    output: Path,
    expected_strategies: list[str],
) -> tuple[int, int]:
    aggregates: dict[str, dict[str, int | set[str]]] = {}
    summary_row_count = 0
    for summaries_path in summary_paths:
        with gzip.open(summaries_path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "gene_id",
                "strategy",
                "status",
                "event_count",
                "aligned_target_bp",
            }
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
        writer = csv.DictWriter(
            handle,
            fieldnames=STRATEGY_SUMMARY_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
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
        writer = csv.DictWriter(
            handle,
            fieldnames=STRATEGY_SUMMARY_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for strategy in sorted(aggregates):
            writer.writerow({"strategy": strategy, **aggregates[strategy]})
    return len(aggregates)


def concatenate_tsv_gz(paths: list[Path], output: Path) -> int:
    """Concatenate compatible small analytics tables under one gzip header."""

    output.parent.mkdir(parents=True, exist_ok=True)
    expected_header: list[str] | None = None
    row_count = 0
    with gzip.open(output, "wt", newline="") as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        for path in paths:
            with gzip.open(path, "rt", newline="") as input_handle:
                reader = csv.reader(input_handle, delimiter="\t")
                header = next(reader, None)
                if header is None:
                    raise ValueError(f"Analytics table has no header: {path}")
                if expected_header is None:
                    expected_header = header
                    writer.writerow(header)
                elif header != expected_header:
                    raise ValueError(
                        f"Header mismatch while merging {path}: "
                        f"expected {expected_header}, observed {header}"
                    )
                for row in reader:
                    writer.writerow(row)
                    row_count += 1
    if expected_header is None:
        raise ValueError("No analytics tables supplied for concatenation")
    return row_count
