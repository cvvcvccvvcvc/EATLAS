"""Aggregate taxonomic SNV evidence into compact report inputs."""

from __future__ import annotations

import bisect
import csv
import gzip
from collections import defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path

from taxonomic_evidence import COUNT_KEYS, SCOPE_ORDER, UNIT_ORDER


CONTEXT_PRIORITY = ("cds", "utr", "exon", "intron")
SUMMARY_FIELDS = [
    "strategy",
    "target_context",
    "taxonomic_scope",
    "evidence_unit",
    "site_aligned_count",
    "alt_support_count",
    "gnomad_found_count",
    "gnomad_not_found_count",
    "gnomad_lookup_failed_count",
]


def _read_context_intervals(
    paths: Iterable[Path],
) -> tuple[dict[str, list[tuple[int, int, str]]], dict[str, list[int]]]:
    features: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    lengths: dict[str, int] = {}
    required = {"gene_id", "feature_type", "target_start0", "target_end0"}
    for path in paths:
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"Target features {path} missing columns: {', '.join(sorted(missing))}"
                )
            for row in reader:
                gene_id = str(row["gene_id"])
                feature = str(row["feature_type"]).lower()
                start = int(row["target_start0"])
                end = int(row["target_end0"])
                if feature == "gene":
                    lengths[gene_id] = max(lengths.get(gene_id, 0), end)
                if feature in CONTEXT_PRIORITY and end > start:
                    features[gene_id][feature].append((start, end))

    contexts: dict[str, list[tuple[int, int, str]]] = {}
    starts: dict[str, list[int]] = {}
    for gene_id, length in lengths.items():
        boundaries = {0, length}
        for intervals in features.get(gene_id, {}).values():
            for start, end in intervals:
                boundaries.add(max(0, start))
                boundaries.add(min(length, end))
        ordered = sorted(boundaries)
        disjoint: list[tuple[int, int, str]] = []
        for start, end in zip(ordered, ordered[1:]):
            if end <= start:
                continue
            context = "other"
            for candidate in CONTEXT_PRIORITY:
                if any(
                    left < end and right > start
                    for left, right in features.get(gene_id, {}).get(candidate, [])
                ):
                    context = "other_exon" if candidate == "exon" else candidate
                    break
            if disjoint and disjoint[-1][1] == start and disjoint[-1][2] == context:
                disjoint[-1] = (disjoint[-1][0], end, context)
            else:
                disjoint.append((start, end, context))
        contexts[gene_id] = disjoint
        starts[gene_id] = [start for start, _end, _context in disjoint]
    return contexts, starts


def _context_at(
    intervals: list[tuple[int, int, str]],
    starts: list[int],
    position: int,
) -> str:
    index = bisect.bisect_right(starts, position) - 1
    if index >= 0 and intervals[index][0] <= position < intervals[index][1]:
        return intervals[index][2]
    return "other"


def _iter_rows(path: Path, required: set[str]) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")
        yield from reader


def _site_key(row: Mapping[str, str]) -> tuple[str, str, int]:
    return str(row["gene_id"]), str(row["strategy"]), int(row["target_start0"])


def write_ortholog_evidence_summary(
    depth_path: Path,
    alt_support_path: Path,
    target_feature_paths: Iterable[Path],
    gnomad_statuses: Mapping[tuple[str, int, str, str], str],
    output: Path,
) -> int:
    """Stream site and ALT counts into a bounded-size histogram."""
    contexts, context_starts = _read_context_intervals(target_feature_paths)
    depth_required = {"gene_id", "strategy", "target_start0", *COUNT_KEYS}
    alt_required = {"gene_id", "strategy", "target_start0", "ref", "alt", *COUNT_KEYS}
    depth_rows = iter(_iter_rows(depth_path, depth_required))
    depth_row = next(depth_rows, None)
    previous_alt_key: tuple[str, str, int] | None = None
    totals: dict[tuple[str, str, str, str, int, int], list[int]] = {}
    status_index = {"found": 0, "not_found": 1, "lookup_failed": 2}

    for alt_row in _iter_rows(alt_support_path, alt_required):
        alt_key = _site_key(alt_row)
        if previous_alt_key is not None and alt_key < previous_alt_key:
            raise ValueError(f"ALT taxonomic support is not sorted at {alt_key}")
        previous_alt_key = alt_key
        while depth_row is not None and _site_key(depth_row) < alt_key:
            depth_row = next(depth_rows, None)
        if depth_row is None or _site_key(depth_row) != alt_key:
            raise ValueError(f"Missing taxonomic site depth for SNV {alt_key}")

        gene_id, strategy, position = alt_key
        context = _context_at(
            contexts.get(gene_id, []),
            context_starts.get(gene_id, []),
            position,
        )
        ref = str(alt_row["ref"]).upper()
        alt = str(alt_row["alt"]).upper()
        status = gnomad_statuses.get((gene_id, position, ref, alt), "lookup_failed")
        if status not in {"found", "not_found", "lookup_failed"}:
            raise ValueError(f"Unknown gnomAD status for {gene_id}:{position}:{ref}>{alt}: {status}")

        for count_key in COUNT_KEYS:
            site_count = int(depth_row[count_key])
            alt_count = int(alt_row[count_key])
            if alt_count < 0 or site_count < 0 or alt_count > site_count:
                raise ValueError(
                    "Invalid taxonomic ortholog evidence for "
                    f"{gene_id}:{position}:{ref}>{alt} / {strategy} / {count_key}: "
                    f"ALT={alt_count}, site={site_count}"
                )
            if site_count == 0:
                continue
            scope, unit = count_key.split("__", 1)
            group = (strategy, context, scope, unit, site_count, alt_count)
            group_counts = totals.setdefault(group, [0, 0, 0])
            group_counts[status_index[status]] += 1

    output.parent.mkdir(parents=True, exist_ok=True)
    scope_index = {scope: index for index, scope in enumerate(SCOPE_ORDER)}
    unit_index = {unit: index for index, unit in enumerate(UNIT_ORDER)}
    row_count = 0
    with gzip.open(output, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for group in sorted(
            totals,
            key=lambda item: (
                item[0],
                item[1],
                scope_index[item[2]],
                unit_index[item[3]],
                item[4],
                item[5],
            ),
        ):
            strategy, context, scope, unit, site_count, alt_count = group
            writer.writerow(
                {
                    "strategy": strategy,
                    "target_context": context,
                    "taxonomic_scope": scope,
                    "evidence_unit": unit,
                    "site_aligned_count": site_count,
                    "alt_support_count": alt_count,
                    "gnomad_found_count": totals[group][0],
                    "gnomad_not_found_count": totals[group][1],
                    "gnomad_lookup_failed_count": totals[group][2],
                }
            )
            row_count += 1
    return row_count
