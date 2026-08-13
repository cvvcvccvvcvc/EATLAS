"""Keep the dominant-strand monotonic backbone of BWA pseudo-read alignments."""

from __future__ import annotations

import json
import logging
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam


logger = logging.getLogger("bam_filtering")


@dataclass
class FilterResult:
    """Output metadata for one gene-level BAM filtering run."""

    input_bam: Path
    output_bam: Path
    output_bai: Path
    per_homologue_stats_json: Path
    overall_stats_json: Path
    filtering_stats: dict[str, dict[str, Any]]
    overall: dict[str, Any]


def _read_key(read: pysam.AlignedSegment) -> tuple[Any, ...]:
    """Stable identity for one BAM alignment record."""
    return (
        read.query_name,
        read.flag,
        read.reference_id,
        read.reference_start,
        read.reference_end,
        read.cigarstring,
        read.mapping_quality,
    )


def parse_read_name(read_name: str) -> tuple[str, int, int, int]:
    """Parse ``<ortholog>_pseudo_<number>_<start>-<end>``."""
    parts = read_name.split("_pseudo_")
    if len(parts) != 2:
        raise ValueError(f"Unexpected read name format: {read_name}")

    homologue_id = parts[0]
    read_num_str, pos_range = parts[1].split("_", 1)
    start_str, end_str = pos_range.split("-")
    return homologue_id, int(read_num_str), int(start_str), int(end_str)


def _build_read_record(
    read: pysam.AlignedSegment,
    parsed: tuple[str, int, int, int],
) -> dict[str, Any]:
    """Convert one alignment row to the compact representation used by LIS."""
    return {
        "read_key": _read_key(read),
        "actual_read_num": parsed[1],
        "alignment_pos": read.reference_start,
        "is_reverse": read.is_reverse,
    }


def collect_homologue_data(
    bam_path: Path,
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, dict[str, Any]]]:
    """Collect mapped reads and dominant-strand statistics in one BAM pass."""
    stats: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "total_reads": 0,
            "forward_reads": 0,
            "reverse_reads": 0,
        }
    )
    reads_by_homologue: dict[str, list[dict[str, Any]]] = defaultdict(list)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            parsed = parse_read_name(read.query_name)
            homologue_id = parsed[0]
            stats[homologue_id]["total_reads"] += 1
            strand_key = "reverse_reads" if read.is_reverse else "forward_reads"
            stats[homologue_id][strand_key] += 1
            reads_by_homologue[homologue_id].append(_build_read_record(read, parsed))

    for homologue_id, homologue_stats in stats.items():
        homologue_stats["dominant_strand"] = (
            "forward"
            if homologue_stats["forward_reads"] >= homologue_stats["reverse_reads"]
            else "reverse"
        )
        reads_by_homologue[homologue_id].sort(key=lambda row: row["alignment_pos"])

    return dict(reads_by_homologue), dict(stats)


def _lis_indices(values: list[int]) -> set[int]:
    """Return indices of one longest strictly increasing subsequence."""
    if not values:
        return set()

    tails_values: list[int] = []
    tails_indices: list[int] = []
    previous = [-1] * len(values)

    for index, value in enumerate(values):
        position = bisect_left(tails_values, value)
        if position == len(tails_values):
            tails_values.append(value)
            tails_indices.append(index)
        else:
            tails_values[position] = value
            tails_indices[position] = index
        if position > 0:
            previous[index] = tails_indices[position - 1]

    selected: set[int] = set()
    current = tails_indices[-1]
    while current != -1:
        selected.add(current)
        current = previous[current]
    return selected


def _backbone_indices(
    reads: list[dict[str, Any]],
    dominant_strand: str,
) -> set[int]:
    """Find the forward LIS or reverse LDS in target-coordinate order."""
    sequence = [int(read["actual_read_num"]) for read in reads]
    if dominant_strand == "reverse":
        sequence = [-value for value in sequence]
    return _lis_indices(sequence)


def _filter_homologue_reads(
    reads: list[dict[str, Any]],
    dominant_strand: str,
) -> tuple[set[tuple[Any, ...]], dict[str, Any]]:
    """Apply the mandatory dominant-strand and LIS/LDS filters."""
    if not reads:
        return set(), {}
    if dominant_strand not in {"forward", "reverse"}:
        raise ValueError(f"Unsupported dominant strand: {dominant_strand}")

    initial_count = len(reads)
    forward_initial = sum(not read["is_reverse"] for read in reads)
    reverse_initial = initial_count - forward_initial
    reads_after_strand = [
        read
        for read in reads
        if read["is_reverse"] == (dominant_strand == "reverse")
    ]
    backbone = _backbone_indices(reads_after_strand, dominant_strand)
    retained_reads = [
        read for index, read in enumerate(reads_after_strand) if index in backbone
    ]

    retained_count = len(retained_reads)
    stats = {
        "initial_count": initial_count,
        "forward_initial": forward_initial,
        "reverse_initial": reverse_initial,
        "dominant_strand": dominant_strand,
        "after_strand_filter": len(reads_after_strand),
        "after_order_filter": retained_count,
        "filtered_by_strand": initial_count - len(reads_after_strand),
        "filtered_by_order": len(reads_after_strand) - retained_count,
        "total_filtered": initial_count - retained_count,
        "pct_kept": 100.0 * retained_count / initial_count,
    }
    return {read["read_key"] for read in retained_reads}, stats


def filter_bam_for_gene(work_dir: Path) -> FilterResult:
    """Write the mandatory dominant-strand LIS/LDS backbone BAM."""
    work_dir = Path(work_dir)
    input_bam = work_dir / "aln.sorted.bam"
    output_bam = work_dir / "aln.filtered.lis.bam"
    per_homologue_stats_json = work_dir / "bam_filtering_stats.json"
    overall_stats_json = work_dir / "bam_filtering_overall.json"
    if not input_bam.exists():
        raise FileNotFoundError(f"Input BAM not found: {input_bam}")

    reads_by_homologue, homologue_stats = collect_homologue_data(input_bam)
    all_keep_keys: set[tuple[Any, ...]] = set()
    filtering_stats: dict[str, dict[str, Any]] = {}
    for homologue_id in sorted(homologue_stats):
        keep_keys, stats = _filter_homologue_reads(
            reads_by_homologue.get(homologue_id, []),
            homologue_stats[homologue_id]["dominant_strand"],
        )
        all_keep_keys.update(keep_keys)
        filtering_stats[homologue_id] = stats

    with pysam.AlignmentFile(str(input_bam), "rb") as input_handle:
        with pysam.AlignmentFile(str(output_bam), "wb", header=input_handle.header) as output_handle:
            for read in input_handle.fetch():
                if _read_key(read) in all_keep_keys:
                    output_handle.write(read)
    pysam.index(str(output_bam))
    output_bai = Path(f"{output_bam}.bai")

    total_initial = sum(stats["initial_count"] for stats in filtering_stats.values())
    total_kept = sum(stats["after_order_filter"] for stats in filtering_stats.values())
    overall = {
        "input_bam": str(input_bam),
        "output_bam": str(output_bam),
        "per_homologue_stats_json": str(per_homologue_stats_json),
        "filters": ["dominant_strand", "lis_backbone"],
        "homologue_count": len(filtering_stats),
        "total_initial_reads": total_initial,
        "total_kept_reads": total_kept,
        "total_filtered_reads": total_initial - total_kept,
        "pct_kept": 100.0 * total_kept / total_initial if total_initial else 0.0,
    }
    per_homologue_stats_json.write_text(
        json.dumps(filtering_stats, indent=2, sort_keys=True) + "\n"
    )
    overall_stats_json.write_text(json.dumps(overall, indent=2, sort_keys=True) + "\n")

    logger.info(
        "BAM filtering complete: output=%s kept=%s/%s (%.2f%%)",
        output_bam,
        total_kept,
        total_initial,
        overall["pct_kept"],
    )
    return FilterResult(
        input_bam=input_bam,
        output_bam=output_bam,
        output_bai=output_bai,
        per_homologue_stats_json=per_homologue_stats_json,
        overall_stats_json=overall_stats_json,
        filtering_stats=filtering_stats,
        overall=overall,
    )
