"""Assign target-sequence positions to exclusive genomic contexts."""

from __future__ import annotations

import bisect
from collections import defaultdict
from pathlib import Path
from typing import Mapping

import pandas as pd


CONTEXT_PRIORITY = ("cds", "utr", "exon", "intron")


def read_disjoint_contexts(
    path: Path,
    gene_lengths: Mapping[str, int],
) -> dict[str, list[tuple[int, int, str]]]:
    """Partition each target gene into exclusive CDS/UTR/exon/intron/other segments."""
    frame = pd.read_csv(
        path,
        sep="\t",
        compression="gzip" if path.suffix == ".gz" else None,
        keep_default_na=False,
        usecols=["gene_id", "feature_type", "target_start0", "target_end0"],
    )
    by_gene: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    for row in frame.itertuples(index=False):
        feature = str(row.feature_type).lower()
        if feature in CONTEXT_PRIORITY:
            by_gene[str(row.gene_id)][feature].append((int(row.target_start0), int(row.target_end0)))

    result = {}
    for gene_id, length_value in gene_lengths.items():
        length = int(length_value)
        feature_intervals = by_gene.get(str(gene_id), {})
        boundaries = {0, length}
        for intervals in feature_intervals.values():
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
                if any(left < end and right > start for left, right in feature_intervals.get(candidate, [])):
                    context = "other_exon" if candidate == "exon" else candidate
                    break
            if disjoint and disjoint[-1][2] == context and disjoint[-1][1] == start:
                disjoint[-1] = (disjoint[-1][0], end, context)
            else:
                disjoint.append((start, end, context))
        result[str(gene_id)] = disjoint
    return result


def context_at(
    intervals: list[tuple[int, int, str]],
    position: int,
    starts: list[int] | None = None,
) -> str:
    """Return the exclusive target context at a zero-based target position."""
    if starts is None:
        starts = [item[0] for item in intervals]
    index = bisect.bisect_right(starts, position) - 1
    if index >= 0 and intervals[index][0] <= position < intervals[index][1]:
        return intervals[index][2]
    return "other"
