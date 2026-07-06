"""Feature coverage summarization shared by alignment and merge steps."""

from __future__ import annotations

import bisect
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path


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
