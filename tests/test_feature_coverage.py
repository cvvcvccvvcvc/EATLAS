from __future__ import annotations

import csv
import gzip
import shutil
import sys
from pathlib import Path

import pytest


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from feature_coverage import (  # noqa: E402
    site_aligned_ortholog_counts,
    summarize_feature_coverage,
    summarize_feature_coverage_rows,
)


pytestmark = pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools is not installed")


def write_tsv_gz(path: Path, rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_row_and_path_feature_coverage_are_equivalent(tmp_path: Path) -> None:
    features = [
        {
            "gene_id": "1",
            "feature_type": "exon",
            "feature_id": "exon_1",
            "genomic_accession": "NC_1",
            "genomic_start1": "101",
            "genomic_end1": "110",
            "target_start0": "0",
            "target_end0": "10",
            "length_bp": "10",
        },
        {
            "gene_id": "1",
            "feature_type": "exon",
            "feature_id": "exon_2",
            "genomic_accession": "NC_1",
            "genomic_start1": "111",
            "genomic_end1": "115",
            "target_start0": "10",
            "target_end0": "15",
            "length_bp": "5",
        },
    ]
    summaries = [
        {"gene_id": "1", "strategy": "test", "ortholog_gene_id": "101"},
        {"gene_id": "1", "strategy": "test", "ortholog_gene_id": "102"},
        {"gene_id": "1", "strategy": "test", "ortholog_gene_id": "103"},
    ]
    segments = [
        {
            "gene_id": "1",
            "strategy": "test",
            "ortholog_gene_id": "101",
            "target_start0": "0",
            "target_end0": "6",
        },
        {
            "gene_id": "1",
            "strategy": "test",
            "ortholog_gene_id": "101",
            "target_start0": "4",
            "target_end0": "8",
        },
        {
            "gene_id": "1",
            "strategy": "test",
            "ortholog_gene_id": "102",
            "target_start0": "5",
            "target_end0": "10",
        },
    ]

    features_path = tmp_path / "features.tsv.gz"
    summaries_path = tmp_path / "summaries.tsv.gz"
    segments_path = tmp_path / "segments.tsv.gz"
    rows_output = tmp_path / "rows.tsv.gz"
    paths_output = tmp_path / "paths.tsv.gz"
    write_tsv_gz(features_path, features)
    write_tsv_gz(summaries_path, summaries)
    write_tsv_gz(segments_path, segments)

    row_count = summarize_feature_coverage_rows(
        features_path,
        summaries,
        segments,
        rows_output,
    )
    path_count = summarize_feature_coverage(
        features_path,
        summaries_path,
        segments_path,
        paths_output,
    )

    assert row_count == path_count == 2
    assert gzip.open(rows_output, "rt").read() == gzip.open(paths_output, "rt").read()
    assert read_tsv_gz(rows_output) == [
        {
            "gene_id": "1",
            "strategy": "test",
            "feature_type": "exon",
            "feature_id": "exon_1",
            "genomic_accession": "NC_1",
            "genomic_start1": "101",
            "genomic_end1": "110",
            "target_start0": "0",
            "target_end0": "10",
            "length_bp": "10",
            "ortholog_count": "3",
            "orthologs_covered": "2",
            "covered_bases": "10",
            "coverage_breadth": "1.000000",
            "depth_bases": "13",
            "mean_depth": "1.300000",
        },
        {
            "gene_id": "1",
            "strategy": "test",
            "feature_type": "exon",
            "feature_id": "exon_2",
            "genomic_accession": "NC_1",
            "genomic_start1": "111",
            "genomic_end1": "115",
            "target_start0": "10",
            "target_end0": "15",
            "length_bp": "5",
            "ortholog_count": "3",
            "orthologs_covered": "0",
            "covered_bases": "0",
            "coverage_breadth": "0.000000",
            "depth_bases": "0",
            "mean_depth": "0.000000",
        },
    ]


def test_summary_without_segments_produces_zero_coverage(tmp_path: Path) -> None:
    features = [
        {
            "gene_id": "2",
            "feature_type": "gene",
            "feature_id": "gene_2",
            "genomic_accession": "NC_2",
            "genomic_start1": "201",
            "genomic_end1": "210",
            "target_start0": "0",
            "target_end0": "10",
            "length_bp": "10",
        }
    ]
    features_path = tmp_path / "features.tsv.gz"
    output = tmp_path / "coverage.tsv.gz"
    write_tsv_gz(features_path, features)

    row_count = summarize_feature_coverage_rows(
        features_path,
        [{"gene_id": "2", "strategy": "test", "ortholog_gene_id": "201"}],
        [],
        output,
    )

    assert row_count == 1
    row = read_tsv_gz(output)[0]
    assert row["ortholog_count"] == "1"
    assert row["orthologs_covered"] == "0"
    assert row["covered_bases"] == "0"
    assert row["depth_bases"] == "0"


def test_site_depth_counts_distinct_primary_orthologs(tmp_path: Path) -> None:
    segments_path = tmp_path / "segments.tsv.gz"
    write_tsv_gz(
        segments_path,
        [
            {
                "gene_id": "1",
                "strategy": "test",
                "ortholog_gene_id": "101",
                "target_start0": "0",
                "target_end0": "6",
                "is_primary": "true",
            },
            {
                "gene_id": "1",
                "strategy": "test",
                "ortholog_gene_id": "101",
                "target_start0": "4",
                "target_end0": "8",
                "is_primary": "true",
            },
            {
                "gene_id": "1",
                "strategy": "test",
                "ortholog_gene_id": "102",
                "target_start0": "5",
                "target_end0": "10",
                "is_primary": "true",
            },
            {
                "gene_id": "1",
                "strategy": "test",
                "ortholog_gene_id": "103",
                "target_start0": "0",
                "target_end0": "10",
                "is_primary": "false",
            },
        ],
    )

    counts = site_aligned_ortholog_counts(
        segments_path,
        [
            {
                "variant_key": "1:5:A>G",
                "gene_id": "1",
                "strategy": "test",
                "target_start0": "4",
            },
            {
                "variant_key": "1:6:A>G",
                "gene_id": "1",
                "strategy": "test",
                "target_start0": "5",
            },
        ],
        tmp_path,
    )

    assert counts == {
        ("1", "test", "1:5:A>G"): 1,
        ("1", "test", "1:6:A>G"): 2,
    }
