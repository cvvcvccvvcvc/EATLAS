from __future__ import annotations

import gzip
import sys
from pathlib import Path


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from feature_coverage import (  # noqa: E402
    read_tsv_gz,
    summarize_feature_coverage,
    summarize_feature_coverage_rows,
    write_tsv_gz,
)


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
        }
    ]
    summaries = [
        {"gene_id": "1", "strategy": "test", "ortholog_gene_id": "101"},
        {"gene_id": "1", "strategy": "test", "ortholog_gene_id": "102"},
    ]
    segments = [
        {
            "gene_id": "1",
            "strategy": "test",
            "ortholog_gene_id": "101",
            "target_start0": "0",
            "target_end0": "5",
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
    write_tsv_gz(features_path, list(features[0]), features)
    write_tsv_gz(summaries_path, list(summaries[0]), summaries)
    write_tsv_gz(segments_path, list(segments[0]), segments)

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

    assert row_count == path_count == 1
    assert gzip.open(rows_output, "rt").read() == gzip.open(paths_output, "rt").read()
    assert read_tsv_gz(rows_output)[0] == {
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
        "ortholog_count": "2",
        "orthologs_covered": "2",
        "covered_bases": "10",
        "coverage_breadth": "1.000000",
        "depth_bases": "10",
        "mean_depth": "1.000000",
    }
