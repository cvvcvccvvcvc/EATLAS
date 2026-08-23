from __future__ import annotations

import csv
import gzip
import shutil
from pathlib import Path

import pytest


from analytics.derivations.feature_coverage import (
    load_snv_site_depth,
    summarize_feature_coverage,
    summarize_feature_coverage_rows,
    write_snv_site_depth,
    write_snv_taxonomic_depth,
)
from analytics.derivations.taxonomy import COUNT_KEYS
from genomics.taxonomy import TAXONOMY_FIELDS


pytestmark = pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools is not installed")


def write_tsv_gz(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str] | None = None,
) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames or list(rows[0]), delimiter="\t")
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


def test_site_depth_counts_distinct_orthologs_across_segments(tmp_path: Path) -> None:
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

    output = tmp_path / "snv_site_depth.tsv.gz"
    write_snv_site_depth(
        [segments_path],
        [
            {
                "gene_id": "1",
                "strategy": "test",
                "target_start0": "4",
            },
            {
                "gene_id": "1",
                "strategy": "test",
                "target_start0": "5",
            },
        ],
        output,
        tmp_path,
    )
    counts = load_snv_site_depth(output)

    assert counts == {
        ("1", "test", 4): 2,
        ("1", "test", 5): 3,
    }


def test_taxonomic_site_depth_collapses_members_by_rank(tmp_path: Path) -> None:
    taxonomy = tmp_path / "taxonomy.tsv.gz"
    write_tsv_gz(
        taxonomy,
        [
            {
                "tax_id": "11",
                "taxonomy_status": "resolved",
                "species_id": "11",
                "genus_id": "10",
                "family_id": "9",
                "order_id": "8",
                "lineage_tax_ids": "2759,33208,7742,32523,32524,40674,11",
            },
            {
                "tax_id": "12",
                "taxonomy_status": "resolved",
                "species_id": "12",
                "genus_id": "10",
                "family_id": "9",
                "order_id": "8",
                "lineage_tax_ids": "2759,33208,7742,32523,32524,40674,12",
            },
        ],
        TAXONOMY_FIELDS,
    )
    segments = tmp_path / "segments.tsv.gz"
    write_tsv_gz(
        segments,
        [
            {
                "gene_id": "1",
                "strategy": "s1",
                "ortholog_gene_id": ortholog,
                "tax_id": tax_id,
                "target_start0": 0,
                "target_end0": 10,
            }
            for ortholog, tax_id in [("101", "11"), ("102", "12")]
        ],
    )
    output = tmp_path / "depth.tsv.gz"

    count = write_snv_taxonomic_depth(
        [segments],
        [{"gene_id": "1", "strategy": "s1", "target_start0": 4}],
        taxonomy,
        output,
        tmp_path,
    )

    assert count == 1
    row = read_tsv_gz(output)[0]
    assert set(COUNT_KEYS) <= set(row)
    assert row["all__ortholog"] == "2"
    assert row["all__species"] == "2"
    assert row["all__genus"] == "1"
    assert row["mammalia__family"] == "1"
    assert row["primates__ortholog"] == "0"


def test_taxonomic_site_depth_sorts_gene_id_prefixes_by_output_columns(
    tmp_path: Path,
) -> None:
    taxonomy = tmp_path / "taxonomy.tsv.gz"
    write_tsv_gz(
        taxonomy,
        [
            {
                "tax_id": "11",
                "taxonomy_status": "resolved",
                "species_id": "11",
                "genus_id": "10",
                "family_id": "9",
                "order_id": "8",
                "lineage_tax_ids": "2759,33208,7742,32523,32524,40674,11",
            }
        ],
        TAXONOMY_FIELDS,
    )
    segments = tmp_path / "segments.tsv.gz"
    write_tsv_gz(
        segments,
        [
            {
                "gene_id": gene_id,
                "strategy": "s1",
                "ortholog_gene_id": f"{gene_id}01",
                "tax_id": "11",
                "target_start0": 0,
                "target_end0": 20,
            }
            for gene_id in ("466", "4665")
        ],
    )
    output = tmp_path / "depth.tsv.gz"

    count = write_snv_taxonomic_depth(
        [segments],
        [
            {"gene_id": "466", "strategy": "s1", "target_start0": 13},
            {"gene_id": "4665", "strategy": "s1", "target_start0": 13},
        ],
        taxonomy,
        output,
        tmp_path,
    )

    assert count == 2
    assert [row["gene_id"] for row in read_tsv_gz(output)] == ["466", "4665"]
