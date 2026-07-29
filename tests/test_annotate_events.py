from __future__ import annotations

import sys
from pathlib import Path

import pytest


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from annotate_events import (  # noqa: E402
    VARIANT_ANNOTATION_FIELDS,
    VARIANT_STRATEGY_SUPPORT_FIELDS,
    add_strategy_support,
    build_variant_strategy_support,
    event_vcf_key,
    iter_variant_strategy_snv_sites,
    variant_aggregate_key,
)
from finalize_annotation_partitions import (  # noqa: E402
    COUNT_FIELDS,
    merge_gnomad_shared_cache,
)


def test_variant_annotation_schema_is_analysis_ready() -> None:
    assert VARIANT_ANNOTATION_FIELDS == [
        "variant_key",
        "gene_id",
        "event_type",
        "ref",
        "alt",
        "lookup_status",
        "support_row_count",
        "support_ortholog_count",
        "strategies",
        "clinvar_sig",
        "clinvar_revstat",
        "clinvar_review_stars",
        "clinvar_review_stars_status",
        "clinvar_id",
        "clinvar_allele_id",
        "clinvar_scv_count",
        "clinvar_hgvs",
        "clinvar_disease",
        "clinvar_variant_type",
        "gnomad_af",
        "gnomad_af_source",
        "gnomad_csq",
    ]


def test_partitioned_manifest_keeps_non_concrete_exclusion_count() -> None:
    assert "excluded_non_concrete_event_count" in COUNT_FIELDS


def test_variant_strategy_support_schema_includes_site_depth() -> None:
    assert VARIANT_STRATEGY_SUPPORT_FIELDS[-1] == "site_aligned_ortholog_count"
    assert "variant_strategy_site_depth_count" in COUNT_FIELDS


def test_partitioned_manifest_aggregates_shared_gnomad_cache_metrics(tmp_path: Path) -> None:
    identity = {
        "enabled": True,
        "directory": "/cache/gnomad",
        "schema_version": 1,
        "dataset": "gnomad_r4",
        "reference_genome": "GRCh38",
        "tile_size_bp": 25_000,
    }
    partitions = [
        (
            tmp_path / "one",
            {
                "gnomad_shared_cache": {
                    **identity,
                    "tile_hit_count": 2,
                    "tile_miss_count": 3,
                    "tile_write_count": 3,
                    "fetch_batch_count": 1,
                }
            },
        ),
        (
            tmp_path / "two",
            {
                "gnomad_shared_cache": {
                    **identity,
                    "tile_hit_count": 5,
                    "tile_miss_count": 0,
                    "tile_write_count": 0,
                    "fetch_batch_count": 0,
                }
            },
        ),
    ]

    merged = merge_gnomad_shared_cache(partitions)

    assert merged is not None
    assert merged["tile_hit_count"] == 7
    assert merged["tile_miss_count"] == 3
    assert merged["tile_write_count"] == 3
    assert merged["fetch_batch_count"] == 1


def test_variant_strategy_support_counts_distinct_orthologs() -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "_support_by_strategy": {},
    }
    for strategy, ortholog in [
        ("s1", "101"),
        ("s1", "101"),
        ("s1", "102"),
        ("s2", "101"),
    ]:
        add_strategy_support(
            aggregate,
            {
                "strategy": strategy,
                "ortholog_gene_id": ortholog,
            },
        )

    rows, missing_key_count = build_variant_strategy_support([aggregate])

    assert missing_key_count == 0
    assert rows == [
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "strategy": "s1",
            "alt_support_row_count": 3,
            "alt_support_ortholog_count": 2,
            "site_aligned_ortholog_count": "",
        },
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "strategy": "s2",
            "alt_support_row_count": 1,
            "alt_support_ortholog_count": 1,
            "site_aligned_ortholog_count": "",
        },
    ]


def test_snv_support_uses_site_aligned_depth() -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "event_type": "snv",
        "target_start0": "9",
        "_support_by_strategy": {},
    }
    add_strategy_support(
        aggregate,
        {"strategy": "s1", "ortholog_gene_id": "101"},
    )

    assert list(iter_variant_strategy_snv_sites([aggregate])) == [
        {
            "gene_id": "1",
            "strategy": "s1",
            "target_start0": "9",
        }
    ]
    rows, _missing_key_count = build_variant_strategy_support(
        [aggregate],
        {("1", "s1", 9): 4},
    )
    assert rows[0]["site_aligned_ortholog_count"] == 4


def test_snv_support_rejects_alt_count_above_site_depth() -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "event_type": "snv",
        "_support_by_strategy": {},
    }
    add_strategy_support(
        aggregate,
        {"strategy": "s1", "support_ortholog_count": "2"},
    )

    with pytest.raises(ValueError, match="exceeds site-aligned"):
        build_variant_strategy_support(
            [aggregate],
            {("1", "s1", 0): 1},
        )


def test_variant_strategy_support_accepts_single_strategy_compact_counts() -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "_support_by_strategy": {},
    }
    add_strategy_support(
        aggregate,
        {
            "strategy": "s1",
            "support_row_count": "10",
            "support_ortholog_count": "3",
        },
    )

    rows, _missing_key_count = build_variant_strategy_support([aggregate])

    assert rows[0]["alt_support_row_count"] == 10
    assert rows[0]["alt_support_ortholog_count"] == 3


def test_variant_strategy_support_rejects_cross_strategy_compact_counts() -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "_support_by_strategy": {},
    }

    with pytest.raises(ValueError, match="aggregated across multiple strategies"):
        add_strategy_support(
            aggregate,
            {
                "strategies": "s1,s2",
                "support_row_count": "10",
                "support_ortholog_count": "3",
            },
        )


def test_canonical_variant_key_collapses_raw_indel_representations() -> None:
    left = {
        "gene_id": "1",
        "event_type": "del",
        "target_start0": "10",
        "target_end0": "11",
        "ref": "A",
        "alt": "",
    }
    right = {
        "gene_id": "1",
        "event_type": "del",
        "target_start0": "11",
        "target_end0": "12",
        "ref": "A",
        "alt": "",
    }

    assert variant_aggregate_key(left, "1:100:AA>A") == variant_aggregate_key(
        right,
        "1:100:AA>A",
    )
    assert variant_aggregate_key(left, "") != variant_aggregate_key(right, "")


@pytest.mark.parametrize(
    ("ref", "alt"),
    [
        ("G", "Y"),
        ("T", "R"),
        ("A", "N"),
        ("", "."),
    ],
)
def test_event_vcf_key_rejects_non_concrete_alleles(ref: str, alt: str) -> None:
    key, status = event_vcf_key(
        {
            "gene_id": "1",
            "event_type": "snv",
            "genomic_accession": "NC_000001.11",
            "genomic_start1": "100",
            "target_start0": "0",
            "ref": ref,
            "alt": alt,
        },
        {},
    )

    assert key is None
    assert status == "non_concrete_allele"


def test_event_vcf_key_keeps_concrete_alleles() -> None:
    key, status = event_vcf_key(
        {
            "gene_id": "1",
            "event_type": "snv",
            "genomic_accession": "NC_000001.11",
            "genomic_start1": "100",
            "target_start0": "0",
            "ref": "A",
            "alt": "G",
        },
        {},
    )

    assert key == ("1", 100, "A", "G")
    assert status == "raw_no_context"
