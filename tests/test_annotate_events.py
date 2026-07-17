from __future__ import annotations

import sys
from pathlib import Path

import pytest


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from annotate_events import (  # noqa: E402
    VARIANT_ANNOTATION_FIELDS,
    add_strategy_support,
    build_variant_strategy_support,
    variant_aggregate_key,
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
        },
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "strategy": "s2",
            "alt_support_row_count": 1,
            "alt_support_ortholog_count": 1,
        },
    ]


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
