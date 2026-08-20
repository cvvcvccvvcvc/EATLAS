from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from analytics.analyses.minimap_concordance import (
    ASM10,
    ASM20,
    MinimapConcordanceAnalysis,
    compute_minimap_candidate_summary,
    minimap_group_eligibility,
    minimap_group_memberships,
)
from analytics.reporting.minimap_concordance import minimap_candidate_view


duckdb = pytest.importorskip("duckdb")


def test_minimap_candidate_summary_separates_union_intersection_and_only(
    tmp_path: Path,
) -> None:
    path = tmp_path / "filter_scores.parquet"
    frame = pd.DataFrame(
        [
            {
                "variant_key": "both",
                "strategy": ASM10,
                "variant_type": "snv",
                "gnomad_status": "found",
            },
            {
                "variant_key": "both",
                "strategy": ASM20,
                "variant_type": "snv",
                "gnomad_status": "found",
            },
            {
                "variant_key": "asm10",
                "strategy": ASM10,
                "variant_type": "snv",
                "gnomad_status": "not_found",
            },
            {
                "variant_key": "asm20",
                "strategy": ASM20,
                "variant_type": "snv",
                "gnomad_status": "lookup_failed",
            },
        ]
    )
    with duckdb.connect() as connection:
        connection.register("scores", frame)
        connection.execute("COPY scores TO ? (FORMAT PARQUET)", [str(path)])

    summary = compute_minimap_candidate_summary(path).set_index("group_key")

    assert summary.loc["minimap2_either", "variant_count"] == 3
    assert summary.loc["minimap2_both", "variant_count"] == 1
    assert summary.loc["minimap2_asm10_only", "variant_count"] == 1
    assert summary.loc["minimap2_asm20_only", "variant_count"] == 1
    assert summary.loc["minimap2_both", "gnomad_found_fraction"] == 1.0
    assert summary.loc["minimap2_either", "allele_fraction"] == 1.0


def test_minimap_membership_and_only_group_eligibility_are_strict() -> None:
    observed = {
        (ASM10, "snv"): {"shared", "ten"},
        (ASM20, "snv"): {"shared", "twenty"},
    }

    memberships = minimap_group_memberships(observed)
    eligibility = minimap_group_eligibility(
        {ASM10: {"1", "2"}, ASM20: {"2", "3"}}
    )

    assert memberships[("minimap2_either", "snv")] == {"shared", "ten", "twenty"}
    assert memberships[("minimap2_both", "snv")] == {"shared"}
    assert memberships[("minimap2_asm10_only", "snv")] == {"ten"}
    assert eligibility["minimap2_either"] == {"1", "2", "3"}
    assert eligibility["minimap2_both"] == {"2"}
    assert eligibility["minimap2_asm10_only"] == {"2"}


def test_minimap_candidate_view_uses_safe_payload(tmp_path: Path) -> None:
    summary = pd.DataFrame(
        [
            {
                "group_key": "minimap2_either",
                "variant_type": "snv",
                "variant_count": 3,
                "gnomad_found_count": 1,
                "gnomad_eligible_count": 2,
                "gnomad_lookup_failed_count": 1,
                "allele_fraction": 1.0,
                "gnomad_found_fraction": 0.5,
            }
        ]
    )
    analysis = MinimapConcordanceAnalysis(True, "", summary, None)

    html = minimap_candidate_view(analysis)

    assert "Either preset" in html
    assert 'data-role="variant-type"' in html
    assert "NaN" not in html
