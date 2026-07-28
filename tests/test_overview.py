from __future__ import annotations

from types import SimpleNamespace

import pandas as pd

from analytics.strategy_report import (
    alignment_summary_for_report,
    build_overview,
    overview_strategy_table,
    target_gene_coverage_for_report,
)


def test_overview_uses_clear_strategy_metrics_and_gene_level_coverage() -> None:
    alignment = pd.DataFrame(
        {
            "strategy": ["s1", "s2"],
            "summary_row_count": [10, 8],
            "aligned_summary_row_count": [8, 4],
            "event_count": [100, 200],
        }
    )
    aligned = alignment_summary_for_report(alignment)
    assert aligned["Orthologs evaluated"].tolist() == [10, 8]
    assert aligned["Orthologs aligned"].tolist() == [8, 4]
    assert aligned["Orthologs aligned %"].tolist() == [0.8, 0.5]

    coverage = pd.DataFrame(
        [
            {"gene_id": "1", "strategy": "s1", "feature_type": "gene", "length_bp": 100, "covered_bases": 50},
            {"gene_id": "2", "strategy": "s1", "feature_type": "gene", "length_bp": 300, "covered_bases": 300},
            {"gene_id": "1", "strategy": "s1", "feature_type": "cds", "length_bp": 10, "covered_bases": 0},
            {"gene_id": "1", "strategy": "s2", "feature_type": "gene", "length_bp": 100, "covered_bases": 25},
        ]
    )
    target_coverage = target_gene_coverage_for_report(coverage).set_index("Strategy")
    assert target_coverage.loc["s1", "Target bases covered %"] == 0.875
    assert target_coverage.loc["s2", "Target bases covered %"] == 0.25

    strategy_stats = pd.DataFrame(
        {
            "Strategy": ["s1", "s2"],
            "Unique Variants": [10, 20],
            "Found in ClinVar": [2, 5],
            "gnomAD Found": [4, 8],
            "gnomAD Eligible": [10, 20],
            "Genes with result": [2, 1],
            "Orthologs evaluated": [10, 8],
            "Orthologs aligned": [8, 4],
        }
    )
    summary = SimpleNamespace(
        unique_contribution=pd.DataFrame(
            {"Strategy": ["s1", "s2"], "Unique To Strategy": [1, 5]}
        ),
        unique_variant_count=25,
        all_strategy_variant_count=5,
        strategies=["s1", "s2"],
        gene_count=2,
    )

    table = overview_strategy_table(summary, coverage, strategy_stats, input_gene_count=2)
    assert table["Strategy"].tolist() == ["s2", "s1"]
    assert table.loc[0, "Only this strategy"] == "5 (25.0%)"
    assert table.loc[1, "Orthologs aligned"] == "8 / 10 (80.0%)"

    html = "".join(
        build_overview(
            summary,
            coverage,
            strategy_stats,
            {"failure_count": 3},
            input_gene_count=2,
        )
    )
    assert "Candidates found by all strategies" in html
    assert "5 (20.0%)" in html
    assert "Raw support events" not in html
    assert "Found in ClinVar" not in html
    assert "Found in gnomAD" not in html
    assert "Orthologs aligned" in html
