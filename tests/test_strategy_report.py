from __future__ import annotations

from pathlib import Path

import pandas as pd

from analytics.core.negative_controls import TargetSpaceNullAnalysis
from analytics.strategy_report import (
    build_target_space_null_sections,
    conservation_selector_view,
    dataframe_records,
    format_table_dataframe,
)


def test_target_space_null_section_reports_consequence_matched_design(tmp_path: Path) -> None:
    analysis = TargetSpaceNullAnalysis(
        summary=pd.DataFrame(
            [
                {
                    "strategy": "s1",
                    "matched_focals": 2,
                    "observed_median": 0.5,
                    "null_median": 1.0,
                    "null_ci_low": 0.8,
                    "null_ci_high": 1.2,
                    "median_difference": -0.5,
                }
            ]
        ),
        consequence_summary=pd.DataFrame(),
        ecdf=pd.DataFrame(),
        gnomad_summary=pd.DataFrame(
            [
                {
                    "strategy": "s1",
                    "metric": "found_fraction",
                    "observed_value": 0.4,
                    "null_value": 0.2,
                    "null_ci_low": 0.1,
                    "null_ci_high": 0.3,
                },
                {
                    "strategy": "s1",
                    "metric": "median_af",
                    "observed_value": 0.001,
                    "null_value": 0.0005,
                    "null_ci_low": 0.0001,
                    "null_ci_high": 0.002,
                },
            ]
        ),
        clinvar_summary=pd.DataFrame(
            [
                {
                    "strategy": "s1",
                    "observed_value": 0.1,
                    "null_value": 0.05,
                    "null_ci_low": 0.02,
                    "null_ci_high": 0.08,
                }
            ]
        ),
        clinvar_class_summary=pd.DataFrame(
            [
                {
                    "strategy": "s1",
                    "clinvar_class": category,
                    "observed_value": observed,
                    "null_value": null,
                    "null_ci_low": max(0.0, null - 0.05),
                    "null_ci_high": min(1.0, null + 0.05),
                }
                for category, observed, null in [
                    ("B/LB", 0.5, 0.4),
                    ("P/LP", 0.1, 0.2),
                    ("VUS", 0.3, 0.3),
                    ("Other", 0.1, 0.1),
                ]
            ]
        ),
        manifest={
            "inputs": {"sample_size_per_strategy": 25_000},
            "sampled_focal_count": 2,
            "vep_annotated_focal_count": 2,
            "matched_focal_count": 2,
            "focal_vep": {"release": "116"},
            "conservation": {"status": "complete"},
            "external_evidence": {"gnomad": {"failed_region_count": 0}},
        },
        manifest_path=tmp_path / "manifest.json",
        matched_path=tmp_path / "target_space_null.snv.tsv.gz",
        conservation_path=tmp_path / "target_space_null.phyloP100way.tsv.gz",
        vep_cache_path=tmp_path / "vep.sqlite",
        external_evidence_path=tmp_path / "target_space_null.external_evidence.tsv.gz",
        external_evidence_manifest_path=tmp_path / "target_space_null.external_evidence.manifest.json",
        resamples=1_000,
    )

    html = "".join(build_target_space_null_sections(analysis, include_plotly=False))

    assert "Target-Space Null" in html
    assert "same genomic REF&gt;ALT substitution" in html
    assert "Matched Callable" not in html
    assert "Same-Position" not in html
    assert "<details><summary>Strategy Summary</summary>" in html
    assert "Exact alleles found in gnomAD" in html
    assert "gnomAD allele frequency among exact hits" in html
    assert "Exact alleles found in ClinVar" in html
    assert "ClinVar class composition" in html


def test_target_space_null_section_reports_disabled_state() -> None:
    html = "".join(
        build_target_space_null_sections(
            None,
            include_plotly=False,
            enabled=False,
        )
    )

    assert "was disabled for this report run" in html
    assert "--target-space-null" in html


def test_phyloP_quantiles_are_not_formatted_as_percentages() -> None:
    frame = pd.DataFrame(
        {
            "Background median Q2.5": [0.095],
            "Comparator rate Q2.5": [0.095],
        }
    )

    shown = format_table_dataframe(frame)

    assert shown.loc[0, "Background median Q2.5"] == "0.095"
    assert shown.loc[0, "Comparator rate Q2.5"] == "9.5%"


def test_conservation_selector_serializes_sparse_results_and_has_all_controls() -> None:
    primary = pd.DataFrame(
        [
            {
                "strategy": "s1",
                "variant_type": "snv",
                "consequence": "missense",
                "odds_ratio_mh": float("inf"),
                "ci_low": float("nan"),
                "ci_high": float("nan"),
                "cmh_p": float("nan"),
                "cmh_q": float("nan"),
                "usable_rows": 10,
                "status": "not_estimable",
                "reason": "Sparse data",
            }
        ]
    )
    detail = pd.DataFrame()

    html = conservation_selector_view(
        view_id="fixed-test",
        strategies=["s1"],
        primary=primary,
        detail=detail,
        mode="fixed",
    )

    assert 'data-role="strategy"' in html
    assert 'data-role="variant-type"' in html
    assert 'data-role="consequence"' in html
    assert "Missense" in html
    assert "Infinity" not in html
    assert "NaN" not in html


def test_dataframe_records_replaces_nonfinite_values() -> None:
    records = dataframe_records(pd.DataFrame({"value": [1.0, float("inf"), float("nan")]}))
    assert records == [{"value": 1.0}, {"value": "inf"}, {"value": None}]
