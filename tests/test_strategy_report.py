from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from analytics.core.negative_controls import TargetSpaceNullAnalysis
from analytics.strategy_report import (
    RunInputs,
    build_variant_sections,
    build_target_space_null_sections,
    build_target_space_null_qc_sections,
    clinvar_association_view,
    dataframe_records,
    format_table_dataframe,
    gnomad_stratification_figure,
    validate_report_inputs,
)


def test_report_preflight_accepts_compact_production_contract(tmp_path: Path) -> None:
    def write_table(name: str, columns: list[str]) -> Path:
        path = tmp_path / name
        pd.DataFrame(columns=columns).to_csv(path, sep="\t", index=False, compression="gzip")
        return path

    annotations = write_table(
        "annotations.tsv.gz",
        [
            "variant_key", "gene_id", "event_type", "ref", "alt", "lookup_status",
            "strategies", "support_row_count", "support_ortholog_count", "clinvar_id",
            "clinvar_sig", "clinvar_review_stars", "clinvar_scv_count", "gnomad_af", "gnomad_csq",
        ],
    )
    genes = write_table("genes.tsv.gz", ["gene_id", "chromosome", "begin", "end", "sequence_length"])
    features = write_table(
        "features.tsv.gz", ["gene_id", "feature_type", "target_start0", "target_end0"]
    )
    coverage = write_table("coverage.tsv.gz", ["gene_id", "strategy", "feature_type"])
    summary = write_table(
        "summary.tsv.gz",
        ["strategy", "gene_count", "summary_row_count", "aligned_summary_row_count", "event_count"],
    )
    targets = tmp_path / "targets"
    targets.mkdir()
    inputs = RunInputs(
        tmp_path,
        genes,
        features,
        targets,
        annotations,
        tmp_path / "annotation_manifest.json",
        tmp_path / "annotation_failures.tsv.gz",
        coverage,
        tmp_path / "alignment_segments.tsv.gz",
        tmp_path / "alignment_manifest.json",
        summary,
    )

    validate_report_inputs(inputs)


def test_single_strategy_candidate_profile_loads_plotly_without_overlap() -> None:
    summary = SimpleNamespace(
        overlap=None,
        event_counts=pd.DataFrame(
            [{"strategy": "s1", "event_type": "snv", "Variant_Count": 10}]
        ),
        target_context_counts=pd.DataFrame(),
    )
    stats = pd.DataFrame([{"Strategy": "s1", "Unique Variants": 10}])

    html = "".join(build_variant_sections(summary, stats, include_plotly=True))

    assert "cdn.plot.ly" in html
    assert "Variant type composition by strategy" in html


def test_gnomad_stratification_places_found_and_not_found_bars_side_by_side() -> None:
    counts = pd.DataFrame(
        [
            {"strategy": "s1", "gnomad_status": "found", "kind": "SNV", "Variant_Count": 3},
            {"strategy": "s1", "gnomad_status": "found", "kind": "INDEL", "Variant_Count": 1},
            {"strategy": "s1", "gnomad_status": "not_found", "kind": "SNV", "Variant_Count": 1},
            {"strategy": "s1", "gnomad_status": "not_found", "kind": "INDEL", "Variant_Count": 3},
        ]
    )

    figure = gnomad_stratification_figure(counts, "kind", ["SNV", "INDEL"], ["s1"], "Test")

    assert figure.layout.barmode == "stack"
    assert list(figure.data[0].x[0]) == ["s1", "s1"]
    assert list(figure.data[0].x[1]) == ["Found", "Not found"]
    assert list(figure.data[0].y) == [0.75, 0.25]


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

    assert "Matched Control" in html
    assert "Sampled / matched focal SNVs" in html
    assert "2 / 2" in html
    assert "Matched Callable" not in html
    assert "Same-Position" not in html
    assert "Strategy Summary" not in html
    assert "VEP release" not in html
    assert "Exact alleles found in gnomAD" in html
    assert "gnomAD allele frequency among exact hits" in html
    assert "Exact alleles found in ClinVar" in html
    assert "ClinVar class composition" in html

    qc_html = "".join(build_target_space_null_qc_sections(analysis))
    assert "Matched-control QC" in qc_html
    assert "Strategy summary" in qc_html
    assert "VEP release" in qc_html


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


def test_clinvar_association_serializes_sparse_results_and_has_all_controls() -> None:
    unadjusted = pd.DataFrame(
        [
            {
                "strategy": "s1",
                "variant_type": "snv",
                "consequence": "missense",
                "odds_ratio": float("inf"),
                "ci_low": float("nan"),
                "ci_high": float("nan"),
                "fisher_p": float("nan"),
                "fisher_q": float("nan"),
                "usable_rows": 10,
                "benign_observed": 5,
                "pathogenic_observed": 0,
                "benign_not_observed": 4,
                "pathogenic_not_observed": 1,
                "status": "not_estimable",
                "reason": "Sparse data",
            }
        ]
    )
    validation = SimpleNamespace(
        unadjusted=unadjusted,
        fixed_adjusted=pd.DataFrame(),
        continuous=pd.DataFrame(),
        fixed_bins=pd.DataFrame(),
        distributions=pd.DataFrame(),
    )
    html = clinvar_association_view(validation)

    assert 'data-role="strategy"' in html
    assert 'data-role="mode"' in html
    assert 'data-role="variant-type"' in html
    assert 'data-role="consequence"' in html
    assert "phyloP fixed bands" in html
    assert "phyloP continuous" in html
    assert "Missense" in html
    assert "Infinity" not in html
    assert "NaN" not in html


def test_dataframe_records_replaces_nonfinite_values() -> None:
    records = dataframe_records(pd.DataFrame({"value": [1.0, float("inf"), float("nan")]}))
    assert records == [{"value": 1.0}, {"value": "inf"}, {"value": None}]
