from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from analytics.analyses.candidate_conservation import CandidateConservation
from analytics.analyses.matched_control import TargetSpaceNullAnalysis
from analytics.annotation.consequences import UNANNOTATED_CONSEQUENCE
from analytics.io.run_inputs import RunInputs, resolve_run_inputs, validate_report_inputs
from analytics.reporting.components import dataframe_records, format_table_dataframe
from analytics.reporting.config import CONSEQUENCE_GROUP_COLORS
from analytics.reporting.conservation import (
    clinvar_association_view,
    hidden_clinvar_association_views,
)
from analytics.reporting.matched_control import (
    build_target_space_null_qc_sections,
    build_target_space_null_sections,
)
from analytics.reporting.document import render_html
from analytics.reporting.ortholog_evidence import (
    build_ortholog_evidence_sections,
    ortholog_evidence_distribution_figure,
    ortholog_evidence_distribution_stats,
    ortholog_evidence_figure,
)
from analytics.reporting.overview import overview_strategy_table
from analytics.reporting.qc import vep_qc_tables
from analytics.reporting.variant_profile import (
    build_gnomad_stratification_sections,
    build_variant_sections,
    gene_variant_distribution_counts,
    gene_variant_distribution_figure,
    gnomad_stratification_figure,
    group_consequence_counts,
    pathogenic_variant_table,
    top_gene_contribution_counts,
    top_gene_contribution_figure,
)
from analytics.strategy_report import _default_phylop_bigwig


def test_local_phylop_default_prefers_environment_or_existing_gaph_root(
    tmp_path: Path,
    monkeypatch,
) -> None:
    explicit = tmp_path / "explicit.bw"
    monkeypatch.setenv("GAPH_PHYLOP_BIGWIG", str(explicit))
    assert _default_phylop_bigwig() == explicit

    monkeypatch.delenv("GAPH_PHYLOP_BIGWIG")
    monkeypatch.setenv("GAPH_ROOT", str(tmp_path))
    partial = tmp_path / "reference" / "ucsc" / "hg38.phyloP100way.bw.partial"
    partial.parent.mkdir(parents=True)
    partial.write_bytes(b"incomplete bigwig")
    assert _default_phylop_bigwig() is None
    discovered = tmp_path / "reference" / "ucsc" / "hg38.phyloP100way.bw"
    discovered.write_bytes(b"bigwig")
    assert _default_phylop_bigwig() == discovered


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
            "clinvar_sig", "clinvar_review_stars", "clinvar_scv_count", "gnomad_af",
            "vep_status", "vep_primary_consequence", "vep_consequence_terms",
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
    support = write_table(
        "support.tsv.gz",
        [
            "variant_key",
            "gene_id",
            "strategy",
            "alt_support_row_count",
            "alt_support_ortholog_count",
        ],
    )
    targets = tmp_path / "targets"
    targets.mkdir()
    inputs = RunInputs(
        run_dir=tmp_path,
        fetch_manifest_json=tmp_path / "fetch_manifest.json",
        genes_tsv=genes,
        target_features_tsv=features,
        target_sequences_dir=targets,
        variant_annotations_tsv=annotations,
        variant_strategy_support_tsv=support,
        ortholog_evidence_summary_tsv=tmp_path / "ortholog_evidence_summary.tsv.gz",
        annotation_manifest_json=tmp_path / "annotation_manifest.json",
        annotation_failures_tsv=tmp_path / "annotation_failures.tsv.gz",
        feature_coverage_tsv=coverage,
        alignment_segments_tsv=tmp_path / "alignment_segments.tsv.gz",
        alignment_manifest_json=tmp_path / "alignment_manifest.json",
        strategy_summary_tsv=summary,
        taxonomy_summary_tsv=tmp_path / "taxonomy_summary.tsv.gz",
    )

    validate_report_inputs(inputs)


def test_report_inputs_use_matching_completed_vep_artifact(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    annotation_dir = run_dir / "annotation"
    artifact_dir = run_dir / "analytics" / "vep_consequences"
    (run_dir / "fetch" / "sequences" / "targets").mkdir(parents=True)
    annotation_dir.mkdir()
    artifact_dir.mkdir(parents=True)
    source = annotation_dir / "variant_annotations.tsv.gz"
    pd.DataFrame(columns=["variant_key"]).to_csv(
        source, sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame(columns=["gene_id"]).to_csv(
        run_dir / "fetch" / "genes.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame(columns=["gene_id"]).to_csv(
        run_dir / "fetch" / "target_features.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    output = artifact_dir / "variant_annotations.vep.tsv.gz"
    pd.DataFrame(columns=["variant_key", "vep_status"]).to_csv(
        output, sep="\t", index=False, compression="gzip"
    )
    source_stat = source.stat()
    output_stat = output.stat()
    (artifact_dir / "manifest.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "source": {
                    "path": str(source.resolve()),
                    "size_bytes": source_stat.st_size,
                    "mtime_ns": source_stat.st_mtime_ns,
                },
                "output": {
                    "size_bytes": output_stat.st_size,
                    "mtime_ns": output_stat.st_mtime_ns,
                },
            }
        )
    )

    inputs = resolve_run_inputs(run_dir)

    assert inputs.variant_annotations_tsv == output


def test_report_inputs_require_finalized_vep_artifact(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    annotation_dir = run_dir / "annotation"
    (run_dir / "fetch" / "sequences" / "targets").mkdir(parents=True)
    annotation_dir.mkdir()
    pd.DataFrame(columns=["variant_key"]).to_csv(
        annotation_dir / "variant_annotations.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    with pytest.raises(FileNotFoundError, match="Missing finalized bulk VEP artifact"):
        resolve_run_inputs(run_dir)


def test_vep_qc_reports_candidate_and_clinvar_statuses() -> None:
    candidate_manifest = {
        "row_count": 100,
        "config": {"backend": "local", "release": "116"},
        "status_counts": {"ok": 99, "invalid_variant_key": 1},
    }
    validation = SimpleNamespace(
        manifest={
            "vep": {
                "allele_count": 20,
                "contract": {"backend": "local", "release": "116"},
                "status_counts": {"ok": 18, "no_target_gene": 2},
            }
        }
    )

    configuration, statuses = vep_qc_tables(candidate_manifest, validation)

    assert configuration["VEP release"].tolist() == ["116", "116"]
    assert statuses.to_dict(orient="records") == [
        {"Dataset": "Candidates", "VEP status": "ok", "Rows": 99, "Fraction": "99.000%"},
        {
            "Dataset": "Candidates",
            "VEP status": "invalid_variant_key",
            "Rows": 1,
            "Fraction": "1.000%",
        },
        {"Dataset": "ClinVar universe", "VEP status": "ok", "Rows": 18, "Fraction": "90.000%"},
        {
            "Dataset": "ClinVar universe",
            "VEP status": "no_target_gene",
            "Rows": 2,
            "Fraction": "10.000%",
        },
    ]


def test_unannotated_vep_rows_are_a_separate_grey_consequence_group() -> None:
    grouped = group_consequence_counts(
        pd.DataFrame(
            [
                {"strategy": "s1", "value": "missense_variant", "Variant_Count": 9},
                {"strategy": "s1", "value": UNANNOTATED_CONSEQUENCE, "Variant_Count": 1},
            ]
        )
    )

    unannotated = grouped[grouped["Consequence group"].astype(str).eq("Not annotated")]
    assert unannotated.iloc[0]["Variant_Count"] == 1
    assert unannotated.iloc[0]["Fraction"] == pytest.approx(0.1)
    assert CONSEQUENCE_GROUP_COLORS["Not annotated"] == "#d9d9d9"


def test_pathogenic_table_exposes_only_vep_consequence() -> None:
    variants = pd.DataFrame(
        [
            {
                "variant_id": "1:100:A>G",
                "gene_id": "1",
                "event_type": "snv",
                "clinvar_category": "P/LP",
                "clinvar_sig": "Pathogenic",
                "clinvar_review_stars": "2",
                "clinvar_revstat": "criteria_provided",
                "clinvar_scv_count": "3",
                "clinvar_id": "VCV1",
                "clinvar_allele_id": "1",
                "clinvar_disease": "Example disease",
                "clinvar_hgvs": "NC_000001.11:g.100A>G",
                "clinvar_variant_type": "single nucleotide variant",
                "gnomad_af": "0.01",
                "vep_primary_consequence": "missense_variant",
                "vep_status": "ok",
                "support_ortholog_mean": 2.0,
                "support_ortholog_min": 1,
                "support_ortholog_max": 3,
                "strategies": "s1",
            }
        ]
    )

    table = pathogenic_variant_table(variants)

    assert "gnomAD consequence" not in table.columns
    assert table.loc[0, "VEP consequence"] == "missense_variant"
    assert table.loc[0, "VEP status"] == "ok"


def test_single_strategy_candidate_profile_loads_plotly_without_overlap() -> None:
    summary = SimpleNamespace(
        overlap=None,
        gene_variant_counts=pd.DataFrame(
            [{"strategy": "s1", "gene_id": "1", "Variant_Count": 10}]
        ),
        event_counts=pd.DataFrame(
            [{"strategy": "s1", "event_type": "snv", "Variant_Count": 10}]
        ),
        target_context_counts=pd.DataFrame(),
    )
    stats = pd.DataFrame(
        [{"Strategy": "s1", "Unique Variants": 10, "Genes with result": 2}]
    )

    html = "".join(build_variant_sections(summary, stats))

    assert "cdn.plot.ly" not in html
    assert "Candidate variants per gene" in html
    assert "Top 5 contributing genes by strategy" in html
    assert "Variant type composition by strategy" in html


def test_gene_variant_distribution_includes_zero_candidate_genes() -> None:
    gene_counts = pd.DataFrame(
        [
            {"strategy": "s1", "gene_id": "1", "Variant_Count": 8},
            {"strategy": "s1", "gene_id": "2", "Variant_Count": 2},
        ]
    )
    stats = pd.DataFrame([{"Strategy": "s1", "Genes with result": 4}])

    distribution = gene_variant_distribution_counts(gene_counts, stats).set_index("Bin")

    assert int(distribution.loc["0", "Gene_Count"]) == 2
    assert int(distribution.loc["2-3", "Gene_Count"]) == 1
    assert int(distribution.loc["8-15", "Gene_Count"]) == 1
    assert distribution["Gene_Fraction"].sum() == 1.0

    figure = gene_variant_distribution_figure(gene_counts, stats)
    assert figure is not None
    assert [trace.type for trace in figure.data] == ["scatter"]


def test_top_gene_contribution_uses_strategy_denominator_and_equal_share() -> None:
    gene_counts = pd.DataFrame(
        [
            {"strategy": "s1", "gene_id": "1", "Variant_Count": 8},
            {"strategy": "s1", "gene_id": "2", "Variant_Count": 2},
        ]
    )
    stats = pd.DataFrame([{"Strategy": "s1", "Genes with result": 4}])

    top = top_gene_contribution_counts(gene_counts, stats, limit=1).iloc[0]

    assert top["gene_id"] == "1"
    assert top["Variant_Fraction"] == 0.8
    assert top["Equal_Share"] == 0.25
    assert top["Equal_Share_Ratio"] == 3.2
    assert top["Top_Share"] == 0.8


def test_top_gene_contribution_orders_bars_by_contribution() -> None:
    gene_counts = pd.DataFrame(
        [
            {"strategy": "s1", "gene_id": "30", "Variant_Count": 2},
            {"strategy": "s1", "gene_id": "10", "Variant_Count": 8},
            {"strategy": "s1", "gene_id": "20", "Variant_Count": 5},
        ]
    )
    stats = pd.DataFrame([{"Strategy": "s1", "Genes with result": 3}])

    figure = top_gene_contribution_figure(gene_counts, stats, limit=3)

    assert figure is not None
    assert list(figure.data[0].x[1]) == ["10", "20", "30"]


def test_report_document_loads_plotly_once() -> None:
    html = render_html([("one", "One", ["<div>plot</div>"]), ("two", "Two", [])])

    assert html.count("cdn.plot.ly") == 1


def test_ortholog_evidence_section_renders_three_context_heatmaps() -> None:
    cells = pd.DataFrame(
        [
            {
                "strategy": "s1",
                "target_context": context,
                "taxonomic_scope": "all",
                "evidence_unit": "ortholog",
                "quantile_count": quantile_count,
                "depth_bin": 0,
                "alt_bin": 0,
                "depth_label": "1-4",
                "alt_label": "1-2",
                "gnomad_found_count": 1,
                "gnomad_eligible_count": 2,
                "gnomad_found_fraction": 0.5,
            }
            for context in ["cds", "utr", "intron"]
            for quantile_count in [2, 4, 10]
        ]
    )
    distributions = pd.DataFrame(
        [
            {
                "strategy": "s1",
                "taxonomic_scope": "all",
                "evidence_unit": "ortholog",
                "metric": metric,
                "value": value,
                "variant_count": count,
            }
            for metric, values in {
                "site_aligned": [(1, 1), (4, 3)],
                "exact_alt": [(0, 2), (2, 2)],
            }.items()
            for value, count in values
        ]
    )
    summary = SimpleNamespace(
        ortholog_evidence_available=True,
        ortholog_evidence_cells=cells,
        ortholog_evidence_distributions=distributions,
        strategies=["s1", "precomputed_ensembl_92_mammals_epo_extended"],
    )
    taxonomy_summary = pd.DataFrame(
        [
            {
                "taxonomic_scope": "all",
                "evidence_unit": "ortholog",
                "ortholog_count": 100,
                "taxon_count": 20,
                "orthologs_per_gene_median": 50.0,
                "units_per_gene_median": 50.0,
                "unit_count": 100,
            }
        ]
    )

    figure = ortholog_evidence_figure(cells, "s1", 2)
    distribution_figure = ortholog_evidence_distribution_figure(distributions, "s1")
    distribution_stats = ortholog_evidence_distribution_stats(
        distributions,
        "s1",
        "all",
        "ortholog",
    )
    html = "".join(
        build_ortholog_evidence_sections(
            summary,
            taxonomy_summary=taxonomy_summary,
        )
    )

    assert len(figure.data) == 3
    assert len(distribution_figure.data) == 2
    assert all(float(trace.y[-1]) == 1.0 for trace in distribution_figure.data)
    assert distribution_stats == [
        {"label": "SNVs", "value": "4"},
        {"label": "Site-aligned median [IQR]", "value": "4 [1-4]"},
        {"label": "Exact-ALT median [IQR]", "value": "0 [0-2]"},
    ]
    assert [annotation.text for annotation in figure.layout.annotations] == ["CDS", "UTR", "Intron"]
    assert "Median" in html
    assert "Quartiles" in html
    assert "Deciles" in html
    assert "Taxonomic scope" in html
    assert "Evidence unit" in html
    assert "Exact-ALT support" in html
    assert "Evidence distributions" in html
    assert "Cumulative SNVs at or below X" in html
    assert "Median selected orthologs/gene: 50.0" in html
    assert "taxonomy unavailable" in html
    assert "Plotly.react('ortholog-evidence-plot'" in html


def test_ortholog_evidence_section_explains_legacy_output() -> None:
    summary = SimpleNamespace(
        ortholog_evidence_available=False,
        ortholog_evidence_cells=pd.DataFrame(),
        ortholog_evidence_distributions=pd.DataFrame(),
    )

    html = "".join(build_ortholog_evidence_sections(summary))

    assert "predates site-aligned ortholog depth" in html


def test_overview_reports_strategy_gene_completeness_and_gnomad_eligible_denominator() -> None:
    summary = SimpleNamespace(
        unique_contribution=pd.DataFrame([{"Strategy": "s1", "Unique To Strategy": 2}])
    )
    coverage = pd.DataFrame(
        [
            {
                "strategy": "s1",
                "feature_type": "gene",
                "length_bp": 100,
                "covered_bases": 80,
            }
        ]
    )
    stats = pd.DataFrame(
        [
            {
                "Strategy": "s1",
                "Unique Variants": 10,
                "gnomAD Found": 4,
                "gnomAD Eligible": 8,
                "Found in ClinVar": 3,
                "Orthologs aligned": 7,
                "Orthologs evaluated": 10,
                "Genes with result": 9,
            }
        ]
    )

    table = overview_strategy_table(summary, coverage, stats, input_gene_count=10)

    assert table.loc[0, "Genes with result"] == "9 / 10 (90.0%)"
    assert table.loc[0, "gnomAD matches"] == "4 (50.0%)"


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


def test_gnomad_stratification_uses_vep_consequence_groups(tmp_path: Path) -> None:
    summary = SimpleNamespace(
        consequence_source="Ensembl VEP",
        gnomad_event_counts=pd.DataFrame(),
        gnomad_context_counts=pd.DataFrame(),
        gnomad_consequence_counts=pd.DataFrame(
            [
                {
                    "strategy": "s1",
                    "gnomad_status": "found",
                    "value": "missense_variant",
                    "Variant_Count": 3,
                },
                {
                    "strategy": "s1",
                    "gnomad_status": "not_found",
                    "value": UNANNOTATED_CONSEQUENCE,
                    "Variant_Count": 1,
                },
            ]
        ),
    )
    candidate = CandidateConservation(
        distributions_path=tmp_path / "distributions.tsv.gz",
        histograms_path=tmp_path / "histograms.tsv.gz",
        manifest_path=tmp_path / "manifest.json",
        distributions=pd.DataFrame(),
        histograms=pd.DataFrame(),
        manifest={},
    )

    sections = build_gnomad_stratification_sections(
        summary,
        pd.DataFrame([{"Strategy": "s1", "gnomAD found %": 0.5}]),
        candidate,
    )
    rendered = "".join(sections)

    assert "Ensembl VEP consequence: gnomAD hits versus non-hits" in rendered
    assert "Not annotated" in rendered
    assert "#d9d9d9" in rendered
    assert "Not computed" not in rendered


def test_target_space_null_section_reports_consequence_matched_design(tmp_path: Path) -> None:
    analysis = TargetSpaceNullAnalysis(
        summary=pd.DataFrame(
            [
                {
                    "strategy": "s1",
                    "matched_focals": 2,
                    "observed_median": 0.5,
                    "observed_ci_low": 0.3,
                    "observed_ci_high": 0.7,
                    "null_median": 1.0,
                    "null_ci_low": 0.8,
                    "null_ci_high": 1.2,
                    "median_difference": -0.5,
                    "difference_ci_low": -0.8,
                    "difference_ci_high": -0.2,
                    "valid_resamples": 1_000,
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
                    "matched_focals": 2,
                    "observed_value": 0.4,
                    "observed_ci_low": 0.2,
                    "observed_ci_high": 0.6,
                    "null_value": 0.2,
                    "null_ci_low": 0.1,
                    "null_ci_high": 0.3,
                    "difference": 0.2,
                    "difference_ci_low": -0.1,
                    "difference_ci_high": 0.5,
                    "valid_resamples": 1_000,
                },
                {
                    "strategy": "s1",
                    "metric": "median_af",
                    "matched_focals": 2,
                    "observed_value": 0.001,
                    "observed_ci_low": 0.0005,
                    "observed_ci_high": 0.002,
                    "null_value": 0.0005,
                    "null_ci_low": 0.0001,
                    "null_ci_high": 0.002,
                    "difference": 0.0005,
                    "difference_ci_low": -0.0002,
                    "difference_ci_high": 0.001,
                    "valid_resamples": 1_000,
                },
            ]
        ),
        clinvar_summary=pd.DataFrame(
            [
                {
                    "strategy": "s1",
                    "matched_focals": 2,
                    "observed_value": 0.1,
                    "observed_ci_low": 0.05,
                    "observed_ci_high": 0.15,
                    "null_value": 0.05,
                    "null_ci_low": 0.02,
                    "null_ci_high": 0.08,
                    "difference": 0.05,
                    "difference_ci_low": -0.02,
                    "difference_ci_high": 0.12,
                    "valid_resamples": 1_000,
                }
            ]
        ),
        clinvar_class_summary=pd.DataFrame(
            [
                {
                    "strategy": "s1",
                    "clinvar_class": category,
                    "matched_focals": 2,
                    "observed_value": observed,
                    "observed_ci_low": max(0.0, observed - 0.05),
                    "observed_ci_high": min(1.0, observed + 0.05),
                    "null_value": null,
                    "null_ci_low": max(0.0, null - 0.05),
                    "null_ci_high": min(1.0, null + 0.05),
                    "difference": observed - null,
                    "difference_ci_low": observed - null - 0.05,
                    "difference_ci_high": observed - null + 0.05,
                    "valid_resamples": 1_000,
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

    html = "".join(build_target_space_null_sections(analysis))

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
    assert "GAPH fraction (95% paired bootstrap interval)" in html
    assert "Matched-control fraction (95% paired bootstrap interval)" in html

    qc_html = "".join(build_target_space_null_qc_sections(analysis))
    assert "Matched-control QC" in qc_html
    assert "Strategy summary" in qc_html
    assert "VEP release" in qc_html
    assert "Paired difference" in qc_html
    assert "Difference Q2.5" in qc_html


def test_target_space_null_section_reports_disabled_state() -> None:
    html = "".join(
        build_target_space_null_sections(
            None,
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
                "target_context": "all",
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
    assert 'data-role="target-context"' in html
    assert 'data-role="consequence"' in html
    assert "phyloP fixed bands" in html
    assert "phyloP continuous" in html
    assert "Missense" in html
    assert "All target contexts" in html
    assert "Infinity" not in html
    assert "NaN" not in html


def test_dataframe_records_replaces_nonfinite_values() -> None:
    records = dataframe_records(pd.DataFrame({"value": [1.0, float("inf"), float("nan")]}))
    assert records == [{"value": 1.0}, {"value": "inf"}, {"value": None}]


def test_hidden_clinvar_association_views_reports_suppressed_selector_counts() -> None:
    frame = pd.DataFrame(
        [
            {
                "strategy": "s1",
                "variant_type": "snv",
                "target_context": "all",
                "consequence": "all",
                "usable_rows": 100,
                "status": "estimated",
                "reason": "",
            },
            {
                "strategy": "s1",
                "variant_type": "snv",
                "target_context": "all",
                "consequence": "splice_region",
                "usable_rows": 0,
                "status": "not_estimable",
                "reason": "No scored alleles in this band.",
            },
        ]
    )
    validation = SimpleNamespace(
        unadjusted=frame,
        fixed_adjusted=frame,
        continuous=frame,
    )

    summary, hidden = hidden_clinvar_association_views(validation)

    assert summary["Displayed selector combinations"].tolist() == [1, 1, 1]
    assert summary["Hidden selector combinations"].tolist() == [1, 1, 1]
    assert hidden["Consequence subset"].tolist() == ["Splice region"] * 3
    assert hidden["N across strategies"].tolist() == ["0-0"] * 3
