from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pandas as pd

from analytics.analyses.basic_filtering import (
    build_or_load_filter_score_store,
    candidate_curves_from_histograms,
    compute_clinvar_filter_curves,
    read_filter_score_histograms,
)
from analytics.analyses.conservation_validation import ConservationCohort


def _write_gzip(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_filter_score_store_and_candidate_curves_are_allele_level(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    annotation_fields = [
        "variant_key",
        "gene_id",
        "event_type",
        "lookup_status",
        "strategies",
        "gnomad_af",
        "vep_status",
        "vep_primary_consequence",
    ]
    _write_gzip(
        annotations,
        [
            {
                "variant_key": "1:10:A>G",
                "gene_id": "1",
                "event_type": "snv",
                "lookup_status": "ok",
                "strategies": "s1,s2",
                "gnomad_af": "0.01",
                "vep_status": "ok",
                "vep_primary_consequence": "missense_variant",
            },
            {
                "variant_key": "1:20:C>T",
                "gene_id": "1",
                "event_type": "snv",
                "lookup_status": "ok",
                "strategies": "s1",
                "gnomad_af": "",
                "vep_status": "ok",
                "vep_primary_consequence": "missense_variant",
            },
        ],
        annotation_fields,
    )
    support = tmp_path / "variant_strategy_support.tsv.gz"
    support_fields = [
        "variant_key",
        "gene_id",
        "strategy",
        "alt_support_ortholog_count",
        "alt_support_family_count",
        "site_aligned_ortholog_count",
    ]
    _write_gzip(
        support,
        [
            {
                "variant_key": "1:10:A>G",
                "gene_id": "1",
                "strategy": "s1",
                "alt_support_ortholog_count": 3,
                "alt_support_family_count": 2,
                "site_aligned_ortholog_count": 10,
            },
            {
                "variant_key": "1:10:A>G",
                "gene_id": "1",
                "strategy": "s2",
                "alt_support_ortholog_count": 2,
                "alt_support_family_count": 1,
                "site_aligned_ortholog_count": 10,
            },
            {
                "variant_key": "1:10:A>G",
                "gene_id": "2",
                "strategy": "s1",
                "alt_support_ortholog_count": 2,
                "alt_support_family_count": 1,
                "site_aligned_ortholog_count": 20,
            },
            {
                "variant_key": "1:20:C>T",
                "gene_id": "1",
                "strategy": "s1",
                "alt_support_ortholog_count": 1,
                "alt_support_family_count": 1,
                "site_aligned_ortholog_count": 10,
            },
        ],
        support_fields,
    )
    failures = tmp_path / "failures.tsv.gz"
    _write_gzip(
        failures,
        [],
        ["source", "scope", "chrom", "start", "end", "failure_type", "message"],
    )

    score_path, _manifest, cache_hit = build_or_load_filter_score_store(
        variant_annotations_source=annotations,
        variant_strategy_support_tsv=support,
        annotation_failures_tsv=failures,
        analytics_dir=tmp_path / "analytics",
        strategies=["s1", "s2"],
    )
    curves = candidate_curves_from_histograms(read_filter_score_histograms(score_path))

    assert cache_hit is False
    strategy_curve = curves[
        curves["strategy"].eq("s1")
        & curves["variant_type"].eq("snv")
        & curves["filter_key"].eq("strategy")
    ].set_index("threshold")
    assert strategy_curve.loc[1, "retained_variant_count"] == 2
    assert strategy_curve.loc[2, "retained_variant_count"] == 1
    assert strategy_curve.loc[2, "gnomad_found_fraction"] == 1.0
    union = curves[curves["strategy"].eq("union") & curves["filter_key"].eq("ortholog")]
    assert union.set_index("threshold").loc[1, "total_variant_count"] == 2
    assert union.set_index("threshold").loc[3, "retained_variant_count"] == 1
    assert union.set_index("threshold").loc[4, "retained_variant_count"] == 0
    minimum = curves[
        curves["strategy"].eq("s1") & curves["filter_key"].eq("aligned_min")
    ].set_index("threshold")
    maximum = curves[
        curves["strategy"].eq("s1") & curves["filter_key"].eq("aligned_max")
    ].set_index("threshold")
    assert minimum.loc[11, "retained_variant_count"] == 1
    assert maximum.loc[10, "retained_variant_count"] == 2


def test_clinvar_filter_curve_uses_support_threshold_as_observation(tmp_path: Path) -> None:
    cohort_frame = pd.DataFrame(
        [
            {
                "variant_key": "b-high",
                "variant_subtype": "snv",
                "target_context": "cds",
                "consequence_mask": 0,
                "label_class": "benign",
                "gene_ids": "1",
                "phyloP100way": 2.0,
            },
            {
                "variant_key": "b-low",
                "variant_subtype": "snv",
                "target_context": "cds",
                "consequence_mask": 0,
                "label_class": "benign",
                "gene_ids": "1",
                "phyloP100way": 0.0,
            },
            {
                "variant_key": "p-low-1",
                "variant_subtype": "snv",
                "target_context": "cds",
                "consequence_mask": 0,
                "label_class": "pathogenic",
                "gene_ids": "1",
                "phyloP100way": -2.0,
            },
            {
                "variant_key": "p-low-2",
                "variant_subtype": "snv",
                "target_context": "cds",
                "consequence_mask": 0,
                "label_class": "pathogenic",
                "gene_ids": "1",
                "phyloP100way": 0.5,
            },
        ]
    )
    scores = pd.DataFrame(
        [
            {
                "variant_key": "b-high",
                "strategy": "s1",
                "ortholog_support": 3,
                "strategy_support": 1,
                "family_support": 2,
                "site_aligned_min": 10,
                "site_aligned_max": 10,
            },
            {
                "variant_key": "b-low",
                "strategy": "s1",
                "ortholog_support": 1,
                "strategy_support": 1,
                "family_support": 1,
                "site_aligned_min": 10,
                "site_aligned_max": 10,
            },
            {
                "variant_key": "p-low-1",
                "strategy": "s1",
                "ortholog_support": 1,
                "strategy_support": 1,
                "family_support": 1,
                "site_aligned_min": 10,
                "site_aligned_max": 10,
            },
            {
                "variant_key": "p-low-2",
                "strategy": "s1",
                "ortholog_support": 1,
                "strategy_support": 1,
                "family_support": 1,
                "site_aligned_min": 10,
                "site_aligned_max": 10,
            },
        ]
    )
    curves = compute_clinvar_filter_curves(
        cohort=ConservationCohort(cohort_frame, {}),
        clinvar_scores=scores,
        strategies=["s1"],
        eligible_gene_ids_by_strategy={"s1": {"1"}},
    )
    row = curves[
        curves["mode"].eq("unadjusted")
        & curves["filter_key"].eq("ortholog")
        & curves["variant_type"].eq("snv")
        & curves["target_context"].eq("all")
        & curves["consequence"].eq("all")
        & curves["threshold"].eq(2)
    ].iloc[0]

    assert row["benign_observed"] == 1
    assert row["pathogenic_observed"] == 0
    assert row["status"] == "estimated"
    assert row["result_or"] == float("inf")


def test_local_thresholds_union_and_upper_bound_do_not_create_calls(tmp_path: Path) -> None:
    from analytics.vep.consequences import VALIDATION_CONSEQUENCE_BITS

    frame = pd.DataFrame(
        [
            {
                "variant_key": key,
                "variant_subtype": "snv",
                "target_context": context,
                "consequence_mask": VALIDATION_CONSEQUENCE_BITS["intronic"],
                "label_class": label,
                "gene_ids": "1",
                "phyloP100way": phylop,
            }
            for key, context, label, phylop in [
                ("b1", "intron", "benign", -2),
                ("b2", "intron", "benign", 2),
                ("p1", "intron", "pathogenic", -2),
                ("p2", "intron", "pathogenic", 2),
                ("missing", "intron", "benign", 0),
                ("cds", "cds", "benign", 2),
            ]
        ]
    )
    scores = pd.DataFrame(
        [
            {
                "variant_key": key,
                "strategy": strategy,
                "ortholog_support": score,
                "strategy_support": 2 if key == "b1" else 1,
                "family_support": 0,
                "site_aligned_min": depth,
                "site_aligned_max": depth,
            }
            for key, strategy, score, depth in [
                ("b1", "s1", 1, 2),
                ("b1", "s2", 3, 8),
                ("b2", "s1", 2, 4),
                ("p1", "s1", 1, 2),
                ("p2", "s1", 3, 8),
                ("cds", "s1", 500, 500),
            ]
        ]
    )
    curves = compute_clinvar_filter_curves(
        cohort=ConservationCohort(frame, {}),
        clinvar_scores=scores,
        strategies=["s1", "s2"],
        eligible_gene_ids_by_strategy={"s1": {"1"}, "s2": {"1"}},
    )
    local = curves[
        curves["variant_type"].eq("snv")
        & curves["strategy"].eq("s1")
        & curves["target_context"].eq("intron")
        & curves["consequence"].eq("intronic")
        & curves["mode"].eq("fixed")
        & curves["filter_key"].eq("ortholog")
    ]
    assert local["threshold"].tolist() == [1, 2, 3, 4]
    assert local.iloc[-1]["status"] == "not_estimable"
    upper = curves[
        curves["strategy"].eq("union")
        & curves["target_context"].eq("intron")
        & curves["consequence"].eq("intronic")
        & curves["mode"].eq("unadjusted")
        & curves["filter_key"].eq("aligned_max")
    ].set_index("threshold")
    assert upper.loc[2, "benign_observed"] == 1
    assert upper.loc[8, "benign_observed"] == 2
    assert upper.loc[8, "benign_not_observed"] == 1
    assert upper.loc[0, "benign_observed"] == 0
    family = curves[
        curves["strategy"].eq("union")
        & curves["target_context"].eq("intron")
        & curves["consequence"].eq("intronic")
        & curves["mode"].eq("unadjusted")
        & curves["filter_key"].eq("family")
    ].set_index("threshold")
    assert family.loc[0, "benign_observed"] == 2
    assert family.loc[1, "benign_observed"] == 0

    from analytics.analyses.basic_filtering import BasicFilteringAnalysis
    from analytics.reporting.basic_filtering import build_basic_filtering_sections
    from analytics.reporting.document import render_html

    histograms = pd.DataFrame(
        [
            {
                "strategy": strategy,
                "variant_type": "snv",
                "filter_key": key,
                "score": score,
                "gnomad_status": status,
                "variant_count": 1,
            }
            for strategy in ("s1", "s2", "union")
            for key in ("ortholog", "family", "strategy", "aligned_min", "aligned_max")
            for score, status in ((1, "found"), (2, "not_found"), (3, "found"))
        ]
    )
    analysis = BasicFilteringAnalysis(
        tmp_path / "scores.parquet",
        tmp_path / "manifest.json",
        candidate_curves_from_histograms(histograms),
        curves,
    )
    (tmp_path / "filtering_smoke.html").write_text(
        render_html(
            [
                ("basic-filtering", "Basic Filtering", build_basic_filtering_sections(analysis)),
            ]
        )
    )
