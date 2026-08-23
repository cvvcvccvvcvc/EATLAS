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
        "alt_support_genus_count",
    ]
    _write_gzip(
        support,
        [
            {
                "variant_key": "1:10:A>G",
                "gene_id": "1",
                "strategy": "s1",
                "alt_support_ortholog_count": 3,
                "alt_support_genus_count": 2,
            },
            {
                "variant_key": "1:10:A>G",
                "gene_id": "1",
                "strategy": "s2",
                "alt_support_ortholog_count": 2,
                "alt_support_genus_count": 1,
            },
            {
                "variant_key": "1:20:C>T",
                "gene_id": "1",
                "strategy": "s1",
                "alt_support_ortholog_count": 1,
                "alt_support_genus_count": 1,
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
                "genus_support": 2,
            },
            {
                "variant_key": "b-low",
                "strategy": "s1",
                "ortholog_support": 1,
                "strategy_support": 1,
                "genus_support": 1,
            },
            {
                "variant_key": "p-low-1",
                "strategy": "s1",
                "ortholog_support": 1,
                "strategy_support": 1,
                "genus_support": 1,
            },
            {
                "variant_key": "p-low-2",
                "strategy": "s1",
                "ortholog_support": 1,
                "strategy_support": 1,
                "genus_support": 1,
            },
        ]
    )
    candidate_curves = pd.DataFrame(
        [
            {
                "strategy": "s1",
                "variant_type": "snv",
                "filter_key": "ortholog",
                "threshold": threshold,
            }
            for threshold in (1, 2, 3)
        ]
    )

    curves = compute_clinvar_filter_curves(
        cohort=ConservationCohort(cohort_frame, {}),
        clinvar_scores=scores,
        candidate_curves=candidate_curves,
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
