from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from analytics.core import conservation as conservation_module
from analytics.core import conservation_validation as validation_module
from analytics.core.clinvar_validation import parse_molecular_consequences
from analytics.core.conservation import DEFAULT_TRACK_NAMES, PositionScores, Track, annotate_track, score_positions
from analytics.core.stats import benjamini_hochberg
from analytics.core.conservation_validation import (
    CONSEQUENCE_OPTIONS,
    SCORE_COLUMN,
    TARGET_CONTEXT_OPTIONS,
    VARIANT_TYPE_OPTIONS,
    assign_phylop_band,
    build_conservation_cohort,
    compute_continuous_firth,
    compute_fixed_band_enrichment,
    compute_unadjusted_enrichment,
    consequence_membership_mask,
    consequence_memberships,
)


def test_phyloP_is_the_only_default_conservation_track() -> None:
    assert DEFAULT_TRACK_NAMES == "phyloP100way"


def test_allele_score_positions_exclude_indel_padding() -> None:
    assert score_positions(100, "A", "G") == ([99], "substituted_base")
    assert score_positions(100, "ATC", "A") == ([100, 101], "deleted_reference_bases_mean")
    assert score_positions(100, "A", "AGG") == ([99, 100], "insertion_flanks_mean")


def test_clinvar_mc_parser_preserves_all_so_terms() -> None:
    assert parse_molecular_consequences(
        "SO:0001583|missense_variant,SO:0001630|splice_region_variant"
    ) == [
        ("SO:0001583", "missense_variant"),
        ("SO:0001630", "splice_region_variant"),
    ]


def test_track_annotation_averages_all_required_indel_bases(monkeypatch: pytest.MonkeyPatch) -> None:
    class FakeBigWig:
        def chroms(self):
            return {"chr1": 1_000}

        def values(self, _chrom: str, start: int, end: int):
            return [float(position) for position in range(start, end)]

        def close(self):
            return None

    class FakePyBigWig:
        @staticmethod
        def open(_url: str):
            return FakeBigWig()

    monkeypatch.setattr(conservation_module, "pyBigWig", FakePyBigWig())
    rows = [
        {"variant_key": "1:100:A>G", "chrom": "1", "pos": "100", "ref": "A", "alt": "G"},
        {"variant_key": "1:100:ATC>A", "chrom": "1", "pos": "100", "ref": "ATC", "alt": "A"},
        {"variant_key": "1:100:A>AG", "chrom": "1", "pos": "100", "ref": "A", "alt": "AG"},
    ]

    summary = annotate_track(
        rows=rows,
        track=Track("test", "memory://test", "ucsc"),
        max_block_bp=1_000,
        max_gap_bp=1_000,
        remote_retries=1,
        retry_sleep_seconds=0,
        precision=6,
    )

    assert [float(row["test"]) for row in rows] == [99.0, 100.5, 99.5]
    assert summary["annotated_variants"] == 3


def test_track_annotation_reuses_precomputed_position_scores() -> None:
    track = Track("test", "memory://test", "ucsc")
    rows = [{"variant_key": "1:100:A>G", "chrom": "1", "pos": "100", "ref": "A", "alt": "G"}]
    scores = PositionScores(
        track,
        {("chr1", 99): 2.5},
        {"status": "complete", "failed_block_count": 0},
    )

    summary = annotate_track(
        rows=rows,
        track=track,
        max_block_bp=1_000,
        max_gap_bp=1_000,
        remote_retries=1,
        retry_sleep_seconds=0,
        precision=6,
        position_scores=scores,
    )

    assert rows[0]["test"] == "2.5"
    assert summary["source"] == "shared_position_read"
    assert summary["annotated_variants"] == 1


def test_consequence_membership_is_nonexclusive() -> None:
    groups = consequence_memberships("missense_variant|splice_region_variant")
    assert groups == {"missense", "splice_region"}
    assert consequence_memberships("") == {"other"}
    assert consequence_memberships("missense_variant|unrecognized_term") == {"missense", "other"}


def test_clinvar_association_selector_contract_keeps_all_consequences() -> None:
    assert [key for key, _label in VARIANT_TYPE_OPTIONS] == ["snv", "indel"]
    assert [key for key, _label in TARGET_CONTEXT_OPTIONS] == ["all", "cds", "utr", "intron"]
    assert [key for key, _label in CONSEQUENCE_OPTIONS] == [
        "all",
        "missense",
        "synonymous",
        "protein_truncating",
        "canonical_splice",
        "inframe_protein_altering",
        "splice_region",
        "intronic",
        "utr_noncoding",
        "other",
    ]


def test_cohort_keeps_snv_and_indel_subtypes_multiple_consequences_and_target_context(
    tmp_path: Path,
) -> None:
    universe = pd.DataFrame(
        [
            universe_row("1:10:A>G", "snv", "A", "G", "missense_variant|splice_region_variant"),
            universe_row("1:20:A>AT", "indel", "A", "AT", "inframe_insertion"),
            universe_row("1:30:AT>A", "indel", "AT", "A", "frameshift_variant"),
        ]
    )
    conservation = pd.DataFrame(
        {"variant_key": universe["variant_key"], SCORE_COLUMN: [-1.0, 0.2, 2.0]}
    )
    genes_path = tmp_path / "genes.tsv.gz"
    pd.DataFrame([{"gene_id": "1", "begin": 1, "sequence_length": 100}]).to_csv(
        genes_path,
        sep="\t",
        index=False,
        compression="gzip",
    )
    features_path = tmp_path / "target_features.tsv.gz"
    pd.DataFrame(
        [
            {"gene_id": "1", "feature_type": "cds", "target_start0": 0, "target_end0": 15},
            {"gene_id": "1", "feature_type": "intron", "target_start0": 15, "target_end0": 100},
        ]
    ).to_csv(features_path, sep="\t", index=False, compression="gzip")

    cohort = build_conservation_cohort(
        universe=universe,
        conservation=conservation,
        genes_tsv=genes_path,
        target_features_tsv=features_path,
    )

    assert cohort.variants["variant_subtype"].tolist() == ["snv", "insertion", "deletion"]
    assert cohort.variants["target_context"].tolist() == ["cds", "intron", "intron"]
    assert cohort.variants.loc[0, "consequence_groups"] == "missense|splice_region"
    assert cohort.summary["multiple_consequence_group_count"] == 1


def test_fixed_bands_use_prespecified_boundaries_and_all_selectors(tmp_path: Path) -> None:
    values = pd.Series([-1.30103, -1.0, 1.301029, 1.30103])
    assert assign_phylop_band(values).astype(str).tolist() == [
        "acceleration",
        "central",
        "central",
        "conservation",
    ]

    cohort = synthetic_cohort(120)
    observed = set(cohort.loc[cohort.index % 3 != 0, "variant_key"])
    bins, adjusted = compute_fixed_band_enrichment(
        cohort=cohort,
        observed_by_strategy_type={("s1", "snv"): observed, ("s1", "indel"): set()},
        strategies=["s1"],
    )

    selected_bins = bins[
        (bins["strategy"] == "s1")
        & (bins["variant_type"] == "snv")
        & (bins["target_context"] == "all")
        & (bins["consequence"] == "missense")
    ]
    selected_adjusted = adjusted[
        (adjusted["strategy"] == "s1")
        & (adjusted["variant_type"] == "snv")
        & (adjusted["target_context"] == "all")
        & (adjusted["consequence"] == "missense")
    ].iloc[0]
    assert len(selected_bins) == 3
    assert selected_adjusted["status"] == "estimated"
    assert np.isfinite(selected_adjusted["cmh_p"])
    assert bins["fisher_q"].notna().any()
    empty_indel_bins = bins[
        (bins["strategy"] == "s1")
        & (bins["variant_type"] == "indel")
        & (bins["target_context"] == "all")
        & (bins["consequence"] == "all")
    ]
    assert set(empty_indel_bins["status"]) == {"not_estimable"}
    assert empty_indel_bins["fisher_p"].isna().all()


def test_unadjusted_enrichment_supports_shared_selectors_and_strategy_level_fdr() -> None:
    cohort = synthetic_cohort(120)
    observed_s1 = set(cohort.loc[cohort.index % 3 != 0, "variant_key"])
    observed_s2 = set(cohort.loc[cohort.index % 4 != 0, "variant_key"])
    results = compute_unadjusted_enrichment(
        cohort=cohort,
        observed_by_strategy_type={
            ("s1", "snv"): observed_s1,
            ("s1", "indel"): set(),
            ("s2", "snv"): observed_s2,
            ("s2", "indel"): set(),
        },
        strategies=["s1", "s2"],
    )

    for target_context in ["all", "cds"]:
        selected = results[
            results["variant_type"].eq("snv")
            & results["target_context"].eq(target_context)
            & results["consequence"].eq("missense")
        ].sort_values("strategy")
        assert len(selected) == 2
        assert selected["usable_rows"].tolist() == [120, 120]
        assert selected["fisher_q"].tolist() == pytest.approx(
            benjamini_hochberg(selected["fisher_p"].tolist())
        )


def test_unadjusted_enrichment_uses_strategy_specific_gene_denominator() -> None:
    cohort = synthetic_cohort(120)
    cohort["gene_ids"] = np.where(cohort.index < 60, "1", "2")
    observed = set(cohort.loc[cohort.index % 3 == 0, "variant_key"])

    results = compute_unadjusted_enrichment(
        cohort=cohort,
        observed_by_strategy_type={
            ("s1", "snv"): observed,
            ("s1", "indel"): set(),
        },
        strategies=["s1"],
        eligible_gene_ids_by_strategy={"s1": {"1"}},
    )

    selected = results[
        results["variant_type"].eq("snv")
        & results["target_context"].eq("all")
        & results["consequence"].eq("missense")
    ].iloc[0]
    assert selected["usable_rows"] == 60


def test_continuous_precheck_rejects_nonoverlapping_score_ranges(tmp_path: Path) -> None:
    cohort = synthetic_cohort(80)
    observed_keys = set(cohort.loc[cohort[SCORE_COLUMN] > 0, "variant_key"])

    results, distributions, versions = compute_continuous_firth(
        cohort=cohort,
        observed_by_strategy_type={("s1", "snv"): observed_keys, ("s1", "indel"): set()},
        strategies=["s1"],
        analytics_dir=tmp_path,
        rscript="/path/not/used",
    )

    selected = results[
        (results["variant_type"] == "snv")
        & (results["target_context"] == "all")
        & (results["consequence"] == "missense")
    ].iloc[0]
    assert selected["status"] == "not_estimable"
    assert "do not overlap" in selected["reason"]
    assert not distributions.empty
    assert {
        "bin_left",
        "bin_right",
        "fraction",
        "q1",
        "median",
        "q3",
        "lower_whisker",
        "upper_whisker",
    }.issubset(distributions.columns)
    assert versions == {}


def test_continuous_model_exports_strategy_eligibility_mask(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cohort = synthetic_cohort(120)
    cohort["gene_ids"] = np.where(cohort.index < 60, "1", "2")
    observed_keys = set(cohort.loc[cohort.index % 3 != 0, "variant_key"])
    captured = {}

    def fake_run_firth_models(*, model_data, specs, **_kwargs):
        captured["model_data"] = model_data
        captured["specs"] = specs
        fitted = specs[["analysis_id"]].copy()
        fitted["odds_ratio"] = 1.0
        fitted["ci_low"] = 0.5
        fitted["ci_high"] = 2.0
        fitted["plr_p"] = 1.0
        fitted["status"] = "estimated"
        fitted["reason"] = ""
        return fitted, {"R": "test", "logistf": "test"}

    monkeypatch.setattr(validation_module, "run_firth_models", fake_run_firth_models)
    results, _distributions, _versions = compute_continuous_firth(
        cohort=cohort,
        observed_by_strategy_type={
            ("s1", "snv"): observed_keys,
            ("s1", "indel"): set(),
        },
        strategies=["s1"],
        analytics_dir=tmp_path,
        eligible_gene_ids_by_strategy={"s1": {"1"}},
    )

    assert captured["model_data"]["eligible_0"].sum() == 60
    assert set(captured["specs"]["eligibility_column"]) == {"eligible_0"}
    assert set(captured["specs"]["target_context"]) == {"all", "cds"}
    selected = results[
        results["variant_type"].eq("snv")
        & results["target_context"].eq("all")
        & results["consequence"].eq("missense")
    ].iloc[0]
    assert selected["usable_rows"] == 60


def test_firth_model_returns_profile_likelihood_result_when_r_is_available(tmp_path: Path) -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not installed")
    package_check = subprocess.run(
        [rscript, "--vanilla", "-e", "quit(status=ifelse(requireNamespace('logistf', quietly=TRUE),0,1))"],
        check=False,
    )
    if package_check.returncode != 0:
        pytest.skip("R package logistf is not installed")

    rng = np.random.default_rng(7)
    count = 250
    score = rng.normal(size=count)
    observed = rng.binomial(1, 0.5, size=count)
    probability = 1 / (1 + np.exp(-(-0.4 + 0.9 * observed - 0.35 * score + 0.1 * score**2)))
    benign = rng.binomial(1, probability)
    cohort = synthetic_cohort(count)
    cohort[SCORE_COLUMN] = score
    cohort["label_class"] = np.where(benign == 1, "benign", "pathogenic")
    observed_keys = set(cohort.loc[observed == 1, "variant_key"])

    results, _distributions, versions = compute_continuous_firth(
        cohort=cohort,
        observed_by_strategy_type={("s1", "snv"): observed_keys, ("s1", "indel"): set()},
        strategies=["s1"],
        analytics_dir=tmp_path,
        rscript=rscript,
    )

    selected = results[
        (results["variant_type"] == "snv")
        & (results["target_context"] == "all")
        & (results["consequence"] == "missense")
    ].iloc[0]
    assert selected["status"] == "estimated"
    assert selected["odds_ratio"] > 1
    assert selected["ci_low"] > 1
    assert np.isfinite(selected["plr_p"])
    assert set(versions) == {"R", "logistf"}


def universe_row(key: str, variant_type: str, ref: str, alt: str, terms: str) -> dict[str, str]:
    return {
        "variant_key": key,
        "variant_type": variant_type,
        "ref": ref,
        "alt": alt,
        "label_class": "benign",
        "clinvar_mc_terms": terms,
        "gene_ids": "1",
    }


def synthetic_cohort(count: int) -> pd.DataFrame:
    index = np.arange(count)
    return pd.DataFrame(
        {
            "variant_key": [f"1:{value}:A>G" for value in index],
            "variant_type": "snv",
            "variant_subtype": "snv",
            "label_class": np.where(index % 4 == 0, "pathogenic", "benign"),
            "consequence_groups": "missense",
            "consequence_mask": consequence_membership_mask("missense_variant"),
            "target_context": "cds",
            SCORE_COLUMN: np.linspace(-3, 3, count),
        }
    )
