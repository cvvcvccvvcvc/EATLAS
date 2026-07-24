from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from analytics.core.conservation import DEFAULT_TRACK_NAMES
from analytics.core.intronic_conservation import (
    PRIMARY_SCOPE,
    assign_conservation_category,
    build_intronic_cohort,
    compute_categorical_enrichment,
    compute_continuous_enrichment,
)


def test_phyloP_is_the_default_conservation_track() -> None:
    assert DEFAULT_TRACK_NAMES == "phyloP100way"


def test_intronic_cohort_uses_gene_aware_intervals_and_splice_distance(tmp_path: Path) -> None:
    features_path = tmp_path / "target_features.tsv.gz"
    features = pd.DataFrame(
        [
            feature("1", "exon", 100, 109),
            feature("1", "intron", 110, 190),
            feature("1", "exon", 191, 200),
            feature("2", "intron", 100, 149),
            feature("2", "exon", 150, 160),
            feature("2", "intron", 161, 200),
        ]
    )
    features.to_csv(features_path, sep="\t", index=False, compression="gzip")
    universe = pd.DataFrame(
        [
            universe_row("1:110:A>G", 110, "1"),
            universe_row("1:117:A>G", 117, "1"),
            universe_row("1:118:A>G", 118, "1"),
            universe_row("1:150:A>G", 150, "1|2"),
            universe_row("1:105:A>G", 105, "1"),
            universe_row("1:250:A>G", 250, "1"),
        ]
    )
    conservation = pd.DataFrame(
        {
            "variant_key": universe["variant_key"],
            "phyloP100way": np.linspace(-2, 2, len(universe)),
        }
    )

    result = build_intronic_cohort(
        universe=universe,
        conservation=conservation,
        target_features_tsv=features_path,
        score_columns=["phyloP100way"],
    )

    cohort = result.variants.set_index("variant_key")
    assert set(cohort.index) == {"1:110:A>G", "1:117:A>G", "1:118:A>G"}
    assert cohort.loc["1:110:A>G", "intron_boundary_distance"] == 1
    assert bool(cohort.loc["1:117:A>G", "splice_proximal"])
    assert not bool(cohort.loc["1:118:A>G", "splice_proximal"])
    assert result.summary["mixed_exon_intron_count"] == 1
    assert result.summary["exonic_count"] == 1
    assert result.summary["unclassified_count"] == 1


def test_prespecified_conservation_category_boundaries() -> None:
    phylo = assign_conservation_category(
        pd.Series([-1.30103, -1.0, 1.301029, 1.30103]), "phyloP100way"
    ).astype(str).tolist()
    assert phylo == [
        "Nominal acceleration band",
        "Central phyloP band",
        "Central phyloP band",
        "Nominal conservation band",
    ]

    phastcons = assign_conservation_category(pd.Series([0.499999, 0.5]), "phastCons100way")
    assert phastcons.astype(str).tolist() == [
        "Lower conservation probability",
        "Higher conservation probability",
    ]

    gerp = assign_conservation_category(
        pd.Series([0.0, 0.0001, 1.9999, 2.0, 3.9999, 4.0]), "GERP_RS_92mammals"
    )
    assert gerp.astype(str).tolist() == [
        "No constraint / substitution surplus",
        "Weak constraint",
        "Weak constraint",
        "Moderate constraint",
        "Moderate constraint",
        "Strong constraint",
    ]


def test_categorical_analysis_reports_fixed_strata_and_adjusted_result() -> None:
    cohort = synthetic_cohort(80)
    observed = set(cohort.loc[cohort.index % 3 != 0, "variant_key"])
    bins, adjusted = compute_categorical_enrichment(
        cohort=cohort,
        observed_by_strategy_type={("s1", "snv"): observed},
        strategies=["s1"],
        score_columns=["phyloP100way", "phastCons100way", "GERP_RS_92mammals"],
    )

    primary_bins = bins[bins["scope"] == PRIMARY_SCOPE]
    assert len(primary_bins[primary_bins["score"] == "phyloP100way"]) == 3
    assert len(primary_bins[primary_bins["score"] == "phastCons100way"]) == 2
    assert len(primary_bins[primary_bins["score"] == "GERP_RS_92mammals"]) == 4
    primary_adjusted = adjusted[adjusted["scope"] == PRIMARY_SCOPE]
    assert set(primary_adjusted["status"]) == {"estimated"}
    assert primary_adjusted["cmh_p"].notna().all()


def test_continuous_analysis_estimates_alt_effect_and_rejects_one_class() -> None:
    rng = np.random.default_rng(7)
    count = 500
    score = rng.normal(size=count)
    observed = rng.binomial(1, 0.5, size=count)
    benign_probability = 1 / (1 + np.exp(-(-0.4 + 0.9 * observed - 0.35 * score + 0.1 * score**2)))
    benign = rng.binomial(1, benign_probability)
    cohort = pd.DataFrame(
        {
            "variant_key": [f"1:{index}:A>G" for index in range(count)],
            "label_class": np.where(benign == 1, "benign", "pathogenic"),
            "splice_proximal": False,
            "phyloP100way": score,
        }
    )
    observed_keys = set(cohort.loc[observed == 1, "variant_key"])

    results = compute_continuous_enrichment(
        cohort=cohort,
        observed_by_strategy_type={("s1", "snv"): observed_keys},
        strategies=["s1"],
        score_columns=["phyloP100way"],
    )
    primary = results[results["scope"] == PRIMARY_SCOPE].iloc[0]
    assert primary["status"] == "estimated"
    assert primary["odds_ratio"] > 1
    assert primary["ci_low"] > 1

    one_class = cohort.copy()
    one_class["label_class"] = "benign"
    sparse = compute_continuous_enrichment(
        cohort=one_class,
        observed_by_strategy_type={("s1", "snv"): observed_keys},
        strategies=["s1"],
        score_columns=["phyloP100way"],
    )
    assert set(sparse["status"]) == {"not_estimable"}
    assert sparse["reason"].str.contains("Both B/LB and P/LP").all()


def feature(gene_id: str, feature_type: str, start: int, end: int) -> dict[str, object]:
    return {
        "gene_id": gene_id,
        "feature_type": feature_type,
        "genomic_start1": start,
        "genomic_end1": end,
    }


def universe_row(key: str, pos: int, gene_ids: str) -> dict[str, object]:
    return {
        "variant_key": key,
        "variant_type": "snv",
        "label_class": "benign",
        "gene_ids": gene_ids,
        "pos": pos,
    }


def synthetic_cohort(count: int) -> pd.DataFrame:
    index = np.arange(count)
    return pd.DataFrame(
        {
            "variant_key": [f"1:{value}:A>G" for value in index],
            "label_class": np.where(index % 4 == 0, "pathogenic", "benign"),
            "splice_proximal": index % 10 == 0,
            "phyloP100way": np.linspace(-3, 3, count),
            "phastCons100way": np.linspace(0, 1, count),
            "GERP_RS_92mammals": np.linspace(-2, 6, count),
        }
    )
