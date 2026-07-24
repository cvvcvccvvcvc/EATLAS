"""Intronic ClinVar validation adjusted for site-level conservation."""

from __future__ import annotations

import math
import warnings
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from patsy import dmatrix
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning

from .stats import benjamini_hochberg, enrichment_result, mantel_haenszel_adjusted


SPLICE_PROXIMAL_BP = 8
PRIMARY_SCOPE = "all_intronic"
SENSITIVITY_SCOPE = "excluding_splice_proximal"
SCOPE_LABELS = {
    PRIMARY_SCOPE: "All intronic SNVs",
    SENSITIVITY_SCOPE: f"Excluding first/last {SPLICE_PROXIMAL_BP} intronic bp",
}


@dataclass(frozen=True)
class ConservationCategory:
    label: str
    range_text: str


CATEGORY_DEFINITIONS = {
    "phyloP100way": [
        ConservationCategory("Nominal acceleration band", "<= -1.30103"),
        ConservationCategory("Central phyloP band", "-1.30103 to 1.30103"),
        ConservationCategory("Nominal conservation band", ">= 1.30103"),
    ],
    "phastCons100way": [
        ConservationCategory("Lower conservation probability", "< 0.5"),
        ConservationCategory("Higher conservation probability", ">= 0.5"),
    ],
    "GERP_RS_92mammals": [
        ConservationCategory("No constraint / substitution surplus", "<= 0"),
        ConservationCategory("Weak constraint", "0 to < 2"),
        ConservationCategory("Moderate constraint", "2 to < 4"),
        ConservationCategory("Strong constraint", ">= 4"),
    ],
}


@dataclass(frozen=True)
class IntervalLookup:
    starts: tuple[int, ...]
    ends: tuple[int, ...]

    def containing(self, position: int) -> tuple[int, int] | None:
        index = bisect_right(self.starts, position) - 1
        if index >= 0 and position <= self.ends[index]:
            return self.starts[index], self.ends[index]
        return None


@dataclass(frozen=True)
class IntronicCohort:
    variants: pd.DataFrame
    summary: dict[str, int]


def build_intronic_cohort(
    *,
    universe: pd.DataFrame,
    conservation: pd.DataFrame,
    target_features_tsv: Path,
    score_columns: list[str],
) -> IntronicCohort:
    """Return usable ClinVar SNVs located in unambiguous target introns."""
    required = {"variant_key", "variant_type", "label_class", "gene_ids", "pos"}
    missing = required - set(universe.columns)
    if missing:
        raise ValueError(f"ClinVar universe missing required columns: {', '.join(sorted(missing))}")

    features = pd.read_csv(
        target_features_tsv,
        sep="\t",
        compression="infer",
        dtype={"gene_id": str, "feature_type": str},
        keep_default_na=False,
        low_memory=False,
    )
    feature_required = {"gene_id", "feature_type", "genomic_start1", "genomic_end1"}
    feature_missing = feature_required - set(features.columns)
    if feature_missing:
        raise ValueError(
            f"Target feature table missing required columns: {', '.join(sorted(feature_missing))}"
        )
    interval_indexes = build_feature_interval_indexes(features)

    base = universe[universe["variant_type"].astype(str) == "snv"].copy()
    base = base[["variant_key", "label_class", "gene_ids", "pos"]].drop_duplicates("variant_key")
    score_frame = conservation[["variant_key", *score_columns]].drop_duplicates("variant_key")
    base = base.merge(score_frame, on="variant_key", how="left", validate="one_to_one")

    region_classes: list[str] = []
    boundary_distances: list[float] = []
    for row in base.itertuples(index=False):
        region_class, boundary_distance = classify_target_position(
            str(row.gene_ids), int(row.pos), interval_indexes
        )
        region_classes.append(region_class)
        boundary_distances.append(boundary_distance)
    base["target_region_class"] = region_classes
    base["intron_boundary_distance"] = boundary_distances
    base["splice_proximal"] = (
        base["intron_boundary_distance"].notna()
        & (base["intron_boundary_distance"] <= SPLICE_PROXIMAL_BP)
    )

    intronic = base[base["target_region_class"] == "intronic"].copy()
    summary = {
        "usable_snv_count": int(len(base)),
        "intronic_snv_count": int(len(intronic)),
        "intronic_benign_count": int((intronic["label_class"].astype(str) == "benign").sum()),
        "intronic_pathogenic_count": int((intronic["label_class"].astype(str) == "pathogenic").sum()),
        "splice_proximal_count": int(intronic["splice_proximal"].sum()),
        "non_splice_proximal_count": int((~intronic["splice_proximal"]).sum()),
        "mixed_exon_intron_count": int((base["target_region_class"] == "mixed_exon_intron").sum()),
        "exonic_count": int((base["target_region_class"] == "exonic").sum()),
        "unclassified_count": int((base["target_region_class"] == "unclassified").sum()),
    }
    return IntronicCohort(intronic, summary)


def build_feature_interval_indexes(
    features: pd.DataFrame,
) -> dict[str, dict[str, IntervalLookup]]:
    selected = features[features["feature_type"].astype(str).isin({"exon", "intron"})].copy()
    selected["genomic_start1"] = pd.to_numeric(selected["genomic_start1"], errors="raise").astype(int)
    selected["genomic_end1"] = pd.to_numeric(selected["genomic_end1"], errors="raise").astype(int)
    indexes: dict[str, dict[str, IntervalLookup]] = {}
    for (gene_id, feature_type), group in selected.groupby(["gene_id", "feature_type"], sort=False):
        intervals = sorted(zip(group["genomic_start1"], group["genomic_end1"], strict=True))
        indexes.setdefault(str(gene_id), {})[str(feature_type)] = IntervalLookup(
            tuple(int(start) for start, _end in intervals),
            tuple(int(end) for _start, end in intervals),
        )
    return indexes


def classify_target_position(
    gene_ids_text: str,
    position: int,
    interval_indexes: dict[str, dict[str, IntervalLookup]],
) -> tuple[str, float]:
    gene_ids = [value.strip() for value in gene_ids_text.split("|") if value.strip()]
    intron_intervals = []
    in_exon = False
    for gene_id in gene_ids:
        indexes = interval_indexes.get(gene_id, {})
        intron_lookup = indexes.get("intron")
        exon_lookup = indexes.get("exon")
        if intron_lookup is not None:
            interval = intron_lookup.containing(position)
            if interval is not None:
                intron_intervals.append(interval)
        if exon_lookup is not None and exon_lookup.containing(position) is not None:
            in_exon = True

    if intron_intervals and in_exon:
        return "mixed_exon_intron", float("nan")
    if intron_intervals:
        distance = min(
            min(position - start + 1, end - position + 1)
            for start, end in intron_intervals
        )
        return "intronic", float(distance)
    if in_exon:
        return "exonic", float("nan")
    return "unclassified", float("nan")


def analysis_scopes(cohort: pd.DataFrame) -> dict[str, pd.DataFrame]:
    return {
        PRIMARY_SCOPE: cohort,
        SENSITIVITY_SCOPE: cohort[~cohort["splice_proximal"].astype(bool)].copy(),
    }


def assign_conservation_category(values: pd.Series, score: str) -> pd.Categorical:
    numeric = pd.to_numeric(values, errors="coerce")
    definitions = CATEGORY_DEFINITIONS.get(score)
    if definitions is None:
        raise ValueError(f"No prespecified categorical thresholds for conservation score {score!r}")
    labels = [definition.label for definition in definitions]

    categories = pd.Series(pd.NA, index=numeric.index, dtype="object")
    if score == "phyloP100way":
        categories.loc[numeric <= -1.30103] = labels[0]
        categories.loc[(numeric > -1.30103) & (numeric < 1.30103)] = labels[1]
        categories.loc[numeric >= 1.30103] = labels[2]
    elif score == "phastCons100way":
        categories.loc[numeric < 0.5] = labels[0]
        categories.loc[numeric >= 0.5] = labels[1]
    elif score == "GERP_RS_92mammals":
        categories.loc[numeric <= 0] = labels[0]
        categories.loc[(numeric > 0) & (numeric < 2)] = labels[1]
        categories.loc[(numeric >= 2) & (numeric < 4)] = labels[2]
        categories.loc[numeric >= 4] = labels[3]
    return pd.Categorical(categories, categories=labels, ordered=True)


def compute_categorical_enrichment(
    *,
    cohort: pd.DataFrame,
    observed_by_strategy_type: dict[tuple[str, str], set[str]],
    strategies: list[str],
    score_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    bin_rows: list[dict[str, object]] = []
    adjusted_rows: list[dict[str, object]] = []
    for scope, scope_frame in analysis_scopes(cohort).items():
        for score in score_columns:
            if score not in CATEGORY_DEFINITIONS:
                continue
            working = scope_frame.copy()
            working["conservation_category"] = assign_conservation_category(working[score], score)
            working = working[working["conservation_category"].notna()].copy()
            definitions = CATEGORY_DEFINITIONS[score]
            for strategy in strategies:
                observed_keys = observed_by_strategy_type.get((strategy, "snv"), set())
                working["ALT_observed"] = working["variant_key"].astype(str).isin(observed_keys)
                results = []
                for index, definition in enumerate(definitions, start=1):
                    subset = working[working["conservation_category"] == definition.label]
                    result = enrichment_for_subset(subset, definition.label)
                    results.append(result)
                    bin_rows.append(
                        {
                            "scope": scope,
                            "score": score,
                            "strategy": strategy,
                            "bin_index": index,
                            "bin_label": definition.label,
                            "bin_range": definition.range_text,
                            "row_count": int(len(subset)),
                            "benign_observed": result.benign_observed,
                            "pathogenic_observed": result.pathogenic_observed,
                            "benign_not_observed": result.benign_not_observed,
                            "pathogenic_not_observed": result.pathogenic_not_observed,
                            "odds_ratio": result.odds_ratio,
                            "ci_low": result.ci_low,
                            "ci_high": result.ci_high,
                            "fisher_p": result.fisher_p,
                        }
                    )

                status, reason = categorical_estimability(working)
                adjusted = mantel_haenszel_adjusted(results) if status == "estimated" else None
                if adjusted is None or math.isnan(adjusted.odds_ratio_mh) or not math.isfinite(adjusted.cmh_p):
                    if not reason:
                        reason = "The category tables do not provide an estimable common odds ratio."
                    adjusted_rows.append(empty_adjusted_row(scope, score, strategy, len(working), reason))
                else:
                    adjusted_rows.append(
                        {
                            "scope": scope,
                            "score": score,
                            "strategy": strategy,
                            "usable_rows": int(len(working)),
                            "bin_count": len(definitions),
                            "odds_ratio_mh": adjusted.odds_ratio_mh,
                            "ci_low": adjusted.ci_low,
                            "ci_high": adjusted.ci_high,
                            "cmh_chi2": adjusted.cmh_chi2,
                            "cmh_p": adjusted.cmh_p,
                            "status": "estimated",
                            "reason": "",
                        }
                    )

    bin_results = pd.DataFrame(bin_rows, columns=bin_result_columns())
    adjusted_results = pd.DataFrame(adjusted_rows, columns=adjusted_result_columns())
    for scope in SCOPE_LABELS:
        bin_mask = bin_results["scope"] == scope
        bin_results.loc[bin_mask, "fisher_q"] = benjamini_hochberg(
            bin_results.loc[bin_mask, "fisher_p"].tolist()
        )
        adjusted_mask = adjusted_results["scope"] == scope
        adjusted_results.loc[adjusted_mask, "cmh_q"] = benjamini_hochberg(
            adjusted_results.loc[adjusted_mask, "cmh_p"].tolist()
        )
    return bin_results, adjusted_results


def compute_continuous_enrichment(
    *,
    cohort: pd.DataFrame,
    observed_by_strategy_type: dict[tuple[str, str], set[str]],
    strategies: list[str],
    score_columns: list[str],
    spline_df: int = 3,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for scope, scope_frame in analysis_scopes(cohort).items():
        for score in score_columns:
            for strategy in strategies:
                observed_keys = observed_by_strategy_type.get((strategy, "snv"), set())
                working = pd.DataFrame(
                    {
                        "label_class": scope_frame["label_class"].astype(str),
                        "ALT_observed": scope_frame["variant_key"].astype(str).isin(observed_keys).astype(int),
                        "score": pd.to_numeric(scope_frame[score], errors="coerce"),
                    }
                )
                working = working[np.isfinite(working["score"])].copy()
                rows.append(fit_continuous_model(working, scope, score, strategy, spline_df))

    results = pd.DataFrame(rows, columns=continuous_result_columns())
    for scope in SCOPE_LABELS:
        mask = results["scope"] == scope
        results.loc[mask, "wald_q"] = benjamini_hochberg(results.loc[mask, "wald_p"].tolist())
    return results


def fit_continuous_model(
    working: pd.DataFrame,
    scope: str,
    score: str,
    strategy: str,
    spline_df: int,
) -> dict[str, object]:
    base = {
        "scope": scope,
        "score": score,
        "strategy": strategy,
        "usable_rows": int(len(working)),
        "benign_rows": int((working["label_class"] == "benign").sum()),
        "pathogenic_rows": int((working["label_class"] == "pathogenic").sum()),
        "observed_rows": int(working["ALT_observed"].sum()),
        "not_observed_rows": int((working["ALT_observed"] == 0).sum()),
        "score_min": float(working["score"].min()) if not working.empty else float("nan"),
        "score_max": float(working["score"].max()) if not working.empty else float("nan"),
        "spline_df": spline_df,
        "odds_ratio": float("nan"),
        "ci_low": float("nan"),
        "ci_high": float("nan"),
        "wald_p": float("nan"),
        "status": "not_estimable",
        "reason": "",
    }
    reason = continuous_estimability_reason(working, spline_df)
    if reason:
        base["reason"] = reason
        return base

    model_data = pd.DataFrame(
        {
            "benign": (working["label_class"] == "benign").astype(int),
            "ALT_observed": working["ALT_observed"].astype(float),
            "score": working["score"].astype(float),
        }
    )
    try:
        design = dmatrix(
            f"ALT_observed + cr(score, df={int(spline_df)}, constraints='center')",
            model_data,
            return_type="dataframe",
        )
        if np.linalg.matrix_rank(design.to_numpy()) < design.shape[1]:
            base["reason"] = "The spline design matrix is rank deficient."
            return base
        with warnings.catch_warnings(record=True) as fit_warnings:
            warnings.simplefilter("always")
            fitted = sm.GLM(model_data["benign"], design, family=sm.families.Binomial()).fit(
                maxiter=100,
                disp=0,
            )
        if any(issubclass(item.category, PerfectSeparationWarning) for item in fit_warnings):
            base["reason"] = "The ALT effect is not estimable because the data are separated."
            return base
        if not bool(getattr(fitted, "converged", True)):
            base["reason"] = "The logistic model did not converge."
            return base
        coefficient = float(fitted.params["ALT_observed"])
        ci = fitted.conf_int(alpha=0.05).loc["ALT_observed"]
        ci_low_log, ci_high_log = float(ci.iloc[0]), float(ci.iloc[1])
        pvalue = float(fitted.pvalues["ALT_observed"])
        if not all(math.isfinite(value) for value in [coefficient, ci_low_log, ci_high_log, pvalue]):
            base["reason"] = "The fitted ALT coefficient or its uncertainty is non-finite."
            return base
    except Exception as exc:
        base["reason"] = f"Model fitting failed: {type(exc).__name__}."
        return base

    base.update(
        {
            "odds_ratio": math.exp(coefficient),
            "ci_low": math.exp(ci_low_log),
            "ci_high": math.exp(ci_high_log),
            "wald_p": pvalue,
            "status": "estimated",
            "reason": "",
        }
    )
    return base


def categorical_estimability(working: pd.DataFrame) -> tuple[str, str]:
    if working.empty:
        return "not_estimable", "No intronic ClinVar SNVs have a score."
    if working["label_class"].nunique() < 2:
        return "not_estimable", "Both B/LB and P/LP alleles are required."
    if working["ALT_observed"].nunique() < 2:
        return "not_estimable", "Both observed and not-observed alleles are required."
    if working["conservation_category"].nunique() < 2:
        return "not_estimable", "At least two conservation categories with data are required."
    return "estimated", ""


def continuous_estimability_reason(working: pd.DataFrame, spline_df: int) -> str:
    if working.empty:
        return "No intronic ClinVar SNVs have a score."
    if working["label_class"].nunique() < 2:
        return "Both B/LB and P/LP alleles are required."
    if working["ALT_observed"].nunique() < 2:
        return "Both observed and not-observed alleles are required."
    if working["score"].nunique() < spline_df + 1:
        return f"At least {spline_df + 1} distinct score values are required for the spline."
    if len(working) <= spline_df + 2:
        return "Too few rows for the spline-adjusted model."
    return ""


def enrichment_for_subset(df: pd.DataFrame, name: str):
    benign = df["label_class"].astype(str) == "benign"
    pathogenic = df["label_class"].astype(str) == "pathogenic"
    observed = df["ALT_observed"].astype(bool)
    return enrichment_result(
        name,
        int((benign & observed).sum()),
        int((pathogenic & observed).sum()),
        int((benign & ~observed).sum()),
        int((pathogenic & ~observed).sum()),
    )


def empty_adjusted_row(
    scope: str,
    score: str,
    strategy: str,
    usable_rows: int,
    reason: str,
) -> dict[str, object]:
    return {
        "scope": scope,
        "score": score,
        "strategy": strategy,
        "usable_rows": int(usable_rows),
        "bin_count": len(CATEGORY_DEFINITIONS[score]),
        "odds_ratio_mh": float("nan"),
        "ci_low": float("nan"),
        "ci_high": float("nan"),
        "cmh_chi2": float("nan"),
        "cmh_p": float("nan"),
        "status": "not_estimable",
        "reason": reason,
    }


def bin_result_columns() -> list[str]:
    return [
        "scope",
        "score",
        "strategy",
        "bin_index",
        "bin_label",
        "bin_range",
        "row_count",
        "benign_observed",
        "pathogenic_observed",
        "benign_not_observed",
        "pathogenic_not_observed",
        "odds_ratio",
        "ci_low",
        "ci_high",
        "fisher_p",
        "fisher_q",
    ]


def adjusted_result_columns() -> list[str]:
    return [
        "scope",
        "score",
        "strategy",
        "usable_rows",
        "bin_count",
        "odds_ratio_mh",
        "ci_low",
        "ci_high",
        "cmh_chi2",
        "cmh_p",
        "cmh_q",
        "status",
        "reason",
    ]


def continuous_result_columns() -> list[str]:
    return [
        "scope",
        "score",
        "strategy",
        "usable_rows",
        "benign_rows",
        "pathogenic_rows",
        "observed_rows",
        "not_observed_rows",
        "score_min",
        "score_max",
        "spline_df",
        "odds_ratio",
        "ci_low",
        "ci_high",
        "wald_p",
        "wald_q",
        "status",
        "reason",
    ]
