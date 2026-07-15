"""Conservation-stratified ClinVar enrichment calculations."""

from __future__ import annotations

import pandas as pd

from .stats import benjamini_hochberg, enrichment_result, mantel_haenszel_adjusted


def compute_conservation_stratified_enrichment(
    *,
    universe: pd.DataFrame,
    conservation: pd.DataFrame,
    observed_by_strategy_type: dict[tuple[str, str], set[str]],
    strategies: list[str],
    score_columns: list[str],
    bins: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if bins < 2:
        raise ValueError("conservation bins must be >= 2")
    base = universe[universe["variant_type"].astype(str) == "snv"].copy()
    base = base[["variant_key", "label_class"]].drop_duplicates("variant_key")
    if base.empty:
        return pd.DataFrame(columns=bin_result_columns()), pd.DataFrame(columns=adjusted_result_columns())

    scores = conservation[["variant_key", *score_columns]].copy()
    merged = base.merge(scores, on="variant_key", how="left")
    bin_rows = []
    adjusted_rows = []
    for score_column in score_columns:
        score_values = pd.to_numeric(merged[score_column], errors="coerce")
        working = merged[score_values.notna()].copy()
        working[score_column] = score_values[score_values.notna()]
        if working.empty or working[score_column].nunique() < 2:
            continue
        q = max(2, min(int(bins), int(working[score_column].nunique())))
        try:
            working["conservation_bin"] = pd.qcut(working[score_column], q=q, duplicates="drop", precision=6)
        except ValueError:
            continue
        categories = list(working["conservation_bin"].cat.categories)
        if not categories:
            continue

        for strategy in strategies:
            observed_keys = observed_by_strategy_type.get((strategy, "snv"), set())
            working["ALT_observed"] = working["variant_key"].astype(str).isin(observed_keys)
            results = []
            for index, interval in enumerate(categories, start=1):
                subset = working[working["conservation_bin"] == interval]
                result = enrichment_for_subset(subset, f"Q{index}")
                results.append(result)
                bin_rows.append(
                    {
                        "score": score_column,
                        "strategy": strategy,
                        "bin_index": index,
                        "bin_label": f"Q{index}",
                        "bin_low": float(interval.left),
                        "bin_high": float(interval.right),
                        "row_count": result.benign_observed
                        + result.pathogenic_observed
                        + result.benign_not_observed
                        + result.pathogenic_not_observed,
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
            adjusted = mantel_haenszel_adjusted(results)
            if adjusted is None:
                continue
            adjusted_rows.append(
                {
                    "score": score_column,
                    "strategy": strategy,
                    "usable_rows": int(len(working)),
                    "bin_count": len(categories),
                    "odds_ratio_mh": adjusted.odds_ratio_mh,
                    "ci_low": adjusted.ci_low,
                    "ci_high": adjusted.ci_high,
                    "cmh_chi2": adjusted.cmh_chi2,
                    "cmh_p": adjusted.cmh_p,
                }
            )

    bin_results = pd.DataFrame(bin_rows)
    adjusted_results = pd.DataFrame(adjusted_rows)
    if not bin_results.empty:
        bin_results["fisher_q"] = benjamini_hochberg(bin_results["fisher_p"].tolist())
    if not adjusted_results.empty:
        adjusted_results["cmh_q"] = benjamini_hochberg(adjusted_results["cmh_p"].tolist())
    return (
        bin_results.reindex(columns=bin_result_columns()),
        adjusted_results.reindex(columns=adjusted_result_columns()),
    )


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


def bin_result_columns() -> list[str]:
    return [
        "score",
        "strategy",
        "bin_index",
        "bin_label",
        "bin_low",
        "bin_high",
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
    ]
