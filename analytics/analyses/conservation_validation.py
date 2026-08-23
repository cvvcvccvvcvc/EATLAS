"""ClinVar enrichment analyses adjusted for phyloP100way."""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from analytics.vep.consequences import (
    VALIDATION_CONSEQUENCE_BITS as CONSEQUENCE_BITS,
    VALIDATION_CONSEQUENCE_OPTIONS as CONSEQUENCE_OPTIONS,
    validation_consequence_membership_mask as consequence_membership_mask,
    validation_consequence_memberships_text as consequence_memberships_text,
)
from analytics.io.performance import PerformanceProfile, profile_stage
from .statistics import benjamini_hochberg, enrichment_result, mantel_haenszel_adjusted
from .target_context import context_at, read_disjoint_contexts
from genomics.variants import changed_target_position, parse_variant_key


SCORE_COLUMN = "phyloP100way"
SPLINE_DF = 3
MAX_DISTRIBUTION_BINS = 40

PHYLOP_BANDS = [
    ("acceleration", "Nominal acceleration", "<= -1.30103"),
    ("central", "Central band", "-1.30103 to 1.30103"),
    ("conservation", "Nominal conservation", ">= 1.30103"),
]

VARIANT_TYPE_OPTIONS = [
    ("snv", "SNV"),
    ("indel", "INDEL"),
]

TARGET_CONTEXT_OPTIONS = [
    ("all", "All target contexts"),
    ("cds", "CDS"),
    ("utr", "UTR"),
    ("intron", "Intron"),
]
TARGET_CONTEXT_VALUES = ("cds", "utr", "other_exon", "intron", "other")

@dataclass(frozen=True)
class ConservationCohort:
    variants: pd.DataFrame
    summary: dict[str, int]


@dataclass(frozen=True)
class ConservationValidation:
    cohort: ConservationCohort
    unadjusted: pd.DataFrame
    fixed_bins: pd.DataFrame
    fixed_adjusted: pd.DataFrame
    continuous: pd.DataFrame
    distributions: pd.DataFrame
    r_versions: dict[str, str]


def build_conservation_cohort(
    *,
    universe: pd.DataFrame,
    conservation: pd.DataFrame,
    genes_tsv: Path,
    target_features_tsv: Path,
    score_column: str = SCORE_COLUMN,
    consequence_column: str = "clinvar_mc_terms",
) -> ConservationCohort:
    required = {
        "variant_key",
        "variant_type",
        "ref",
        "alt",
        "label_class",
        consequence_column,
        "gene_ids",
    }
    missing = required - set(universe.columns)
    if missing:
        raise ValueError(f"ClinVar universe missing required columns: {', '.join(sorted(missing))}")
    if score_column not in conservation.columns:
        raise ValueError(f"Conservation cache missing required score column: {score_column}")

    columns = [
        "variant_key",
        "variant_type",
        "ref",
        "alt",
        "label_class",
        consequence_column,
    ]
    columns.append("gene_ids")
    base = universe[columns].drop_duplicates("variant_key").copy()
    scores = conservation[["variant_key", score_column]].drop_duplicates("variant_key").copy()
    base = base.merge(scores, on="variant_key", how="left", validate="one_to_one")
    base[score_column] = pd.to_numeric(base[score_column], errors="coerce")
    base["variant_subtype"] = "other"
    snv_mask = base["variant_type"].astype(str) == "snv"
    ref_lengths = base["ref"].astype(str).str.len()
    alt_lengths = base["alt"].astype(str).str.len()
    base.loc[snv_mask, "variant_subtype"] = "snv"
    base.loc[~snv_mask & (alt_lengths > ref_lengths), "variant_subtype"] = "insertion"
    base.loc[~snv_mask & (ref_lengths > alt_lengths), "variant_subtype"] = "deletion"
    base["consequence_groups"] = base[consequence_column].map(consequence_memberships_text)
    base["consequence_mask"] = base[consequence_column].map(consequence_membership_mask)
    base["target_context"] = assign_target_contexts(
        base,
        genes_tsv=genes_tsv,
        target_features_tsv=target_features_tsv,
    )

    summary = {
        "allele_count": int(len(base)),
        "scored_allele_count": int(base[score_column].notna().sum()),
        "missing_score_count": int(base[score_column].isna().sum()),
        "snv_count": int((base["variant_subtype"] == "snv").sum()),
        "insertion_count": int((base["variant_subtype"] == "insertion").sum()),
        "deletion_count": int((base["variant_subtype"] == "deletion").sum()),
        "missing_consequence_count": int((base[consequence_column].astype(str) == "").sum()),
        "multiple_consequence_group_count": int(
            base["consequence_groups"].map(lambda value: len(split_memberships(value)) > 1).sum()
        ),
        **{
            f"target_context_{context}_count": int((base["target_context"] == context).sum())
            for context in TARGET_CONTEXT_VALUES
        },
    }
    return ConservationCohort(base, summary)


def assign_target_contexts(
    variants: pd.DataFrame,
    *,
    genes_tsv: Path,
    target_features_tsv: Path,
) -> pd.Series:
    genes = pd.read_csv(
        genes_tsv,
        sep="\t",
        compression="gzip" if genes_tsv.suffix == ".gz" else None,
        keep_default_na=False,
        usecols=["gene_id", "begin", "sequence_length"],
        dtype={"gene_id": str},
    )
    gene_begins = dict(zip(genes["gene_id"], genes["begin"].astype(int)))
    intervals = read_disjoint_contexts(
        target_features_tsv,
        dict(zip(genes["gene_id"], genes["sequence_length"].astype(int))),
    )
    starts = {
        gene_id: [start for start, _end, _context in values]
        for gene_id, values in intervals.items()
    }
    priority = {context: index for index, context in enumerate(TARGET_CONTEXT_VALUES)}

    def classify(variant_key: object, gene_ids: object) -> str:
        key = parse_variant_key(variant_key)
        if key is None:
            raise ValueError(f"Invalid normalized ClinVar variant_key: {variant_key}")
        memberships = []
        for gene_id in str(gene_ids).split("|"):
            if gene_id not in gene_begins or gene_id not in intervals:
                raise ValueError(
                    f"ClinVar allele {variant_key} references unknown target gene {gene_id}"
                )
            target_position = changed_target_position(key, gene_begins[gene_id])
            memberships.append(context_at(intervals[gene_id], target_position, starts[gene_id]))
        if not memberships:
            raise ValueError(f"ClinVar allele {variant_key} has no target gene membership")
        return min(set(memberships), key=lambda context: priority.get(context, len(priority)))

    values = [
        classify(row.variant_key, row.gene_ids)
        for row in variants[["variant_key", "gene_ids"]].itertuples(index=False)
    ]
    return pd.Series(values, index=variants.index, dtype="object")


def split_memberships(value: str) -> set[str]:
    return {item for item in str(value or "").split("|") if item}


def assign_phylop_band(values: pd.Series) -> pd.Categorical:
    numeric = pd.to_numeric(values, errors="coerce")
    labels = [key for key, _label, _range in PHYLOP_BANDS]
    categories = pd.Series(pd.NA, index=numeric.index, dtype="object")
    categories.loc[numeric <= -1.30103] = labels[0]
    categories.loc[(numeric > -1.30103) & (numeric < 1.30103)] = labels[1]
    categories.loc[numeric >= 1.30103] = labels[2]
    return pd.Categorical(categories, categories=labels, ordered=True)


def compute_conservation_validation(
    *,
    cohort: ConservationCohort,
    observed_by_strategy_type: dict[tuple[str, str], set[str]],
    strategies: list[str],
    analytics_dir: Path,
    eligible_gene_ids_by_strategy: dict[str, set[str]] | None = None,
    rscript: str | None = None,
    firth_workers: int = 1,
    performance_profile: PerformanceProfile | None = None,
) -> ConservationValidation:
    with profile_stage(performance_profile, "Unadjusted ClinVar association"):
        unadjusted = compute_unadjusted_enrichment(
            cohort=cohort.variants,
            observed_by_strategy_type=observed_by_strategy_type,
            strategies=strategies,
            eligible_gene_ids_by_strategy=eligible_gene_ids_by_strategy,
        )
    with profile_stage(performance_profile, "Fixed-band ClinVar association"):
        fixed_bins, fixed_adjusted = compute_fixed_band_enrichment(
            cohort=cohort.variants,
            observed_by_strategy_type=observed_by_strategy_type,
            strategies=strategies,
            eligible_gene_ids_by_strategy=eligible_gene_ids_by_strategy,
        )
    with profile_stage(performance_profile, "Continuous ClinVar association"):
        continuous, distributions, versions = compute_continuous_firth(
            cohort=cohort.variants,
            observed_by_strategy_type=observed_by_strategy_type,
            strategies=strategies,
            analytics_dir=analytics_dir,
            eligible_gene_ids_by_strategy=eligible_gene_ids_by_strategy,
            rscript=rscript,
            firth_workers=firth_workers,
            performance_profile=performance_profile,
        )
    return ConservationValidation(
        cohort,
        unadjusted,
        fixed_bins,
        fixed_adjusted,
        continuous,
        distributions,
        versions,
    )


def compute_unadjusted_enrichment(
    *,
    cohort: pd.DataFrame,
    observed_by_strategy_type: dict[tuple[str, str], set[str]],
    strategies: list[str],
    eligible_gene_ids_by_strategy: dict[str, set[str]] | None = None,
) -> pd.DataFrame:
    rows = []
    for strategy, variant_key, target_context_key, consequence_key, working in selector_frames(
        cohort, observed_by_strategy_type, strategies, eligible_gene_ids_by_strategy
    ):
        result = enrichment_for_subset(working, strategy)
        reason = two_by_two_estimability_reason(working)
        rows.append(
            {
                "strategy": strategy,
                "variant_type": variant_key,
                "target_context": target_context_key,
                "consequence": consequence_key,
                "usable_rows": int(len(working)),
                "benign_observed": result.benign_observed,
                "pathogenic_observed": result.pathogenic_observed,
                "benign_not_observed": result.benign_not_observed,
                "pathogenic_not_observed": result.pathogenic_not_observed,
                "odds_ratio": result.odds_ratio if not reason else float("nan"),
                "ci_low": result.ci_low if not reason else float("nan"),
                "ci_high": result.ci_high if not reason else float("nan"),
                "fisher_p": result.fisher_p if not reason else float("nan"),
                "status": "not_estimable" if reason else "estimated",
                "reason": reason,
            }
        )
    results = pd.DataFrame(rows)
    add_grouped_bh(
        results,
        "fisher_p",
        "fisher_q",
        ["variant_type", "target_context", "consequence"],
    )
    return results


def compute_fixed_band_enrichment(
    *,
    cohort: pd.DataFrame,
    observed_by_strategy_type: dict[tuple[str, str], set[str]],
    strategies: list[str],
    eligible_gene_ids_by_strategy: dict[str, set[str]] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    bin_rows: list[dict[str, object]] = []
    adjusted_rows: list[dict[str, object]] = []
    for strategy, variant_key, target_context_key, consequence_key, working in selector_frames(
        cohort, observed_by_strategy_type, strategies, eligible_gene_ids_by_strategy
    ):
        working = working[np.isfinite(working[SCORE_COLUMN])].copy()
        working["band"] = assign_phylop_band(working[SCORE_COLUMN])
        strata = []
        for index, (band_key, band_label, range_text) in enumerate(PHYLOP_BANDS, start=1):
            subset = working[working["band"] == band_key]
            result = enrichment_for_subset(subset, band_label)
            strata.append(result)
            table_reason = two_by_two_estimability_reason(subset)
            bin_rows.append(
                {
                    "strategy": strategy,
                    "variant_type": variant_key,
                    "target_context": target_context_key,
                    "consequence": consequence_key,
                    "band_index": index,
                    "band": band_key,
                    "band_label": band_label,
                    "band_range": range_text,
                    "row_count": int(len(subset)),
                    "benign_observed": result.benign_observed,
                    "pathogenic_observed": result.pathogenic_observed,
                    "benign_not_observed": result.benign_not_observed,
                    "pathogenic_not_observed": result.pathogenic_not_observed,
                    "odds_ratio": result.odds_ratio if not table_reason else float("nan"),
                    "ci_low": result.ci_low if not table_reason else float("nan"),
                    "ci_high": result.ci_high if not table_reason else float("nan"),
                    "fisher_p": result.fisher_p if not table_reason else float("nan"),
                    "status": "not_estimable" if table_reason else "estimated",
                    "reason": table_reason,
                }
            )

        reason = fixed_estimability_reason(working)
        adjusted = mantel_haenszel_adjusted(strata) if not reason else None
        if adjusted is None or not math.isfinite(adjusted.cmh_p):
            reason = reason or "The fixed-band tables do not provide an estimable common odds ratio."
            adjusted_rows.append(
                fixed_adjusted_row(
                    strategy,
                    variant_key,
                    target_context_key,
                    consequence_key,
                    working,
                    reason=reason,
                )
            )
        else:
            finite_effect = (
                math.isfinite(adjusted.odds_ratio_mh)
                and math.isfinite(adjusted.ci_low)
                and math.isfinite(adjusted.ci_high)
            )
            adjusted_rows.append(
                fixed_adjusted_row(
                    strategy,
                    variant_key,
                    target_context_key,
                    consequence_key,
                    working,
                    odds_ratio=adjusted.odds_ratio_mh,
                    ci_low=adjusted.ci_low,
                    ci_high=adjusted.ci_high,
                    cmh_chi2=adjusted.cmh_chi2,
                    cmh_p=adjusted.cmh_p,
                    status="estimated" if finite_effect else "test_only",
                    reason=(
                        "CMH p-value is available, but the common OR and CI are not finite because of separation."
                        if not finite_effect
                        else ""
                    ),
                )
            )

    bins = pd.DataFrame(bin_rows)
    adjusted = pd.DataFrame(adjusted_rows)
    add_grouped_bh(
        bins,
        "fisher_p",
        "fisher_q",
        ["variant_type", "target_context", "consequence", "band"],
    )
    add_grouped_bh(
        adjusted,
        "cmh_p",
        "cmh_q",
        ["variant_type", "target_context", "consequence"],
    )
    return bins, adjusted


def selector_frames(
    cohort: pd.DataFrame,
    observed_by_strategy_type: dict[tuple[str, str], set[str]],
    strategies: Iterable[str],
    eligible_gene_ids_by_strategy: dict[str, set[str]] | None = None,
):
    for strategy in strategies:
        observed_keys = set(observed_by_strategy_type.get((strategy, "snv"), set()))
        observed_keys.update(observed_by_strategy_type.get((strategy, "indel"), set()))
        strategy_frame = cohort
        if eligible_gene_ids_by_strategy is not None:
            eligible = eligible_gene_ids_by_strategy.get(strategy, set())
            strategy_frame = strategy_frame[
                strategy_frame["gene_ids"].map(
                    lambda value: bool(eligible.intersection(str(value).split("|")))
                )
            ]
        strategy_frame = strategy_frame.copy()
        strategy_frame["ALT_observed"] = strategy_frame["variant_key"].astype(str).isin(observed_keys).astype(int)
        for variant_key, _variant_label in VARIANT_TYPE_OPTIONS:
            type_frame = filter_variant_type(strategy_frame, variant_key)
            for target_context_key, _target_context_label in TARGET_CONTEXT_OPTIONS:
                context_frame = filter_target_context(type_frame, target_context_key)
                for consequence_key, _consequence_label in CONSEQUENCE_OPTIONS:
                    yield (
                        strategy,
                        variant_key,
                        target_context_key,
                        consequence_key,
                        filter_consequence(context_frame, consequence_key),
                    )


def filter_variant_type(frame: pd.DataFrame, variant_key: str) -> pd.DataFrame:
    if variant_key == "indel":
        return frame[frame["variant_subtype"].isin({"insertion", "deletion"})]
    return frame[frame["variant_subtype"] == variant_key]


def filter_target_context(frame: pd.DataFrame, target_context_key: str) -> pd.DataFrame:
    if target_context_key == "all":
        return frame
    return frame[frame["target_context"] == target_context_key]


def filter_consequence(frame: pd.DataFrame, consequence_key: str) -> pd.DataFrame:
    if consequence_key == "all":
        return frame
    bit = CONSEQUENCE_BITS[consequence_key]
    mask = (frame["consequence_mask"].astype(int) & bit) != 0
    return frame.loc[mask]


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


def two_by_two_estimability_reason(working: pd.DataFrame) -> str:
    if working.empty:
        return "No scored alleles in this band."
    if working["label_class"].nunique() < 2:
        return "Both B/LB and P/LP alleles are required."
    if working["ALT_observed"].nunique() < 2:
        return "Both ALT-observed and ALT-not-observed alleles are required."
    return ""


def fixed_estimability_reason(working: pd.DataFrame) -> str:
    if working.empty:
        return "No ClinVar B/LB or P/LP alleles in this view have phyloP100way scores."
    if working["label_class"].nunique() < 2:
        return "Both B/LB and P/LP alleles are required."
    if working["ALT_observed"].nunique() < 2:
        return "Both ALT-observed and ALT-not-observed alleles are required."
    if working["band"].nunique() < 2:
        return "At least two populated phyloP bands are required for a pooled adjusted estimate."
    return ""


def fixed_adjusted_row(
    strategy: str,
    variant_type: str,
    target_context: str,
    consequence: str,
    working: pd.DataFrame,
    *,
    reason: str = "",
    odds_ratio: float = float("nan"),
    ci_low: float = float("nan"),
    ci_high: float = float("nan"),
    cmh_chi2: float = float("nan"),
    cmh_p: float = float("nan"),
    status: str = "",
) -> dict[str, object]:
    return {
        "strategy": strategy,
        "variant_type": variant_type,
        "target_context": target_context,
        "consequence": consequence,
        "usable_rows": int(len(working)),
        "populated_bands": int(working["band"].nunique()) if not working.empty else 0,
        "odds_ratio_mh": odds_ratio,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "cmh_chi2": cmh_chi2,
        "cmh_p": cmh_p,
        "status": status or ("not_estimable" if reason else "estimated"),
        "reason": reason,
    }


def compute_continuous_firth(
    *,
    cohort: pd.DataFrame,
    observed_by_strategy_type: dict[tuple[str, str], set[str]],
    strategies: list[str],
    analytics_dir: Path,
    eligible_gene_ids_by_strategy: dict[str, set[str]] | None = None,
    rscript: str | None = None,
    firth_workers: int = 1,
    performance_profile: PerformanceProfile | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    if firth_workers < 1:
        raise ValueError("firth_workers must be >= 1")
    model_rows: list[dict[str, object]] = []
    distribution_rows: list[dict[str, object]] = []
    fit_specs: list[dict[str, str]] = []
    model_data = cohort[np.isfinite(cohort[SCORE_COLUMN])].copy()
    observation_columns: dict[str, str] = {}
    eligibility_columns: dict[str, str] = {}
    for index, strategy in enumerate(strategies):
        observation_column = f"observed_{index}"
        eligibility_column = f"eligible_{index}"
        observation_columns[strategy] = observation_column
        eligibility_columns[strategy] = eligibility_column
        keys = set(observed_by_strategy_type.get((strategy, "snv"), set()))
        keys.update(observed_by_strategy_type.get((strategy, "indel"), set()))
        model_data[observation_column] = model_data["variant_key"].astype(str).isin(keys).astype(int)
        if eligible_gene_ids_by_strategy is None:
            model_data[eligibility_column] = 1
        else:
            eligible = eligible_gene_ids_by_strategy.get(strategy, set())
            model_data[eligibility_column] = model_data["gene_ids"].map(
                lambda value: int(bool(eligible.intersection(str(value).split("|"))))
            )

    with profile_stage(performance_profile, "Firth selector preparation") as timing:
        for analysis_index, (
            strategy,
            variant_key,
            target_context_key,
            consequence_key,
            working,
        ) in enumerate(
            selector_frames(
                cohort,
                observed_by_strategy_type,
                strategies,
                eligible_gene_ids_by_strategy,
            )
        ):
            working = working[np.isfinite(working[SCORE_COLUMN])].copy()
            analysis_id = f"model_{analysis_index}"
            metrics = continuous_metrics(working)
            reason = continuous_estimability_reason(working)
            model_rows.append(
                {
                    "analysis_id": analysis_id,
                    "strategy": strategy,
                    "variant_type": variant_key,
                    "target_context": target_context_key,
                    "consequence": consequence_key,
                    **metrics,
                    "spline_df": SPLINE_DF,
                    "odds_ratio": float("nan"),
                    "ci_low": float("nan"),
                    "ci_high": float("nan"),
                    "plr_p": float("nan"),
                    "status": "not_estimable" if reason else "pending",
                    "reason": reason,
                }
            )
            distribution_rows.extend(
                distribution_detail_rows(
                    working,
                    strategy,
                    variant_key,
                    target_context_key,
                    consequence_key,
                )
            )
            if not reason:
                fit_specs.append(
                    {
                        "analysis_id": analysis_id,
                        "strategy": strategy,
                        "observation_column": observation_columns[strategy],
                        "eligibility_column": eligibility_columns[strategy],
                        "variant_type": variant_key,
                        "target_context": target_context_key,
                        "consequence": consequence_key,
                    }
                )
        timing["metrics"] = {
            "selectors": int(len(model_rows)),
            "estimable_models": int(len(fit_specs)),
        }

    results = pd.DataFrame(model_rows)
    versions: dict[str, str] = {}
    if fit_specs:
        with profile_stage(performance_profile, "Firth model fitting") as timing:
            fitted, versions = run_firth_models(
                model_data=model_data,
                specs=pd.DataFrame(fit_specs),
                analytics_dir=analytics_dir,
                rscript=rscript,
                workers=firth_workers,
            )
            timing["metrics"] = {
                "models": int(len(fit_specs)),
                "workers": int(min(firth_workers, len(fit_specs))),
            }
        fitted_by_id = fitted.set_index("analysis_id").to_dict("index")
        for index, row in results.iterrows():
            fitted_row = fitted_by_id.get(str(row["analysis_id"]))
            if fitted_row is None:
                continue
            for column in ["odds_ratio", "ci_low", "ci_high", "plr_p", "status", "reason"]:
                results.at[index, column] = fitted_row.get(column, results.at[index, column])

    add_grouped_bh(
        results,
        "plr_p",
        "plr_q",
        ["variant_type", "target_context", "consequence"],
    )
    return results, pd.DataFrame(distribution_rows), versions


def continuous_metrics(working: pd.DataFrame) -> dict[str, object]:
    observed_scores = working.loc[working["ALT_observed"] == 1, SCORE_COLUMN]
    not_observed_scores = working.loc[working["ALT_observed"] == 0, SCORE_COLUMN]
    overlap_low = max(observed_scores.min(), not_observed_scores.min()) if len(observed_scores) and len(not_observed_scores) else float("nan")
    overlap_high = min(observed_scores.max(), not_observed_scores.max()) if len(observed_scores) and len(not_observed_scores) else float("nan")
    total_low = working[SCORE_COLUMN].min() if not working.empty else float("nan")
    total_high = working[SCORE_COLUMN].max() if not working.empty else float("nan")
    total_span = total_high - total_low if math.isfinite(total_low) and math.isfinite(total_high) else float("nan")
    overlap_span = max(0.0, overlap_high - overlap_low) if math.isfinite(overlap_low) and math.isfinite(overlap_high) else float("nan")
    overlap_fraction = overlap_span / total_span if math.isfinite(total_span) and total_span > 0 else float("nan")
    return {
        "usable_rows": int(len(working)),
        "benign_rows": int((working["label_class"] == "benign").sum()),
        "pathogenic_rows": int((working["label_class"] == "pathogenic").sum()),
        "observed_rows": int((working["ALT_observed"] == 1).sum()),
        "not_observed_rows": int((working["ALT_observed"] == 0).sum()),
        "benign_observed": int(((working["label_class"] == "benign") & (working["ALT_observed"] == 1)).sum()),
        "pathogenic_observed": int(
            ((working["label_class"] == "pathogenic") & (working["ALT_observed"] == 1)).sum()
        ),
        "benign_not_observed": int(
            ((working["label_class"] == "benign") & (working["ALT_observed"] == 0)).sum()
        ),
        "pathogenic_not_observed": int(
            ((working["label_class"] == "pathogenic") & (working["ALT_observed"] == 0)).sum()
        ),
        "score_min": float(total_low) if not working.empty else float("nan"),
        "score_max": float(total_high) if not working.empty else float("nan"),
        "overlap_low": float(overlap_low),
        "overlap_high": float(overlap_high),
        "overlap_fraction": float(overlap_fraction),
        "overlap_warning": bool(math.isfinite(overlap_fraction) and overlap_fraction < 0.1),
    }


def continuous_estimability_reason(working: pd.DataFrame) -> str:
    if working.empty:
        return "No ClinVar B/LB or P/LP alleles in this view have phyloP100way scores."
    if working["label_class"].nunique() < 2:
        return "Both B/LB and P/LP alleles are required."
    if working["ALT_observed"].nunique() < 2:
        return "Both ALT-observed and ALT-not-observed alleles are required."
    if working[SCORE_COLUMN].nunique() < SPLINE_DF + 1:
        return f"At least {SPLINE_DF + 1} distinct phyloP values are required for the spline."
    observed = working.loc[working["ALT_observed"] == 1, SCORE_COLUMN]
    not_observed = working.loc[working["ALT_observed"] == 0, SCORE_COLUMN]
    if min(observed.max(), not_observed.max()) < max(observed.min(), not_observed.min()):
        return "ALT-observed and ALT-not-observed phyloP ranges do not overlap."
    return ""


def distribution_detail_rows(
    working: pd.DataFrame,
    strategy: str,
    variant_type: str,
    target_context: str,
    consequence: str,
) -> list[dict[str, object]]:
    groups = [(1, "ALT observed"), (0, "ALT not observed")]
    values_by_group = {
        label: working.loc[working["ALT_observed"] == observed_value, SCORE_COLUMN].to_numpy(float)
        for observed_value, label in groups
    }
    nonempty = [values for values in values_by_group.values() if len(values)]
    if not nonempty:
        return []
    combined = np.concatenate(nonempty)
    if float(np.min(combined)) == float(np.max(combined)):
        center = float(combined[0])
        padding = max(abs(center) * 0.05, 0.5)
        edges = np.asarray([center - padding, center + padding])
    else:
        edges = np.histogram_bin_edges(combined, bins="fd")
        if len(edges) - 1 > MAX_DISTRIBUTION_BINS:
            edges = np.linspace(float(np.min(combined)), float(np.max(combined)), MAX_DISTRIBUTION_BINS + 1)

    rows = []
    for _observed_value, label in groups:
        values = values_by_group[label]
        if not len(values):
            continue
        q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
        iqr = q3 - q1
        lower_whisker = float(np.min(values[values >= q1 - 1.5 * iqr]))
        upper_whisker = float(np.max(values[values <= q3 + 1.5 * iqr]))
        counts, _ = np.histogram(values, bins=edges)
        rows.extend(
            {
                "strategy": strategy,
                "variant_type": variant_type,
                "target_context": target_context,
                "consequence": consequence,
                "group": label,
                "bin_left": float(left),
                "bin_right": float(right),
                "count": int(count),
                "fraction": float(count / len(values)),
                "group_count": int(len(values)),
                "q1": float(q1),
                "median": float(median),
                "q3": float(q3),
                "lower_whisker": lower_whisker,
                "upper_whisker": upper_whisker,
            }
            for left, right, count in zip(edges[:-1], edges[1:], counts)
        )
    return rows


def add_grouped_bh(
    frame: pd.DataFrame,
    p_column: str,
    q_column: str,
    group_columns: list[str],
) -> None:
    if frame.empty:
        frame[q_column] = pd.Series(dtype=float)
        return
    frame[q_column] = float("nan")
    for _group, indices in frame.groupby(group_columns, sort=False).groups.items():
        frame.loc[indices, q_column] = benjamini_hochberg(frame.loc[indices, p_column].tolist())


def run_firth_models(
    *,
    model_data: pd.DataFrame,
    specs: pd.DataFrame,
    analytics_dir: Path,
    rscript: str | None = None,
    workers: int = 1,
) -> tuple[pd.DataFrame, dict[str, str]]:
    if workers < 1:
        raise ValueError("workers must be >= 1")
    executable = rscript or shutil.which("Rscript")
    if not executable:
        raise FileNotFoundError(
            "Rscript was not found. Run the report in envs/analytics.yml; ordinary logistic regression is not used as a fallback."
        )
    script_path = Path(__file__).with_name("firth_logistic.R")
    if not script_path.exists():
        raise FileNotFoundError(script_path)
    analytics_dir.mkdir(parents=True, exist_ok=True)
    columns = [
        "variant_key",
        "label_class",
        "variant_subtype",
        "target_context",
        "consequence_groups",
        SCORE_COLUMN,
        *specs["observation_column"].drop_duplicates().tolist(),
        *specs["eligibility_column"].drop_duplicates().tolist(),
    ]
    export = model_data[columns].copy()
    export = export.rename(columns={SCORE_COLUMN: "score"})
    export["benign"] = (export.pop("label_class") == "benign").astype(int)

    with tempfile.TemporaryDirectory(prefix=".firth_", dir=analytics_dir) as temporary:
        temp_dir = Path(temporary)
        data_path = temp_dir / "cohort.tsv"
        specs_path = temp_dir / "specs.tsv"
        output_path = temp_dir / "results.tsv"
        versions_path = temp_dir / "versions.tsv"
        export.to_csv(data_path, sep="\t", index=False)
        specs.to_csv(specs_path, sep="\t", index=False)
        environment = os.environ.copy()
        for variable in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            environment[variable] = "1"
        proc = subprocess.run(
            [
                str(executable),
                "--vanilla",
                str(script_path),
                str(data_path),
                str(specs_path),
                str(output_path),
                str(versions_path),
                str(min(workers, len(specs))),
            ],
            text=True,
            capture_output=True,
            check=False,
            env=environment,
        )
        if proc.returncode != 0:
            message = proc.stderr.strip() or proc.stdout.strip() or "unknown R error"
            raise RuntimeError(f"Firth logistic regression failed: {message}")
        fitted = pd.read_csv(output_path, sep="\t", keep_default_na=False)
        for column in ["odds_ratio", "ci_low", "ci_high", "plr_p"]:
            fitted[column] = pd.to_numeric(fitted[column], errors="coerce")
        versions_frame = pd.read_csv(versions_path, sep="\t", keep_default_na=False)
        versions = dict(zip(versions_frame["component"], versions_frame["version"]))
    return fitted, versions


def validate_firth_runtime(rscript: str | None = None) -> dict[str, str]:
    """Fail before expensive conservation work when the required R package is absent."""

    executable = rscript or shutil.which("Rscript")
    if not executable:
        raise FileNotFoundError(
            "Rscript was not found. Run the report in envs/analytics.yml; "
            "ordinary logistic regression is not used as a fallback."
        )
    expression = (
        "suppressPackageStartupMessages(library(logistf)); "
        "cat(as.character(getRversion()), '\\t', "
        "as.character(packageVersion('logistf')), '\\n', sep='')"
    )
    proc = subprocess.run(
        [str(executable), "--vanilla", "-e", expression],
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        message = proc.stderr.strip() or proc.stdout.strip() or "unknown R error"
        raise RuntimeError(f"Firth runtime preflight failed: {message}")
    fields = proc.stdout.strip().split("\t")
    if len(fields) != 2 or not all(fields):
        raise RuntimeError(
            "Firth runtime preflight returned an invalid R/logistf version response"
        )
    return {"R": fields[0], "logistf": fields[1]}
