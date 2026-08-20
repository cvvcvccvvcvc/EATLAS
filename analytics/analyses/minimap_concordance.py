"""Concordance analysis for the fixed minimap2 asm10 and asm20 presets."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from analytics.io.performance import PerformanceProfile
from analytics.io.variant_source import sql_string
from .conservation_validation import (
    ConservationValidation,
    add_grouped_bh,
    compute_conservation_validation,
)


ASM10 = "minimap2_asm10"
ASM20 = "minimap2_asm20"
GROUP_OPTIONS = [
    (ASM10, "minimap2 asm10"),
    (ASM20, "minimap2 asm20"),
    ("minimap2_either", "Either preset"),
    ("minimap2_both", "Both presets"),
    ("minimap2_asm10_only", "asm10 only"),
    ("minimap2_asm20_only", "asm20 only"),
]


@dataclass(frozen=True)
class MinimapConcordanceAnalysis:
    available: bool
    reason: str
    candidate_summary: pd.DataFrame
    validation: ConservationValidation | None


def build_minimap_concordance_analysis(
    *,
    score_path: Path,
    observed_by_strategy_type: dict[tuple[str, str], set[str]],
    base_validation: ConservationValidation,
    strategies: list[str],
    eligible_gene_ids_by_strategy: dict[str, set[str]],
    analytics_dir: Path,
    firth_workers: int = 1,
    performance_profile: PerformanceProfile | None = None,
) -> MinimapConcordanceAnalysis:
    if not {ASM10, ASM20}.issubset(strategies):
        return MinimapConcordanceAnalysis(
            False,
            "Both minimap2_asm10 and minimap2_asm20 are required for this comparison.",
            pd.DataFrame(),
            None,
        )

    candidate_summary = compute_minimap_candidate_summary(score_path)
    synthetic_memberships = minimap_group_memberships(observed_by_strategy_type)
    synthetic_eligibility = minimap_group_eligibility(eligible_gene_ids_by_strategy)
    synthetic_groups = [key for key, _label in GROUP_OPTIONS if key not in {ASM10, ASM20}]
    synthetic_validation = compute_conservation_validation(
        cohort=base_validation.cohort,
        observed_by_strategy_type=synthetic_memberships,
        strategies=synthetic_groups,
        analytics_dir=analytics_dir,
        eligible_gene_ids_by_strategy=synthetic_eligibility,
        firth_workers=firth_workers,
        performance_profile=performance_profile,
    )
    combined = combine_minimap_validation(base_validation, synthetic_validation)
    return MinimapConcordanceAnalysis(True, "", candidate_summary, combined)


def compute_minimap_candidate_summary(score_path: Path) -> pd.DataFrame:
    duckdb = _import_duckdb()
    predicates = {
        ASM10: "has_asm10",
        ASM20: "has_asm20",
        "minimap2_either": "has_asm10 OR has_asm20",
        "minimap2_both": "has_asm10 AND has_asm20",
        "minimap2_asm10_only": "has_asm10 AND NOT has_asm20",
        "minimap2_asm20_only": "has_asm20 AND NOT has_asm10",
    }
    unions = []
    for group_key, predicate in predicates.items():
        unions.append(
            f"SELECT {sql_string(group_key)} AS group_key, variant_type, "
            "count(*) AS variant_count, count_if(gnomad_status = 'found') AS gnomad_found_count, "
            "count_if(gnomad_status IN ('found','not_found')) AS gnomad_eligible_count, "
            "count_if(gnomad_status = 'lookup_failed') AS gnomad_lookup_failed_count "
            f"FROM flags WHERE {predicate} GROUP BY variant_type"
        )
    query = (
        "WITH flags AS (SELECT variant_key, first(variant_type) AS variant_type, "
        "first(gnomad_status) AS gnomad_status, "
        f"bool_or(strategy = {sql_string(ASM10)}) AS has_asm10, "
        f"bool_or(strategy = {sql_string(ASM20)}) AS has_asm20 "
        "FROM read_parquet(?) WHERE strategy IN (?, ?) GROUP BY variant_key) "
        + " UNION ALL ".join(unions)
    )
    with duckdb.connect() as connection:
        frame = connection.execute(query, [str(score_path), ASM10, ASM20]).fetchdf()
    if frame.empty:
        return frame
    variant_types = frame.loc[
        frame["group_key"].eq("minimap2_either"), "variant_type"
    ].astype(str).tolist()
    complete_index = pd.MultiIndex.from_product(
        [[key for key, _label in GROUP_OPTIONS], variant_types],
        names=["group_key", "variant_type"],
    )
    frame = (
        frame.set_index(["group_key", "variant_type"])
        .reindex(complete_index)
        .reset_index()
    )
    count_columns = [
        "variant_count",
        "gnomad_found_count",
        "gnomad_eligible_count",
        "gnomad_lookup_failed_count",
    ]
    frame[count_columns] = frame[count_columns].fillna(0).astype(int)
    denominators = (
        frame[frame["group_key"].eq("minimap2_either")]
        .set_index("variant_type")["variant_count"]
        .astype(int)
        .to_dict()
    )
    frame["allele_fraction"] = [
        int(row.variant_count) / denominators.get(str(row.variant_type), 1)
        for row in frame.itertuples(index=False)
    ]
    frame["gnomad_found_fraction"] = [
        int(row.gnomad_found_count) / int(row.gnomad_eligible_count)
        if int(row.gnomad_eligible_count)
        else float("nan")
        for row in frame.itertuples(index=False)
    ]
    order = {key: index for index, (key, _label) in enumerate(GROUP_OPTIONS)}
    frame["group_order"] = frame["group_key"].map(order)
    return frame.sort_values(["variant_type", "group_order"]).drop(columns="group_order")


def minimap_group_memberships(
    observed: dict[tuple[str, str], set[str]],
) -> dict[tuple[str, str], set[str]]:
    memberships: dict[tuple[str, str], set[str]] = {}
    for variant_type in ("snv", "indel"):
        asm10 = set(observed.get((ASM10, variant_type), set()))
        asm20 = set(observed.get((ASM20, variant_type), set()))
        memberships.update(
            {
                ("minimap2_either", variant_type): asm10 | asm20,
                ("minimap2_both", variant_type): asm10 & asm20,
                ("minimap2_asm10_only", variant_type): asm10 - asm20,
                ("minimap2_asm20_only", variant_type): asm20 - asm10,
            }
        )
    return memberships


def minimap_group_eligibility(
    eligible: dict[str, set[str]],
) -> dict[str, set[str]]:
    asm10 = set(eligible.get(ASM10, set()))
    asm20 = set(eligible.get(ASM20, set()))
    shared = asm10 & asm20
    return {
        "minimap2_either": asm10 | asm20,
        "minimap2_both": shared,
        "minimap2_asm10_only": shared,
        "minimap2_asm20_only": shared,
    }


def combine_minimap_validation(
    base: ConservationValidation,
    synthetic: ConservationValidation,
) -> ConservationValidation:
    base_strategies = {ASM10, ASM20}

    def combine(frame: pd.DataFrame, synthetic_frame: pd.DataFrame) -> pd.DataFrame:
        selected = frame[frame["strategy"].astype(str).isin(base_strategies)].copy()
        return pd.concat([selected, synthetic_frame.copy()], ignore_index=True)

    unadjusted = combine(base.unadjusted, synthetic.unadjusted)
    fixed_bins = combine(base.fixed_bins, synthetic.fixed_bins)
    fixed_adjusted = combine(base.fixed_adjusted, synthetic.fixed_adjusted)
    continuous = combine(base.continuous, synthetic.continuous)
    distributions = combine(base.distributions, synthetic.distributions)
    add_grouped_bh(
        unadjusted,
        "fisher_p",
        "fisher_q",
        ["variant_type", "target_context", "consequence"],
    )
    add_grouped_bh(
        fixed_bins,
        "fisher_p",
        "fisher_q",
        ["variant_type", "target_context", "consequence", "band"],
    )
    add_grouped_bh(
        fixed_adjusted,
        "cmh_p",
        "cmh_q",
        ["variant_type", "target_context", "consequence"],
    )
    add_grouped_bh(
        continuous,
        "plr_p",
        "plr_q",
        ["variant_type", "target_context", "consequence"],
    )
    return ConservationValidation(
        base.cohort,
        unadjusted,
        fixed_bins,
        fixed_adjusted,
        continuous,
        distributions,
        {**base.r_versions, **synthetic.r_versions},
    )


def _import_duckdb():
    try:
        import duckdb
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Minimap2 concordance requires python-duckdb") from exc
    return duckdb
