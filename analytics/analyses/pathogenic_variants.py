"""Focused characterization of candidate variants classified as P/LP by ClinVar."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from analytics.analyses.variant_summary import VariantSummary
from .clinvar_conditions import (
    condition_rows,
    compare_condition_distributions,
    global_condition_distribution,
)
from analytics.io.artifacts import write_tsv_atomic
from analytics.vep.consequences import display_consequence_group
from genomics.clinvar import pathogenic_subtype


PATHOGENIC_VARIANT_COLUMNS = [
    "variant_key",
    "gene_ids",
    "event_type",
    "pathogenic_subtype",
    "low_penetrance",
    "clinvar_review_stars",
    "clinvar_review_status",
    "clinvar_scv_count",
    "clinvar_ids",
    "clinvar_allele_id",
    "conditions",
    "condition_ids",
    "clinvar_hgvs",
    "clinvar_variant_type",
    "vep_primary_consequence",
    "vep_consequence_group",
    "vep_consequence_terms",
    "vep_transcript_id",
    "vep_mane_select",
    "vep_status",
    "gnomad_af",
    "support_ortholog_mean",
    "support_ortholog_min",
    "support_ortholog_max",
    "strategies",
]


@dataclass(frozen=True)
class PathogenicVariantAnalysis:
    variants_path: Path
    variants: pd.DataFrame
    star_counts: pd.DataFrame
    consequence_counts: pd.DataFrame
    condition_counts: pd.DataFrame
    support_rows: pd.DataFrame


def build_pathogenic_variant_analysis(
    *,
    summary: VariantSummary,
    clinvar_universe: pd.DataFrame,
    clinvar_vcf: Path,
    condition_cache_dir: Path,
    eligible_gene_ids_by_strategy: dict[str, set[str]],
    analytics_dir: Path,
) -> PathogenicVariantAnalysis:
    """Build one compact P/LP dataset and the plotting relations derived from it."""

    variants_path = analytics_dir / "pathogenic_clinvar_hits.tsv.gz"
    source = summary.pathogenic_rows.copy()
    if source.empty:
        variants = pd.DataFrame(columns=PATHOGENIC_VARIANT_COLUMNS)
        write_tsv_atomic(variants_path, variants)
        return PathogenicVariantAnalysis(
            variants_path=variants_path,
            variants=variants,
            star_counts=pd.DataFrame(),
            consequence_counts=summary.pathogenic_consequence_counts.copy(),
            condition_counts=pd.DataFrame(),
            support_rows=pd.DataFrame(),
        )

    source = source[source["clinvar_category"].astype(str).eq("P/LP")].copy()
    if source["variant_id"].duplicated().any():
        raise ValueError("P/LP detail rows repeat a normalized variant key")
    if len(source) != summary.pathogenic_variant_count:
        raise ValueError(
            "P/LP detail row count does not match the unique P/LP variant count: "
            f"{len(source)} != {summary.pathogenic_variant_count}"
        )

    subtype_values = source["clinvar_sig"].map(pathogenic_subtype)
    source["pathogenic_subtype"] = subtype_values.map(lambda value: value[0] or "P/LP")
    source["low_penetrance"] = subtype_values.map(lambda value: bool(value[1]))
    source["variant_key"] = source["variant_key"].fillna("").astype(str)
    if source["variant_key"].eq("").any():
        raise ValueError("P/LP detail rows require a normalized variant_key")

    universe_columns = (
        "variant_key",
        "clinvar_disease_names",
        "clinvar_disease_ids",
    )
    missing_universe_columns = set(universe_columns) - set(clinvar_universe.columns)
    if missing_universe_columns:
        raise ValueError(
            "ClinVar universe is missing P/LP condition columns: "
            + ", ".join(sorted(missing_universe_columns))
        )
    universe = clinvar_universe[list(universe_columns)].copy()
    if universe["variant_key"].duplicated().any():
        raise ValueError("ClinVar universe repeats a normalized variant key")
    source = source.merge(universe, on="variant_key", how="left", validate="one_to_one")

    memberships = condition_rows(source)
    condition_summary = (
        memberships.groupby("variant_key", sort=False)
        .agg(
            conditions=("condition", lambda values: "; ".join(dict.fromkeys(values))),
            condition_ids=("condition_ids", _join_condition_ids),
        )
        .reset_index()
        if not memberships.empty
        else pd.DataFrame(columns=["variant_key", "conditions", "condition_ids"])
    )
    source = source.merge(condition_summary, on="variant_key", how="left", validate="one_to_one")

    variants = pd.DataFrame(
        {
            "variant_key": source["variant_key"],
            "gene_ids": source["gene_id"],
            "event_type": source["event_type"],
            "pathogenic_subtype": source["pathogenic_subtype"],
            "low_penetrance": source["low_penetrance"],
            "clinvar_review_stars": source["clinvar_review_stars"],
            "clinvar_review_status": source["clinvar_revstat"],
            "clinvar_scv_count": source["clinvar_scv_count"],
            "clinvar_ids": source["clinvar_id"],
            "clinvar_allele_id": source["clinvar_allele_id"],
            "conditions": source.get("conditions", ""),
            "condition_ids": source.get("condition_ids", ""),
            "clinvar_hgvs": source["clinvar_hgvs"],
            "clinvar_variant_type": source["clinvar_variant_type"],
            "vep_primary_consequence": source["vep_primary_consequence"],
            "vep_consequence_group": [
                _consequence_groups(consequence, status)
                for consequence, status in zip(
                    source["vep_primary_consequence"], source["vep_status"]
                )
            ],
            "vep_consequence_terms": source["vep_consequence_terms"],
            "vep_transcript_id": source["vep_transcript_id"],
            "vep_mane_select": source["vep_mane_select"],
            "vep_status": source["vep_status"],
            "gnomad_af": pd.to_numeric(source["gnomad_af"], errors="coerce"),
            "support_ortholog_mean": pd.to_numeric(
                source["support_ortholog_mean"], errors="coerce"
            ),
            "support_ortholog_min": pd.to_numeric(
                source["support_ortholog_min"], errors="coerce"
            ),
            "support_ortholog_max": pd.to_numeric(
                source["support_ortholog_max"], errors="coerce"
            ),
            "strategies": source["strategies"],
        }
    )
    variants["_star_sort"] = pd.to_numeric(
        variants["clinvar_review_stars"], errors="coerce"
    )
    variants = (
        variants.sort_values(
            ["_star_sort", "support_ortholog_max", "variant_key"],
            ascending=[False, False, True],
            na_position="last",
            kind="mergesort",
        )
        .drop(columns="_star_sort")
        .reset_index(drop=True)
    )
    write_tsv_atomic(variants_path, variants[PATHOGENIC_VARIANT_COLUMNS])

    exploded = _explode_strategies(variants)
    star_counts = (
        exploded.groupby(
            ["strategy", "clinvar_review_stars"],
            as_index=False,
            dropna=False,
        )["variant_key"]
        .nunique()
        .rename(columns={"variant_key": "variant_count"})
    )

    condition_counts = compare_condition_distributions(
        variants,
        clinvar_universe,
        eligible_gene_ids_by_strategy,
        global_condition_distribution(clinvar_vcf, condition_cache_dir),
    )

    support_rows = _support_rows(summary.pathogenic_support_rows, variants)
    return PathogenicVariantAnalysis(
        variants_path=variants_path,
        variants=variants[PATHOGENIC_VARIANT_COLUMNS],
        star_counts=star_counts,
        consequence_counts=summary.pathogenic_consequence_counts.copy(),
        condition_counts=condition_counts,
        support_rows=support_rows,
    )


def _explode_strategies(variants: pd.DataFrame) -> pd.DataFrame:
    rows = variants[
        ["variant_key", "pathogenic_subtype", "clinvar_review_stars", "strategies"]
    ].copy()
    rows["strategy"] = rows["strategies"].fillna("").astype(str).str.split(",")
    rows = rows.explode("strategy")
    rows["strategy"] = rows["strategy"].astype(str).str.strip()
    return rows[rows["strategy"].ne("")].drop(columns="strategies")


def _consequence_groups(consequences: object, status: object) -> str:
    if pd.isna(consequences) or not str(consequences).strip() or "ok" not in str(status).split("|"):
        return "Not annotated"
    groups = [
        display_consequence_group(term)
        for term in str(consequences).split("|")
        if term
    ]
    return "|".join(dict.fromkeys(groups))


def _join_condition_ids(values: pd.Series) -> str:
    identifiers = []
    for value in values:
        identifiers.extend(item for item in str(value or "").split(",") if item)
    return "; ".join(dict.fromkeys(identifiers))


def _support_rows(support: pd.DataFrame, variants: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "variant_key",
        "gene_id",
        "strategy",
        "alt_support_ortholog_count",
        "alt_support_family_count",
        "site_aligned_ortholog_count",
    ]
    rows = support[
        support["variant_key"].isin(variants.loc[variants["event_type"].eq("snv"), "variant_key"])
    ].copy()
    if rows.empty:
        return pd.DataFrame(columns=columns)
    for column in columns[3:]:
        rows[column] = pd.to_numeric(rows[column], errors="raise")
    invalid = (
        rows[columns[3:]].isna().any(axis=1)
        | rows["alt_support_ortholog_count"].lt(1)
        | rows["alt_support_family_count"].lt(0)
        | rows["alt_support_family_count"].gt(rows["alt_support_ortholog_count"])
        | rows["alt_support_ortholog_count"].gt(rows["site_aligned_ortholog_count"])
    )
    if invalid.any():
        raise ValueError("Invalid exact support counts for P/LP SNVs")
    # Select a complete row, preserving its gene/depth/family context on ties.
    return (
        rows[columns]
        .sort_values(
            ["strategy", "variant_key", "alt_support_ortholog_count", "gene_id"],
            ascending=[True, True, False, True],
            kind="mergesort",
        )
        .drop_duplicates(["strategy", "variant_key"])
        .reset_index(drop=True)
    )
