"""Focused characterization of candidate variants classified as P/LP by ClinVar."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

from analytics.analyses.variant_summary import VariantSummary
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
    "phylop100way",
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
    evolution_rows: pd.DataFrame


def build_pathogenic_variant_analysis(
    *,
    summary: VariantSummary,
    clinvar_universe: pd.DataFrame,
    conservation_cohort: pd.DataFrame,
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
            evolution_rows=pd.DataFrame(),
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

    conservation = conservation_cohort[
        [
            column
            for column in conservation_cohort.columns
            if column in {"variant_key", "phyloP100way"}
        ]
    ].copy()
    missing_conservation_columns = {"variant_key", "phyloP100way"} - set(
        conservation.columns
    )
    if missing_conservation_columns:
        raise ValueError(
            "ClinVar conservation annotations are missing columns: "
            + ", ".join(sorted(missing_conservation_columns))
        )
    conservation = conservation[["variant_key", "phyloP100way"]]
    if conservation["variant_key"].duplicated().any():
        raise ValueError("ClinVar conservation annotations repeat a variant key")
    source = source.merge(conservation, on="variant_key", how="left", validate="one_to_one")

    condition_rows = _condition_rows(source)
    condition_summary = (
        condition_rows.groupby("variant_key", sort=False)
        .agg(
            conditions=("condition", lambda values: "; ".join(dict.fromkeys(values))),
            condition_ids=("condition_ids", _join_condition_ids),
        )
        .reset_index()
        if not condition_rows.empty
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
            "phylop100way": pd.to_numeric(source["phyloP100way"], errors="coerce"),
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
            ["strategy", "pathogenic_subtype", "clinvar_review_stars"],
            as_index=False,
            dropna=False,
        )["variant_key"]
        .nunique()
        .rename(columns={"variant_key": "variant_count"})
    )

    if condition_rows.empty:
        condition_counts = pd.DataFrame(
            columns=["strategy", "pathogenic_subtype", "condition", "variant_count"]
        )
    else:
        condition_memberships = condition_rows.merge(
            exploded[["variant_key", "strategy", "pathogenic_subtype"]],
            on="variant_key",
            how="inner",
            validate="many_to_many",
        )
        condition_counts = (
            condition_memberships.groupby(
                ["strategy", "pathogenic_subtype", "condition"], as_index=False
            )["variant_key"]
            .nunique()
            .rename(columns={"variant_key": "variant_count"})
        )

    evolution_rows = _evolution_rows(summary.pathogenic_support_rows, variants)
    return PathogenicVariantAnalysis(
        variants_path=variants_path,
        variants=variants[PATHOGENIC_VARIANT_COLUMNS],
        star_counts=star_counts,
        consequence_counts=summary.pathogenic_consequence_counts.copy(),
        condition_counts=condition_counts,
        evolution_rows=evolution_rows,
    )


def _explode_strategies(variants: pd.DataFrame) -> pd.DataFrame:
    rows = variants[
        ["variant_key", "pathogenic_subtype", "clinvar_review_stars", "strategies"]
    ].copy()
    rows["strategy"] = rows["strategies"].fillna("").astype(str).str.split(",")
    rows = rows.explode("strategy")
    rows["strategy"] = rows["strategy"].astype(str).str.strip()
    return rows[rows["strategy"].ne("")].drop(columns="strategies")


def _condition_rows(source: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for row in source.itertuples(index=False):
        found_named_condition = False
        names_text = str(getattr(row, "clinvar_disease_names", "") or "")
        ids_text = str(getattr(row, "clinvar_disease_ids", "") or "")
        for names_record, ids_record in zip_longest(
            names_text.split(";") if names_text else [],
            ids_text.split(";") if ids_text else [],
            fillvalue="",
        ):
            names = str(names_record).split("|")
            identifiers = str(ids_record).split("|")
            for name, condition_ids in zip_longest(names, identifiers, fillvalue=""):
                cleaned = _clean_condition_name(name)
                if cleaned is None:
                    continue
                found_named_condition = True
                rows.append(
                    {
                        "variant_key": str(row.variant_key),
                        "condition": cleaned,
                        "condition_ids": "" if condition_ids == "." else condition_ids,
                    }
                )
        if not found_named_condition:
            for name in str(getattr(row, "clinvar_disease", "") or "").split("|"):
                cleaned = _clean_condition_name(name)
                if cleaned is not None:
                    rows.append(
                        {
                            "variant_key": str(row.variant_key),
                            "condition": cleaned,
                            "condition_ids": "",
                        }
                    )
    if not rows:
        return pd.DataFrame(columns=["variant_key", "condition", "condition_ids"])
    return pd.DataFrame(rows).drop_duplicates(
        ["variant_key", "condition", "condition_ids"]
    )


def _clean_condition_name(value: object) -> str | None:
    raw = unquote(str(value or "")).strip()
    normalized = raw.lower().replace(" ", "_")
    if normalized in {"", ".", "not_provided", "not_specified", "not_applicable"}:
        return None
    return raw.replace("_", " ")


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


def _evolution_rows(support: pd.DataFrame, variants: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "variant_key",
        "gene_id",
        "strategy",
        "pathogenic_subtype",
        "phylop100way",
        "alt_support_ortholog_count",
        "alt_support_genus_count",
        "site_aligned_ortholog_count",
        "alt_support_fraction",
    ]
    if support.empty or variants.empty:
        return pd.DataFrame(columns=columns)
    rows = support.merge(
        variants[
            ["variant_key", "event_type", "pathogenic_subtype", "phylop100way"]
        ],
        on="variant_key",
        how="inner",
        validate="many_to_one",
    )
    rows = rows[rows["event_type"].astype(str).eq("snv")].copy()
    for column in (
        "phylop100way",
        "alt_support_ortholog_count",
        "alt_support_genus_count",
        "site_aligned_ortholog_count",
    ):
        rows[column] = pd.to_numeric(rows[column], errors="coerce")
    rows = rows[
        rows["phylop100way"].notna()
        & rows["site_aligned_ortholog_count"].gt(0)
        & rows["alt_support_ortholog_count"].notna()
    ].copy()
    rows["alt_support_fraction"] = (
        rows["alt_support_ortholog_count"] / rows["site_aligned_ortholog_count"]
    )
    invalid_fraction = ~rows["alt_support_fraction"].between(0, 1)
    if invalid_fraction.any():
        example = rows.loc[invalid_fraction, "variant_key"].iloc[0]
        raise ValueError(
            "Exact ALT-support fraction falls outside [0, 1] for P/LP SNV "
            f"{example}"
        )
    return rows[columns].sort_values(
        ["strategy", "variant_key", "gene_id"], kind="mergesort"
    )
