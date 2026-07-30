"""Overview table and run-level report metrics."""

from __future__ import annotations

import numpy as np
import pandas as pd

from analytics.analyses.variant_summary import VariantSummary
from .components import (
    format_count_ratio,
    format_count_share,
    format_int,
    metric_cards,
    strategy_label,
    table_html,
)

def alignment_summary_for_report(summary: pd.DataFrame) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame()
    report = pd.DataFrame({"Strategy": summary["strategy"].map(strategy_label)})
    if "gene_count" in summary.columns:
        report["Genes with result"] = summary["gene_count"]
    report["Orthologs aligned %"] = (
        summary["aligned_summary_row_count"] / summary["summary_row_count"].replace(0, np.nan)
    )
    report["Orthologs evaluated"] = summary["summary_row_count"]
    report["Orthologs aligned"] = summary["aligned_summary_row_count"]
    report["Raw support events"] = summary["event_count"]
    if "aligned_target_bp" in summary.columns:
        report["Aligned target bp"] = summary["aligned_target_bp"]
    return report


def merge_alignment_summary(strategy_stats: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    report_summary = alignment_summary_for_report(summary)
    if report_summary.empty:
        return strategy_stats
    return strategy_stats.merge(report_summary, on="Strategy", how="left")


def target_gene_coverage_for_report(cov: pd.DataFrame) -> pd.DataFrame:
    required = {"strategy", "feature_type", "length_bp", "covered_bases"}
    if cov.empty or not required.issubset(cov.columns):
        return pd.DataFrame(columns=["Strategy", "Target bases covered %"])
    genes = cov[cov["feature_type"].astype(str).str.lower().eq("gene")]
    if genes.empty:
        return pd.DataFrame(columns=["Strategy", "Target bases covered %"])
    coverage = (
        genes.groupby("strategy", as_index=False)
        .agg(Target_Length_bp=("length_bp", "sum"), Covered_Bases=("covered_bases", "sum"))
    )
    coverage["Target bases covered %"] = (
        coverage["Covered_Bases"] / coverage["Target_Length_bp"].replace(0, np.nan)
    )
    coverage["Strategy"] = coverage["strategy"].map(strategy_label)
    return coverage[["Strategy", "Target bases covered %"]]


def overview_strategy_table(
    variant_summary: VariantSummary,
    cov: pd.DataFrame,
    strategy_stats: pd.DataFrame,
    input_gene_count: int,
) -> pd.DataFrame:
    stats = strategy_stats.merge(variant_summary.unique_contribution, on="Strategy", how="left")
    stats = stats.merge(target_gene_coverage_for_report(cov), on="Strategy", how="left")
    stats["Unique To Strategy"] = stats["Unique To Strategy"].fillna(0)
    stats = stats.sort_values("Unique Variants", ascending=False, kind="mergesort")

    table = pd.DataFrame(
        {
            "Strategy": stats["Strategy"],
            "Genes with result": [
                format_count_ratio(count, input_gene_count)
                for count in stats["Genes with result"]
            ],
            "Candidate variants": stats["Unique Variants"],
            "Only this strategy": [
                format_count_share(count, total)
                for count, total in zip(stats["Unique To Strategy"], stats["Unique Variants"])
            ],
            "gnomAD matches": [
                format_count_share(count, total)
                for count, total in zip(stats["gnomAD Found"], stats["gnomAD Eligible"])
            ],
            "ClinVar matches": [
                format_count_share(count, total)
                for count, total in zip(stats["Found in ClinVar"], stats["Unique Variants"])
            ],
            "Orthologs aligned": [
                format_count_ratio(count, total)
                for count, total in zip(stats["Orthologs aligned"], stats["Orthologs evaluated"])
            ],
            "Target bases covered %": stats["Target bases covered %"],
        }
    )
    return table.reset_index(drop=True)


def build_overview(
    variant_summary: VariantSummary,
    cov: pd.DataFrame,
    strategy_stats: pd.DataFrame,
    annotation_manifest: dict,
    input_gene_count: int,
) -> list[str]:
    unique_variant_count = variant_summary.unique_variant_count
    all_strategy_count = variant_summary.all_strategy_variant_count
    annotation_warnings = int(annotation_manifest.get("failure_count", 0) or 0)
    cards = [
        ("Unique candidate variants", format_int(unique_variant_count)),
        ("Strategies", format_int(len(variant_summary.strategies))),
        ("Input genes", format_int(input_gene_count)),
        ("Candidates found by all strategies", format_count_share(all_strategy_count, unique_variant_count)),
        ("Annotation warnings", format_int(annotation_warnings)),
    ]
    sections = [metric_cards(cards)]
    sections.append("<h2>Strategies</h2>")
    sections.append(
        table_html(
            overview_strategy_table(variant_summary, cov, strategy_stats, input_gene_count),
            classes="table overview-table",
        )
    )
    return sections
