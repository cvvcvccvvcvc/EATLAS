"""Candidate-variant, external-evidence, and feature-profile sections."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from analytics.analyses.candidate_conservation import CandidateConservation
from analytics.analyses.variant_summary import StrategyOverlap, VariantSummary
from analytics.annotation.consequences import display_consequence_group as consequence_group
from .components import compact_figure, fig_html, sort_by_metric, strategy_label
from .config import (
    CLINVAR_COLORS,
    CLINVAR_ORDER,
    CONSEQUENCE_GROUP_COLORS,
    CONSEQUENCE_GROUP_ORDER,
    CONSEQUENCE_GROUP_TERMS,
    FEATURE_ORDER,
    PROFILE_FEATURE_ORDER,
    REVIEW_STAR_COLORS,
    REVIEW_STAR_ORDER,
    TARGET_CONTEXT_COLORS,
    TARGET_CONTEXT_LABELS,
    TARGET_CONTEXT_ORDER,
)
from .conservation import candidate_phylop_figure, candidate_phylop_summary_figure

def strategy_overlap_figure(overlap: StrategyOverlap | None):
    if overlap is None:
        return None
    order = np.arange(len(overlap.strategies))
    if len(order) > 2:
        distance = np.clip(1.0 - overlap.jaccard, 0.0, 1.0)
        np.fill_diagonal(distance, 0.0)
        order = leaves_list(linkage(squareform(distance, checks=False), method="average", optimal_ordering=True))
        left, right = int(order[0]), int(order[-1])
        mean_distance = distance.mean(axis=1)
        if (mean_distance[left], overlap.strategies[left]) < (mean_distance[right], overlap.strategies[right]):
            order = order[::-1]
    labels = [strategy_label(overlap.strategies[index]) for index in order]
    jaccard = overlap.jaccard[np.ix_(order, order)]
    intersections = overlap.intersections[np.ix_(order, order)]
    unions = overlap.unions[np.ix_(order, order)]

    fig = go.Figure(
        data=go.Heatmap(
            z=jaccard,
            x=labels,
            y=labels,
            text=np.vectorize(lambda value: f"{value:.0%}")(jaccard),
            customdata=np.dstack([intersections, unions]),
            colorscale="Blues",
            zmin=0,
            zmax=1,
            colorbar={"title": "Jaccard"},
            hovertemplate=(
                "%{y} vs %{x}<br>"
                "Jaccard: %{z:.1%}<br>"
                "Shared variants: %{customdata[0]:,}<br>"
                "Union variants: %{customdata[1]:,}<extra></extra>"
            ),
        )
    )
    fig.update_traces(texttemplate="%{text}", textfont_size=12)
    fig.update_layout(
        title="Pairwise strategy overlap",
        height=520,
        margin={"l": 135, "r": 30, "t": 95, "b": 35},
        template="plotly_white",
    )
    fig.update_xaxes(side="top", tickangle=-35, title_text=None, automargin=True)
    fig.update_yaxes(title_text=None, automargin=True)
    return fig


def review_star_category(row: pd.Series) -> str:
    stars = str(row.get("clinvar_review_stars", "") or "").strip()
    if stars in {"0", "1", "2", "3", "4"}:
        return stars
    return "Unmapped"


def group_consequence_counts(raw_counts: pd.DataFrame) -> pd.DataFrame:
    if raw_counts.empty:
        return pd.DataFrame(columns=["Strategy", "Consequence group", "Variant_Count", "Fraction"])
    counts = raw_counts.rename(columns={"strategy": "Strategy"}).copy()
    counts["Consequence group"] = counts["value"].map(consequence_group)
    counts = (
        counts.groupby(["Strategy", "Consequence group"], observed=True)["Variant_Count"]
        .sum()
        .reset_index()
    )
    totals = counts.groupby("Strategy", observed=True)["Variant_Count"].transform("sum")
    counts["Fraction"] = counts["Variant_Count"] / totals.replace(0, np.nan)
    counts["Consequence group"] = pd.Categorical(
        counts["Consequence group"], categories=CONSEQUENCE_GROUP_ORDER, ordered=True
    )
    return counts.sort_values(["Strategy", "Consequence group"])


def consequence_strategy_order(counts: pd.DataFrame) -> list[str]:
    if counts.empty:
        return []
    pivot = counts.pivot_table(
        index="Strategy",
        columns="Consequence group",
        values="Fraction",
        aggfunc="sum",
        fill_value=0,
        observed=True,
    )
    for column in CONSEQUENCE_GROUP_ORDER:
        if column not in pivot.columns:
            pivot[column] = 0.0
    pivot["impact_fraction"] = pivot["LoF/splice"] + pivot["Missense/inframe"]
    pivot["total_count"] = counts.groupby("Strategy", observed=True)["Variant_Count"].sum()
    return pivot.sort_values(["impact_fraction", "total_count"], ascending=False).index.tolist()


def gene_variant_distribution_counts(
    gene_counts: pd.DataFrame,
    strategy_stats: pd.DataFrame,
) -> pd.DataFrame:
    columns = [
        "Strategy",
        "Bin",
        "Bin_Order",
        "Gene_Count",
        "Gene_Fraction",
        "Genes_With_Result",
    ]
    if gene_counts.empty:
        return pd.DataFrame(columns=columns)

    counts = gene_counts.rename(columns={"strategy": "Strategy"}).copy()
    counts["Variant_Count"] = pd.to_numeric(counts["Variant_Count"], errors="raise").astype(int)

    eligible = {}
    if {"Strategy", "Genes with result"}.issubset(strategy_stats.columns):
        eligible = {
            str(strategy): int(value)
            for strategy, value in zip(
                strategy_stats["Strategy"],
                strategy_stats["Genes with result"],
            )
            if not pd.isna(value)
        }

    strategy_order = [str(value) for value in strategy_stats.get("Strategy", pd.Series(dtype=str))]
    strategy_order.extend(
        strategy for strategy in counts["Strategy"].astype(str).unique() if strategy not in strategy_order
    )

    def bin_order(value: int) -> int:
        return 0 if value == 0 else value.bit_length()

    def bin_label(order: int) -> str:
        if order == 0:
            return "0"
        if order == 1:
            return "1"
        lower = 1 << (order - 1)
        return f"{lower}-{(lower << 1) - 1}"

    counts["Bin_Order"] = counts["Variant_Count"].map(bin_order)
    grouped = (
        counts.groupby(["Strategy", "Bin_Order"], as_index=False, observed=True)
        .agg(Gene_Count=("gene_id", "nunique"))
    )
    zero_rows = []
    for strategy in strategy_order:
        observed = int(counts.loc[counts["Strategy"].astype(str).eq(strategy), "gene_id"].nunique())
        total = max(eligible.get(strategy, observed), observed)
        if total > observed:
            zero_rows.append({"Strategy": strategy, "Bin_Order": 0, "Gene_Count": total - observed})
    if zero_rows:
        grouped = pd.concat([grouped, pd.DataFrame(zero_rows)], ignore_index=True)

    max_order = int(grouped["Bin_Order"].max())
    complete = pd.MultiIndex.from_product(
        [strategy_order, range(max_order + 1)], names=["Strategy", "Bin_Order"]
    )
    grouped = (
        grouped.groupby(["Strategy", "Bin_Order"], as_index=False, observed=True)["Gene_Count"]
        .sum()
        .set_index(["Strategy", "Bin_Order"])
        .reindex(complete, fill_value=0)
        .reset_index()
    )
    observed_totals = grouped.groupby("Strategy", observed=True)["Gene_Count"].sum().astype(int)
    grouped["Genes_With_Result"] = grouped["Strategy"].map(
        lambda strategy: max(eligible.get(str(strategy), 0), int(observed_totals[strategy]))
    )
    grouped["Gene_Fraction"] = grouped["Gene_Count"] / grouped["Genes_With_Result"].replace(0, np.nan)
    grouped["Bin"] = grouped["Bin_Order"].map(bin_label)
    return grouped[columns]


def gene_variant_distribution_figure(
    gene_counts: pd.DataFrame,
    strategy_stats: pd.DataFrame,
):
    distribution = gene_variant_distribution_counts(gene_counts, strategy_stats)
    if distribution.empty:
        return None
    fig = go.Figure()
    for strategy, values in distribution.groupby("Strategy", sort=False, observed=True):
        fig.add_trace(
            go.Bar(
                x=values["Bin"],
                y=values["Gene_Count"],
                name=str(strategy),
                customdata=values[["Gene_Fraction", "Genes_With_Result"]],
                hovertemplate=(
                    "%{fullData.name}<br>Candidates per gene: %{x}<br>"
                    "Genes: %{y:,} (%{customdata[0]:.1%})<br>"
                    "Genes with result: %{customdata[1]:,}<extra></extra>"
                ),
            )
        )
    fig.update_layout(
        title="Candidate variants per gene",
        xaxis_title="Unique candidate variants per gene",
        yaxis_title="Genes",
        barmode="group",
        bargap=0.12,
        hovermode="closest",
    )
    compact_figure(fig, height=400, show_x_title=True)
    return fig


def top_gene_contribution_counts(
    gene_counts: pd.DataFrame,
    strategy_stats: pd.DataFrame,
    limit: int = 5,
) -> pd.DataFrame:
    columns = [
        "Strategy",
        "gene_id",
        "Rank",
        "Variant_Count",
        "Variant_Fraction",
        "Equal_Share",
        "Equal_Share_Ratio",
        "Top_Share",
    ]
    if gene_counts.empty:
        return pd.DataFrame(columns=columns)
    counts = gene_counts.rename(columns={"strategy": "Strategy"}).copy()
    counts["Variant_Count"] = pd.to_numeric(counts["Variant_Count"], errors="raise").astype(int)
    counts = counts.sort_values(
        ["Strategy", "Variant_Count", "gene_id"],
        ascending=[True, False, True],
        kind="mergesort",
    )
    counts["Rank"] = counts.groupby("Strategy", sort=False).cumcount() + 1
    totals = counts.groupby("Strategy", observed=True)["Variant_Count"].transform("sum")
    counts["Variant_Fraction"] = counts["Variant_Count"] / totals.replace(0, np.nan)

    eligible = strategy_stats[["Strategy", "Genes with result"]].copy()
    eligible["Equal_Share"] = 1.0 / pd.to_numeric(
        eligible["Genes with result"], errors="coerce"
    ).replace(0, np.nan)
    counts = counts.merge(eligible[["Strategy", "Equal_Share"]], on="Strategy", how="left")
    observed_gene_counts = counts.groupby("Strategy", observed=True)["gene_id"].transform("nunique")
    counts["Equal_Share"] = counts["Equal_Share"].fillna(1.0 / observed_gene_counts.replace(0, np.nan))
    counts["Equal_Share_Ratio"] = counts["Variant_Fraction"] / counts["Equal_Share"]
    counts = counts[counts["Rank"] <= limit].copy()
    counts["Top_Share"] = counts.groupby("Strategy", observed=True)["Variant_Fraction"].transform("sum")
    return counts[columns]


def top_gene_contribution_figure(
    gene_counts: pd.DataFrame,
    strategy_stats: pd.DataFrame,
    limit: int = 5,
):
    top = top_gene_contribution_counts(gene_counts, strategy_stats, limit=limit)
    if top.empty:
        return None
    fig = go.Figure()
    for strategy, values in top.groupby("Strategy", sort=False, observed=True):
        values = values.sort_values("Rank", kind="mergesort")
        rank_labels = [
            f"#{int(rank)}<br>{gene_id}"
            for rank, gene_id in zip(values["Rank"], values["gene_id"], strict=True)
        ]
        fig.add_trace(
            go.Bar(
                x=[[str(strategy)] * len(values), rank_labels],
                y=values["Variant_Fraction"],
                name=str(strategy),
                customdata=values[
                    [
                        "gene_id",
                        "Rank",
                        "Variant_Count",
                        "Equal_Share",
                        "Equal_Share_Ratio",
                        "Top_Share",
                    ]
                ],
                hovertemplate=(
                    "%{fullData.name}<br>Gene: %{customdata[0]}<br>Rank: %{customdata[1]}<br>"
                    "Candidates: %{customdata[2]:,}<br>Strategy share: %{y:.2%}<br>"
                    "Equal-share reference: %{customdata[3]:.2%}<br>"
                    "Observed / equal share: %{customdata[4]:.1f}x<br>"
                    f"Top-{limit} cumulative share: %{{customdata[5]:.2%}}<extra></extra>"
                ),
            )
        )
    fig.update_layout(
        title=f"Top {limit} contributing genes by strategy",
        yaxis_title="Share of strategy candidates",
        bargap=0.16,
    )
    fig.update_yaxes(tickformat=".0%")
    fig.update_xaxes(tickangle=-45)
    compact_figure(fig, height=430)
    return fig


def compact_list_text(value: str, max_items: int = 2, max_chars: int = 90) -> str:
    items = [item for item in re.split(r"[|,]", str(value or "")) if item and item != "."]
    if not items:
        return ""
    shown = "; ".join(items[:max_items])
    if len(items) > max_items:
        shown += f"; +{len(items) - max_items}"
    if len(shown) > max_chars:
        shown = shown[: max_chars - 1].rstrip() + "..."
    return shown


def format_strategy_list(value: str) -> str:
    strategies = [strategy_label(item.strip()) for item in str(value or "").split(",") if item.strip()]
    return ", ".join(strategies)


def pathogenic_variant_table(variants: pd.DataFrame) -> pd.DataFrame:
    pathogenic = variants[variants["clinvar_category"].astype(str) == "P/LP"].copy()
    if pathogenic.empty:
        return pd.DataFrame()

    pathogenic["Stars"] = pathogenic.apply(review_star_category, axis=1)
    pathogenic["Strategies"] = pathogenic["strategies"].map(format_strategy_list)
    pathogenic["Disease"] = pathogenic["clinvar_disease"].map(compact_list_text)
    pathogenic["HGVS"] = pathogenic["clinvar_hgvs"].map(lambda value: compact_list_text(value, max_items=1, max_chars=70))
    pathogenic["gnomAD AF"] = pathogenic["gnomad_af"]
    pathogenic["Ortholog support / strategy"] = pathogenic.apply(
        lambda row: (
            ""
            if pd.isna(row.get("support_ortholog_mean"))
            else (
                f"{float(row['support_ortholog_mean']):.1f} "
                f"({int(row['support_ortholog_min'])}-{int(row['support_ortholog_max'])})"
            )
        ),
        axis=1,
    )
    table = pd.DataFrame(
        {
            "Key": pathogenic["variant_id"],
            "Gene": pathogenic["gene_id"],
            "Event": pathogenic["event_type"],
            "ClinVar sig": pathogenic["clinvar_sig"],
            "Stars": pathogenic["Stars"],
            "Review status": pathogenic["clinvar_revstat"],
            "SCVs": pathogenic["clinvar_scv_count"],
            "ClinVar ID": pathogenic["clinvar_id"],
            "Allele ID": pathogenic["clinvar_allele_id"],
            "Disease": pathogenic["Disease"],
            "HGVS": pathogenic["HGVS"],
            "ClinVar type": pathogenic["clinvar_variant_type"],
            "gnomAD AF": pathogenic["gnomAD AF"],
            "gnomAD consequence": pathogenic["gnomad_csq"],
            **(
                {"VEP consequence": pathogenic["vep_primary_consequence"]}
                if "vep_primary_consequence" in pathogenic.columns
                else {}
            ),
            "Ortholog support / strategy": pathogenic["Ortholog support / strategy"],
            "Strategies": pathogenic["Strategies"],
        }
    )
    star_rank = table["Stars"].map({star: index for index, star in enumerate(REVIEW_STAR_ORDER[::-1])}).fillna(-1)
    table["_star_rank"] = star_rank
    table = table.sort_values(
        ["_star_rank", "SCVs", "Key"],
        ascending=[False, False, True],
        kind="mergesort",
    ).drop(columns=["_star_rank"])
    return table


def coverage_summary(cov: pd.DataFrame, feature_types: list[str] | None = None) -> pd.DataFrame:
    if cov.empty:
        return pd.DataFrame()
    if feature_types:
        cov = cov[cov["feature_type"].isin(feature_types)].copy()
    summary = (
        cov.groupby(["strategy", "feature_type"], as_index=False)
        .agg(
            Gene_Count=("gene_id", "nunique"),
            Feature_Count=("feature_id", "count"),
            Total_Length_bp=("length_bp", "sum"),
            Covered_Bases=("covered_bases", "sum"),
            Depth_Bases=("depth_bases", "sum"),
            Median_Breadth=("coverage_breadth", "median"),
            Median_Mean_Depth=("mean_depth", "median"),
            Median_Orthologs_Covered=("orthologs_covered", "median"),
        )
    )
    summary["Breadth_Weighted"] = summary["Covered_Bases"] / summary["Total_Length_bp"].replace(0, np.nan)
    summary["Mean_Depth_Weighted"] = summary["Depth_Bases"] / summary["Total_Length_bp"].replace(0, np.nan)
    summary["feature_type"] = pd.Categorical(summary["feature_type"], categories=FEATURE_ORDER, ordered=True)
    summary["strategy"] = summary["strategy"].astype(str).map(strategy_label)
    return summary.sort_values(["strategy", "feature_type"])


def build_variant_sections(
    variant_summary: VariantSummary,
    strategy_stats: pd.DataFrame,
) -> list[str]:
    sections = ["<h2>Strategy Concordance</h2>"]
    fig_overlap = strategy_overlap_figure(variant_summary.overlap)
    if fig_overlap is not None:
        sections.append(fig_html(fig_overlap))

    gene_distribution = gene_variant_distribution_figure(
        variant_summary.gene_variant_counts,
        strategy_stats,
    )
    top_genes = top_gene_contribution_figure(
        variant_summary.gene_variant_counts,
        strategy_stats,
    )
    if gene_distribution is not None or top_genes is not None:
        sections.append("<h2>Gene Concentration</h2>")
        if gene_distribution is not None:
            sections.append(fig_html(gene_distribution))
        if top_genes is not None:
            sections.append(fig_html(top_genes))

    sections.append("<h2>Variant Composition</h2>")
    counts = variant_summary.event_counts.copy()
    totals = counts.groupby("strategy", observed=True)["Variant_Count"].transform("sum")
    counts["Fraction"] = counts["Variant_Count"] / totals.replace(0, np.nan)
    snv_order = (
        counts[counts["event_type"].astype(str).str.lower() == "snv"]
        .sort_values("Fraction", ascending=False)
        ["strategy"]
        .tolist()
    )
    all_strategies = strategy_stats["Strategy"].tolist()
    order = snv_order + [strategy for strategy in all_strategies if strategy not in snv_order]
    fig_events = px.bar(
        counts,
        x="strategy",
        y="Fraction",
        color="event_type",
        barmode="stack",
        title="Variant type composition by strategy",
        category_orders={"strategy": order},
        labels={"strategy": "", "Fraction": "Variant fraction", "event_type": "Variant type"},
    )
    fig_events.update_layout(yaxis_tickformat=".0%")
    compact_figure(fig_events, height=360)
    sections.append(fig_html(fig_events))

    contexts = variant_summary.target_context_counts.copy()
    if not contexts.empty:
        totals = contexts.groupby("strategy", observed=True)["Variant_Count"].transform("sum")
        contexts["Fraction"] = contexts["Variant_Count"] / totals.replace(0, np.nan)
        contexts["Target context"] = contexts["target_context"].map(TARGET_CONTEXT_LABELS).fillna("Other")
        cds_order = (
            contexts[contexts["target_context"].astype(str).eq("cds")]
            .sort_values("Fraction", ascending=False)["strategy"]
            .tolist()
        )
        context_order = cds_order + [strategy for strategy in all_strategies if strategy not in cds_order]
        fig_context = px.bar(
            contexts,
            x="strategy",
            y="Fraction",
            color="Target context",
            barmode="stack",
            title="Target context composition by strategy",
            category_orders={
                "strategy": context_order,
                "Target context": [TARGET_CONTEXT_LABELS[item] for item in TARGET_CONTEXT_ORDER],
            },
            color_discrete_map=TARGET_CONTEXT_COLORS,
            labels={"strategy": "", "Fraction": "Variant fraction"},
            custom_data=["Variant_Count"],
        )
        fig_context.update_layout(yaxis_tickformat=".0%")
        fig_context.update_traces(
            hovertemplate="%{x}<br>%{fullData.name}: %{customdata[0]:,} (%{y:.1%})<extra></extra>"
        )
        compact_figure(fig_context, height=360)
        sections.append(fig_html(fig_context))
    return sections


def build_clinvar_gnomad_sections(
    variant_summary: VariantSummary,
    strategy_stats: pd.DataFrame,
) -> list[str]:
    sections = ["<h2>Population Evidence</h2>"]
    gnomad_rate = sort_by_metric(strategy_stats[["Strategy", "gnomAD found %"]], "gnomAD found %")
    gnomad_rate = gnomad_rate.merge(
        strategy_stats[["Strategy", "gnomAD Found", "gnomAD Eligible"]],
        on="Strategy",
        how="left",
    )
    fig_gnomad_rate = px.bar(
        gnomad_rate,
        x="Strategy",
        y="gnomAD found %",
        title="gnomAD hit rate by strategy",
        category_orders={"Strategy": gnomad_rate["Strategy"].tolist()},
        custom_data=["gnomAD Found", "gnomAD Eligible"],
    )
    fig_gnomad_rate.update_layout(yaxis_tickformat=".1%")
    fig_gnomad_rate.update_traces(
        hovertemplate=(
            "%{x}<br>Found in gnomAD: %{customdata[0]:,}<br>"
            "Completed lookups: %{customdata[1]:,}<br>Hit rate: %{y:.2%}<extra></extra>"
        )
    )
    compact_figure(fig_gnomad_rate)
    sections.append(fig_html(fig_gnomad_rate))

    af_summary = variant_summary.gnomad_af_summary.sort_values("Median", ascending=False)
    if not af_summary.empty:
        fig_af = go.Figure()
        for width, low, high, color, name in [
            (3, "Q05", "Q95", "#9ecae1", "5-95% interval"),
            (10, "Q25", "Q75", "#3182bd", "Interquartile interval"),
        ]:
            x_values, y_values = [], []
            for row in af_summary.itertuples(index=False):
                x_values.extend([getattr(row, low), getattr(row, high), None])
                y_values.extend([row.Strategy, row.Strategy, None])
            fig_af.add_trace(go.Scatter(
                x=x_values, y=y_values, mode="lines", line={"width": width, "color": color},
                name=name, hoverinfo="skip",
            ))
        fig_af.add_trace(go.Scatter(
            x=af_summary["Median"],
            y=af_summary["Strategy"],
            mode="markers",
            marker={"size": 9, "color": "#08306b"},
            name="Median",
            customdata=af_summary[["Count", "Q05", "Q25", "Q75", "Q95"]],
            hovertemplate=(
                "%{y}<br>Variants with AF &gt; 0: %{customdata[0]:,}<br>"
                "5th percentile: %{customdata[1]:.3f}<br>Q1: %{customdata[2]:.3f}<br>"
                "Median: %{x:.3f}<br>Q3: %{customdata[3]:.3f}<br>"
                "95th percentile: %{customdata[4]:.3f}<extra></extra>"
            ),
        ))
        fig_af.update_layout(title="gnomAD allele frequency among exact hits", xaxis_title="log10 gnomAD AF")
        compact_figure(fig_af, height=380, show_x_title=True)
        sections.append(fig_html(fig_af))
    else:
        sections.append("<p>No non-zero gnomAD AF values were found.</p>")

    consequence_counts = group_consequence_counts(variant_summary.consequence_counts)
    if not consequence_counts.empty:
        order = consequence_strategy_order(consequence_counts)
        fig_conseq = px.bar(
            consequence_counts,
            x="Strategy",
            y="Fraction",
            color="Consequence group",
            barmode="stack",
            title=f"{variant_summary.consequence_source} consequence mix among candidates",
            category_orders={"Strategy": order, "Consequence group": CONSEQUENCE_GROUP_ORDER},
            color_discrete_map=CONSEQUENCE_GROUP_COLORS,
            labels={"Strategy": "", "Fraction": "Within-strategy fraction", "Consequence group": "Consequence group"},
        )
        fig_conseq.update_layout(yaxis_tickformat=".0%")
        compact_figure(fig_conseq, height=360)
        sections.append(fig_html(fig_conseq))

    sections.append("<h2>Clinical Evidence</h2>")
    clinvar_rate = sort_by_metric(strategy_stats[["Strategy", "ClinVar found %"]], "ClinVar found %")
    clinvar_rate = clinvar_rate.merge(
        strategy_stats[["Strategy", "Found in ClinVar", "Unique Variants"]], on="Strategy", how="left"
    )
    fig_clin_rate = px.bar(
        clinvar_rate,
        x="Strategy",
        y="ClinVar found %",
        title="ClinVar hit rate by strategy",
        category_orders={"Strategy": clinvar_rate["Strategy"].tolist()},
        custom_data=["Found in ClinVar", "Unique Variants"],
    )
    fig_clin_rate.update_layout(yaxis_tickformat=".2%")
    fig_clin_rate.update_traces(
        hovertemplate=(
            "%{x}<br>Found in ClinVar: %{customdata[0]:,}<br>"
            "Candidate variants: %{customdata[1]:,}<br>Hit rate: %{y:.3%}<extra></extra>"
        )
    )
    compact_figure(fig_clin_rate)
    sections.append(fig_html(fig_clin_rate))

    clin_counts = variant_summary.clinvar_counts.copy()
    clin_plot = clin_counts[
        clin_counts["clinvar_category"].astype(str).isin(["P/LP", "B/LB", "VUS", "Other"])
    ].copy()
    totals = clin_plot.groupby("strategy", observed=True)["Variant_Count"].transform("sum")
    clin_plot["Fraction"] = clin_plot["Variant_Count"] / totals.replace(0, np.nan)
    clin_order = (
        clin_plot[clin_plot["clinvar_category"].astype(str) == "B/LB"]
        .sort_values("Fraction", ascending=False)
        ["strategy"]
        .tolist()
    )
    clin_order += [strategy for strategy in clinvar_rate["Strategy"].tolist() if strategy not in clin_order]
    fig_clin = px.bar(
        clin_plot,
        x="strategy",
        y="Fraction",
        color="clinvar_category",
        barmode="stack",
        title="ClinVar classification mix among classified variants",
        category_orders={"strategy": clin_order, "clinvar_category": CLINVAR_ORDER},
        color_discrete_map=CLINVAR_COLORS,
        labels={"strategy": "", "Fraction": "ClinVar class fraction", "clinvar_category": "ClinVar class"},
    )
    fig_clin.update_layout(yaxis_tickformat=".0%")
    compact_figure(fig_clin, height=360)
    sections.append(fig_html(fig_clin))

    star_counts = variant_summary.pathogenic_star_counts.copy()
    if not star_counts.empty:
        present_stars = [star for star in REVIEW_STAR_ORDER if star != "Unmapped"]
        if "Unmapped" in set(star_counts["Review stars"].astype(str)):
            present_stars.append("Unmapped")
        complete_index = pd.MultiIndex.from_product(
            [strategy_stats["Strategy"].tolist(), present_stars], names=["Strategy", "Review stars"]
        )
        star_counts = (
            star_counts.set_index(["Strategy", "Review stars"])
            .reindex(complete_index, fill_value=0)
            .reset_index()
        )
        totals = star_counts.groupby("Strategy", observed=True)["Variant_Count"].sum()
        high_conf = star_counts[star_counts["Review stars"].astype(str).isin(["4", "3", "2"])]
        high_conf_totals = high_conf.groupby("Strategy", observed=True)["Variant_Count"].sum()
        star_order = (
            pd.DataFrame({"total": totals, "high_conf": high_conf_totals})
            .fillna(0)
            .sort_values(["high_conf", "total"], ascending=False)
            .index.tolist()
        )
        fig_stars = px.bar(
            star_counts,
            x="Strategy",
            y="Variant_Count",
            color="Review stars",
            barmode="stack",
            title="Pathogenic ClinVar hits by review stars",
            category_orders={"Strategy": star_order, "Review stars": present_stars},
            color_discrete_map=REVIEW_STAR_COLORS,
            labels={"Strategy": "", "Variant_Count": "P/LP ClinVar variants", "Review stars": "Review stars"},
        )
        compact_figure(fig_stars, height=340)
        sections.append("<h3>Pathogenic ClinVar Evidence</h3>")
        sections.append(fig_html(fig_stars))
    else:
        sections.append("<h3>Pathogenic ClinVar Evidence</h3>")
        sections.append("<p>No P/LP ClinVar variants were found in the candidate set.</p>")

    pathogenic_consequence_counts = group_consequence_counts(variant_summary.pathogenic_consequence_counts)
    if not pathogenic_consequence_counts.empty:
        pathogenic_order = (
            pathogenic_consequence_counts.groupby("Strategy", observed=True)["Variant_Count"]
            .sum()
            .sort_values(ascending=False)
            .index.tolist()
        )
        fig_path_conseq = px.bar(
            pathogenic_consequence_counts,
            x="Strategy",
            y="Variant_Count",
            color="Consequence group",
            barmode="stack",
            title=f"{variant_summary.consequence_source} consequence groups for pathogenic ClinVar hits",
            category_orders={"Strategy": pathogenic_order, "Consequence group": CONSEQUENCE_GROUP_ORDER},
            color_discrete_map=CONSEQUENCE_GROUP_COLORS,
            labels={"Strategy": "", "Variant_Count": "P/LP ClinVar variants", "Consequence group": "Consequence group"},
        )
        compact_figure(fig_path_conseq, height=320)
        sections.append(fig_html(fig_path_conseq))

    return sections


def build_feature_sections(cov: pd.DataFrame) -> list[str]:
    sections = ["<h2>Alignment Coverage</h2>"]
    if cov.empty:
        sections.append("<p>No feature coverage table was found.</p>")
        return sections

    profile_summary = coverage_summary(cov, PROFILE_FEATURE_ORDER)
    if profile_summary.empty:
        sections.append("<p>No CDS/UTR/intron coverage rows were found.</p>")
        return sections

    cds_depth = (
        profile_summary[profile_summary["feature_type"].astype(str) == "cds"]
        .sort_values("Mean_Depth_Weighted", ascending=False)
    )
    strategy_order = cds_depth["strategy"].tolist() or sorted(profile_summary["strategy"].unique())

    fig_breadth = px.bar(
        profile_summary,
        x="strategy",
        y="Breadth_Weighted",
        color="feature_type",
        barmode="group",
        title="Target bases covered by one or more orthologs",
        category_orders={"strategy": strategy_order, "feature_type": PROFILE_FEATURE_ORDER},
        labels={"strategy": "", "Breadth_Weighted": "Target bases covered", "feature_type": "Feature"},
    )
    fig_breadth.update_layout(yaxis_tickformat=".0%")
    compact_figure(fig_breadth, height=360)
    sections.append(fig_html(fig_breadth))

    fig_depth = px.bar(
        profile_summary,
        x="strategy",
        y="Mean_Depth_Weighted",
        color="feature_type",
        barmode="group",
        title="Weighted mean ortholog depth by target feature",
        category_orders={"strategy": strategy_order, "feature_type": PROFILE_FEATURE_ORDER},
        labels={
            "strategy": "",
            "Mean_Depth_Weighted": "Weighted mean ortholog depth",
            "feature_type": "Feature",
        },
    )
    compact_figure(fig_depth, height=360)
    sections.append(fig_html(fig_depth))
    return sections


def gnomad_stratification_figure(
    counts: pd.DataFrame,
    category_column: str,
    category_order: list[str],
    strategy_order: list[str],
    title: str,
    color_map: dict[str, str] | None = None,
):
    if counts.empty:
        return None
    plot = counts.copy()
    plot["strategy"] = plot["strategy"].astype(str)
    plot["gnomad_status"] = plot["gnomad_status"].astype(str)
    plot[category_column] = plot[category_column].astype(str)
    present_categories = set(plot[category_column])
    categories = [category for category in category_order if category in present_categories]
    categories += sorted(present_categories - set(categories))
    statuses = [("found", "Found"), ("not_found", "Not found")]
    combinations = [
        (strategy, status, label)
        for strategy in strategy_order
        for status, label in statuses
    ]
    indexed = plot.set_index(["strategy", "gnomad_status", category_column])["Variant_Count"]
    totals = plot.groupby(["strategy", "gnomad_status"], observed=True)["Variant_Count"].sum()
    x_strategy = [strategy for strategy, _status, _label in combinations]
    x_status = [label for _strategy, _status, label in combinations]
    fig = go.Figure()
    for category in categories:
        counts_for_category = np.asarray(
            [int(indexed.get((strategy, status, category), 0)) for strategy, status, _label in combinations]
        )
        totals_for_group = np.asarray(
            [int(totals.get((strategy, status), 0)) for strategy, status, _label in combinations]
        )
        fractions = np.divide(
            counts_for_category,
            totals_for_group,
            out=np.zeros(len(combinations), dtype=float),
            where=totals_for_group > 0,
        )
        fig.add_trace(
            go.Bar(
                x=[x_strategy, x_status],
                y=fractions,
                name=category,
                marker_color=(color_map or {}).get(category),
                customdata=np.column_stack(
                    [x_strategy, x_status, counts_for_category, totals_for_group]
                ),
                hovertemplate=(
                    "%{customdata[0]}<br>%{customdata[1]} in gnomAD<br>"
                    + category + ": %{customdata[2]:,} / %{customdata[3]:,} "
                    "(%{y:.1%})<extra></extra>"
                ),
            )
        )
    fig.update_layout(title=title, barmode="stack", height=440, bargap=0.18)
    fig.update_yaxes(tickformat=".0%", title_text="Within-stratum fraction", range=[0, 1])
    compact_figure(fig, height=440)
    return fig


def build_gnomad_stratification_sections(
    variant_summary: VariantSummary,
    strategy_stats: pd.DataFrame,
    candidate_conservation: CandidateConservation,
) -> list[str]:
    sections = [
        "<h2>gnomAD Stratification</h2>",
        "<p class=\"lead\">Descriptive comparison of candidate alleles found and not found in gnomAD. "
        "This is not a matched-control analysis.</p>",
    ]
    strategy_order = sort_by_metric(
        strategy_stats[["Strategy", "gnomAD found %"]], "gnomAD found %"
    )["Strategy"].tolist()

    event_counts = variant_summary.gnomad_event_counts.copy()
    event_order = [
        item for item in ["snv", "ins", "del", "mnv", "complex"]
        if item in set(event_counts.get("event_type", pd.Series(dtype=str)).astype(str))
    ]
    event_order += sorted(set(event_counts.get("event_type", pd.Series(dtype=str)).astype(str)) - set(event_order))
    event_fig = gnomad_stratification_figure(
        event_counts,
        "event_type",
        event_order,
        strategy_order,
        "Variant type: gnomAD hits versus non-hits",
    )
    if event_fig is not None:
        sections.append(fig_html(event_fig))

    context_counts = variant_summary.gnomad_context_counts.copy()
    context_counts["Target context"] = context_counts.get("target_context", pd.Series(dtype=str)).map(
        TARGET_CONTEXT_LABELS
    ).fillna("Other")
    context_fig = gnomad_stratification_figure(
        context_counts,
        "Target context",
        [TARGET_CONTEXT_LABELS[item] for item in TARGET_CONTEXT_ORDER],
        strategy_order,
        "Target context: gnomAD hits versus non-hits",
        TARGET_CONTEXT_COLORS,
    )
    if context_fig is not None:
        sections.append(fig_html(context_fig))

    sections.append("<h3>Conservation</h3>")
    phylop_fig = candidate_phylop_figure(candidate_conservation, strategy_order)
    if phylop_fig is not None:
        sections.append(
            "<p class=\"lead\">Candidate phyloP100way distributions are shown separately for exact gnomAD hits "
            "and alleles without a gnomAD hit. Select one strategy to compare the two strata.</p>"
        )
        sections.append(fig_html(phylop_fig))
        phylop_summary_fig = candidate_phylop_summary_figure(candidate_conservation, strategy_order)
        if phylop_summary_fig is not None:
            sections.append(fig_html(phylop_summary_fig))
    else:
        sections.append("<p>No candidate-wide phyloP100way scores were available.</p>")

    sections.extend(
        [
            "<h3>Functional Consequence</h3>",
            "<p class=\"analysis-note\">Not computed. A defensible comparison requires the same VEP release, "
            "transcript set, and consequence-selection rule for both gnomAD strata.</p>",
        ]
    )
    return sections
