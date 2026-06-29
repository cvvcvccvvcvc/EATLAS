#!/usr/bin/env python3
"""Build an interactive HTML report for alignment strategy analytics."""

from __future__ import annotations

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import plotly.express as px


warnings.filterwarnings("ignore", r"All-NaN (slice|axis) encountered")
warnings.filterwarnings("ignore", r"Mean of empty slice")


FEATURE_ORDER = ["gene", "exon", "cds", "utr", "intron"]
DISJOINT_FEATURE_ORDER = ["cds", "utr", "intron"]
CLINVAR_ORDER = ["P/LP", "B/LB", "VUS", "Other", "Not Found"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events-tsv", required=True, help="Annotated or raw alignment events TSV")
    parser.add_argument("--feature-coverage-tsv", help="Feature coverage TSV produced by alignment merge")
    parser.add_argument("--out-html", required=True, help="Path to output HTML report")
    return parser.parse_args()


def read_events(path: str) -> pd.DataFrame:
    print(f"Reading {path}...")
    df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    for col in ["clinvar_sig", "clinvar_revstat", "clinvar_id", "gnomad_af", "gnomad_af_source", "gnomad_csq"]:
        if col not in df.columns:
            df[col] = ""
    df["gnomad_af"] = pd.to_numeric(df["gnomad_af"], errors="coerce")
    df["variant_id"] = (
        df["genomic_accession"].astype(str)
        + ":"
        + df["genomic_start1"].astype(str)
        + ":"
        + df["ref"].fillna("").astype(str)
        + ">"
        + df["alt"].fillna("").astype(str)
    )
    df["clinvar_category"] = df["clinvar_sig"].apply(categorize_clinvar)
    return df


def categorize_clinvar(value) -> str:
    if pd.isna(value) or value == "":
        return "Not Found"
    text = str(value).lower()
    if "conflicting" in text:
        return "Other"
    if "uncertain" in text or "vus" in text:
        return "VUS"
    if "pathogenic" in text:
        return "P/LP"
    if "benign" in text:
        return "B/LB"
    return "Other"


def calc_titv(df: pd.DataFrame) -> float:
    snvs = df[df["event_type"] == "snv"]
    if snvs.empty:
        return np.nan
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    ti = 0
    tv = 0
    for _, row in snvs.iterrows():
        ref = str(row["ref"])
        alt = str(row["alt"])
        if len(ref) != 1 or len(alt) != 1:
            continue
        if (ref, alt) in transitions:
            ti += 1
        else:
            tv += 1
    if tv == 0:
        return np.nan if ti == 0 else float("inf")
    return round(ti / tv, 3)


def unique_variants(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop_duplicates(["strategy", "variant_id"]).copy()


def fig_html(fig, include_plotlyjs: bool = False) -> str:
    return fig.to_html(full_html=False, include_plotlyjs="cdn" if include_plotlyjs else False)


def table_html(df: pd.DataFrame, classes: str = "table table-striped table-bordered", max_rows: int | None = None) -> str:
    shown = df if max_rows is None else df.head(max_rows)
    return shown.to_html(index=False, classes=classes, float_format="%.5g")


def metric_cards(items: list[tuple[str, object]]) -> str:
    cards = []
    for label, value in items:
        cards.append(
            f"""
            <div class="metric-card">
                <div class="metric-label">{label}</div>
                <div class="metric-value">{value}</div>
            </div>
            """
        )
    return f"<div class=\"metric-grid\">{''.join(cards)}</div>"


def strategy_variant_table(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    rows = []
    strategy_variant_sets: dict[str, set[str]] = {}
    for strategy in sorted(df["strategy"].dropna().unique()):
        s_df = df[df["strategy"] == strategy]
        unique_vars = s_df.drop_duplicates("variant_id").copy()
        support_counts = s_df.groupby("variant_id")["ortholog_gene_id"].nunique()
        unique_vars["ortholog_count"] = unique_vars["variant_id"].map(support_counts)
        strategy_variant_sets[strategy] = set(unique_vars["variant_id"])
        clinvar_counts = unique_vars["clinvar_category"].value_counts()
        rows.append(
            {
                "Strategy": strategy,
                "Unique Variants": len(unique_vars),
                "Ti/Tv": calc_titv(unique_vars),
                "ClinVar Found": int((unique_vars["clinvar_category"] != "Not Found").sum()),
                "gnomAD Found": int(unique_vars["gnomad_af"].notna().sum()),
                "P/LP": int(clinvar_counts.get("P/LP", 0)),
                "B/LB": int(clinvar_counts.get("B/LB", 0)),
                "VUS": int(clinvar_counts.get("VUS", 0)),
                "Other ClinVar": int(clinvar_counts.get("Other", 0)),
                "Median gnomAD AF": unique_vars["gnomad_af"].median(),
            }
        )
    return pd.DataFrame(rows), strategy_variant_sets


def unique_contribution_table(df: pd.DataFrame, strategy_variant_sets: dict[str, set[str]]) -> pd.DataFrame:
    rows = []
    for strategy in sorted(strategy_variant_sets):
        union_others: set[str] = set()
        for other_strategy, variants in strategy_variant_sets.items():
            if other_strategy != strategy:
                union_others.update(variants)
        unique_to_strategy = strategy_variant_sets[strategy] - union_others
        s_df = df[(df["strategy"] == strategy) & (df["variant_id"].isin(unique_to_strategy))].drop_duplicates(
            "variant_id"
        )
        clinvar_counts = s_df["clinvar_category"].value_counts()
        rows.append(
            {
                "Strategy": strategy,
                "Unique To Strategy": len(unique_to_strategy),
                "Ti/Tv": calc_titv(s_df),
                "ClinVar Found": int((s_df["clinvar_category"] != "Not Found").sum()),
                "gnomAD Found": int(s_df["gnomad_af"].notna().sum()),
                "P/LP": int(clinvar_counts.get("P/LP", 0)),
                "B/LB": int(clinvar_counts.get("B/LB", 0)),
                "VUS": int(clinvar_counts.get("VUS", 0)),
                "Other ClinVar": int(clinvar_counts.get("Other", 0)),
                "Median gnomAD AF": s_df["gnomad_af"].median(),
            }
        )
    return pd.DataFrame(rows)


def support_breakdown_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for strategy in sorted(df["strategy"].dropna().unique()):
        s_df = df[df["strategy"] == strategy]
        unique_vars = s_df.drop_duplicates("variant_id").copy()
        support_counts = s_df.groupby("variant_id")["ortholog_gene_id"].nunique()
        unique_vars["ortholog_count"] = unique_vars["variant_id"].map(support_counts)
        unique_vars["Support Bucket"] = unique_vars["ortholog_count"].apply(lambda x: str(x) if x < 5 else "5+")
        for bucket in sorted(unique_vars["Support Bucket"].unique(), key=lambda x: int(x.replace("+", ""))):
            b_df = unique_vars[unique_vars["Support Bucket"] == bucket]
            clinvar_counts = b_df["clinvar_category"].value_counts()
            rows.append(
                {
                    "Strategy": strategy,
                    "Ortholog Support": bucket,
                    "Variant Count": len(b_df),
                    "Ti/Tv": calc_titv(b_df),
                    "ClinVar Found": int((b_df["clinvar_category"] != "Not Found").sum()),
                    "gnomAD Found": int(b_df["gnomad_af"].notna().sum()),
                    "P/LP": int(clinvar_counts.get("P/LP", 0)),
                    "B/LB": int(clinvar_counts.get("B/LB", 0)),
                    "VUS": int(clinvar_counts.get("VUS", 0)),
                    "Other ClinVar": int(clinvar_counts.get("Other", 0)),
                    "Median gnomAD AF": b_df["gnomad_af"].median(),
                }
            )
    return pd.DataFrame(rows)


def build_overview(df: pd.DataFrame, cov: pd.DataFrame, strategy_stats: pd.DataFrame) -> list[str]:
    unique_all = df.drop_duplicates("variant_id")
    cards = [
        ("Event Rows", f"{len(df):,}"),
        ("Unique Variants", f"{len(unique_all):,}"),
        ("Strategies", f"{df['strategy'].nunique():,}"),
        ("Genes", f"{df['gene_id'].nunique():,}"),
        ("ClinVar Variants", f"{int((unique_all['clinvar_category'] != 'Not Found').sum()):,}"),
        ("gnomAD Variants", f"{int(unique_all['gnomad_af'].notna().sum()):,}"),
    ]
    if not cov.empty:
        cards.append(("Feature Coverage Rows", f"{len(cov):,}"))
    sections = [metric_cards(cards)]
    sections.append("<h2>Strategy Summary</h2>")
    sections.append(table_html(strategy_stats))
    return sections


def build_variant_sections(
    df: pd.DataFrame,
    strategy_variant_sets: dict[str, set[str]],
    include_plotly: bool,
) -> list[str]:
    sections = ["<h2>Variant Evidence</h2>"]
    support = support_breakdown_table(df)
    unique_contrib = unique_contribution_table(df, strategy_variant_sets)
    sections.append("<h3>Ortholog Support Buckets</h3>")
    sections.append(table_html(support, classes="table table-sm table-striped"))
    sections.append("<h3>Unique Contributions</h3>")
    sections.append(table_html(unique_contrib))

    plot_df = df.drop_duplicates(["strategy", "variant_id"])
    fig_titv = px.bar(
        strategy_variant_table(df)[0],
        x="Strategy",
        y="Ti/Tv",
        title="Ti/Tv by Strategy",
    )
    sections.append(fig_html(fig_titv, include_plotlyjs=include_plotly))

    event_type_counts = (
        plot_df.groupby(["strategy", "event_type"], as_index=False)
        .agg(Variant_Count=("variant_id", "count"))
        .sort_values(["strategy", "event_type"])
    )
    fig_events = px.bar(
        event_type_counts,
        x="strategy",
        y="Variant_Count",
        color="event_type",
        barmode="group",
        title="Unique Variant Event Types by Strategy",
    )
    sections.append(fig_html(fig_events))
    return sections


def build_clinvar_gnomad_sections(df: pd.DataFrame, include_plotly: bool) -> list[str]:
    sections = ["<h2>ClinVar and gnomAD</h2>"]
    unique_all = df.drop_duplicates("variant_id")
    clinvar_found = unique_all[unique_all["clinvar_category"] != "Not Found"]
    gnomad_found = unique_all[unique_all["gnomad_af"].notna()]
    sections.append(
        metric_cards(
            [
                ("ClinVar Unique Variants", f"{len(clinvar_found):,}"),
                ("gnomAD Unique Variants", f"{len(gnomad_found):,}"),
                ("ClinVar Event Rows", f"{int((df['clinvar_category'] != 'Not Found').sum()):,}"),
                ("gnomAD Event Rows", f"{int(df['gnomad_af'].notna().sum()):,}"),
            ]
        )
    )

    plot_df = df.drop_duplicates(["strategy", "variant_id"])
    clin_counts = (
        plot_df.groupby(["strategy", "clinvar_category"], as_index=False)
        .agg(Variant_Count=("variant_id", "count"))
    )
    clin_counts["clinvar_category"] = pd.Categorical(
        clin_counts["clinvar_category"], categories=CLINVAR_ORDER, ordered=True
    )
    clin_counts = clin_counts.sort_values(["strategy", "clinvar_category"])
    sections.append("<h3>ClinVar Categories</h3>")
    sections.append(table_html(clin_counts))
    fig_clin = px.bar(
        clin_counts[clin_counts["clinvar_category"] != "Not Found"],
        x="strategy",
        y="Variant_Count",
        color="clinvar_category",
        barmode="group",
        title="ClinVar Variants by Strategy",
        category_orders={"clinvar_category": CLINVAR_ORDER},
    )
    sections.append(fig_html(fig_clin, include_plotlyjs=include_plotly))

    gnomad = plot_df[plot_df["gnomad_af"].notna() & (plot_df["gnomad_af"] > 0)].copy()
    if not gnomad.empty:
        gnomad["log10_gnomad_af"] = np.log10(gnomad["gnomad_af"])
        fig_af = px.histogram(
            gnomad,
            x="log10_gnomad_af",
            color="strategy",
            nbins=80,
            barmode="overlay",
            opacity=0.55,
            histnorm="probability density",
            title="gnomAD AF Distribution by Strategy",
        )
        fig_af.update_layout(yaxis_title="Density", xaxis_title="log10 gnomAD AF")
        sections.append("<h3>gnomAD AF Distribution</h3>")
        sections.append(fig_html(fig_af))
    else:
        sections.append("<p>No non-zero gnomAD AF values were found.</p>")

    gnomad_counts = (
        plot_df.assign(gnomad_found=plot_df["gnomad_af"].notna())
        .groupby(["strategy", "gnomad_found"], as_index=False)
        .agg(Variant_Count=("variant_id", "count"))
    )
    fig_gnomad_found = px.bar(
        gnomad_counts,
        x="strategy",
        y="Variant_Count",
        color="gnomad_found",
        barmode="group",
        title="Variants Found in gnomAD by Strategy",
    )
    sections.append(fig_html(fig_gnomad_found))
    return sections


def read_feature_coverage(path: str | None) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    print(f"Reading {path}...")
    cov = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    numeric_cols = [
        "length_bp",
        "ortholog_count",
        "orthologs_covered",
        "covered_bases",
        "coverage_breadth",
        "depth_bases",
        "mean_depth",
    ]
    for col in numeric_cols:
        cov[col] = pd.to_numeric(cov[col], errors="coerce")
    return cov


def coverage_summary(cov: pd.DataFrame, feature_types: list[str] | None = None) -> pd.DataFrame:
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
    return summary.sort_values(["strategy", "feature_type"])


def build_feature_sections(cov: pd.DataFrame, include_plotly: bool) -> list[str]:
    sections = ["<h2>Target Feature Coverage</h2>"]
    if cov.empty:
        sections.append("<p>No feature coverage table was provided.</p>")
        return sections

    summary = coverage_summary(cov)
    disjoint_summary = coverage_summary(cov, DISJOINT_FEATURE_ORDER)

    sections.append("<h3>Coverage Summary</h3>")
    sections.append(table_html(summary))

    fig_strategy_breadth = px.bar(
        disjoint_summary,
        x="feature_type",
        y="Breadth_Weighted",
        color="strategy",
        barmode="group",
        title="Weighted Coverage Breadth: Toggle Strategies",
        category_orders={"feature_type": DISJOINT_FEATURE_ORDER},
    )
    sections.append(fig_html(fig_strategy_breadth, include_plotlyjs=include_plotly))

    fig_feature_breadth = px.bar(
        disjoint_summary,
        x="strategy",
        y="Breadth_Weighted",
        color="feature_type",
        barmode="group",
        title="Weighted Coverage Breadth: Toggle Feature Types",
        category_orders={"feature_type": DISJOINT_FEATURE_ORDER},
    )
    sections.append(fig_html(fig_feature_breadth))

    fig_strategy_depth = px.bar(
        disjoint_summary,
        x="feature_type",
        y="Mean_Depth_Weighted",
        color="strategy",
        barmode="group",
        title="Weighted Mean Ortholog Depth: Toggle Strategies",
        category_orders={"feature_type": DISJOINT_FEATURE_ORDER},
    )
    sections.append(fig_html(fig_strategy_depth))

    fig_feature_depth = px.bar(
        disjoint_summary,
        x="strategy",
        y="Mean_Depth_Weighted",
        color="feature_type",
        barmode="group",
        title="Weighted Mean Ortholog Depth: Toggle Feature Types",
        category_orders={"feature_type": DISJOINT_FEATURE_ORDER},
    )
    sections.append(fig_html(fig_feature_depth))

    low_coverage = cov[cov["feature_type"].isin(DISJOINT_FEATURE_ORDER)].copy()
    low_coverage = low_coverage.sort_values(["coverage_breadth", "mean_depth"]).head(50)
    low_coverage = low_coverage[
        [
            "gene_id",
            "strategy",
            "feature_type",
            "feature_id",
            "length_bp",
            "coverage_breadth",
            "mean_depth",
            "orthologs_covered",
            "ortholog_count",
        ]
    ]
    sections.append("<h3>Lowest-Coverage Feature Examples</h3>")
    sections.append(table_html(low_coverage, classes="table table-sm table-striped"))
    return sections


def build_methods_sections(args: argparse.Namespace, df: pd.DataFrame, cov: pd.DataFrame) -> list[str]:
    rows = [
        {"Key": "Events TSV", "Value": args.events_tsv},
        {"Key": "Feature Coverage TSV", "Value": args.feature_coverage_tsv or ""},
        {"Key": "Event Rows Loaded", "Value": f"{len(df):,}"},
        {"Key": "Feature Coverage Rows Loaded", "Value": str(len(cov)) if not cov.empty else "0"},
    ]
    return [
        "<h2>Files and Methods</h2>",
        table_html(pd.DataFrame(rows), classes="table table-sm table-striped"),
        """
        <p>
        Feature coverage breadth is the fraction of bases in a structural interval covered by at least one
        ortholog alignment segment. Mean depth is computed after merging overlapping segments within each
        ortholog, so one ortholog cannot inflate depth by overlapping itself.
        </p>
        """,
    ]


def render_tabs(sections: list[tuple[str, str, list[str]]]) -> str:
    buttons = []
    pages = []
    for index, (tab_id, title, html_parts) in enumerate(sections):
        active = " active" if index == 0 else ""
        buttons.append(f'<button class="tab-button{active}" data-tab="{tab_id}">{title}</button>')
        pages.append(f'<section id="tab-{tab_id}" class="tab-page{active}">{"".join(html_parts)}</section>')
    return f"""
    <div class="tab-bar">{''.join(buttons)}</div>
    {''.join(pages)}
    <script>
    document.querySelectorAll('.tab-button').forEach(button => {{
        button.addEventListener('click', () => {{
            const tab = button.dataset.tab;
            document.querySelectorAll('.tab-button').forEach(item => item.classList.remove('active'));
            document.querySelectorAll('.tab-page').forEach(item => item.classList.remove('active'));
            button.classList.add('active');
            document.getElementById('tab-' + tab).classList.add('active');
            window.dispatchEvent(new Event('resize'));
        }});
    }});
    </script>
    """


def main() -> None:
    args = parse_args()
    df = read_events(args.events_tsv)
    cov = read_feature_coverage(args.feature_coverage_tsv)

    print("Computing variant metrics...")
    strategy_stats, strategy_variant_sets = strategy_variant_table(df)

    sections = [
        ("overview", "Overview", build_overview(df, cov, strategy_stats)),
        ("variants", "Variants", build_variant_sections(df, strategy_variant_sets, include_plotly=True)),
        ("clinvar-gnomad", "ClinVar & gnomAD", build_clinvar_gnomad_sections(df, include_plotly=True)),
        ("coverage", "Feature Coverage", build_feature_sections(cov, include_plotly=True)),
        ("methods", "Files", build_methods_sections(args, df, cov)),
    ]

    print(f"Writing report to {args.out_html}...")
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>Alignment Strategies Report</title>
        <style>
            body {{
                padding: 20px;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
                color: #1f2933;
            }}
            h1 {{ margin-bottom: 4px; }}
            h2 {{ margin-top: 28px; border-bottom: 1px solid #d5d9df; padding-bottom: 8px; }}
            h3 {{ margin-top: 22px; }}
            .lead {{ margin-top: 0; color: #52606d; }}
            .tab-bar {{
                display: flex;
                flex-wrap: wrap;
                gap: 8px;
                margin: 20px 0;
                border-bottom: 1px solid #d5d9df;
            }}
            .tab-button {{
                border: 1px solid #cbd2d9;
                border-bottom: none;
                background: #f5f7fa;
                color: #1f2933;
                padding: 8px 12px;
                cursor: pointer;
                border-radius: 6px 6px 0 0;
                font-size: 14px;
            }}
            .tab-button.active {{
                background: white;
                font-weight: 600;
            }}
            .tab-page {{ display: none; }}
            .tab-page.active {{ display: block; }}
            .metric-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
                gap: 12px;
                margin: 16px 0 24px 0;
            }}
            .metric-card {{
                border: 1px solid #d5d9df;
                border-radius: 6px;
                padding: 12px;
                background: #fff;
            }}
            .metric-label {{ color: #52606d; font-size: 13px; }}
            .metric-value {{ font-size: 24px; font-weight: 650; margin-top: 4px; }}
            table {{
                border-collapse: collapse;
                width: auto;
                max-width: 100%;
                margin-bottom: 18px;
                font-size: 13px;
            }}
            th, td {{ border: 1px solid #d5d9df; padding: 6px 8px; text-align: right; }}
            th {{ background: #f5f7fa; }}
            td:first-child, th:first-child {{ text-align: left; }}
            .plotly-graph-div {{ min-height: 420px; }}
        </style>
    </head>
    <body>
        <h1>Alignment Strategies Report</h1>
        <p class="lead">Variant, annotation, and target-feature coverage analytics.</p>
        {render_tabs(sections)}
    </body>
    </html>
    """
    with open(args.out_html, "w") as handle:
        handle.write(html)
    print("Done!")


if __name__ == "__main__":
    main()
