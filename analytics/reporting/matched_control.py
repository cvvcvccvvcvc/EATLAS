"""Matched-control visualizations and quality-control details."""

from __future__ import annotations

import math

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from analytics.analyses.matched_control import TargetSpaceNullAnalysis
from .components import compact_figure, fig_html, format_int, metric_cards, strategy_label, table_html
from .config import CLINVAR_COLORS

def target_null_interval_figure(
    frame: pd.DataFrame,
    *,
    title: str,
    xaxis_title: str,
    observed_name: str,
    null_name: str,
    tickformat: str | None = None,
    log_x: bool = False,
) -> go.Figure:
    plot = frame.copy()
    plot["Strategy"] = plot["strategy"].map(strategy_label)
    plot["difference"] = plot["observed_value"] - plot["null_value"]
    plot = plot.sort_values("difference", ascending=False, kind="mergesort")
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=plot["null_value"],
            y=plot["Strategy"],
            mode="markers",
            marker={"color": "#8c8c8c", "size": 8},
            error_x={
                "type": "data",
                "symmetric": False,
                "array": plot["null_ci_high"] - plot["null_value"],
                "arrayminus": plot["null_value"] - plot["null_ci_low"],
                "color": "#8c8c8c",
            },
            name=null_name,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=plot["observed_value"],
            y=plot["Strategy"],
            mode="markers",
            marker={"color": "#2166ac", "size": 10, "symbol": "diamond"},
            error_x={
                "type": "data",
                "symmetric": False,
                "array": plot["observed_ci_high"] - plot["observed_value"],
                "arrayminus": plot["observed_value"] - plot["observed_ci_low"],
                "color": "#2166ac",
            },
            name=observed_name,
        )
    )
    xaxis: dict[str, object] = {"title": xaxis_title}
    if tickformat:
        xaxis["tickformat"] = tickformat
    if log_x:
        xaxis["type"] = "log"
    fig.update_layout(
        title=title,
        xaxis=xaxis,
        yaxis={"categoryorder": "array", "categoryarray": plot["Strategy"].tolist()[::-1]},
    )
    compact_figure(fig, height=max(360, 52 * len(plot) + 120), show_x_title=True)
    return fig


def clinvar_class_null_figure(frame: pd.DataFrame) -> go.Figure:
    class_order = ["B/LB", "P/LP", "VUS", "Other"]
    colors = [CLINVAR_COLORS[category] for category in class_order]
    benign = frame[frame["clinvar_class"].eq("B/LB")].copy()
    benign["difference"] = benign["observed_value"] - benign["null_value"]
    strategies = benign.sort_values("difference", ascending=False, kind="mergesort")["strategy"].tolist()
    if not strategies:
        strategies = sorted(frame["strategy"].astype(str).unique())
    columns = 2 if len(strategies) > 1 else 1
    rows = math.ceil(len(strategies) / columns)
    fig = make_subplots(
        rows=rows,
        cols=columns,
        subplot_titles=[strategy_label(strategy) for strategy in strategies],
        shared_yaxes=True,
        vertical_spacing=min(0.16, 0.35 / max(rows, 1)),
    )
    for index, strategy in enumerate(strategies):
        row = index // columns + 1
        column = index % columns + 1
        values = frame[frame["strategy"].eq(strategy)].set_index("clinvar_class").reindex(class_order)
        fig.add_trace(
            go.Bar(
                x=class_order,
                y=values["observed_value"],
                marker_color=colors,
                error_y={
                    "type": "data",
                    "symmetric": False,
                    "array": values["observed_ci_high"] - values["observed_value"],
                    "arrayminus": values["observed_value"] - values["observed_ci_low"],
                    "color": "#2166ac",
                },
                name="GAPH (95% paired bootstrap interval)",
                showlegend=index == 0,
                hovertemplate="%{x}<br>GAPH: %{y:.1%}<extra></extra>",
            ),
            row=row,
            col=column,
        )
        fig.add_trace(
            go.Scatter(
                x=class_order,
                y=values["null_value"],
                mode="markers",
                marker={"color": "#4d4d4d", "size": 8},
                error_y={
                    "type": "data",
                    "symmetric": False,
                    "array": values["null_ci_high"] - values["null_value"],
                    "arrayminus": values["null_value"] - values["null_ci_low"],
                    "color": "#4d4d4d",
                },
                name="Matched control (95% paired bootstrap interval)",
                showlegend=index == 0,
                hovertemplate="%{x}<br>Matched null: %{y:.1%}<extra></extra>",
            ),
            row=row,
            col=column,
        )
    fig.update_yaxes(tickformat=".0%", rangemode="tozero", title_text="Classified ClinVar hits")
    fig.update_layout(
        title="ClinVar class composition",
        barmode="group",
        height=max(390, 265 * rows),
        margin={"l": 65, "r": 25, "t": 75, "b": 90},
        legend={"orientation": "h", "x": 0.5, "xanchor": "center", "y": -0.14, "yanchor": "top"},
    )
    return fig


def build_target_space_null_sections(
    analysis: TargetSpaceNullAnalysis | None,
    *,
    enabled: bool = True,
) -> list[str]:
    sections = ["<h2>Matched Control</h2>"]
    if not enabled:
        sections.append(
            "<p>Matched Control was disabled for this report run. Enable it with "
            "<code>--target-space-null</code>; this analysis uses Ensembl VEP and may take hours.</p>"
        )
        return sections
    if analysis is None:
        sections.append("<p>No target-space-null analysis is available for this run.</p>")
        return sections
    summary = analysis.summary.copy()
    if summary.empty:
        sections.append("<p>No consequence-matched target-space controls could be constructed.</p>")
        return sections

    sampled = analysis.manifest.get("sampled_focal_count", 0)
    matched = analysis.manifest.get("matched_focal_count", 0)
    sections.append(
        metric_cards(
            [
                (
                    "Sample cap per strategy (input)",
                    format_int(analysis.manifest.get("inputs", {}).get("sample_size_per_strategy", 0)),
                ),
                ("Sampled / matched focal SNVs", f"{format_int(sampled)} / {format_int(matched)}"),
                ("Matched-set bootstrap resamples (input)", format_int(analysis.resamples)),
            ]
        )
    )
    conservation_status = analysis.manifest.get("conservation", {}).get("status", "")
    if conservation_status != "complete":
        error = analysis.manifest.get("conservation", {}).get("error", "")
        sections.append(f"<p>phyloP annotation was incomplete: {error or conservation_status}</p>")

    sections.append("<h3>Conservation</h3>")
    ecdf = analysis.ecdf.copy()
    if not ecdf.empty:
        ecdf["Strategy"] = ecdf["strategy"].map(strategy_label)
        ecdf["set"] = ecdf["set"].replace({"Matched target-space null": "Matched control"})
        fig_ecdf = px.line(
            ecdf,
            x="phyloP100way",
            y="fraction_leq",
            color="set",
            facet_col="Strategy",
            facet_col_wrap=2,
            title="phyloP100way distributions: GAPH and matched control",
            labels={"fraction_leq": "Cumulative fraction", "set": ""},
            color_discrete_map={"GAPH": "#2166ac", "Matched control": "#8c8c8c"},
        )
        fig_ecdf.for_each_annotation(lambda item: item.update(text=item.text.split("=")[-1]))
        fig_ecdf.update_yaxes(tickformat=".0%")
        compact_figure(fig_ecdf, height=max(420, 260 * math.ceil(ecdf["Strategy"].nunique() / 2)))
        sections.append(fig_html(fig_ecdf))

    plot = summary.copy()
    plot["Strategy"] = plot["strategy"].map(strategy_label)
    plot = plot.sort_values("median_difference", ascending=False, kind="mergesort")
    fig = go.Figure()
    for row in plot.itertuples(index=False):
        fig.add_trace(
            go.Scatter(
                x=[row.null_median, row.observed_median],
                y=[row.Strategy, row.Strategy],
                mode="lines",
                line={"color": "#c7c7c7", "width": 2},
                hoverinfo="skip",
                showlegend=False,
            )
        )
    fig.add_trace(
        go.Scatter(
            x=plot["null_median"],
            y=plot["Strategy"],
            mode="markers",
            marker={"color": "#8c8c8c", "size": 8},
            error_x={
                "type": "data",
                "symmetric": False,
                "array": plot["null_ci_high"] - plot["null_median"],
                "arrayminus": plot["null_median"] - plot["null_ci_low"],
                "color": "#8c8c8c",
            },
            name="Matched-control median (95% paired bootstrap interval)",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=plot["observed_median"],
            y=plot["Strategy"],
            mode="markers",
            marker={"color": "#2166ac", "size": 10, "symbol": "diamond"},
            error_x={
                "type": "data",
                "symmetric": False,
                "array": plot["observed_ci_high"] - plot["observed_median"],
                "arrayminus": plot["observed_median"] - plot["observed_ci_low"],
                "color": "#2166ac",
            },
            name="GAPH median (95% paired bootstrap interval)",
        )
    )
    fig.update_layout(
        title="phyloP100way median: GAPH vs matched control",
        xaxis_title="phyloP100way",
        yaxis={"categoryorder": "array", "categoryarray": plot["Strategy"].tolist()[::-1]},
    )
    compact_figure(fig, height=max(360, 52 * len(plot) + 120), show_x_title=True)
    sections.append(fig_html(fig))

    gnomad = analysis.gnomad_summary.copy()
    if not gnomad.empty:
        sections.append("<h3>gnomAD</h3>")
        found = gnomad[gnomad["metric"].eq("found_fraction")]
        if not found.empty:
            fig_gnomad_found = target_null_interval_figure(
                found,
                title="Exact alleles found in gnomAD",
                xaxis_title="Fraction found",
                observed_name="GAPH fraction (95% paired bootstrap interval)",
                null_name="Matched-control fraction (95% paired bootstrap interval)",
                tickformat=".0%",
            )
            sections.append(fig_html(fig_gnomad_found))
        af = gnomad[gnomad["metric"].eq("median_af")]
        if not af.empty and (af[["observed_value", "null_value"]].gt(0).any(axis=None)):
            fig_gnomad_af = target_null_interval_figure(
                af,
                title="gnomAD allele frequency among exact hits",
                xaxis_title="Median allele frequency (log scale)",
                observed_name="GAPH median AF (95% paired bootstrap interval)",
                null_name="Matched-control median AF (95% paired bootstrap interval)",
                log_x=True,
            )
            sections.append(fig_html(fig_gnomad_af))
        gnomad_manifest = analysis.manifest.get("external_evidence", {}).get("gnomad", {})
        if gnomad_manifest.get("failed_region_count", 0):
            sections.append(
                f"<p>gnomAD evidence is incomplete: {format_int(gnomad_manifest['failed_region_count'])} "
                "region request(s) failed. Failed regions were treated as missing, not as absence.</p>"
            )

    clinvar = analysis.clinvar_summary.copy()
    if not clinvar.empty:
        sections.append("<h3>ClinVar</h3>")
        fig_clinvar_found = target_null_interval_figure(
            clinvar,
            title="Exact alleles found in ClinVar",
            xaxis_title="Fraction found",
            observed_name="GAPH fraction (95% paired bootstrap interval)",
            null_name="Matched-control fraction (95% paired bootstrap interval)",
            tickformat=".0%",
        )
        sections.append(fig_html(fig_clinvar_found))
        if not analysis.clinvar_class_summary.empty:
            fig_clinvar_class = clinvar_class_null_figure(analysis.clinvar_class_summary)
            sections.append(fig_html(fig_clinvar_class))

    return sections


def build_target_space_null_qc_sections(analysis: TargetSpaceNullAnalysis) -> list[str]:
    sections = ["<details><summary>Matched-control QC</summary>"]
    focal_vep = analysis.manifest.get("focal_vep", {})
    sections.append(
        table_html(
            pd.DataFrame(
                [
                    {"Metric": "VEP backend", "Value": focal_vep.get("backend", "")},
                    {"Metric": "VEP release", "Value": focal_vep.get("release", "")},
                    {
                        "Metric": "VEP-annotated focal SNVs",
                        "Value": analysis.manifest.get("vep_annotated_focal_count", 0),
                    },
                ]
            ),
            classes="table table-sm table-striped",
        )
    )

    table = analysis.summary.rename(
        columns={
            "strategy": "Strategy",
            "matched_focals": "Matched SNVs",
            "observed_median": "GAPH median",
            "observed_ci_low": "GAPH median Q2.5",
            "observed_ci_high": "GAPH median Q97.5",
            "null_median": "Matched-control median",
            "null_ci_low": "Matched-control median Q2.5",
            "null_ci_high": "Matched-control median Q97.5",
            "median_difference": "Median difference",
            "difference_ci_low": "Difference Q2.5",
            "difference_ci_high": "Difference Q97.5",
            "valid_resamples": "Valid bootstrap resamples",
        }
    )
    if not table.empty:
        table["Strategy"] = table["Strategy"].map(strategy_label)
        sections.append("<h4>Strategy summary</h4>")
        sections.append(table_html(table, classes="table table-sm table-striped"))

    matching = pd.DataFrame(analysis.manifest.get("matching_by_consequence", []))
    if not matching.empty:
        matching = matching.rename(
            columns={
                "strategy": "Strategy",
                "primary_consequence": "Primary VEP consequence",
                "eligible_focals": "VEP-annotated focal SNVs",
                "matched_focals": "Matched focal SNVs",
                "match_rate": "Matched %",
            }
        )
        matching["Strategy"] = matching["Strategy"].map(strategy_label)
        sections.append("<h4>Matching yield by consequence</h4>")
        sections.append(table_html(matching, classes="table table-sm table-striped"))

    consequence = analysis.consequence_summary.copy()
    if not consequence.empty:
        consequence = consequence.rename(
            columns={
                "strategy": "Strategy",
                "primary_consequence": "Primary VEP consequence",
                "matched_focals": "Matched SNVs",
                "observed_median": "GAPH median",
                "observed_ci_low": "GAPH median Q2.5",
                "observed_ci_high": "GAPH median Q97.5",
                "null_median": "Matched-control median",
                "null_ci_low": "Matched-control median Q2.5",
                "null_ci_high": "Matched-control median Q97.5",
                "median_difference": "Median difference",
                "difference_ci_low": "Difference Q2.5",
                "difference_ci_high": "Difference Q97.5",
                "valid_resamples": "Valid bootstrap resamples",
            }
        )
        consequence["Strategy"] = consequence["Strategy"].map(strategy_label)
        consequence = consequence[
            [
                "Strategy",
                "Primary VEP consequence",
                "Matched SNVs",
                "GAPH median",
                "GAPH median Q2.5",
                "GAPH median Q97.5",
                "Matched-control median",
                "Matched-control median Q2.5",
                "Matched-control median Q97.5",
                "Median difference",
                "Difference Q2.5",
                "Difference Q97.5",
                "Valid bootstrap resamples",
            ]
        ]
        sections.append("<h4>Results by primary VEP consequence</h4>")
        sections.append(table_html(consequence, classes="table table-sm table-striped"))

    outcome_frames = []
    if not analysis.gnomad_summary.empty:
        gnomad = analysis.gnomad_summary.copy()
        gnomad["Outcome"] = gnomad["metric"].map(
            {
                "found_fraction": "gnomAD found fraction",
                "median_af": "gnomAD median AF among exact hits",
            }
        )
        outcome_frames.append(gnomad)
    if not analysis.clinvar_summary.empty:
        clinvar = analysis.clinvar_summary.copy()
        clinvar["Outcome"] = "ClinVar found fraction"
        outcome_frames.append(clinvar)
    if not analysis.clinvar_class_summary.empty:
        clinvar_class = analysis.clinvar_class_summary.copy()
        clinvar_class["Outcome"] = "ClinVar class: " + clinvar_class["clinvar_class"].astype(str)
        outcome_frames.append(clinvar_class)
    if outcome_frames:
        outcomes = pd.concat(outcome_frames, ignore_index=True)
        outcomes = outcomes.rename(
            columns={
                "strategy": "Strategy",
                "matched_focals": "Matched SNVs",
                "observed_value": "GAPH statistic",
                "observed_ci_low": "GAPH Q2.5",
                "observed_ci_high": "GAPH Q97.5",
                "null_value": "Matched-control statistic",
                "null_ci_low": "Matched-control Q2.5",
                "null_ci_high": "Matched-control Q97.5",
                "difference": "Paired difference",
                "difference_ci_low": "Difference Q2.5",
                "difference_ci_high": "Difference Q97.5",
                "valid_resamples": "Valid bootstrap resamples",
            }
        )
        outcomes["Strategy"] = outcomes["Strategy"].map(strategy_label)
        outcomes = outcomes[
            [
                "Strategy",
                "Outcome",
                "Matched SNVs",
                "GAPH statistic",
                "GAPH Q2.5",
                "GAPH Q97.5",
                "Matched-control statistic",
                "Matched-control Q2.5",
                "Matched-control Q97.5",
                "Paired difference",
                "Difference Q2.5",
                "Difference Q97.5",
                "Valid bootstrap resamples",
            ]
        ]
        sections.append("<h4>External-evidence bootstrap summary</h4>")
        sections.append(table_html(outcomes, classes="table table-sm table-striped"))
    sections.append("</details>")
    return sections
