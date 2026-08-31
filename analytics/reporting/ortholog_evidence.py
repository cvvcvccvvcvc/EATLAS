"""Ortholog-support visualizations and report sections."""

from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from analytics.analyses.variant_summary import VariantSummary
from .components import format_int, strategy_label
from .config import (
    EVIDENCE_UNIT_LABELS,
    EVIDENCE_UNIT_ORDER,
    TAXONOMIC_SCOPE_LABELS,
    TAXONOMIC_SCOPE_ORDER,
)

def ortholog_evidence_figure(
    cells: pd.DataFrame,
    strategy: str,
    quantile_count: int,
    taxonomic_scope: str = "all",
    evidence_unit: str = "ortholog",
):
    contexts = [("cds", "CDS"), ("utr", "UTR"), ("intron", "Intron")]
    figure = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=[label for _context, label in contexts],
        horizontal_spacing=0.08,
    )
    selected = cells[
        cells["strategy"].astype(str).eq(strategy)
        & cells["quantile_count"].astype(int).eq(quantile_count)
        & cells["taxonomic_scope"].astype(str).eq(taxonomic_scope)
        & cells["evidence_unit"].astype(str).eq(evidence_unit)
    ]
    for column, (context, _label) in enumerate(contexts, start=1):
        subset = selected[selected["target_context"].astype(str).eq(context)]
        depth_labels = {int(row.depth_bin): str(row.depth_label) for row in subset.itertuples()}
        alt_labels = {
            int(row.alt_bin): str(row.alt_label)
            for row in subset.itertuples()
        }
        x = [depth_labels.get(index, f"Q{index + 1} (empty)") for index in range(quantile_count)]
        y = [
            alt_labels.get(index, f"Q{index + 1} (empty)")
            for index in range(quantile_count)
        ]
        values = {
            (int(row.depth_bin), int(row.alt_bin)): row
            for row in subset.itertuples()
        }
        z = []
        customdata = []
        for y_index in range(quantile_count):
            z_row = []
            custom_row = []
            for x_index in range(quantile_count):
                row = values.get((x_index, y_index))
                z_row.append(None if row is None else float(row.gnomad_found_fraction))
                custom_row.append(
                    [0, 0]
                    if row is None
                    else [int(row.gnomad_found_count), int(row.gnomad_eligible_count)]
                )
            z.append(z_row)
            customdata.append(custom_row)
        figure.add_trace(
            go.Heatmap(
                x=x,
                y=y,
                z=z,
                customdata=customdata,
                zmin=0,
                zmax=1,
                colorscale="Viridis",
                colorbar={"title": "gnomAD found", "tickformat": ".0%"},
                showscale=column == len(contexts),
                hoverongaps=False,
                hovertemplate=(
                    "Site-aligned orthologs: %{x}<br>"
                    "Exact-ALT support: %{y}<br>"
                    "gnomAD found: %{customdata[0]:,} / %{customdata[1]:,} "
                    "(%{z:.1%})<extra></extra>"
                ),
            ),
            row=1,
            col=column,
        )
        figure.update_xaxes(title_text="Site-aligned evidence units", row=1, col=column)
        if column == 1:
            figure.update_yaxes(title_text="Exact-ALT evidence units", row=1, col=column)
    figure.update_layout(
        height=500,
        margin={"l": 70, "r": 90, "t": 55, "b": 100},
        template="plotly_white",
    )
    figure.update_xaxes(automargin=True)
    figure.update_yaxes(automargin=True)
    return figure


def _distribution_ticks(maximum: int) -> tuple[list[float], list[str]]:
    candidates = [0, 1, 2, 5]
    scale = 10
    while scale <= maximum:
        candidates.extend([scale, 2 * scale, 5 * scale])
        scale *= 10
    values = sorted({value for value in candidates if value <= maximum})
    if maximum not in values:
        values.append(maximum)
    return [math.log1p(value) for value in values], [format_int(value) for value in values]


def _weighted_quantile_from_distribution(distribution: pd.DataFrame, quantile: float) -> int:
    ordered = distribution.sort_values("value", kind="mergesort")
    cumulative = ordered["variant_count"].astype("int64").cumsum().to_numpy()
    threshold = int(ordered["variant_count"].sum()) * quantile
    index = int(np.searchsorted(cumulative, threshold, side="left"))
    return int(ordered.iloc[min(index, len(ordered) - 1)]["value"])


def ortholog_evidence_distribution_figure(
    distributions: pd.DataFrame,
    strategy: str,
    taxonomic_scope: str = "all",
    evidence_unit: str = "ortholog",
):
    selected = distributions[
        distributions["strategy"].astype(str).eq(strategy)
        & distributions["taxonomic_scope"].astype(str).eq(taxonomic_scope)
        & distributions["evidence_unit"].astype(str).eq(evidence_unit)
    ]
    metrics = [
        ("site_aligned", "Site-aligned evidence units", "#2166ac"),
        ("exact_alt", "Exact-ALT evidence units", "#2ca25f"),
    ]
    figure = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=[label for _metric, label, _color in metrics],
        horizontal_spacing=0.12,
    )
    for column, (metric, _label, color) in enumerate(metrics, start=1):
        subset = selected[selected["metric"].astype(str).eq(metric)].copy()
        if subset.empty:
            continue
        subset["value"] = pd.to_numeric(subset["value"], errors="raise").astype("int64")
        subset["variant_count"] = (
            pd.to_numeric(subset["variant_count"], errors="raise").astype("int64")
        )
        subset = subset.groupby("value", as_index=False, sort=True)["variant_count"].sum()
        total = int(subset["variant_count"].sum())
        subset["cumulative_count"] = subset["variant_count"].cumsum()
        subset["exact_fraction"] = subset["variant_count"] / total
        subset["cumulative_fraction"] = subset["cumulative_count"] / total
        customdata = np.column_stack(
            [
                subset["value"],
                subset["variant_count"],
                subset["exact_fraction"],
                subset["cumulative_count"],
                np.full(len(subset), total),
            ]
        )
        figure.add_trace(
            go.Scatter(
                x=np.log1p(subset["value"]),
                y=subset["cumulative_fraction"],
                customdata=customdata,
                mode="lines+markers",
                line={"color": color, "shape": "hv", "width": 2.5},
                marker={"color": color, "size": 5},
                showlegend=False,
                hovertemplate=(
                    "Evidence units: %{customdata[0]:,.0f}<br>"
                    "Exactly this value: %{customdata[1]:,.0f} (%{customdata[2]:.1%})<br>"
                    "Cumulative: %{customdata[3]:,.0f} / %{customdata[4]:,.0f} "
                    "(%{y:.1%})<extra></extra>"
                ),
            ),
            row=1,
            col=column,
        )
        tick_values, tick_labels = _distribution_ticks(int(subset["value"].max()))
        figure.update_xaxes(
            title_text="Evidence units (log1p scale)",
            tickmode="array",
            tickvals=tick_values,
            ticktext=tick_labels,
            row=1,
            col=column,
        )
        figure.add_hline(
            y=0.5,
            line={"color": "#9e9e9e", "dash": "dot", "width": 1},
            row=1,
            col=column,
        )
    figure.update_yaxes(
        title_text="Cumulative SNVs at or below X",
        tickformat=".0%",
        range=[0, 1.02],
        row=1,
        col=1,
    )
    figure.update_yaxes(tickformat=".0%", range=[0, 1.02], row=1, col=2)
    figure.update_layout(
        height=410,
        margin={"l": 70, "r": 30, "t": 55, "b": 70},
        template="plotly_white",
    )
    figure.update_xaxes(automargin=True)
    figure.update_yaxes(automargin=True)
    return figure


def ortholog_evidence_distribution_stats(
    distributions: pd.DataFrame,
    strategy: str,
    taxonomic_scope: str,
    evidence_unit: str,
) -> list[dict[str, str]]:
    selected = distributions[
        distributions["strategy"].astype(str).eq(strategy)
        & distributions["taxonomic_scope"].astype(str).eq(taxonomic_scope)
        & distributions["evidence_unit"].astype(str).eq(evidence_unit)
    ]
    stats = []
    totals = []
    for metric, label in (
        ("site_aligned", "Site-aligned median [IQR]"),
        ("exact_alt", "Exact-ALT median [IQR]"),
    ):
        subset = selected[selected["metric"].astype(str).eq(metric)].copy()
        if subset.empty:
            continue
        subset["value"] = pd.to_numeric(subset["value"], errors="raise").astype("int64")
        subset["variant_count"] = (
            pd.to_numeric(subset["variant_count"], errors="raise").astype("int64")
        )
        totals.append(int(subset["variant_count"].sum()))
        q1 = _weighted_quantile_from_distribution(subset, 0.25)
        median = _weighted_quantile_from_distribution(subset, 0.5)
        q3 = _weighted_quantile_from_distribution(subset, 0.75)
        stats.append({"label": label, "value": f"{format_int(median)} [{format_int(q1)}-{format_int(q3)}]"})
    if len(set(totals)) > 1:
        raise ValueError("Ortholog evidence distributions use inconsistent SNV totals")
    return [{"label": "SNVs", "value": format_int(totals[0] if totals else 0)}, *stats]


def build_ortholog_evidence_sections(
    variant_summary: VariantSummary,
    taxonomy_summary: pd.DataFrame,
) -> list[str]:
    sections = [
        "<h2>Ortholog Evidence</h2>",
    ]
    cells = variant_summary.ortholog_evidence_cells
    distributions = variant_summary.ortholog_evidence_distributions
    if cells.empty:
        sections.append("<p>No eligible SNVs with ortholog evidence and successful gnomAD lookup.</p>")
        return sections

    supported_strategies = [
        strategy
        for strategy in variant_summary.strategies
        if strategy in set(cells["strategy"].astype(str))
    ]
    if not supported_strategies:
        sections.append("<p>No strategies expose taxonomically identified ortholog evidence.</p>")
        return sections
    default_strategy = supported_strategies[0]
    quantile_options = {2: "Median", 4: "Quartiles", 10: "Deciles"}
    available_scopes = set(cells["taxonomic_scope"].astype(str))
    visible_scopes = []
    seen_scope_signatures = set()
    for scope in TAXONOMIC_SCOPE_ORDER:
        if scope not in available_scopes:
            continue
        signature = None
        if not taxonomy_summary.empty:
            row = taxonomy_summary[
                taxonomy_summary["taxonomic_scope"].astype(str).eq(scope)
                & taxonomy_summary["evidence_unit"].astype(str).eq("ortholog")
            ]
            if not row.empty:
                signature = (
                    int(row.iloc[0]["ortholog_count"]),
                    int(row.iloc[0]["taxon_count"]),
                    float(row.iloc[0]["orthologs_per_gene_median"]),
                )
        if signature is not None and signature in seen_scope_signatures:
            continue
        visible_scopes.append(scope)
        if signature is not None:
            seen_scope_signatures.add(signature)
    visible_scopes.extend(sorted(available_scopes - set(visible_scopes) - set(TAXONOMIC_SCOPE_ORDER)))
    available_units = [
        unit
        for unit in EVIDENCE_UNIT_ORDER
        if unit in set(cells["evidence_unit"].astype(str))
    ]
    default_scope = "all" if "all" in visible_scopes else visible_scopes[0]
    default_unit = "ortholog" if "ortholog" in available_units else available_units[0]

    figures = {}
    distribution_figures = {}
    distribution_stats = {}
    for strategy in supported_strategies:
        figures[strategy] = {}
        distribution_figures[strategy] = {}
        distribution_stats[strategy] = {}
        for scope in visible_scopes:
            scoped = cells[
                cells["strategy"].astype(str).eq(strategy)
                & cells["taxonomic_scope"].astype(str).eq(scope)
            ]
            if scoped.empty:
                continue
            figures[strategy][scope] = {}
            distribution_figures[strategy][scope] = {}
            distribution_stats[strategy][scope] = {}
            for unit in available_units:
                if scoped["evidence_unit"].astype(str).eq(unit).sum() == 0:
                    continue
                figures[strategy][scope][unit] = {}
                for quantile_count in quantile_options:
                    figure = ortholog_evidence_figure(
                        cells,
                        strategy,
                        quantile_count,
                        scope,
                        unit,
                    )
                    figures[strategy][scope][unit][str(quantile_count)] = json.loads(
                        figure.to_json()
                    )
                if not distributions.empty:
                    distribution_subset = distributions[
                        distributions["strategy"].astype(str).eq(strategy)
                        & distributions["taxonomic_scope"].astype(str).eq(scope)
                        & distributions["evidence_unit"].astype(str).eq(unit)
                    ]
                    if not distribution_subset.empty:
                        distribution_figure = ortholog_evidence_distribution_figure(
                            distributions,
                            strategy,
                            scope,
                            unit,
                        )
                        distribution_figures[strategy][scope][unit] = json.loads(
                            distribution_figure.to_json()
                        )
                        distribution_stats[strategy][scope][unit] = (
                            ortholog_evidence_distribution_stats(
                                distributions,
                                strategy,
                                scope,
                                unit,
                            )
                        )

    initial = ortholog_evidence_figure(
        cells,
        default_strategy,
        4,
        default_scope,
        default_unit,
    )
    initial_html = initial.to_html(
        full_html=False,
        include_plotlyjs=False,
        div_id="ortholog-evidence-plot",
    )
    strategy_options = "".join(
        f'<option value="{strategy}">{strategy_label(strategy)}</option>'
        for strategy in supported_strategies
    )
    unsupported_options = "".join(
        f'<option disabled>{strategy_label(strategy)} (taxonomy unavailable)</option>'
        for strategy in variant_summary.strategies
        if strategy not in supported_strategies
    )
    scope_options = "".join(
        f'<option value="{scope}"{" selected" if scope == default_scope else ""}>'
        f'{TAXONOMIC_SCOPE_LABELS.get(scope, scope)}</option>'
        for scope in visible_scopes
    )
    unit_options = "".join(
        f'<option value="{unit}"{" selected" if unit == default_unit else ""}>'
        f'{EVIDENCE_UNIT_LABELS.get(unit, unit)}</option>'
        for unit in available_units
    )
    quantile_html = "".join(
        f'<option value="{count}"{" selected" if count == 4 else ""}>{label}</option>'
        for count, label in quantile_options.items()
    )
    payload = json.dumps(figures, separators=(",", ":"))
    distribution_payload = json.dumps(distribution_figures, separators=(",", ":"))
    distribution_stats_payload = json.dumps(distribution_stats, separators=(",", ":"))
    taxonomy_stats = {}
    if not taxonomy_summary.empty:
        for scope in visible_scopes:
            ortholog_row = taxonomy_summary[
                taxonomy_summary["taxonomic_scope"].astype(str).eq(scope)
                & taxonomy_summary["evidence_unit"].astype(str).eq("ortholog")
            ]
            if ortholog_row.empty:
                continue
            taxonomy_stats[scope] = {}
            ortholog_median = float(ortholog_row.iloc[0]["orthologs_per_gene_median"])
            for unit in available_units:
                unit_row = taxonomy_summary[
                    taxonomy_summary["taxonomic_scope"].astype(str).eq(scope)
                    & taxonomy_summary["evidence_unit"].astype(str).eq(unit)
                ]
                if unit_row.empty:
                    continue
                row = unit_row.iloc[0]
                taxonomy_stats[scope][unit] = (
                    f"Median selected orthologs/gene: {ortholog_median:,.1f}; "
                    f"median {EVIDENCE_UNIT_LABELS.get(unit, unit).lower()} units/gene: "
                    f"{float(row['units_per_gene_median']):,.1f}; "
                    f"distinct units in run: {int(row['unit_count']):,}."
                )
    taxonomy_stats_payload = json.dumps(taxonomy_stats, separators=(",", ":"))
    distribution_html = ""
    if (
        default_strategy in distribution_figures
        and default_scope in distribution_figures[default_strategy]
        and default_unit in distribution_figures[default_strategy][default_scope]
    ):
        initial_distribution = ortholog_evidence_distribution_figure(
            distributions,
            default_strategy,
            default_scope,
            default_unit,
        )
        distribution_html = (
            '<h3>Evidence distributions</h3>'
            '<div class="metric-grid" id="ortholog-evidence-distribution-stats"></div>'
            + initial_distribution.to_html(
                full_html=False,
                include_plotlyjs=False,
                div_id="ortholog-evidence-distribution-plot",
            )
        )
    sections.append(
        f"""
        <div class="analysis-controls" id="ortholog-evidence-controls">
            <label>Strategy<select id="ortholog-evidence-strategy">{strategy_options}{unsupported_options}</select></label>
            <label>Taxonomic scope<select id="ortholog-evidence-scope">{scope_options}</select></label>
            <label>Evidence unit<select id="ortholog-evidence-unit">{unit_options}</select></label>
            <label>Groups<select id="ortholog-evidence-quantiles">{quantile_html}</select></label>
        </div>
        <p class="analysis-note" id="ortholog-evidence-stats"></p>
        {initial_html}
        {distribution_html}
        <script>
        (() => {{
            const figures = {payload};
            const taxonomyStats = {taxonomy_stats_payload};
            const distributionFigures = {distribution_payload};
            const distributionStats = {distribution_stats_payload};
            const strategy = document.getElementById('ortholog-evidence-strategy');
            const scope = document.getElementById('ortholog-evidence-scope');
            const unit = document.getElementById('ortholog-evidence-unit');
            const quantiles = document.getElementById('ortholog-evidence-quantiles');
            const summary = document.getElementById('ortholog-evidence-stats');
            const distributionCards = document.getElementById('ortholog-evidence-distribution-stats');
            const firstKey = value => Object.keys(value)[0];
            const render = () => {{
                const strategyFigures = figures[strategy.value];
                if (!strategyFigures[scope.value]) scope.value = firstKey(strategyFigures);
                const scopeFigures = strategyFigures[scope.value];
                if (!scopeFigures[unit.value]) unit.value = firstKey(scopeFigures);
                const figure = scopeFigures[unit.value][quantiles.value];
                summary.textContent = taxonomyStats[scope.value]?.[unit.value] || '';
                Plotly.react('ortholog-evidence-plot', figure.data, figure.layout, {{responsive: true}});
                const distribution = distributionFigures[strategy.value]?.[scope.value]?.[unit.value];
                if (distribution && distributionCards) {{
                    const items = distributionStats[strategy.value][scope.value][unit.value];
                    distributionCards.innerHTML = items.map(item =>
                        `<div class="metric-card"><div class="metric-label">${{item.label}}</div>` +
                        `<div class="metric-value">${{item.value}}</div></div>`
                    ).join('');
                    Plotly.react(
                        'ortholog-evidence-distribution-plot',
                        distribution.data,
                        distribution.layout,
                        {{responsive: true}}
                    );
                }}
            }};
            strategy.addEventListener('change', render);
            scope.addEventListener('change', render);
            unit.addEventListener('change', render);
            quantiles.addEventListener('change', render);
            render();
        }})();
        </script>
        """
    )
    return sections
