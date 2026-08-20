"""Conservation-adjusted validation analysis and report sections."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from analytics.analyses.candidate_conservation import CandidateConservation
from analytics.analyses.conservation_analysis import ConservationAnalysis
from analytics.analyses.conservation_validation import (
    TARGET_CONTEXT_OPTIONS,
    VARIANT_TYPE_OPTIONS,
    ConservationValidation,
)
from analytics.annotation.consequences import (
    VALIDATION_CONSEQUENCE_OPTIONS as CONSEQUENCE_OPTIONS,
)
from .components import compact_figure, dataframe_records, format_int, strategy_label


def hidden_clinvar_association_views(validation: ConservationValidation) -> tuple[pd.DataFrame, pd.DataFrame]:
    mode_frames = [
        ("Unadjusted", validation.unadjusted),
        ("phyloP fixed bands", validation.fixed_adjusted),
        ("phyloP continuous", validation.continuous),
    ]
    variant_labels = dict(VARIANT_TYPE_OPTIONS)
    context_labels = dict(TARGET_CONTEXT_OPTIONS)
    consequence_labels = dict(CONSEQUENCE_OPTIONS)
    hidden_rows = []
    summary_rows = []
    group_columns = ["variant_type", "target_context", "consequence"]
    for mode_label, frame in mode_frames:
        hidden_count = 0
        visible_count = 0
        for keys, group in frame.groupby(group_columns, sort=False, dropna=False):
            visible = group["status"].astype(str).ne("not_estimable").any()
            if visible:
                visible_count += 1
                continue
            hidden_count += 1
            usable = pd.to_numeric(group["usable_rows"], errors="coerce").dropna().astype(int)
            reasons = sorted({str(value) for value in group["reason"] if str(value)})
            variant_type, target_context, consequence = keys
            hidden_rows.append(
                {
                    "Analysis": mode_label,
                    "Variant type": variant_labels.get(str(variant_type), str(variant_type)),
                    "Target context": context_labels.get(str(target_context), str(target_context)),
                    "Consequence subset": consequence_labels.get(str(consequence), str(consequence)),
                    "N across strategies": (
                        f"{format_int(usable.min())}-{format_int(usable.max())}" if not usable.empty else "0"
                    ),
                    "Reason": "; ".join(reasons) or "No estimable strategy result.",
                }
            )
        summary_rows.append(
            {
                "Analysis": mode_label,
                "Displayed selector combinations": visible_count,
                "Hidden selector combinations": hidden_count,
            }
        )
    return pd.DataFrame(summary_rows), pd.DataFrame(hidden_rows)


def candidate_phylop_figure(
    analysis: CandidateConservation,
    strategy_order: list[str],
):
    distributions = analysis.distributions.copy()
    if distributions.empty:
        return None
    distributions["Strategy"] = distributions["strategy"].astype(str).map(strategy_label)
    available = set(distributions["Strategy"])
    ordered = [strategy for strategy in strategy_order if strategy in available]
    ordered += sorted(available - set(ordered))
    status_styles = {
        "found": ("Found in gnomAD", "#2166ac"),
        "not_found": ("Not found in gnomAD", "#b2182b"),
    }
    fig = go.Figure()
    trace_strategies = []
    for strategy_index, strategy in enumerate(ordered):
        for status, (label, color) in status_styles.items():
            subset = distributions[
                distributions["Strategy"].eq(strategy)
                & distributions["gnomad_status"].astype(str).eq(status)
            ].sort_values("quantile")
            if subset.empty:
                continue
            coverage = subset["scored_count"].iloc[0] / subset["variant_count"].iloc[0]
            fig.add_trace(
                go.Scatter(
                    x=subset["phyloP100way"],
                    y=subset["quantile"],
                    mode="lines",
                    name=label,
                    line={"color": color, "width": 3},
                    visible=strategy_index == 0,
                    customdata=np.column_stack(
                        [
                            np.repeat(subset["scored_count"].iloc[0], len(subset)),
                            np.repeat(subset["variant_count"].iloc[0], len(subset)),
                            np.repeat(coverage, len(subset)),
                        ]
                    ),
                    hovertemplate=(
                        label + "<br>phyloP100way: %{x:.3f}<br>Percentile: %{y:.0%}<br>"
                        "Scored variants: %{customdata[0]:,} / %{customdata[1]:,} "
                        "(%{customdata[2]:.1%})<extra></extra>"
                    ),
                )
            )
            trace_strategies.append(strategy)
    if not fig.data:
        return None
    buttons = []
    for strategy in ordered:
        visible = [trace_strategy == strategy for trace_strategy in trace_strategies]
        if any(visible):
            buttons.append(
                {
                    "label": strategy,
                    "method": "update",
                    "args": [
                        {"visible": visible},
                        {"title": f"Candidate phyloP100way distribution: {strategy}"},
                    ],
                }
            )
    first_strategy = buttons[0]["label"] if buttons else ""
    fig.update_layout(
        title=f"Candidate phyloP100way distribution: {first_strategy}",
        xaxis_title="phyloP100way",
        yaxis_title="Cumulative fraction",
        yaxis_tickformat=".0%",
        updatemenus=[
            {
                "buttons": buttons,
                "direction": "down",
                "showactive": True,
                "x": 1.0,
                "xanchor": "right",
                "y": 1.18,
                "yanchor": "top",
            }
        ],
    )
    fig.add_vline(x=0.0, line_dash="dot", line_color="#8c8c8c")
    compact_figure(fig, height=420, show_x_title=True)
    return fig


def candidate_phylop_summary_figure(
    analysis: CandidateConservation,
    strategy_order: list[str],
):
    histograms = analysis.histograms.copy()
    groups = pd.DataFrame(analysis.manifest.get("groups", []))
    if histograms.empty or groups.empty:
        return None
    histograms["Strategy"] = histograms["strategy"].astype(str).map(strategy_label)
    groups["Strategy"] = groups["strategy"].astype(str).map(strategy_label)
    available = set(histograms["Strategy"])
    ordered = [strategy for strategy in strategy_order if strategy in available]
    ordered += sorted(available - set(ordered))
    status_styles = {
        "found": ("Found in gnomAD", "#2166ac"),
        "not_found": ("Not found in gnomAD", "#b2182b"),
    }
    fig = make_subplots(
        rows=1,
        cols=2,
        horizontal_spacing=0.12,
        subplot_titles=["Relative-frequency histogram", "Box plot"],
    )
    trace_strategies = []
    for strategy_index, strategy in enumerate(ordered):
        for status, (label, color) in status_styles.items():
            histogram = histograms[
                histograms["Strategy"].eq(strategy)
                & histograms["gnomad_status"].astype(str).eq(status)
            ].sort_values("bin_left")
            group = groups[
                groups["Strategy"].eq(strategy)
                & groups["gnomad_status"].astype(str).eq(status)
            ]
            if histogram.empty or group.empty:
                continue
            centers = (histogram["bin_left"] + histogram["bin_right"]) / 2
            widths = histogram["bin_right"] - histogram["bin_left"]
            visible = strategy_index == 0
            fig.add_trace(
                go.Bar(
                    x=centers,
                    y=histogram["fraction"],
                    width=widths,
                    name=label,
                    legendgroup=status,
                    marker_color=color,
                    opacity=0.58,
                    visible=visible,
                    customdata=np.column_stack(
                        [histogram["count"], histogram["bin_left"], histogram["bin_right"]]
                    ),
                    hovertemplate=(
                        label + "<br>phyloP: %{customdata[1]:.3f} to %{customdata[2]:.3f}<br>"
                        "Variants: %{customdata[0]:,}<br>Fraction: %{y:.2%}<extra></extra>"
                    ),
                ),
                row=1,
                col=1,
            )
            trace_strategies.append(strategy)
            summary = group.iloc[0]
            fig.add_trace(
                go.Box(
                    q1=[summary["q1"]],
                    median=[summary["median"]],
                    q3=[summary["q3"]],
                    lowerfence=[summary["lower_whisker"]],
                    upperfence=[summary["upper_whisker"]],
                    name=label,
                    legendgroup=status,
                    marker_color=color,
                    boxpoints=False,
                    showlegend=False,
                    visible=visible,
                    hovertemplate=(
                        label + "<br>Q1: %{q1:.3f}<br>Median: %{median:.3f}<br>"
                        "Q3: %{q3:.3f}<extra></extra>"
                    ),
                ),
                row=1,
                col=2,
            )
            trace_strategies.append(strategy)
    if not fig.data:
        return None
    buttons = []
    for strategy in ordered:
        visible = [trace_strategy == strategy for trace_strategy in trace_strategies]
        if any(visible):
            buttons.append(
                {
                    "label": strategy,
                    "method": "update",
                    "args": [
                        {"visible": visible},
                        {"title": f"Candidate phyloP100way: {strategy}"},
                    ],
                }
            )
    first_strategy = buttons[0]["label"] if buttons else ""
    fig.update_layout(
        title=f"Candidate phyloP100way: {first_strategy}",
        barmode="overlay",
        boxmode="group",
        updatemenus=[
            {
                "buttons": buttons,
                "direction": "down",
                "showactive": True,
                "x": 1.0,
                "xanchor": "right",
                "y": 1.2,
                "yanchor": "top",
            }
        ],
    )
    fig.update_xaxes(title_text="phyloP100way", row=1, col=1)
    fig.update_yaxes(title_text="Fraction per bin", tickformat=".0%", row=1, col=1)
    fig.update_yaxes(title_text="phyloP100way", row=1, col=2)
    fig.add_vline(x=0.0, line_dash="dot", line_color="#8c8c8c", row=1, col=1)
    compact_figure(fig, height=430, show_x_title=True)
    return fig


def build_clinvar_association_sections(analysis: ConservationAnalysis) -> list[str]:
    return [
        "<h2>ClinVar Association</h2>",
        clinvar_association_view(analysis.validation),
    ]


def clinvar_association_view(
    validation: ConservationValidation,
    *,
    view_id: str = "clinvar-association",
    strategy_labels: dict[str, str] | None = None,
) -> str:
    display_labels = strategy_labels or {}
    primary_frames = []
    mode_specs = [
        ("unadjusted", validation.unadjusted, "odds_ratio", "fisher_p", "fisher_q"),
        ("fixed", validation.fixed_adjusted, "odds_ratio_mh", "cmh_p", "cmh_q"),
        ("continuous", validation.continuous, "odds_ratio", "plr_p", "plr_q"),
    ]
    for mode, source, odds_ratio, p_value, q_value in mode_specs:
        frame = source.copy()
        if frame.empty:
            continue
        frame["mode"] = mode
        frame["result_or"] = frame[odds_ratio]
        frame["result_p"] = frame[p_value]
        frame["result_q"] = frame[q_value]
        primary_frames.append(frame)
    primary = pd.concat(primary_frames, ignore_index=True) if primary_frames else pd.DataFrame()
    strategies = validation.unadjusted["strategy"].drop_duplicates().astype(str).tolist()
    payload = {
        "viewId": view_id,
        "modes": [
            {"key": "unadjusted", "label": "Unadjusted"},
            {"key": "fixed", "label": "phyloP fixed bands"},
            {"key": "continuous", "label": "phyloP continuous"},
        ],
        "strategies": [
            {
                "key": value,
                "label": display_labels.get(value, strategy_label(value)),
            }
            for value in strategies
        ],
        "variantTypes": [{"key": key, "label": label} for key, label in VARIANT_TYPE_OPTIONS],
        "targetContexts": [{"key": key, "label": label} for key, label in TARGET_CONTEXT_OPTIONS],
        "consequences": [{"key": key, "label": label} for key, label in CONSEQUENCE_OPTIONS],
        "primary": dataframe_records(primary),
        "fixedDetail": dataframe_records(validation.fixed_bins),
        "continuousDetail": dataframe_records(validation.distributions),
    }
    payload_json = json.dumps(payload, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    return f"""
    <div class="analysis-controls" id="{view_id}-controls">
      <label>Analysis<select data-role="mode"></select></label>
      <label>Variant type<select data-role="variant-type"></select></label>
      <label>Target context<select data-role="target-context"></select></label>
      <label>Consequence subset<select data-role="consequence"></select></label>
    </div>
    <div id="{view_id}-status" class="analysis-note" hidden></div>
    <div id="{view_id}-forest" class="analysis-plot"></div>
    <div id="{view_id}-results"></div>
    <div class="analysis-controls analysis-controls-single" id="{view_id}-strategy-control">
      <label>Inspect strategy<select data-role="strategy"></select></label>
    </div>
    <div id="{view_id}-detail-plot" class="analysis-plot"></div>
    <div id="{view_id}-detail-table"></div>
    <script>
    (() => {{
      const config = {payload_json};
      const controls = document.getElementById(config.viewId + '-controls');
      const modeSelect = controls.querySelector('[data-role="mode"]');
      const variantSelect = controls.querySelector('[data-role="variant-type"]');
      const targetContextSelect = controls.querySelector('[data-role="target-context"]');
      const consequenceSelect = controls.querySelector('[data-role="consequence"]');
      const strategySelect = document.querySelector('#' + config.viewId + '-strategy-control [data-role="strategy"]');
      const optionMap = values => Object.fromEntries(values.map(value => [value.key, value.label]));
      const strategyLabels = optionMap(config.strategies);
      const modeLabels = optionMap(config.modes);
      const variantLabels = optionMap(config.variantTypes);
      const targetContextLabels = optionMap(config.targetContexts);
      const consequenceLabels = optionMap(config.consequences);
      const addOptions = (select, values) => values.forEach(value => {{
        const option = document.createElement('option');
        option.value = value.key; option.textContent = value.label; select.appendChild(option);
      }});
      addOptions(modeSelect, config.modes);
      addOptions(variantSelect, config.variantTypes);
      addOptions(targetContextSelect, config.targetContexts);
      addOptions(strategySelect, config.strategies);
      modeSelect.value = 'unadjusted';
      variantSelect.value = 'snv';
      targetContextSelect.value = 'all';
      consequenceSelect.value = 'all';

      const finite = value => value !== null && value !== 'inf' && value !== '-inf' && Number.isFinite(Number(value));
      const number = value => finite(value) ? Number(value) : null;
      const count = value => number(value) === null ? '0' : Math.round(Number(value)).toLocaleString('en-US').replaceAll(',', ' ');
      const fmt = value => {{
        if (value === 'inf') return '∞';
        if (value === '-inf') return '-∞';
        const item = number(value);
        if (item === null) return 'NA';
        if (item === 0) return '0';
        if (Math.abs(item) < 0.001 || Math.abs(item) >= 1000) return item.toExponential(2);
        return item.toPrecision(3);
      }};
      const ci = row => finite(row?.ci_low) && finite(row?.ci_high) ? `${{fmt(row.ci_low)}}–${{fmt(row.ci_high)}}` : 'NA';
      const effect = row => `${{fmt(row?.result_or)}} [${{ci(row)}}]`;
      const cell = value => `<td>${{value}}</td>`;
      const statusText = row => {{
        if (!row) return 'Not available';
        if (row.status === 'estimated') return 'Estimated';
        if (row.status === 'test_only') return 'Test only';
        return row.reason || 'Not estimable';
      }};
      const matchesSelection = row => row.mode === modeSelect.value
        && row.variant_type === variantSelect.value
        && row.target_context === targetContextSelect.value
        && row.consequence === consequenceSelect.value;
      const plotValue = row => {{
        const raw = number(row.result_or);
        if (raw !== null && raw > 0) return raw;
        if (row.result_or === 'inf' && finite(row.ci_low) && finite(row.ci_high)) {{
          return Math.sqrt(Number(row.ci_low) * Number(row.ci_high));
        }}
        return null;
      }};
      const currentRows = () => config.primary.filter(matchesSelection);

      function refreshConsequences() {{
        const available = new Set(config.primary.filter(row =>
          row.mode === modeSelect.value
          && row.variant_type === variantSelect.value
          && row.target_context === targetContextSelect.value
          && row.status !== 'not_estimable'
        ).map(row => row.consequence));
        const previous = consequenceSelect.value || 'all';
        consequenceSelect.replaceChildren();
        const options = config.consequences.filter(value => available.has(value.key));
        addOptions(consequenceSelect, options);
        if (available.has(previous)) consequenceSelect.value = previous;
        else if (available.has('all')) consequenceSelect.value = 'all';
        else if (options.length) consequenceSelect.value = options[0].key;
        else {{
          const option = document.createElement('option');
          option.value = ''; option.textContent = 'No estimable subsets'; option.disabled = true;
          consequenceSelect.appendChild(option);
          consequenceSelect.value = '';
        }}
      }}

      function renderForest(rows) {{
        const plotted = rows.map(row => ({{row, x: plotValue(row)}}))
          .filter(item => item.x !== null && finite(item.row.ci_low) && finite(item.row.ci_high))
          .sort((left, right) => right.x - left.x);
        const consequenceLabel = consequenceLabels[consequenceSelect.value] || 'No estimable subset';
        const title = `${{modeLabels[modeSelect.value]}}: ${{variantLabels[variantSelect.value]}}, ${{targetContextLabels[targetContextSelect.value]}}, ${{consequenceLabel}}`;
        const trace = {{
          type: 'scatter', mode: 'markers',
          x: plotted.map(item => item.x),
          y: plotted.map(item => strategyLabels[item.row.strategy] || item.row.strategy),
          marker: {{size: 10, color: '#356d8f'}},
          error_x: {{
            type: 'data', symmetric: false,
            array: plotted.map(item => Number(item.row.ci_high) - item.x),
            arrayminus: plotted.map(item => item.x - Number(item.row.ci_low)),
          }},
          customdata: plotted.map(item => [
            fmt(item.row.result_or), ci(item.row), fmt(item.row.result_p), fmt(item.row.result_q),
            count(item.row.usable_rows), statusText(item.row),
          ]),
          hovertemplate: '%{{y}}<br>OR: %{{customdata[0]}}<br>95% CI: %{{customdata[1]}}<br>p: %{{customdata[2]}}<br>FDR q: %{{customdata[3]}}<br>N: %{{customdata[4]}}<br>%{{customdata[5]}}<extra></extra>',
        }};
        const annotations = plotted.length ? [] : [{{
          text: 'No finite odds ratio and confidence interval for this selection.',
          x: 0.5, y: 0.5, xref: 'paper', yref: 'paper', showarrow: false,
        }}];
        Plotly.react(config.viewId + '-forest', plotted.length ? [trace] : [], {{
          title, template: 'plotly_white', height: 370,
          margin: {{l: 170, r: 30, t: 52, b: 58}},
          xaxis: {{title: 'Odds ratio (log scale)', type: 'log', dtick: 1}}, yaxis: {{title: ''}},
          shapes: [{{type: 'line', x0: 1, x1: 1, y0: 0, y1: 1, yref: 'paper', line: {{dash: 'dash', color: '#8c8c8c'}}}}],
          annotations,
        }}, {{responsive: true}});
      }}

      function renderResultsTable(rows) {{
        const ordered = [...rows].sort((left, right) => (plotValue(right) || -1) - (plotValue(left) || -1));
        const body = ordered.map(row => `<tr>${{cell(strategyLabels[row.strategy] || row.strategy)}}${{cell(effect(row))}}${{cell(fmt(row.result_p))}}${{cell(fmt(row.result_q))}}${{cell(count(row.usable_rows))}}${{cell(statusText(row))}}</tr>`).join('');
        document.getElementById(config.viewId + '-results').innerHTML = `<table><thead><tr><th>Strategy</th><th>OR [95% CI]</th><th>p</th><th>FDR q</th><th>N</th><th>Status</th></tr></thead><tbody>${{body}}</tbody></table>`;
      }}

      function twoByTwoTable(row) {{
        if (!row) return '';
        const observedTotal = Number(row.benign_observed || 0) + Number(row.pathogenic_observed || 0);
        const notObservedTotal = Number(row.benign_not_observed || 0) + Number(row.pathogenic_not_observed || 0);
        return `<table><thead><tr><th>ALT status</th><th>B/LB</th><th>P/LP</th><th>Total</th></tr></thead><tbody>`
          + `<tr>${{cell('Observed')}}${{cell(count(row.benign_observed))}}${{cell(count(row.pathogenic_observed))}}${{cell(count(observedTotal))}}</tr>`
          + `<tr>${{cell('Not observed')}}${{cell(count(row.benign_not_observed))}}${{cell(count(row.pathogenic_not_observed))}}${{cell(count(notObservedTotal))}}</tr>`
          + `</tbody></table>`;
      }}

      function renderUnadjusted(row) {{
        document.getElementById(config.viewId + '-detail-plot').hidden = true;
        document.getElementById(config.viewId + '-detail-table').innerHTML = twoByTwoTable(row);
      }}

      function renderFixed(row) {{
        const details = config.fixedDetail.filter(item => item.strategy === strategySelect.value
          && item.variant_type === variantSelect.value
          && item.target_context === targetContextSelect.value
          && item.consequence === consequenceSelect.value);
        const groups = [
          ['ALT observed', 'benign_observed', 'pathogenic_observed', '#2166ac'],
          ['ALT not observed', 'benign_not_observed', 'pathogenic_not_observed', '#8c8c8c'],
        ];
        const traces = groups.map(([label, benignKey, pathogenicKey, color]) => {{
          const fractions = details.map(item => {{
            const denominator = Number(item[benignKey] || 0) + Number(item[pathogenicKey] || 0);
            return denominator ? Number(item[benignKey]) / denominator : 0;
          }});
          return {{
            type: 'bar', name: label, x: details.map(item => item.band_label), y: fractions,
            marker: {{color}},
            customdata: details.map(item => [count(item[benignKey]), count(item[pathogenicKey])]),
            hovertemplate: '%{{x}}<br>' + label + '<br>B/LB: %{{customdata[0]}}<br>P/LP: %{{customdata[1]}}<br>B/LB fraction: %{{y:.1%}}<extra></extra>',
          }};
        }});
        const plot = document.getElementById(config.viewId + '-detail-plot'); plot.hidden = false;
        Plotly.react(plot, traces, {{
          title: 'B/LB fraction within phyloP bands', template: 'plotly_white', barmode: 'group', height: 350,
          margin: {{l: 65, r: 25, t: 50, b: 60}}, yaxis: {{title: 'B/LB fraction', tickformat: '.0%', range: [0, 1]}},
          xaxis: {{title: ''}}, legend: {{orientation: 'h', y: 1.12}},
        }}, {{responsive: true}});
        const body = details.map(item => {{
          const label = `${{item.band_label}} (${{item.band_range}})`;
          const observed = `${{count(item.benign_observed)}} / ${{count(item.pathogenic_observed)}}`;
          const notObserved = `${{count(item.benign_not_observed)}} / ${{count(item.pathogenic_not_observed)}}`;
          return `<tr>${{cell(label)}}${{cell(count(item.row_count))}}${{cell(observed)}}${{cell(notObserved)}}${{cell(`${{fmt(item.odds_ratio)}} [${{ci(item)}}]`)}}${{cell(`${{fmt(item.fisher_p)}} / ${{fmt(item.fisher_q)}}`)}}${{cell(statusText(item))}}</tr>`;
        }}).join('');
        document.getElementById(config.viewId + '-detail-table').innerHTML = `<table><thead><tr><th>Band</th><th>N</th><th>Observed B/LB / P/LP</th><th>Not observed B/LB / P/LP</th><th>OR [95% CI]</th><th>p / FDR q</th><th>Status</th></tr></thead><tbody>${{body}}</tbody></table>`;
      }}

      function renderContinuous(row) {{
        const details = config.continuousDetail.filter(item => item.strategy === strategySelect.value
          && item.variant_type === variantSelect.value
          && item.target_context === targetContextSelect.value
          && item.consequence === consequenceSelect.value);
        const styles = {{
          'ALT observed': '#2166ac',
          'ALT not observed': '#8c8c8c',
        }};
        const traces = [];
        Object.entries(styles).forEach(([group, color]) => {{
          const values = details.filter(item => item.group === group);
          if (!values.length) return;
          traces.push({{
            type: 'bar', name: group, legendgroup: group,
            x: values.map(item => (Number(item.bin_left) + Number(item.bin_right)) / 2),
            y: values.map(item => Number(item.fraction)),
            width: values.map(item => Number(item.bin_right) - Number(item.bin_left)),
            marker: {{color}}, opacity: 0.58,
            customdata: values.map(item => [item.bin_left, item.bin_right, count(item.count)]),
            hovertemplate: group + '<br>phyloP: %{{customdata[0]:.3f}} to %{{customdata[1]:.3f}}<br>N: %{{customdata[2]}}<br>Fraction: %{{y:.2%}}<extra></extra>',
          }});
          const summary = values[0];
          traces.push({{
            type: 'box', name: group, legendgroup: group, showlegend: false,
            q1: [summary.q1], median: [summary.median], q3: [summary.q3],
            lowerfence: [summary.lower_whisker], upperfence: [summary.upper_whisker],
            marker: {{color}}, boxpoints: false, xaxis: 'x2', yaxis: 'y2',
          }});
        }});
        const plot = document.getElementById(config.viewId + '-detail-plot'); plot.hidden = false;
        Plotly.react(plot, traces, {{
          title: 'phyloP100way by ALT-observation status', template: 'plotly_white', height: 380,
          margin: {{l: 65, r: 25, t: 50, b: 58}}, barmode: 'overlay', boxmode: 'group',
          xaxis: {{title: 'phyloP100way', domain: [0, 0.68]}}, yaxis: {{title: 'Fraction per bin', tickformat: '.0%'}},
          xaxis2: {{domain: [0.78, 1], anchor: 'y2'}}, yaxis2: {{title: 'phyloP100way', anchor: 'x2'}},
          legend: {{orientation: 'h', y: 1.12}},
        }}, {{responsive: true}});
        document.getElementById(config.viewId + '-detail-table').innerHTML = twoByTwoTable(row);
      }}

      function render() {{
        const rows = currentRows();
        renderForest(rows);
        renderResultsTable(rows);
        if (!rows.some(row => row.strategy === strategySelect.value) && rows.length) strategySelect.value = rows[0].strategy;
        const selected = rows.find(row => row.strategy === strategySelect.value);
        const status = document.getElementById(config.viewId + '-status');
        status.innerHTML = selected && selected.status !== 'estimated' ? `<strong>${{strategyLabels[selected.strategy]}}:</strong> ${{statusText(selected)}}` : '';
        status.hidden = !status.innerHTML;
        if (modeSelect.value === 'unadjusted') renderUnadjusted(selected);
        else if (modeSelect.value === 'fixed') renderFixed(selected);
        else renderContinuous(selected);
      }}
      [modeSelect, variantSelect, targetContextSelect].forEach(select => select.addEventListener('change', () => {{
        refreshConsequences();
        render();
      }}));
      [consequenceSelect, strategySelect].forEach(select => select.addEventListener('change', render));
      refreshConsequences();
      render();
    }})();
    </script>
    """
