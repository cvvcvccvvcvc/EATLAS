"""Independent descriptive and ClinVar controls for candidate filtering."""

from __future__ import annotations

import json

from analytics.analyses.basic_filtering import (
    BasicFilteringAnalysis,
    FILTER_OPTIONS,
    UNION_STRATEGY,
)
from analytics.analyses.conservation_validation import TARGET_CONTEXT_OPTIONS, VARIANT_TYPE_OPTIONS
from analytics.vep.consequences import VALIDATION_CONSEQUENCE_OPTIONS as CONSEQUENCE_OPTIONS
from .components import dataframe_records, strategy_label


def build_basic_filtering_sections(analysis: BasicFilteringAnalysis) -> list[str]:
    return ["<h2>Basic Filtering</h2>", basic_filtering_view(analysis)]


def basic_filtering_view(analysis: BasicFilteringAnalysis) -> str:
    strategies = (
        analysis.candidate_curves["strategy"].drop_duplicates().astype(str).tolist()
        if not analysis.candidate_curves.empty
        else []
    )
    payload = {
        "filters": [{"key": key, "label": label} for key, label, _column in FILTER_OPTIONS],
        "strategies": [
            {"key": key, "label": strategy_label(key)}
            for key in strategies
            if key != UNION_STRATEGY
        ],
        "variantTypes": [{"key": key, "label": label} for key, label in VARIANT_TYPE_OPTIONS],
        "modes": [
            {"key": "unadjusted", "label": "Unadjusted"},
            {"key": "fixed", "label": "phyloP-adjusted"},
        ],
        "contexts": [{"key": key, "label": label} for key, label in TARGET_CONTEXT_OPTIONS],
        "consequences": [{"key": key, "label": label} for key, label in CONSEQUENCE_OPTIONS],
        "candidate": dataframe_records(analysis.candidate_curves),
        "clinvar": dataframe_records(analysis.clinvar_curves),
    }
    payload_json = json.dumps(payload, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    return (
        """
    <div class="analysis-controls" id="basic-filtering-global-controls">
      <label>Filter<select data-role="filter"></select></label>
      <label>Variant type<select data-role="variant-type"></select></label>
    </div>
    <div class="analysis-controls analysis-controls-single" id="basic-filtering-candidate-controls">
      <label>Strategies<select data-role="candidate-strategy"></select></label>
    </div>
    <div id="basic-filtering-candidate-status" class="analysis-note" hidden></div>
    <div id="basic-filtering-retention" class="analysis-plot"></div>
    <div id="basic-filtering-gnomad" class="analysis-plot"></div>
    <h3>ClinVar association</h3>
    <div class="analysis-controls" id="basic-filtering-clinvar-controls">
      <label>Strategy<select data-role="clinvar-strategy"></select></label>
      <label>Analysis<select data-role="mode"></select></label>
      <label>Target context<select data-role="target-context"></select></label>
      <label>Consequence<select data-role="consequence"></select></label>
    </div>
    <div id="basic-filtering-status" class="analysis-note" hidden></div>
    <div id="basic-filtering-clinvar" class="analysis-plot"></div>
    <script>(() => {
      const config = """
        + payload_json
        + """;
      const roles = ['filter', 'variant-type', 'candidate-strategy', 'clinvar-strategy', 'mode', 'target-context', 'consequence'];
      const selects = Object.fromEntries(roles.map(role => [role,
        document.querySelector('#tab-basic-filtering [data-role="' + role + '"]')]));
      const addOptions = (select, options) => options.forEach(item => select.add(new Option(item.label, item.key)));
      const union = {key: 'union', label: 'Union (any strategy)'};
      addOptions(selects.filter, config.filters);
      addOptions(selects['variant-type'], config.variantTypes);
      addOptions(selects['candidate-strategy'], [{key: 'compare', label: 'Compare all strategies'}, ...config.strategies, union]);
      addOptions(selects['clinvar-strategy'], [...config.strategies, union]);
      addOptions(selects.mode, config.modes);
      addOptions(selects['target-context'], config.contexts);
      addOptions(selects.consequence, config.consequences);
      const colors = ['#356d8f', '#e1812c', '#499c68', '#9b5c98', '#c24d4d', '#79706e'];
      const labels = Object.fromEntries([...config.strategies, union].map(item => [item.key, item.label]));
      const color = strategy => strategy === 'union' ? '#242424' : colors[config.strategies.findIndex(item => item.key === strategy) % colors.length];
      const finite = value => value !== null && value !== '' && Number.isFinite(Number(value));
      const count = value => Number(value || 0).toLocaleString('en-US');
      const fmt = value => value === 'inf' ? '∞' : value === '-inf' ? '−∞' : finite(value) ? Number(value).toPrecision(3) : 'NA';
      const atMost = () => selects.filter.value === 'aligned_max';
      const symbol = () => atMost() ? '≤' : '≥';
      const supported = () => selects['variant-type'].value === 'snv' || ['ortholog', 'strategy'].includes(selects.filter.value);
      const commonRows = row => row.filter_key === selects.filter.value && row.variant_type === selects['variant-type'].value;
      const xaxis = () => ({title: {text: 'Threshold (' + symbol() + ' N)'}, type: 'linear', rangemode: 'tozero'});
      function note(id, text) {
        const element = document.getElementById(id);
        element.textContent = text;
        element.hidden = !text;
      }
      function renderCandidate() {
        const selection = selects['candidate-strategy'].value;
        const strategies = selection === 'compare' ? config.strategies.map(item => item.key) : [selection];
        const groups = strategies.map(strategy => ({strategy, rows: config.candidate.filter(row =>
          supported() && commonRows(row) && row.strategy === strategy).sort((a, b) => a.threshold - b.threshold)}));
        const traces = (gnomad) => groups.filter(group => group.rows.length).map(({strategy, rows}) => ({
          type: 'scatter', mode: 'lines', name: labels[strategy], line: {color: color(strategy), width: 2, shape: 'hv'},
          connectgaps: false, x: rows.map(row => row.threshold),
          y: rows.map(row => gnomad ? row.gnomad_found_fraction : row.retained_fraction),
          customdata: rows.map(row => gnomad
            ? [count(row.gnomad_found_count), count(row.gnomad_eligible_count), count(row.gnomad_lookup_failed_count)]
            : [count(row.retained_variant_count), count(row.total_variant_count)]),
          hovertemplate: 'Threshold ' + symbol() + ' %{x}<br>' + (gnomad
            ? 'gnomAD: %{customdata[0]} / %{customdata[1]}<br>Found: %{y:.2%}<br>Failed lookups excluded: %{customdata[2]}'
            : 'Retained: %{customdata[0]} / %{customdata[1]}<br>Fraction: %{y:.2%}') + '<extra>%{fullData.name}</extra>',
        }));
        for (const [id, gnomad, title, ylabel] of [
          ['retention', false, 'Candidate alleles retained', 'Fraction retained'],
          ['gnomad', true, 'Exact alleles found in gnomAD', 'gnomAD found fraction'],
        ]) {
          Plotly.react('basic-filtering-' + id, traces(gnomad), {
            title: {text: title}, template: 'plotly_white', height: 370, margin: {l: 72, r: 25, t: 55, b: 60},
            xaxis: xaxis(), yaxis: {title: {text: ylabel}, tickformat: '.0%', range: [0, 1]},
            legend: {orientation: 'h'},
          }, {responsive: true});
        }
        note('basic-filtering-candidate-status', !supported() ? 'This filter is available for SNVs only.'
          : groups.every(group => !group.rows.length) ? 'No candidate alleles for this selection.' : '');
      }
      function renderClinVar() {
        const rows = config.clinvar.filter(row => supported() && commonRows(row)
          && row.strategy === selects['clinvar-strategy'].value && row.mode === selects.mode.value
          && row.target_context === selects['target-context'].value && row.consequence === selects.consequence.value)
          .sort((a, b) => a.threshold - b.threshold);
        const estimable = row => row.status === 'estimated' && finite(row.result_or) && Number(row.result_or) > 0
          && finite(row.ci_low) && finite(row.ci_high);
        const hover = row => [fmt(row.result_or), fmt(row.ci_low), fmt(row.ci_high), fmt(row.result_p), fmt(row.result_q),
          count(row.benign_observed), count(row.pathogenic_observed), count(row.benign_not_observed), count(row.pathogenic_not_observed)];
        const hovertemplate = 'Threshold ' + symbol() + ' %{x}<br>OR: %{customdata[0]} [%{customdata[1]}–%{customdata[2]}]'
          + '<br>p / FDR q: %{customdata[3]} / %{customdata[4]}<br>Retained B/LB / P/LP: %{customdata[5]} / %{customdata[6]}'
          + '<br>Not retained B/LB / P/LP: %{customdata[7]} / %{customdata[8]}<extra></extra>';
        const traces = [{
          type: 'scatter', mode: 'lines+markers', line: {color: '#7a3e9d', width: 2, shape: 'hv'}, marker: {size: 6},
          connectgaps: false, x: rows.map(row => row.threshold), y: rows.map(row => estimable(row) ? row.result_or : null),
          error_y: {type: 'data', symmetric: false,
            array: rows.map(row => estimable(row) ? row.ci_high - row.result_or : null),
            arrayminus: rows.map(row => estimable(row) ? row.result_or - row.ci_low : null)},
          customdata: rows.map(hover), hovertemplate,
        }];
        const boundary = rows.filter(row => row.status === 'estimated' && (row.result_or === 'inf' || row.result_or === 0));
        const bounds = rows.flatMap(row => [row.ci_low, row.ci_high]).filter(value => finite(value) && Number(value) > 0).map(Number);
        if (boundary.length && bounds.length) traces.push({
          type: 'scatter', mode: 'markers',
          x: boundary.map(row => row.threshold),
          y: boundary.map(row => row.result_or === 'inf' ? Math.max(...bounds) * 1.5 : Math.min(...bounds) / 1.5),
          marker: {color: '#7a3e9d', size: 11, symbol: boundary.map(row => row.result_or === 'inf' ? 'triangle-up' : 'triangle-down')},
          customdata: boundary.map(hover), hovertemplate,
        });
        const any = rows.some(estimable) || boundary.length;
        Plotly.react('basic-filtering-clinvar', traces, {
          title: {text: 'ClinVar B/LB versus P/LP association'}, template: 'plotly_white', height: 410,
          margin: {l: 78, r: 25, t: 55, b: 60}, showlegend: false, xaxis: xaxis(),
          yaxis: {title: {text: 'Odds ratio (log scale)'}, type: 'log'},
          shapes: [{type: 'line', x0: 0, x1: 1, xref: 'paper', y0: 1, y1: 1, line: {dash: 'dash', color: '#8c8c8c'}}],
          annotations: any ? [] : [{text: 'No estimable odds ratio for this selection.', x: 0.5, y: 0.5, xref: 'paper', yref: 'paper', showarrow: false}],
        }, {responsive: true});
        const unavailable = rows.filter(row => row.status !== 'estimated');
        const reasons = [...new Set(unavailable.map(row => row.reason).filter(Boolean))];
        note('basic-filtering-status', !supported() ? 'This filter is available for SNVs only.'
          : [unavailable.length ? unavailable.length + ' threshold(s) without an OR estimate. ' + reasons.join(' ') : '', boundary.length ? 'Triangles mark OR = 0 or ∞ at the plot boundary.' : ''].filter(Boolean).join(' '));
      }
      ['filter', 'variant-type'].forEach(role => selects[role].addEventListener('change', () => {
        renderCandidate(); renderClinVar();
      }));
      selects['candidate-strategy'].addEventListener('change', renderCandidate);
      ['clinvar-strategy', 'mode', 'target-context', 'consequence'].forEach(role => selects[role].addEventListener('change', renderClinVar));
      renderCandidate(); renderClinVar();
    })();</script>
    """
    )
