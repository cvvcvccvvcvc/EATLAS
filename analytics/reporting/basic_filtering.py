"""Interactive report view for simple candidate support thresholds."""

from __future__ import annotations

import json

from analytics.analyses.basic_filtering import BasicFilteringAnalysis, FILTER_OPTIONS
from analytics.analyses.conservation_validation import (
    TARGET_CONTEXT_OPTIONS,
    VARIANT_TYPE_OPTIONS,
)
from analytics.vep.consequences import (
    VALIDATION_CONSEQUENCE_OPTIONS as CONSEQUENCE_OPTIONS,
)
from .components import dataframe_records, strategy_label


def build_basic_filtering_sections(analysis: BasicFilteringAnalysis) -> list[str]:
    return [
        "<h2>Basic Filtering</h2>",
        "<p class='lead'>Support thresholds provide a deliberately simple first-pass "
        "filter. Candidate retention and gnomAD overlap use unique normalized alleles; "
        "failed gnomAD lookups are excluded only from the overlap denominator. ClinVar "
        "uses at most 20 retention-spanning thresholds per curve; non-estimable OR "
        "points are omitted.</p>",
        "<p>Ortholog and genus scores are computed within the selected strategy and "
        "use the maximum across overlapping target-gene contexts for one normalized "
        "allele. Strategy support counts distinct registered strategies; the two "
        "minimap2 presets are separate in this view. ClinVar OR &gt; 1 indicates "
        "relative enrichment of B/LB over P/LP alleles among retained calls.</p>",
        basic_filtering_view(analysis),
    ]


def basic_filtering_view(analysis: BasicFilteringAnalysis) -> str:
    strategies = (
        analysis.candidate_curves["strategy"].drop_duplicates().astype(str).tolist()
        if not analysis.candidate_curves.empty
        else []
    )
    payload = {
        "viewId": "basic-filtering",
        "filters": [
            {"key": key, "label": label}
            for key, label, _column in FILTER_OPTIONS
        ],
        "strategies": [
            {"key": value, "label": strategy_label(value)} for value in strategies
        ],
        "variantTypes": [
            {"key": key, "label": label} for key, label in VARIANT_TYPE_OPTIONS
        ],
        "modes": [
            {"key": "unadjusted", "label": "Unadjusted"},
            {"key": "fixed", "label": "phyloP fixed bands"},
        ],
        "targetContexts": [
            {"key": key, "label": label} for key, label in TARGET_CONTEXT_OPTIONS
        ],
        "consequences": [
            {"key": key, "label": label} for key, label in CONSEQUENCE_OPTIONS
        ],
        "candidate": dataframe_records(analysis.candidate_curves),
        "clinvar": dataframe_records(analysis.clinvar_curves),
    }
    payload_json = json.dumps(payload, separators=(",", ":"), allow_nan=False).replace(
        "</", "<\\/"
    )
    return f"""
    <div class="analysis-controls" id="basic-filtering-controls">
      <label>Filter<select data-role="filter"></select></label>
      <label>Strategy<select data-role="strategy"></select></label>
      <label>Variant type<select data-role="variant-type"></select></label>
      <label>ClinVar analysis<select data-role="mode"></select></label>
      <label>Target context<select data-role="target-context"></select></label>
      <label>Consequence subset<select data-role="consequence"></select></label>
    </div>
    <div id="basic-filtering-status" class="analysis-note" hidden></div>
    <div id="basic-filtering-retention" class="analysis-plot"></div>
    <div id="basic-filtering-gnomad" class="analysis-plot"></div>
    <div id="basic-filtering-clinvar" class="analysis-plot"></div>
    <script>
    (() => {{
      const config = {payload_json};
      const controls = document.getElementById(config.viewId + '-controls');
      const selects = Object.fromEntries(
        ['filter', 'strategy', 'variant-type', 'mode', 'target-context', 'consequence']
          .map(role => [role, controls.querySelector(`[data-role="${{role}}"]`)])
      );
      const optionMap = values => Object.fromEntries(values.map(value => [value.key, value.label]));
      const labels = {{
        filter: optionMap(config.filters), strategy: optionMap(config.strategies),
        variant: optionMap(config.variantTypes), mode: optionMap(config.modes),
        context: optionMap(config.targetContexts), consequence: optionMap(config.consequences),
      }};
      const addOptions = (select, values) => values.forEach(value => {{
        const option = document.createElement('option');
        option.value = value.key; option.textContent = value.label; select.appendChild(option);
      }});
      addOptions(selects.filter, config.filters);
      addOptions(selects.strategy, config.strategies);
      addOptions(selects['variant-type'], config.variantTypes);
      addOptions(selects.mode, config.modes);
      addOptions(selects['target-context'], config.targetContexts);
      addOptions(selects.consequence, config.consequences);
      selects.filter.value = 'ortholog';
      selects['variant-type'].value = 'snv';
      selects.mode.value = 'unadjusted';
      selects['target-context'].value = 'all';
      selects.consequence.value = 'all';

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
      const candidateRows = () => config.candidate.filter(row =>
        row.filter_key === selects.filter.value
        && row.strategy === selects.strategy.value
        && row.variant_type === selects['variant-type'].value
      ).sort((left, right) => Number(left.threshold) - Number(right.threshold));
      const clinvarRows = () => config.clinvar.filter(row =>
        row.filter_key === selects.filter.value
        && row.strategy === selects.strategy.value
        && row.variant_type === selects['variant-type'].value
        && row.mode === selects.mode.value
        && row.target_context === selects['target-context'].value
        && row.consequence === selects.consequence.value
      ).sort((left, right) => Number(left.threshold) - Number(right.threshold));

      function renderCandidate(rows) {{
        const common = {{
          template: 'plotly_white', height: 350,
          margin: {{l: 72, r: 25, t: 52, b: 58}},
          xaxis: {{title: 'Minimum support threshold', type: 'linear', rangemode: 'tozero'}},
        }};
        Plotly.react(config.viewId + '-retention', [{{
          type: 'scatter', mode: 'lines', line: {{width: 3, color: '#356d8f'}},
          x: rows.map(row => row.threshold), y: rows.map(row => row.retained_fraction),
          customdata: rows.map(row => [count(row.retained_variant_count), count(row.total_variant_count)]),
          hovertemplate: 'Threshold ≥ %{{x}}<br>Retained: %{{customdata[0]}} / %{{customdata[1]}}<br>Fraction: %{{y:.2%}}<extra></extra>',
        }}], {{...common, title: 'Candidate alleles retained', yaxis: {{title: 'Fraction retained', tickformat: '.0%', range: [0, 1]}}}}, {{responsive: true}});
        Plotly.react(config.viewId + '-gnomad', [{{
          type: 'scatter', mode: 'lines+markers', line: {{width: 3, color: '#2ca25f'}},
          marker: {{size: 5}}, x: rows.map(row => row.threshold),
          y: rows.map(row => row.gnomad_found_fraction),
          customdata: rows.map(row => [count(row.gnomad_found_count), count(row.gnomad_eligible_count), count(row.gnomad_lookup_failed_count)]),
          hovertemplate: 'Threshold ≥ %{{x}}<br>gnomAD: %{{customdata[0]}} / %{{customdata[1]}}<br>Found fraction: %{{y:.2%}}<br>Failed lookups excluded: %{{customdata[2]}}<extra></extra>',
        }}], {{...common, title: 'Exact alleles found in gnomAD', yaxis: {{title: 'gnomAD found fraction', tickformat: '.0%', range: [0, 1]}}}}, {{responsive: true}});
      }}

      function renderClinVar(rows) {{
        const plotted = rows.map(row => {{
          let value = number(row.result_or);
          if (value === null && row.result_or === 'inf' && finite(row.ci_low) && finite(row.ci_high))
            value = Math.sqrt(Number(row.ci_low) * Number(row.ci_high));
          return {{row, value}};
        }}).filter(item => item.value !== null && item.value > 0 && finite(item.row.ci_low) && finite(item.row.ci_high));
        const trace = {{
          type: 'scatter', mode: 'lines+markers', line: {{width: 2, color: '#7a3e9d'}},
          marker: {{size: 7}}, x: plotted.map(item => item.row.threshold),
          y: plotted.map(item => item.value),
          error_y: {{type: 'data', symmetric: false,
            array: plotted.map(item => Number(item.row.ci_high) - item.value),
            arrayminus: plotted.map(item => item.value - Number(item.row.ci_low))}},
          customdata: plotted.map(item => [
            fmt(item.row.result_or), fmt(item.row.ci_low), fmt(item.row.ci_high),
            fmt(item.row.result_p), fmt(item.row.result_q), count(item.row.usable_rows),
            count(item.row.benign_observed), count(item.row.pathogenic_observed),
          ]),
          hovertemplate: 'Threshold ≥ %{{x}}<br>OR: %{{customdata[0]}} [%{{customdata[1]}}–%{{customdata[2]}}]<br>p / FDR q: %{{customdata[3]}} / %{{customdata[4]}}<br>ClinVar N: %{{customdata[5]}}<br>Observed B/LB / P/LP: %{{customdata[6]}} / %{{customdata[7]}}<extra></extra>',
        }};
        const annotations = plotted.length ? [] : [{{
          text: 'No estimable finite confidence interval for this selection.',
          x: 0.5, y: 0.5, xref: 'paper', yref: 'paper', showarrow: false,
        }}];
        Plotly.react(config.viewId + '-clinvar', plotted.length ? [trace] : [], {{
          title: 'ClinVar B/LB versus P/LP association along the filter threshold',
          template: 'plotly_white', height: 390, margin: {{l: 76, r: 25, t: 55, b: 58}},
          xaxis: {{title: 'Minimum support threshold', type: 'linear', rangemode: 'tozero'}},
          yaxis: {{title: 'Odds ratio (log scale)', type: 'log'}},
          shapes: [{{type: 'line', x0: 0, x1: 1, xref: 'paper', y0: 1, y1: 1, line: {{dash: 'dash', color: '#8c8c8c'}}}}],
          annotations,
        }}, {{responsive: true}});
        const hidden = rows.filter(row => row.status === 'not_estimable').length;
        const status = document.getElementById(config.viewId + '-status');
        status.textContent = hidden ? `${{hidden}} threshold point(s) are not estimable and are omitted from the OR line.` : '';
        status.hidden = !status.textContent;
      }}

      function refreshConsequences() {{
        const available = new Set(config.clinvar.filter(row =>
          row.filter_key === selects.filter.value
          && row.strategy === selects.strategy.value
          && row.variant_type === selects['variant-type'].value
          && row.mode === selects.mode.value
          && row.target_context === selects['target-context'].value
          && row.status !== 'not_estimable'
        ).map(row => row.consequence));
        const previous = selects.consequence.value || 'all';
        selects.consequence.replaceChildren();
        const options = config.consequences.filter(value => available.has(value.key));
        addOptions(selects.consequence, options);
        if (available.has(previous)) selects.consequence.value = previous;
        else if (available.has('all')) selects.consequence.value = 'all';
        else if (options.length) selects.consequence.value = options[0].key;
      }}

      function ensureAvailableFilter() {{
        const available = new Set(config.candidate.filter(row =>
          row.strategy === selects.strategy.value
          && row.variant_type === selects['variant-type'].value
        ).map(row => row.filter_key));
        [...selects.filter.options].forEach(option => option.disabled = !available.has(option.value));
        if (!available.has(selects.filter.value) && available.size)
          selects.filter.value = config.filters.find(value => available.has(value.key)).key;
      }}

      function render() {{
        renderCandidate(candidateRows());
        renderClinVar(clinvarRows());
      }}
      Object.entries(selects).forEach(([role, select]) => select.addEventListener('change', () => {{
        ensureAvailableFilter();
        if (role !== 'consequence') refreshConsequences();
        render();
      }}));
      ensureAvailableFilter();
      refreshConsequences();
      render();
    }})();
    </script>
    """
