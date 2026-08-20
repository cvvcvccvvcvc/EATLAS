"""Report sections for minimap2 asm10/asm20 concordance."""

from __future__ import annotations

import json

from analytics.analyses.conservation_validation import VARIANT_TYPE_OPTIONS
from analytics.analyses.minimap_concordance import (
    GROUP_OPTIONS,
    MinimapConcordanceAnalysis,
)
from .components import dataframe_records
from .conservation import clinvar_association_view


def build_minimap_concordance_sections(
    analysis: MinimapConcordanceAnalysis,
) -> list[str]:
    if not analysis.available or analysis.validation is None:
        return [
            "<h2>Minimap2 Concordance</h2>",
            f"<p class='analysis-note'>{analysis.reason}</p>",
        ]
    labels = dict(GROUP_OPTIONS)
    return [
        "<h2>Minimap2 Concordance</h2>",
        "<p class='lead'>This diagnostic separates union, intersection, and "
        "preset-specific calls. In the ClinVar models, ‘only’ groups use genes "
        "eligible for both presets, so absence from the other preset is not confused "
        "with a strategy that was not evaluable for that gene.</p>",
        minimap_candidate_view(analysis),
        "<h3>ClinVar association</h3>",
        "<p>OR &gt; 1 indicates relative enrichment of B/LB over P/LP alleles "
        "among calls in the selected group.</p>",
        clinvar_association_view(
            analysis.validation,
            view_id="minimap-concordance-clinvar",
            strategy_labels=labels,
        ),
    ]


def minimap_candidate_view(analysis: MinimapConcordanceAnalysis) -> str:
    available_variant_types = (
        set(analysis.candidate_summary["variant_type"].astype(str))
        if "variant_type" in analysis.candidate_summary
        else set()
    )
    payload = {
        "viewId": "minimap-concordance-candidates",
        "groups": [{"key": key, "label": label} for key, label in GROUP_OPTIONS],
        "variantTypes": [
            {"key": key, "label": label} for key, label in VARIANT_TYPE_OPTIONS
            if key in available_variant_types
        ],
        "rows": dataframe_records(analysis.candidate_summary),
    }
    payload_json = json.dumps(payload, separators=(",", ":"), allow_nan=False).replace(
        "</", "<\\/"
    )
    return f"""
    <div class="analysis-controls analysis-controls-single" id="minimap-concordance-candidates-controls">
      <label>Variant type<select data-role="variant-type"></select></label>
    </div>
    <div id="minimap-concordance-candidates-plot" class="analysis-plot"></div>
    <div id="minimap-concordance-candidates-table"></div>
    <script>
    (() => {{
      const config = {payload_json};
      const select = document.querySelector('#' + config.viewId + '-controls [data-role="variant-type"]');
      const groupLabels = Object.fromEntries(config.groups.map(value => [value.key, value.label]));
      config.variantTypes.forEach(value => {{
        const option = document.createElement('option'); option.value = value.key;
        option.textContent = value.label; select.appendChild(option);
      }});
      select.value = config.variantTypes.some(value => value.key === 'snv')
        ? 'snv' : config.variantTypes[0]?.key || '';
      const count = value => Math.round(Number(value || 0)).toLocaleString('en-US').replaceAll(',', ' ');
      const percent = value => value === null ? 'NA' : (Number(value) * 100).toFixed(1) + '%';
      const cell = value => `<td>${{value}}</td>`;
      function render() {{
        const rows = config.rows.filter(row => row.variant_type === select.value);
        const labels = rows.map(row => groupLabels[row.group_key] || row.group_key);
        const traces = [
          {{type: 'bar', name: 'Alleles / union', x: labels, y: rows.map(row => row.allele_fraction), marker: {{color: '#356d8f'}},
            customdata: rows.map(row => [count(row.variant_count)]),
            hovertemplate: '%{{x}}<br>Alleles: %{{customdata[0]}}<br>Fraction of union: %{{y:.2%}}<extra></extra>'}},
          {{type: 'bar', name: 'Found in gnomAD', x: labels, y: rows.map(row => row.gnomad_found_fraction), marker: {{color: '#2ca25f'}},
            customdata: rows.map(row => [count(row.gnomad_found_count), count(row.gnomad_eligible_count), count(row.gnomad_lookup_failed_count)]),
            hovertemplate: '%{{x}}<br>gnomAD: %{{customdata[0]}} / %{{customdata[1]}}<br>Found fraction: %{{y:.2%}}<br>Failed lookups: %{{customdata[2]}}<extra></extra>'}},
        ];
        Plotly.react(config.viewId + '-plot', traces, {{
          title: 'minimap2 preset concordance and gnomAD overlap', template: 'plotly_white',
          height: 410, barmode: 'group', margin: {{l: 70, r: 25, t: 55, b: 95}},
          yaxis: {{title: 'Fraction', tickformat: '.0%', range: [0, 1]}},
          xaxis: {{title: '', tickangle: -20}}, legend: {{orientation: 'h', y: 1.12}},
        }}, {{responsive: true}});
        const body = rows.map(row => `<tr>${{cell(groupLabels[row.group_key] || row.group_key)}}${{cell(count(row.variant_count))}}${{cell(percent(row.allele_fraction))}}${{cell(count(row.gnomad_found_count) + ' / ' + count(row.gnomad_eligible_count))}}${{cell(percent(row.gnomad_found_fraction))}}${{cell(count(row.gnomad_lookup_failed_count))}}</tr>`).join('');
        document.getElementById(config.viewId + '-table').innerHTML = `<table><thead><tr><th>Call group</th><th>Alleles</th><th>Of union</th><th>gnomAD found / eligible</th><th>gnomAD found</th><th>Lookup failed</th></tr></thead><tbody>${{body}}</tbody></table>`;
      }}
      select.addEventListener('change', render); render();
    }})();
    </script>
    """
