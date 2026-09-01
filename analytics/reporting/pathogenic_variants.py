"""Report sections for focused characterization of P/LP ClinVar hits."""

from __future__ import annotations

import html
import json
from urllib.parse import quote

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from analytics.analyses.pathogenic_variants import PathogenicVariantAnalysis
from .components import (
    compact_figure,
    dataframe_records,
    fig_html,
    format_int,
    metric_cards,
    strategy_label,
)
from .config import (
    CONSEQUENCE_GROUP_COLORS,
    CONSEQUENCE_GROUP_ORDER,
    REVIEW_STAR_COLORS,
    REVIEW_STAR_ORDER,
)
from .variant_profile import group_consequence_counts


def build_pathogenic_variant_sections(analysis: PathogenicVariantAnalysis) -> list[str]:
    variants = analysis.variants
    if variants.empty:
        return ["<h2>Pathogenic ClinVar Hits</h2>", "<p>No P/LP candidate alleles were found.</p>"]
    stars = pd.to_numeric(variants["clinvar_review_stars"], errors="coerce")
    named = analysis.condition_counts
    named = named[named["cohort"].eq("gaph") & named["variant_count"].gt(0)]
    sections = [
        "<h2>Pathogenic ClinVar Hits</h2>",
        metric_cards(
            [
                ("Unique P/LP alleles", format_int(len(variants))),
                ("With ≥2 review stars", format_int(stars.ge(2).sum())),
                ("Named conditions", format_int(named["condition_key"].nunique())),
            ]
        ),
    ]
    for title, figure in (
        ("ClinVar review stars", pathogenic_star_figure(analysis.star_counts)),
        ("Molecular effect", pathogenic_consequence_figure(analysis.consequence_counts)),
        ("Supporting orthologs among P/LP hits", pathogenic_support_figure(analysis.support_rows)),
    ):
        sections.extend(
            [
                f"<h3>{title}</h3>",
                fig_html(figure) if figure is not None else "<p>No eligible observations.</p>",
            ]
        )
    sections.extend(
        [
            "<h3>Associated conditions</h3>",
            pathogenic_condition_view(analysis.condition_counts),
            "<h3>Variant details</h3>",
            f'<p><a download href="../derived/{quote(analysis.variants_path.name)}">Download complete P/LP TSV</a></p>',
            pathogenic_variant_table_html(variants),
        ]
    )
    return sections


def pathogenic_star_figure(counts: pd.DataFrame):
    if counts.empty:
        return None
    shown = counts.copy()
    shown["Strategy"] = shown["strategy"].map(strategy_label)
    shown["Review stars"] = shown["clinvar_review_stars"].astype(str)
    shown.loc[~shown["Review stars"].isin({"0", "1", "2", "3", "4"}), "Review stars"] = (
        "Unmapped"
    )
    present_stars = [star for star in REVIEW_STAR_ORDER if star in set(shown["Review stars"])]
    strategy_order = (
        shown.groupby("Strategy")["variant_count"].sum().sort_values(ascending=False).index.tolist()
    )
    figure = px.bar(
        shown,
        x="Strategy",
        y="variant_count",
        color="Review stars",
        barmode="stack",
        title="P/LP alleles by ClinVar review stars",
        category_orders={
            "Strategy": strategy_order,
            "Review stars": present_stars,
        },
        color_discrete_map=REVIEW_STAR_COLORS,
        labels={
            "variant_count": "Unique P/LP alleles",
        },
    )
    compact_figure(figure, height=390)
    return figure


def pathogenic_consequence_figure(raw_counts: pd.DataFrame):
    counts = group_consequence_counts(raw_counts)
    if counts.empty:
        return None
    strategy_order = (
        counts.groupby("Strategy", observed=True)["Variant_Count"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    )
    figure = px.bar(
        counts,
        x="Strategy",
        y="Variant_Count",
        color="Consequence group",
        barmode="stack",
        title="VEP consequence groups among P/LP candidate hits",
        category_orders={
            "Strategy": strategy_order,
            "Consequence group": CONSEQUENCE_GROUP_ORDER,
        },
        color_discrete_map=CONSEQUENCE_GROUP_COLORS,
        labels={"Strategy": "", "Variant_Count": "Allele–target-gene observations"},
    )
    compact_figure(figure, height=340)
    return figure


def pathogenic_support_figure(rows: pd.DataFrame):
    if rows.empty:
        return None
    figure = go.Figure()
    for strategy in sorted(rows["strategy"].unique(), key=strategy_label):
        values = rows[rows["strategy"].eq(strategy)]
        counts = values["alt_support_ortholog_count"].to_numpy(dtype=int)
        common = dict(
            y=np.log10(counts),
            name=strategy_label(strategy),
            customdata=values[
                [
                    "variant_key",
                    "gene_id",
                    "alt_support_ortholog_count",
                    "site_aligned_ortholog_count",
                    "alt_support_family_count",
                ]
            ],
            marker_color="#356d8f",
            line_color="#356d8f",
            hovertemplate=(
                "%{customdata[0]}<br>Gene: %{customdata[1]}<br>Supporting orthologs: "
                "%{customdata[2]}<br>Site-aligned orthologs: %{customdata[3]}<br>"
                "Supporting families: %{customdata[4]}<extra>%{fullData.name}</extra>"
            ),
        )
        if len(set(counts)) >= 3:
            figure.add_trace(
                go.Violin(
                    **common,
                    points="all",
                    jitter=0.25,
                    pointpos=0,
                    box_visible=True,
                    spanmode="hard",
                    hoveron="points",
                )
            )
        else:
            figure.add_trace(go.Box(**common, boxpoints="all", jitter=0.25, pointpos=0))
    maximum = int(rows["alt_support_ortholog_count"].max())
    ticks = [
        multiplier * 10**power
        for power in range(len(str(maximum)))
        for multiplier in (1, 2, 5)
        if multiplier * 10**power <= maximum
    ]
    figure.update_layout(
        showlegend=False,
        yaxis=dict(
            title="Supporting orthologs (log scale)",
            tickmode="array",
            tickvals=np.log10(ticks).tolist(),
            ticktext=[str(value) for value in ticks],
        ),
    )
    compact_figure(figure, height=440)
    return figure


def pathogenic_condition_view(counts: pd.DataFrame) -> str:
    if counts.empty:
        return "<p>No ClinVar condition data.</p>"
    strategies = counts.loc[counts["cohort"].eq("gaph"), "strategy"].drop_duplicates().tolist()
    payload = json.dumps(
        {
            "rows": dataframe_records(counts),
            "strategies": [{"key": key, "label": strategy_label(key)} for key in strategies],
        },
        separators=(",", ":"),
        allow_nan=False,
    ).replace("</", "<\\/")
    return (
        """
    <div class="analysis-controls" id="pathogenic-conditions-controls">
      <label>Strategy<select data-role="strategy"></select></label>
      <label>Variant type<select data-role="type"><option value="all">SNV + INDEL</option><option value="snv">SNV</option><option value="indel">INDEL</option></select></label>
      <label>ClinVar background<select data-role="background"><option value="target">Matching target regions</option><option value="global">Whole GRCh38 VCF</option></select></label>
      <label>Find condition<input data-role="search" type="search" placeholder="Condition name"></label>
    </div>
    <div id="pathogenic-conditions-plot" class="analysis-plot"></div>
    <h4>ClinVar condition distribution</h4>
    <div id="pathogenic-clinvar-distribution-plot" class="analysis-plot"></div>
    <script>(() => {
      const config = """
        + payload
        + """;
      const controls = document.getElementById('pathogenic-conditions-controls');
      const select = role => controls.querySelector('[data-role="' + role + '"]');
      config.strategies.forEach(item => select('strategy').add(new Option(item.label, item.key)));
      function render() {
        const rows = config.rows.filter(row => row.variant_type === select('type').value);
        const gaph = rows.filter(row => row.cohort === 'gaph' && row.strategy === select('strategy').value);
        const background = rows.filter(row => row.cohort === select('background').value &&
          (row.cohort === 'global' || row.strategy === select('strategy').value));
        const a = new Map(gaph.map(row => [row.condition_key, row]));
        const b = new Map(background.map(row => [row.condition_key, row]));
        const fraction = row => row && row.total_variant_count ? row.variant_count / row.total_variant_count : 0;
        const query = select('search').value.trim().toLowerCase();
        const matches = (mapping, key) => key && mapping.get(key).condition.toLowerCase().includes(query);
        const keys = [...new Set([...a.keys(), ...b.keys()])].filter(key => key &&
          (a.get(key) || b.get(key)).condition.toLowerCase().includes(query))
          .sort((x, y) => Math.max(fraction(a.get(y)), fraction(b.get(y))) - Math.max(fraction(a.get(x)), fraction(b.get(x))))
          .slice(0, 15).reverse();
        const trace = (mapping, group, name, color, shownKeys) => ({
          type: 'bar', orientation: 'h', name, marker: {color},
          y: shownKeys.map(key => mapping.get(key)?.condition || a.get(key)?.condition || b.get(key)?.condition),
          x: shownKeys.map(key => fraction(mapping.get(key))),
          customdata: shownKeys.map(key => [mapping.get(key)?.variant_count || 0, group[0]?.total_variant_count || 0,
                                      group[0]?.named_variant_count || 0]),
          hovertemplate: '%{y}<br>Alleles: %{customdata[0]} / %{customdata[1]}<br>Fraction: %{x:.2%}<br>Alleles with named conditions: %{customdata[2]}<extra>%{fullData.name}</extra>',
        });
        const traces = [[a, gaph, 'GAPH', '#356d8f'], [b, background, 'ClinVar', '#9ca3af']]
          .map(([mapping, group, name, color]) => trace(mapping, group, name, color, keys));
        Plotly.react('pathogenic-conditions-plot', traces, {
          template: 'plotly_white', barmode: 'group', height: 600,
          margin: {l: 300, r: 25, t: 25, b: 60},
          xaxis: {title: {text: 'Fraction of P/LP alleles'}, tickformat: '.0%', rangemode: 'tozero'},
          yaxis: {automargin: true},
          annotations: keys.length ? [] : [{text: 'No matching named conditions.', x: 0.5, y: 0.5, xref: 'paper', yref: 'paper', showarrow: false}],
        }, {responsive: true});
        const clinvarKeys = [...b.keys()].filter(key => matches(b, key))
          .sort((x, y) => fraction(b.get(y)) - fraction(b.get(x)))
          .slice(0, 10).reverse();
        Plotly.react(
          'pathogenic-clinvar-distribution-plot',
          [trace(b, background, 'ClinVar', '#9ca3af', clinvarKeys)],
          {
            template: 'plotly_white', height: 460, showlegend: false,
            margin: {l: 300, r: 25, t: 25, b: 60},
            xaxis: {title: {text: 'Fraction of P/LP alleles'}, tickformat: '.0%', rangemode: 'tozero'},
            yaxis: {automargin: true},
            annotations: clinvarKeys.length ? [] : [{text: 'No matching named conditions.', x: 0.5, y: 0.5, xref: 'paper', yref: 'paper', showarrow: false}],
          },
          {responsive: true},
        );
      }
      controls.querySelectorAll('select').forEach(item => item.addEventListener('change', render));
      select('search').addEventListener('input', render);
      render();
    })();</script>
    """
    )


def pathogenic_variant_table_html(variants: pd.DataFrame) -> str:
    table = pd.DataFrame(
        {
            "Key": variants["variant_key"],
            "Gene": variants["gene_ids"],
            "Subtype": variants["pathogenic_subtype"],
            "Low penetrance": variants["low_penetrance"].map(lambda value: "Yes" if value else ""),
            "Stars": variants["clinvar_review_stars"].map(
                lambda value: str(value) if str(value) in {"0", "1", "2", "3", "4"} else "Unmapped"
            ),
            "Review status": variants["clinvar_review_status"],
            "SCVs": pd.to_numeric(variants["clinvar_scv_count"], errors="coerce"),
            "ClinVar ID": variants["clinvar_ids"],
            "Allele ID": variants["clinvar_allele_id"],
            "Conditions": variants["conditions"],
            "Condition IDs": variants["condition_ids"],
            "HGVS": variants["clinvar_hgvs"],
            "Event": variants["event_type"],
            "VEP consequence": variants["vep_primary_consequence"],
            "Consequence group": variants["vep_consequence_group"],
            "VEP terms": variants["vep_consequence_terms"],
            "VEP transcript": variants["vep_transcript_id"],
            "MANE Select": variants["vep_mane_select"],
            "VEP status": variants["vep_status"],
            "gnomAD AF": variants["gnomad_af"],
            "Mean ortholog support": variants["support_ortholog_mean"],
            "Min ortholog support": variants["support_ortholog_min"],
            "Max ortholog support": variants["support_ortholog_max"],
            "Strategies": variants["strategies"].map(
                lambda value: ", ".join(
                    strategy_label(item.strip())
                    for item in str(value or "").split(",")
                    if item.strip()
                )
            ),
        }
    )
    numeric_columns = {
        "SCVs",
        "gnomAD AF",
        "Mean ortholog support",
        "Min ortholog support",
        "Max ortholog support",
    }
    columns = [
        {
            "key": column,
            "type": (
                "stars"
                if column == "Stars"
                else "number" if column in numeric_columns else "string"
            ),
        }
        for column in table.columns
    ]
    data_json = json.dumps(
        dataframe_records(table), ensure_ascii=False, separators=(",", ":")
    ).replace("</", "<\\/")
    columns_json = json.dumps(columns, ensure_ascii=False, separators=(",", ":"))
    header = "".join(f"<th>{html.escape(column)}</th>" for column in table.columns)
    options = "".join(
        f'<option value="{html.escape(column)}">{html.escape(column)}</option>'
        for column in table.columns
    )
    return f"""
    <div class="pathogenic-table-wrap">
      <table id="pathogenic-variant-table" class="pathogenic-table">
        <thead><tr>{header}</tr></thead><tbody></tbody>
      </table>
    </div>
    <div class="pathogenic-table-footer">
      <div class="pathogenic-sort-controls">
        <label>Primary sort<select id="pathogenic-primary-sort">{options}</select></label>
        <label>Direction<select id="pathogenic-primary-direction">
          <option value="desc">Descending</option><option value="asc">Ascending</option>
        </select></label>
        <label>Secondary sort<select id="pathogenic-secondary-sort">
          <option value="">None</option>{options}
        </select></label>
        <label>Direction<select id="pathogenic-secondary-direction">
          <option value="desc">Descending</option><option value="asc">Ascending</option>
        </select></label>
      </div>
      <div class="pathogenic-pagination">
        <button type="button" id="pathogenic-prev">Previous</button>
        <span id="pathogenic-page-status"></span>
        <button type="button" id="pathogenic-next">Next</button>
      </div>
    </div>
    <script>
    (() => {{
      const rows = {data_json};
      const columns = {columns_json};
      const pageSize = 100;
      let page = 0;
      const tableBody = document.querySelector('#pathogenic-variant-table tbody');
      const primary = document.getElementById('pathogenic-primary-sort');
      const primaryDirection = document.getElementById('pathogenic-primary-direction');
      const secondary = document.getElementById('pathogenic-secondary-sort');
      const secondaryDirection = document.getElementById('pathogenic-secondary-direction');
      const previous = document.getElementById('pathogenic-prev');
      const next = document.getElementById('pathogenic-next');
      const status = document.getElementById('pathogenic-page-status');
      const columnTypes = Object.fromEntries(columns.map(column => [column.key, column.type]));

      primary.value = 'Stars';
      secondary.value = 'Max ortholog support';

      function compareValues(left, right, key, direction) {{
        let leftValue = left[key];
        let rightValue = right[key];
        if (columnTypes[key] === 'stars') {{
          leftValue = Number.isFinite(Number(leftValue)) ? Number(leftValue) : null;
          rightValue = Number.isFinite(Number(rightValue)) ? Number(rightValue) : null;
        }}
        const leftMissing = leftValue === null || leftValue === '';
        const rightMissing = rightValue === null || rightValue === '';
        if (leftMissing || rightMissing) {{
          if (leftMissing && rightMissing) return 0;
          return leftMissing ? 1 : -1;
        }}
        let result;
        if (columnTypes[key] === 'number' || columnTypes[key] === 'stars') {{
          result = Number(leftValue) - Number(rightValue);
        }} else {{
          result = String(leftValue).localeCompare(String(rightValue), undefined, {{numeric: true}});
        }}
        return direction === 'asc' ? result : -result;
      }}

      function sortedRows() {{
        return rows.map((row, index) => ({{row, index}})).sort((left, right) => {{
          let result = compareValues(left.row, right.row, primary.value, primaryDirection.value);
          if (!result && secondary.value && secondary.value !== primary.value) {{
            result = compareValues(
              left.row, right.row, secondary.value, secondaryDirection.value
            );
          }}
          if (!result) result = compareValues(left.row, right.row, 'Key', 'asc');
          return result || left.index - right.index;
        }}).map(item => item.row);
      }}

      function displayValue(value, key) {{
        if (value === null || value === '') return '';
        if (columnTypes[key] !== 'number') return String(value);
        const numeric = Number(value);
        if (!Number.isFinite(numeric)) return String(value);
        if (key === 'gnomAD AF' && numeric !== 0) return numeric.toExponential(3);
        return numeric.toLocaleString(undefined, {{maximumFractionDigits: 4}});
      }}

      function render() {{
        const ordered = sortedRows();
        const pages = Math.max(1, Math.ceil(ordered.length / pageSize));
        page = Math.min(page, pages - 1);
        tableBody.replaceChildren();
        ordered.slice(page * pageSize, (page + 1) * pageSize).forEach(row => {{
          const tr = document.createElement('tr');
          columns.forEach(column => {{
            const td = document.createElement('td');
            td.textContent = displayValue(row[column.key], column.key);
            tr.appendChild(td);
          }});
          tableBody.appendChild(tr);
        }});
        const first = ordered.length ? page * pageSize + 1 : 0;
        const last = Math.min((page + 1) * pageSize, ordered.length);
        status.textContent = `${{first.toLocaleString()}}–${{last.toLocaleString()}} of ${{ordered.length.toLocaleString()}}`;
        previous.disabled = page === 0;
        next.disabled = page >= pages - 1;
      }}

      [primary, primaryDirection, secondary, secondaryDirection].forEach(control => {{
        control.addEventListener('change', () => {{ page = 0; render(); }});
      }});
      previous.addEventListener('click', () => {{ if (page > 0) {{ page -= 1; render(); }} }});
      next.addEventListener('click', () => {{ page += 1; render(); }});
      render();
    }})();
    </script>
    """
