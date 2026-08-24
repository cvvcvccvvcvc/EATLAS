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
from genomics.clinvar import PATHOGENIC_SUBTYPE_ORDER
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


SUBTYPE_COLORS = {
    "Pathogenic": "#b2182b",
    "Likely pathogenic": "#ef8a62",
    "Pathogenic / likely pathogenic": "#7b3294",
    "P/LP": "#666666",
}


def build_pathogenic_variant_sections(
    analysis: PathogenicVariantAnalysis,
) -> list[str]:
    variants = analysis.variants
    if variants.empty:
        return [
            "<h2>Pathogenic ClinVar Hits</h2>",
            "<p>No candidate variants with an unambiguous pathogenic or likely "
            "pathogenic ClinVar classification were found.</p>",
        ]

    stars = pd.to_numeric(variants["clinvar_review_stars"], errors="coerce")
    sections = [
        "<h2>Pathogenic ClinVar Hits</h2>",
        "<p class=\"lead\">This tab characterizes exact normalized candidate alleles "
        "classified by ClinVar as pathogenic or likely pathogenic. Counts in the review-star "
        "and disease panels are unique alleles per strategy; consequence counts are "
        "allele–target-gene observations.</p>",
        metric_cards(
            [
                ("Unique P/LP alleles", format_int(len(variants))),
                ("With ≥2 review stars", format_int(stars.ge(2).sum())),
                (
                    "Named conditions",
                    format_int(
                        analysis.condition_counts["condition"].nunique()
                        if not analysis.condition_counts.empty
                        else 0
                    ),
                ),
                ("SNV support rows plotted", format_int(len(analysis.evolution_rows))),
            ]
        ),
    ]

    star_figure = pathogenic_star_figure(analysis.star_counts)
    if star_figure is not None:
        sections.extend(["<h3>ClinVar assertion strength</h3>", fig_html(star_figure)])

    consequence_figure = pathogenic_consequence_figure(analysis.consequence_counts)
    if consequence_figure is not None:
        sections.extend(["<h3>Molecular effect</h3>", fig_html(consequence_figure)])

    sections.append("<h3>Why could a P/LP allele pass the evolutionary filter?</h3>")
    evolution_figure = pathogenic_evolution_figure(analysis.evolution_rows)
    if evolution_figure is None:
        sections.append(
            "<p>No P/LP SNV had both a phyloP100way score and an aligned-site denominator. "
            "Indels remain in the consequence analysis and detail table, but no denominator "
            "is inferred for them.</p>"
        )
    else:
        sections.extend(
            [
                "<p class=\"lead\">Each point is one P/LP SNV × target gene × strategy. "
                "The y-axis is the exact ALT-supporting ortholog fraction among orthologs "
                "aligned at that site; point size is the number of supporting genera. "
                "This separates weak evolutionary constraint from sparse or concentrated "
                "ortholog support.</p>",
                fig_html(evolution_figure),
            ]
        )

    sections.append("<h3>Associated conditions</h3>")
    condition_figure = pathogenic_condition_figure(analysis.condition_counts)
    if condition_figure is None:
        sections.append(
            "<p>No named ClinVar conditions were available. <code>not provided</code> and "
            "<code>not specified</code> are intentionally excluded.</p>"
        )
    else:
        sections.extend(
            [
                "<p class=\"lead\">Top 15 ClinVar conditions for the selected strategy, "
                "ranked by unique P/LP alleles. Broad disease categories are not inferred. "
                "Bars retain the P versus LP distinction.</p>",
                fig_html(condition_figure),
            ]
        )

    sections.extend(
        [
            "<h3>Variant details</h3>",
            "<p class=\"lead\">All P/LP alleles are available below and in the compressed "
            "TSV artifact. Ortholog-support mean/min/max summarize target-gene × strategy "
            "rows. Sorting is stable: after the selected primary and secondary columns, "
            "the normalized variant key is used as the final tie-breaker.</p>",
            "<p><a download href=\"../derived/"
            f"{quote(analysis.variants_path.name)}\">Download complete P/LP TSV</a> "
            f"(<code>{html.escape(str(analysis.variants_path))}</code>)</p>",
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
    subtype_order = [
        subtype
        for subtype in PATHOGENIC_SUBTYPE_ORDER
        if subtype in set(shown["pathogenic_subtype"])
    ]
    subtype_order.extend(
        subtype
        for subtype in shown["pathogenic_subtype"].unique()
        if subtype not in subtype_order
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
        facet_col="pathogenic_subtype",
        barmode="stack",
        title="P/LP alleles by ClinVar review stars and assertion subtype",
        category_orders={
            "Strategy": strategy_order,
            "Review stars": present_stars,
            "pathogenic_subtype": subtype_order,
        },
        color_discrete_map=REVIEW_STAR_COLORS,
        labels={
            "variant_count": "Unique P/LP alleles",
            "pathogenic_subtype": "ClinVar subtype",
        },
    )
    figure.for_each_annotation(lambda annotation: annotation.update(text=annotation.text.split("=")[-1]))
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


def pathogenic_evolution_figure(rows: pd.DataFrame):
    if rows.empty:
        return None
    strategies = sorted(rows["strategy"].astype(str).unique(), key=strategy_label)
    first_strategy = strategies[0]
    figure = go.Figure()
    trace_strategies: list[str] = []
    max_genera = pd.to_numeric(
        rows["alt_support_genus_count"], errors="coerce"
    ).max()
    maximum_size = max(1.0, float(max_genera)) if pd.notna(max_genera) else 1.0
    size_ref = 2.0 * maximum_size / 22.0**2
    for strategy in strategies:
        strategy_rows = rows[rows["strategy"].astype(str).eq(strategy)]
        subtypes = [
            subtype
            for subtype in PATHOGENIC_SUBTYPE_ORDER
            if subtype in set(strategy_rows["pathogenic_subtype"])
        ]
        subtypes.extend(
            subtype
            for subtype in strategy_rows["pathogenic_subtype"].unique()
            if subtype not in subtypes
        )
        for subtype in subtypes:
            values = strategy_rows[strategy_rows["pathogenic_subtype"].eq(subtype)]
            genera = pd.to_numeric(values["alt_support_genus_count"], errors="coerce").fillna(1)
            figure.add_trace(
                go.Scatter(
                    x=values["phylop100way"],
                    y=values["alt_support_fraction"],
                    mode="markers",
                    name=str(subtype),
                    legendgroup=str(subtype),
                    visible=strategy == first_strategy,
                    marker={
                        "size": np.maximum(genera.to_numpy(dtype=float), 1),
                        "sizemode": "area",
                        "sizeref": size_ref,
                        "sizemin": 5,
                        "color": SUBTYPE_COLORS.get(str(subtype), "#666666"),
                        "opacity": 0.72,
                    },
                    customdata=values[
                        [
                            "variant_key",
                            "gene_id",
                            "alt_support_ortholog_count",
                            "site_aligned_ortholog_count",
                            "alt_support_genus_count",
                        ]
                    ],
                    hovertemplate=(
                        "%{customdata[0]}<br>Gene: %{customdata[1]}<br>"
                        "phyloP100way: %{x:.3f}<br>Exact ALT support: "
                        "%{customdata[2]:.0f}/%{customdata[3]:.0f} (%{y:.1%})<br>"
                        "Supporting genera: %{customdata[4]:.0f}<extra>%{fullData.name}</extra>"
                    ),
                )
            )
            trace_strategies.append(strategy)
    buttons = [
        {
            "label": strategy_label(strategy),
            "method": "update",
            "args": [
                {"visible": [owner == strategy for owner in trace_strategies]},
                {"title": f"Conservation and exact ortholog support — {strategy_label(strategy)}"},
            ],
        }
        for strategy in strategies
    ]
    figure.update_layout(
        title=f"Conservation and exact ortholog support — {strategy_label(first_strategy)}",
        xaxis_title="phyloP100way score",
        yaxis_title="Exact ALT-supporting ortholog fraction",
        yaxis_tickformat=".0%",
        updatemenus=[
            {
                "buttons": buttons,
                "direction": "down",
                "x": 1.0,
                "xanchor": "right",
                "y": 1.18,
                "yanchor": "top",
            }
        ],
    )
    compact_figure(figure, height=440, show_x_title=True)
    return figure


def pathogenic_condition_figure(counts: pd.DataFrame):
    if counts.empty:
        return None
    strategies = sorted(counts["strategy"].astype(str).unique(), key=strategy_label)
    first_strategy = strategies[0]
    figure = go.Figure()
    trace_strategies: list[str] = []
    condition_orders: dict[str, list[str]] = {}
    for strategy in strategies:
        values = counts[counts["strategy"].astype(str).eq(strategy)]
        top = (
            values.groupby("condition")["variant_count"]
            .sum()
            .sort_values(ascending=False, kind="mergesort")
            .head(15)
        )
        conditions = top.sort_values().index.tolist()
        condition_orders[strategy] = conditions
        for subtype in PATHOGENIC_SUBTYPE_ORDER:
            subtype_values = (
                values[values["pathogenic_subtype"].eq(subtype)]
                .set_index("condition")["variant_count"]
                .reindex(conditions, fill_value=0)
            )
            if not subtype_values.any():
                continue
            figure.add_trace(
                go.Bar(
                    x=subtype_values.to_numpy(),
                    y=conditions,
                    orientation="h",
                    name=subtype,
                    legendgroup=subtype,
                    marker_color=SUBTYPE_COLORS[subtype],
                    visible=strategy == first_strategy,
                    hovertemplate="%{y}<br>Unique alleles: %{x:,}<extra>%{fullData.name}</extra>",
                )
            )
            trace_strategies.append(strategy)
    if not figure.data:
        return None
    buttons = [
        {
            "label": strategy_label(strategy),
            "method": "update",
            "args": [
                {"visible": [owner == strategy for owner in trace_strategies]},
                {
                    "title": f"Top ClinVar conditions — {strategy_label(strategy)}",
                    "yaxis.categoryorder": "array",
                    "yaxis.categoryarray": condition_orders[strategy],
                },
            ],
        }
        for strategy in strategies
    ]
    figure.update_layout(
        barmode="stack",
        title=f"Top ClinVar conditions — {strategy_label(first_strategy)}",
        xaxis_title="Unique P/LP alleles",
        yaxis={"categoryorder": "array", "categoryarray": condition_orders[first_strategy]},
        updatemenus=[
            {
                "buttons": buttons,
                "direction": "down",
                "x": 1.0,
                "xanchor": "right",
                "y": 1.18,
                "yanchor": "top",
            }
        ],
    )
    compact_figure(figure, height=520, show_x_title=True)
    return figure


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
            "phyloP100way": variants["phylop100way"],
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
        "phyloP100way",
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
