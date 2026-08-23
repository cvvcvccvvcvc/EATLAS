"""Reusable formatting and figure components for analytics reports."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from .config import STRATEGY_LABELS

def strategy_label(value: str) -> str:
    return STRATEGY_LABELS.get(str(value), str(value))


def sort_by_metric(df: pd.DataFrame, column: str, ascending: bool = False) -> pd.DataFrame:
    if df.empty or column not in df.columns:
        return df
    return df.sort_values(column, ascending=ascending, kind="mergesort")


def format_int(value) -> str:
    if pd.isna(value):
        return ""
    return f"{int(round(float(value))):,}".replace(",", " ")


def format_float(value, digits: int = 3) -> str:
    if pd.isna(value):
        return ""
    return f"{float(value):.{digits}f}".rstrip("0").rstrip(".")


def format_percent(value, digits: int = 1) -> str:
    if pd.isna(value):
        return ""
    return f"{float(value) * 100:.{digits}f}%"


def format_table_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    shown = df.copy()
    for column in shown.columns:
        if column == "Strategy":
            continue
        if column.endswith("%") or " rate" in column.lower() or "breadth" in column.lower():
            shown[column] = shown[column].map(format_percent)
        elif any(token in column.lower() for token in ["variant", "found", "event", "ortholog", "gene", "row", "bp"]):
            numeric = pd.to_numeric(shown[column], errors="coerce")
            nonempty = shown[column].notna() & shown[column].astype(str).ne("")
            if bool(nonempty.any()) and numeric[nonempty].notna().all():
                shown[column] = numeric.map(format_int)
        elif pd.api.types.is_integer_dtype(shown[column]):
            shown[column] = shown[column].map(format_int)
        elif pd.api.types.is_float_dtype(shown[column]):
            shown[column] = shown[column].map(lambda value: format_float(value, 3))
    return shown


def fig_html(fig) -> str:
    return fig.to_html(full_html=False, include_plotlyjs=False)


def compact_figure(fig, height: int = 340, show_x_title: bool = False):
    fig.update_layout(
        height=height,
        margin={"l": 55, "r": 20, "t": 52, "b": 58},
        template="plotly_white",
        legend_title_text="",
    )
    fig.update_xaxes(automargin=True)
    fig.update_yaxes(automargin=True)
    if not show_x_title:
        fig.update_xaxes(title_text=None)
    return fig


def table_html(df: pd.DataFrame, classes: str = "table table-striped table-bordered", max_rows: int | None = None) -> str:
    shown = df if max_rows is None else df.head(max_rows)
    shown = format_table_dataframe(shown)
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


def format_count_share(count: object, total: object) -> str:
    if pd.isna(count) or pd.isna(total):
        return "n/a"
    count_value = int(count)
    total_value = int(total)
    if total_value == 0:
        return f"{format_int(count_value)} (n/a)"
    fraction = count_value / total_value
    percent_digits = 1 if fraction == 0 else 3 if fraction < 0.001 else 2 if fraction < 0.01 else 1
    return f"{format_int(count_value)} ({format_percent(fraction, percent_digits)})"


def format_count_ratio(count: object, total: object) -> str:
    if pd.isna(count) or pd.isna(total):
        return "n/a"
    count_value = int(count)
    total_value = int(total)
    if total_value == 0:
        return f"{format_int(count_value)} / 0 (n/a)"
    fraction = count_value / total_value
    return f"{format_int(count_value)} / {format_int(total_value)} ({format_percent(fraction)})"


def dataframe_records(frame: pd.DataFrame) -> list[dict[str, object]]:
    if frame.empty:
        return []
    records = frame.to_dict(orient="records")
    for row in records:
        for key, value in row.items():
            if value is None or pd.isna(value):
                row[key] = None
            elif isinstance(value, (float, np.floating)) and math.isinf(float(value)):
                row[key] = "inf" if float(value) > 0 else "-inf"
    return records
