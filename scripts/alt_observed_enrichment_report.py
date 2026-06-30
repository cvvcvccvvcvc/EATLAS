#!/usr/bin/env python3
"""Build a scientific ALT-observed enrichment HTML report.

The script expects a variant-level GAPH feature table. Conservation scores
should already be present in the input table, for example from:

    scripts/annotate_variant_conservation.py

It performs two analyses:

1. Global enrichment of ALT_observed in ClinVar benign versus pathogenic rows.
2. The same enrichment within conservation-score bins, with a pooled
   Mantel-Haenszel adjusted odds ratio across bins.
"""

from __future__ import annotations

import argparse
import html
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


BENIGN_LABELS = {"0", "false", "benign", "likely_benign", "b/lb", "lb/b", "b", "lb"}
PATHOGENIC_LABELS = {"1", "true", "pathogenic", "likely_pathogenic", "p/lp", "lp/p", "p", "lp"}
AUTO_CONSERVATION_COLUMNS = [
    "phyloP100way",
    "phastCons100way",
    "GERP_RS_92mammals",
    "phyloP",
    "phastCons",
    "GERP",
    "gerp",
]


@dataclass(frozen=True)
class EnrichmentResult:
    name: str
    row_count: int
    benign_alt: int
    pathogenic_alt: int
    benign_no_alt: int
    pathogenic_no_alt: int
    odds_ratio: float
    ci_low: float
    ci_high: float
    fisher_p: float


@dataclass(frozen=True)
class AdjustedResult:
    odds_ratio_mh: float
    ci_low: float
    ci_high: float
    cmh_chi2: float
    cmh_p: float


@dataclass(frozen=True)
class ConservationAnalysis:
    column: str
    usable_rows: int
    bin_count: int
    bins: list[EnrichmentResult]
    adjusted: AdjustedResult | None
    warning: str


@dataclass(frozen=True)
class Report:
    df: pd.DataFrame
    strategy: str
    global_result: EnrichmentResult
    conservation: list[ConservationAnalysis]
    missing_conservation_columns: list[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features-tsv", required=True, type=Path, help="Input GAPH feature TSV/TSV.GZ.")
    parser.add_argument(
        "--out-html",
        type=Path,
        help="Output HTML path. Default: reports/<features-tsv-name>_alt_observed_enrichment.html",
    )
    parser.add_argument(
        "--report-name",
        help="Short report file name inside the project reports/ directory. '.html' is added if omitted.",
    )
    parser.add_argument("--label-column", default="label")
    parser.add_argument("--alt-count-column", default="gaph_all_alt_count")
    parser.add_argument("--strategy", help="Strategy to analyze. Required if the input has multiple strategies.")
    parser.add_argument(
        "--conservation-columns",
        default="auto",
        help=(
            "Comma-separated conservation columns. Use 'auto' to analyze known columns "
            "present in the table, or 'none' to skip conservation-stratified analysis."
        ),
    )
    parser.add_argument(
        "--conservation-column",
        action="append",
        default=[],
        help="Single conservation column to analyze. Can be repeated. Overrides --conservation-columns.",
    )
    parser.add_argument("--conservation-bins", type=int, default=4)
    return parser.parse_args()


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def reports_dir() -> Path:
    return project_root() / "reports"


def strip_table_suffix(path: Path) -> str:
    name = path.name
    for suffix in [".tsv.gz", ".txt.gz", ".csv.gz", ".tsv", ".txt", ".csv", ".gz"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def safe_report_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    cleaned = cleaned.strip("._-")
    return cleaned or "report"


def resolve_out_html(args: argparse.Namespace) -> Path:
    if args.out_html:
        return args.out_html
    if args.report_name:
        name = safe_report_name(Path(args.report_name).name)
        if not name.endswith(".html"):
            name += ".html"
        return reports_dir() / name
    base = safe_report_name(strip_table_suffix(args.features_tsv))
    return reports_dir() / f"{base}_alt_observed_enrichment.html"


def normalize_label(value: object) -> str | None:
    text = str(value).strip().lower().replace(" ", "_")
    if text in BENIGN_LABELS:
        return "benign"
    if text in PATHOGENIC_LABELS:
        return "pathogenic"
    return None


def require_columns(df: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(missing)}")


def select_strategy(df: pd.DataFrame, strategy: str | None) -> tuple[pd.DataFrame, str]:
    if "strategy" not in df.columns:
        return df.copy(), ""
    strategies = sorted(str(value) for value in df["strategy"].dropna().unique())
    if strategy:
        if strategy not in strategies:
            raise ValueError(f"Requested strategy {strategy!r} not found. Available: {', '.join(strategies)}")
        return df[df["strategy"].astype(str) == strategy].copy(), strategy
    if len(strategies) == 1:
        return df.copy(), strategies[0]
    raise ValueError(
        "Input contains multiple strategies; pass --strategy. "
        f"Available strategies: {', '.join(strategies)}"
    )


def log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_log_prob(k: int, row1: int, row2: int, col1: int, total: int) -> float:
    return log_comb(row1, k) + log_comb(row2, col1 - k) - log_comb(total, col1)


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact test for [[a, b], [c, d]]."""

    row1 = a + b
    row2 = c + d
    col1 = a + c
    total = row1 + row2
    if total == 0:
        return float("nan")
    low = max(0, col1 - row2)
    high = min(row1, col1)
    observed_log = hypergeom_log_prob(a, row1, row2, col1, total)
    included_logs: list[float] = []
    for k in range(low, high + 1):
        log_prob = hypergeom_log_prob(k, row1, row2, col1, total)
        if log_prob <= observed_log + 1e-12:
            included_logs.append(log_prob)
    if not included_logs:
        return 0.0
    max_log = max(included_logs)
    log_p = max_log + math.log(sum(math.exp(value - max_log) for value in included_logs))
    return min(1.0, math.exp(log_p))


def odds_ratio_and_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    """Odds ratio and approximate 95% CI for [[a, b], [c, d]].

    The raw OR is a*d/(b*c). The CI uses the normal approximation on log(OR)
    with a 0.5 Haldane-Anscombe correction when any cell is zero.
    """

    numerator = a * d
    denominator = b * c
    if denominator == 0:
        raw_or = float("inf") if numerator > 0 else float("nan")
    else:
        raw_or = numerator / denominator

    aa, bb, cc, dd = float(a), float(b), float(c), float(d)
    if min(aa, bb, cc, dd) == 0.0:
        aa += 0.5
        bb += 0.5
        cc += 0.5
        dd += 0.5
    corrected_or = (aa * dd) / (bb * cc)
    se = math.sqrt((1.0 / aa) + (1.0 / bb) + (1.0 / cc) + (1.0 / dd))
    log_or = math.log(corrected_or)
    return raw_or, math.exp(log_or - 1.96 * se), math.exp(log_or + 1.96 * se)


def enrichment_for_subset(df: pd.DataFrame, name: str) -> EnrichmentResult:
    benign = df["label_class"] == "benign"
    pathogenic = df["label_class"] == "pathogenic"
    alt = df["ALT_observed"]
    a = int((benign & alt).sum())
    b = int((pathogenic & alt).sum())
    c = int((benign & ~alt).sum())
    d = int((pathogenic & ~alt).sum())
    odds_ratio, ci_low, ci_high = odds_ratio_and_ci(a, b, c, d)
    return EnrichmentResult(
        name=name,
        row_count=int(len(df)),
        benign_alt=a,
        pathogenic_alt=b,
        benign_no_alt=c,
        pathogenic_no_alt=d,
        odds_ratio=odds_ratio,
        ci_low=ci_low,
        ci_high=ci_high,
        fisher_p=fisher_exact_two_sided(a, b, c, d),
    )


def resolve_conservation_columns(df: pd.DataFrame, args: argparse.Namespace) -> tuple[list[str], list[str]]:
    requested: list[str]
    if args.conservation_column:
        requested = []
        for item in args.conservation_column:
            requested.extend(part.strip() for part in item.split(",") if part.strip())
    elif args.conservation_columns.strip().lower() == "auto":
        requested = [column for column in AUTO_CONSERVATION_COLUMNS if column in df.columns]
    elif args.conservation_columns.strip().lower() in {"", "none", "skip"}:
        requested = []
    else:
        requested = [part.strip() for part in args.conservation_columns.split(",") if part.strip()]

    seen = set()
    columns = []
    missing = []
    for column in requested:
        if column in seen:
            continue
        seen.add(column)
        if column in df.columns:
            columns.append(column)
        else:
            missing.append(column)
    return columns, missing


def conservation_bin_results(df: pd.DataFrame, column: str, requested_bins: int) -> ConservationAnalysis:
    values = pd.to_numeric(df[column], errors="coerce")
    working = df[values.notna()].copy()
    working[column] = values[values.notna()]
    if working.empty:
        return ConservationAnalysis(column, 0, 0, [], None, "No numeric conservation values.")

    unique_values = int(working[column].nunique())
    if unique_values < 2:
        return ConservationAnalysis(column, len(working), 0, [], None, "Fewer than two unique conservation values.")

    q = max(2, min(requested_bins, unique_values))
    try:
        working["conservation_bin"] = pd.qcut(working[column], q=q, duplicates="drop")
    except ValueError as exc:
        return ConservationAnalysis(column, len(working), 0, [], None, f"Could not form bins: {exc}")

    categories = list(working["conservation_bin"].cat.categories)
    if len(categories) < 2:
        return ConservationAnalysis(column, len(working), len(categories), [], None, "Only one bin after qcut.")

    bins = []
    for index, interval in enumerate(categories, start=1):
        subset = working[working["conservation_bin"] == interval]
        name = f"Q{index}: {fmt_number(interval.left)} to {fmt_number(interval.right)}"
        bins.append(enrichment_for_subset(subset, name))

    adjusted = mantel_haenszel_adjusted(bins)
    return ConservationAnalysis(column, len(working), len(bins), bins, adjusted, "")


def mantel_haenszel_adjusted(results: list[EnrichmentResult]) -> AdjustedResult | None:
    """Pooled stratified association summary across conservation bins.

    The reported adjusted OR is the Mantel-Haenszel common odds ratio. The
    confidence interval is a transparent fixed-effect log-OR approximation using
    stratum-level Haldane-corrected variances. The p-value is the CMH
    chi-square test with 1 degree of freedom.
    """

    numerator = 0.0
    denominator = 0.0
    observed_minus_expected = 0.0
    variance_sum = 0.0
    weighted_log_or = 0.0
    weight_sum = 0.0

    for result in results:
        a = result.benign_alt
        b = result.pathogenic_alt
        c = result.benign_no_alt
        d = result.pathogenic_no_alt
        n = a + b + c + d
        if n <= 1:
            continue
        numerator += (a * d) / n
        denominator += (b * c) / n

        alt_total = a + b
        no_alt_total = c + d
        benign_total = a + c
        pathogenic_total = b + d
        expected_a = alt_total * benign_total / n
        variance_a = alt_total * no_alt_total * benign_total * pathogenic_total / (n * n * (n - 1))
        observed_minus_expected += a - expected_a
        variance_sum += variance_a

        aa, bb, cc, dd = float(a), float(b), float(c), float(d)
        if min(aa, bb, cc, dd) == 0.0:
            aa += 0.5
            bb += 0.5
            cc += 0.5
            dd += 0.5
        var_log_or = (1.0 / aa) + (1.0 / bb) + (1.0 / cc) + (1.0 / dd)
        if var_log_or > 0:
            weight = 1.0 / var_log_or
            weighted_log_or += weight * math.log((aa * dd) / (bb * cc))
            weight_sum += weight

    if denominator == 0.0 and numerator == 0.0:
        odds_ratio_mh = float("nan")
    elif denominator == 0.0:
        odds_ratio_mh = float("inf")
    else:
        odds_ratio_mh = numerator / denominator

    if weight_sum > 0:
        pooled_log_or = weighted_log_or / weight_sum
        se = math.sqrt(1.0 / weight_sum)
        ci_low = math.exp(pooled_log_or - 1.96 * se)
        ci_high = math.exp(pooled_log_or + 1.96 * se)
    else:
        ci_low = float("nan")
        ci_high = float("nan")

    if variance_sum > 0:
        cmh_chi2 = (observed_minus_expected * observed_minus_expected) / variance_sum
        cmh_p = math.erfc(math.sqrt(cmh_chi2 / 2.0))
    else:
        cmh_chi2 = float("nan")
        cmh_p = float("nan")

    return AdjustedResult(odds_ratio_mh, ci_low, ci_high, cmh_chi2, cmh_p)


def prepare_dataset(args: argparse.Namespace) -> tuple[pd.DataFrame, str, list[str], list[str]]:
    df = pd.read_csv(args.features_tsv, sep="\t", low_memory=False)
    df, strategy = select_strategy(df, args.strategy)
    require_columns(df, [args.label_column, args.alt_count_column])

    labels = df[args.label_column].map(normalize_label)
    df = df[labels.notna()].copy()
    df["label_class"] = labels[labels.notna()].to_numpy()
    df[args.alt_count_column] = pd.to_numeric(df[args.alt_count_column], errors="coerce").fillna(0)
    df["ALT_observed"] = df[args.alt_count_column] > 0

    if df.empty:
        raise ValueError("No rows with usable benign/pathogenic labels")
    if set(df["label_class"]) != {"benign", "pathogenic"}:
        raise ValueError("Both benign and pathogenic rows are required")

    conservation_columns, missing_columns = resolve_conservation_columns(df, args)
    return df, strategy, conservation_columns, missing_columns


def compute_report(args: argparse.Namespace) -> Report:
    df, strategy, conservation_columns, missing_columns = prepare_dataset(args)
    global_result = enrichment_for_subset(df, "all variants")
    conservation = [
        conservation_bin_results(df, column, args.conservation_bins)
        for column in conservation_columns
    ]
    return Report(df, strategy, global_result, conservation, missing_columns)


def fmt_number(value: float, digits: int = 4) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    if isinstance(value, float) and math.isinf(value):
        return "inf"
    if abs(float(value)) < 0.0001 and float(value) != 0.0:
        return f"{float(value):.3e}"
    return f"{float(value):.{digits}g}"


def escape(value: object) -> str:
    return html.escape(str(value), quote=True)


def table_html(headers: list[str], rows: list[list[object]]) -> str:
    return (
        "<table><thead><tr>"
        + "".join(f"<th>{escape(header)}</th>" for header in headers)
        + "</tr></thead><tbody>"
        + "".join(
            "<tr>" + "".join(f"<td>{escape(cell)}</td>" for cell in row) + "</tr>"
            for row in rows
        )
        + "</tbody></table>"
    )


def enrichment_table(results: list[EnrichmentResult], include_name: bool = True) -> str:
    headers = []
    if include_name:
        headers.append("Subset")
    headers.extend(
        [
            "N",
            "Benign ALT+",
            "Pathogenic ALT+",
            "Benign ALT-",
            "Pathogenic ALT-",
            "OR",
            "95% CI",
            "Fisher p",
        ]
    )
    rows = []
    for result in results:
        row = []
        if include_name:
            row.append(result.name)
        row.extend(
            [
                result.row_count,
                result.benign_alt,
                result.pathogenic_alt,
                result.benign_no_alt,
                result.pathogenic_no_alt,
                fmt_number(result.odds_ratio),
                f"{fmt_number(result.ci_low)} to {fmt_number(result.ci_high)}",
                fmt_number(result.fisher_p),
            ]
        )
        rows.append(row)
    return table_html(headers, rows)


def metric_cards(items: list[tuple[str, object]]) -> str:
    cards = []
    for label, value in items:
        cards.append(
            f"""
            <div class="metric-card">
                <div class="metric-label">{escape(label)}</div>
                <div class="metric-value">{escape(value)}</div>
            </div>
            """
        )
    return f"<div class=\"metric-grid\">{''.join(cards)}</div>"


def build_overview(report: Report, args: argparse.Namespace) -> list[str]:
    df = report.df
    label_counts = df["label_class"].value_counts().to_dict()
    alt_count = int(df["ALT_observed"].sum())
    no_alt_count = int((~df["ALT_observed"]).sum())
    conservation_names = [analysis.column for analysis in report.conservation]
    if report.missing_conservation_columns:
        conservation_names.extend(f"{column} (missing)" for column in report.missing_conservation_columns)

    return [
        "<h2>Overview</h2>",
        metric_cards(
            [
                ("Rows analyzed", f"{len(df):,}"),
                ("Benign rows", f"{label_counts.get('benign', 0):,}"),
                ("Pathogenic rows", f"{label_counts.get('pathogenic', 0):,}"),
                ("ALT observed", f"{alt_count:,}"),
                ("ALT not observed", f"{no_alt_count:,}"),
                ("Global OR", fmt_number(report.global_result.odds_ratio)),
                ("Global Fisher p", fmt_number(report.global_result.fisher_p)),
                ("Conservation scores", ", ".join(conservation_names) if conservation_names else "none"),
            ]
        ),
        "<h3>Input</h3>",
        table_html(
            ["Field", "Value"],
            [
                ["Feature TSV", args.features_tsv],
                ["Strategy", report.strategy or "not provided"],
                ["ALT observed definition", f"{args.alt_count_column} > 0"],
                ["Conservation bins", args.conservation_bins],
            ],
        ),
    ]


def build_global_section(report: Report) -> list[str]:
    return [
        "<h2>Experiment 1: Global ALT-Observed Enrichment</h2>",
        """
        <p>
        This test asks whether ClinVar benign SNVs are more likely than
        pathogenic SNVs to have the exact human ALT allele observed in at least
        one ortholog.
        </p>
        """,
        enrichment_table([report.global_result], include_name=False),
    ]


def adjusted_table(analysis: ConservationAnalysis) -> str:
    if analysis.adjusted is None:
        return "<p class=\"warning\">No adjusted stratified estimate could be computed.</p>"
    adjusted = analysis.adjusted
    return table_html(
        ["Conservation score", "Usable rows", "Bins", "MH adjusted OR", "Approx 95% CI", "CMH chi2", "CMH p"],
        [
            [
                analysis.column,
                analysis.usable_rows,
                analysis.bin_count,
                fmt_number(adjusted.odds_ratio_mh),
                f"{fmt_number(adjusted.ci_low)} to {fmt_number(adjusted.ci_high)}",
                fmt_number(adjusted.cmh_chi2),
                fmt_number(adjusted.cmh_p),
            ]
        ],
    )


def build_conservation_section(report: Report) -> list[str]:
    parts = [
        "<h2>Experiment 2: Conservation-Stratified Enrichment</h2>",
        """
        <p>
        Each conservation score is split into quantile bins. Within every bin,
        the report repeats the same 2x2 ALT-observed enrichment test. The pooled
        summary reports a Mantel-Haenszel adjusted odds ratio across bins.
        </p>
        """,
    ]
    if report.missing_conservation_columns:
        parts.append(
            "<p class=\"warning\">Missing requested conservation columns: "
            + escape(", ".join(report.missing_conservation_columns))
            + ".</p>"
        )
    if not report.conservation:
        parts.append(
            """
            <p class="warning">
            No conservation columns were analyzed. Run
            <code>scripts/annotate_variant_conservation.py</code> first, or pass
            <code>--conservation-columns</code> with existing numeric columns.
            </p>
            """
        )
        return parts

    for analysis in report.conservation:
        parts.append(f"<h3>{escape(analysis.column)}</h3>")
        if analysis.warning:
            parts.append(f"<p class=\"warning\">{escape(analysis.warning)}</p>")
            continue
        parts.append(adjusted_table(analysis))
        parts.append(enrichment_table(analysis.bins, include_name=True))
    return parts


def build_methods_section(args: argparse.Namespace) -> list[str]:
    return [
        "<h2>Methods</h2>",
        """
        <p>
        The primary binary feature is <code>ALT_observed</code>, defined as
        <code>gaph_all_alt_count &gt; 0</code> by default. Labels are collapsed
        to two classes: benign / likely benign and pathogenic / likely
        pathogenic.
        </p>
        <p>
        For each 2x2 table, the odds ratio is
        <code>(benign_ALT+ / benign_ALT-) / (pathogenic_ALT+ / pathogenic_ALT-)</code>.
        Fisher p-values are two-sided. The per-table 95% confidence interval is
        the normal approximation on log(OR), using a 0.5 Haldane-Anscombe
        correction when any cell is zero.
        </p>
        <p>
        Conservation-stratified results use quantile bins. The adjusted odds
        ratio is the Mantel-Haenszel common odds ratio across bins. The adjusted
        confidence interval shown in the table is a fixed-effect log-OR
        approximation from stratum-level corrected variances. The stratified
        p-value is the Cochran-Mantel-Haenszel chi-square test with 1 degree of
        freedom.
        </p>
        """,
        "<h3>Files</h3>",
        table_html(
            ["Field", "Value"],
            [
                ["Input feature TSV", args.features_tsv],
                ["Output HTML", args.out_html],
                ["Label column", args.label_column],
                ["ALT count column", args.alt_count_column],
                ["Conservation columns argument", args.conservation_columns],
                ["Repeated --conservation-column", ", ".join(args.conservation_column) or ""],
            ],
        ),
    ]


def render_tabs(sections: list[tuple[str, str, list[str]]]) -> str:
    buttons = []
    pages = []
    for index, (tab_id, title, html_parts) in enumerate(sections):
        active = " active" if index == 0 else ""
        buttons.append(f'<button class="tab-button{active}" data-tab="{tab_id}">{escape(title)}</button>')
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
        }});
    }});
    </script>
    """


def build_html(report: Report, args: argparse.Namespace) -> str:
    sections = [
        ("overview", "Overview", build_overview(report, args)),
        ("global", "Global Enrichment", build_global_section(report)),
        ("conservation", "Conservation Strata", build_conservation_section(report)),
        ("methods", "Methods", build_methods_section(args)),
    ]
    return f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>GAPH ALT-Observed Enrichment Report</title>
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
            code {{
                background: #f5f7fa;
                border: 1px solid #d5d9df;
                border-radius: 4px;
                padding: 1px 4px;
            }}
            .warning {{
                border-left: 4px solid #d64545;
                background: #fff5f5;
                padding: 10px 12px;
            }}
        </style>
    </head>
    <body>
        <h1>GAPH ALT-Observed Enrichment Report</h1>
        <p class="lead">ClinVar benign/pathogenic enrichment and conservation-stratified controls.</p>
        {render_tabs(sections)}
    </body>
    </html>
    """


def main() -> None:
    args = parse_args()
    args.out_html = resolve_out_html(args)
    report = compute_report(args)
    html_text = build_html(report, args)
    out_html = os.path.abspath(args.out_html)
    os.makedirs(os.path.dirname(out_html), exist_ok=True)
    with open(out_html, "w") as handle:
        handle.write(html_text)
    print(f"Wrote {out_html}")


if __name__ == "__main__":
    main()
