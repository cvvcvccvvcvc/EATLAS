#!/usr/bin/env python3
"""Build an HTML report for one completed GAPH run."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px


warnings.filterwarnings("ignore", r"All-NaN (slice|axis) encountered")
warnings.filterwarnings("ignore", r"Mean of empty slice")


FEATURE_ORDER = ["gene", "exon", "cds", "utr", "intron"]
DISJOINT_FEATURE_ORDER = ["cds", "utr", "intron"]
CLINVAR_ORDER = ["P/LP", "B/LB", "VUS", "Other", "Not Found"]
CLINVAR_COLORS = {
    "B/LB": "#2ca25f",
    "P/LP": "#de2d26",
    "VUS": "#f1c40f",
    "Other": "#8c8c8c",
}
STRATEGY_LABELS = {
    "bwa_pseudoreads": "BWA pseudo",
    "bwa_pseudoreads_varscan": "BWA VarScan",
    "minimap2_asm10": "minimap2 asm10",
    "minimap2_asm20": "minimap2 asm20",
    "minimap2_taxonomy_adaptive": "minimap2 adaptive",
    "nucmer": "nucmer",
    "precomputed_ensembl_92_mammals_epo_extended": "Ensembl EPO",
}

VARIANT_USECOLS = [
    "variant_key",
    "gene_id",
    "event_type",
    "ref",
    "alt",
    "lookup_ref",
    "lookup_alt",
    "strategies",
    "support_row_count",
    "support_ortholog_count",
    "support_strategy_count",
    "clinvar_id",
    "clinvar_sig",
    "clinvar_review_stars",
    "gnomad_af",
    "gnomad_af_source",
    "gnomad_csq",
]
VARIANT_REQUIRED = {"variant_key", "gene_id", "event_type", "strategies"}


@dataclass(frozen=True)
class RunInputs:
    run_dir: Path
    variant_annotations_tsv: Path
    annotation_manifest_json: Path
    annotation_failures_tsv: Path
    feature_coverage_tsv: Path
    alignment_manifest_json: Path
    strategy_quick_summary_tsv: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path, help="Completed GAPH run directory.")
    parser.add_argument(
        "--out-html",
        type=Path,
        help="Output HTML path. Default: <run-dir>/reports/strategy_compare.html",
    )
    parser.add_argument(
        "--report-name",
        help="Short report file name inside <run-dir>/reports. '.html' is added if omitted.",
    )
    return parser.parse_args()


def safe_report_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    cleaned = cleaned.strip("._-")
    return cleaned or "strategy_compare"


def resolve_run_inputs(run_dir: Path) -> RunInputs:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory does not exist: {run_dir}")
    if not run_dir.is_dir():
        raise NotADirectoryError(f"--run-dir is not a directory: {run_dir}")

    inputs = RunInputs(
        run_dir=run_dir,
        variant_annotations_tsv=run_dir / "annotation" / "variant_annotations.tsv.gz",
        annotation_manifest_json=run_dir / "annotation" / "manifest.json",
        annotation_failures_tsv=run_dir / "annotation" / "failures.tsv.gz",
        feature_coverage_tsv=run_dir / "alignment" / "feature_coverage.tsv.gz",
        alignment_manifest_json=run_dir / "alignment" / "manifest.json",
        strategy_quick_summary_tsv=run_dir / "reports" / "strategy_quick_summary.tsv",
    )
    if not inputs.variant_annotations_tsv.exists():
        raise FileNotFoundError(
            "Missing annotation/variant_annotations.tsv.gz under --run-dir. "
            "Run the annotation stage before building this report."
        )
    return inputs


def resolve_out_html(args: argparse.Namespace, run_dir: Path) -> Path:
    if args.out_html:
        return args.out_html.expanduser().resolve()
    report_dir = run_dir / "reports"
    if args.report_name:
        name = safe_report_name(Path(args.report_name).name)
        if not name.endswith(".html"):
            name += ".html"
        return report_dir / name
    return report_dir / "strategy_compare.html"


def open_text(path: Path):
    return gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else path.open(newline="")


def tsv_header(path: Path) -> list[str]:
    with open_text(path) as handle:
        header = handle.readline().rstrip("\n")
    return header.split("\t") if header else []


def file_size_label(path: Path) -> str:
    if not path.exists():
        return ""
    size = path.stat().st_size
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    value = float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024
    return str(size)


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
            shown[column] = shown[column].map(format_int)
        elif pd.api.types.is_integer_dtype(shown[column]):
            shown[column] = shown[column].map(format_int)
        elif pd.api.types.is_float_dtype(shown[column]):
            shown[column] = shown[column].map(lambda value: format_float(value, 3))
    return shown


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def categorize_clinvar_series(values: pd.Series) -> pd.Series:
    text = values.fillna("").astype(str).str.lower()
    category = pd.Series("Other", index=values.index, dtype="object")
    category[text.eq("")] = "Not Found"
    category[text.str.contains("conflicting", na=False)] = "Other"
    category[text.str.contains("uncertain|vus", regex=True, na=False)] = "VUS"
    category[text.str.contains("pathogenic", na=False)] = "P/LP"
    category[text.str.contains("benign", na=False)] = "B/LB"
    return pd.Categorical(category, categories=CLINVAR_ORDER, ordered=True)


def add_titv_kind(df: pd.DataFrame) -> None:
    ref = df["lookup_ref"].where(df["lookup_ref"].astype(str) != "", df["ref"]).astype(str).str.upper()
    alt = df["lookup_alt"].where(df["lookup_alt"].astype(str) != "", df["alt"]).astype(str).str.upper()
    event_type = df["event_type"].astype(str)
    valid = event_type.eq("snv") & ref.str.len().eq(1) & alt.str.len().eq(1)
    transitions = ref.str.cat(alt, sep=">").isin(["A>G", "G>A", "C>T", "T>C"])
    kind = pd.Series("", index=df.index, dtype="object")
    kind[valid & transitions] = "ti"
    kind[valid & ~transitions] = "tv"
    df["titv_kind"] = pd.Categorical(kind, categories=["", "ti", "tv"], ordered=False)


def read_variant_annotations(path: Path) -> pd.DataFrame:
    print(f"Reading {path}...")
    header = tsv_header(path)
    missing = VARIANT_REQUIRED - set(header)
    if missing:
        raise ValueError(f"Variant annotations missing required columns: {', '.join(sorted(missing))}")

    usecols = [column for column in VARIANT_USECOLS if column in header]
    df = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=usecols,
        keep_default_na=False,
        low_memory=False,
    )
    for column in VARIANT_USECOLS:
        if column not in df.columns:
            df[column] = ""

    df["variant_id"] = df["variant_key"].astype(str)
    if df["variant_id"].eq("").any():
        missing_key = df["variant_id"].eq("")
        fallback = (
            df.loc[missing_key, "gene_id"].astype(str)
            + ":"
            + df.loc[missing_key, "event_type"].astype(str)
            + ":"
            + df.loc[missing_key, "ref"].astype(str)
            + ">"
            + df.loc[missing_key, "alt"].astype(str)
        )
        df.loc[missing_key, "variant_id"] = fallback

    df["gnomad_af"] = pd.to_numeric(df["gnomad_af"], errors="coerce")
    df["support_row_count"] = pd.to_numeric(df["support_row_count"], errors="coerce").fillna(0).astype("int64")
    df["support_ortholog_count"] = (
        pd.to_numeric(df["support_ortholog_count"], errors="coerce").fillna(0).astype("int64")
    )
    df["support_strategy_count"] = (
        pd.to_numeric(df["support_strategy_count"], errors="coerce").fillna(0).astype("int64")
    )
    df["clinvar_found"] = df["clinvar_id"].astype(str) != ""
    df["clinvar_classified"] = df["clinvar_sig"].astype(str) != ""
    df["clinvar_category"] = categorize_clinvar_series(df["clinvar_sig"])
    add_titv_kind(df)
    return df


def explode_strategy_variants(variants: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "variant_id",
        "gene_id",
        "event_type",
        "strategies",
        "clinvar_found",
        "clinvar_classified",
        "clinvar_category",
        "gnomad_af",
        "titv_kind",
    ]
    long = variants[columns].copy()
    long["strategy"] = long["strategies"].astype(str).str.split(",")
    long = long.explode("strategy")
    long["strategy"] = long["strategy"].fillna("").astype(str).str.strip()
    long = long[long["strategy"] != ""].drop(columns=["strategies"])
    long = long.drop_duplicates(["strategy", "variant_id"])
    long = long.reset_index(drop=True)
    for column in ["strategy", "event_type", "gene_id"]:
        long[column] = long[column].astype("category")
    return long


def read_feature_coverage(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    print(f"Reading {path}...")
    cov = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    numeric_cols = [
        "length_bp",
        "ortholog_count",
        "orthologs_covered",
        "covered_bases",
        "coverage_breadth",
        "depth_bases",
        "mean_depth",
    ]
    for col in numeric_cols:
        if col in cov.columns:
            cov[col] = pd.to_numeric(cov[col], errors="coerce")
    return cov


def read_strategy_quick_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    print(f"Reading {path}...")
    quick = pd.read_csv(path, sep="\t", low_memory=False)
    for column in quick.columns:
        if column != "strategy":
            quick[column] = pd.to_numeric(quick[column], errors="coerce")
    return quick


def read_failures(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    failures = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
    return failures


def calc_titv_ratio(ti: int, tv: int) -> float:
    if tv == 0:
        return np.nan if ti == 0 else float("inf")
    return round(ti / tv, 3)


def titv_by_strategy(df: pd.DataFrame) -> pd.Series:
    counts = (
        df[df["titv_kind"].isin(["ti", "tv"])]
        .groupby(["strategy", "titv_kind"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    values = {}
    for strategy in df["strategy"].cat.categories if hasattr(df["strategy"], "cat") else sorted(df["strategy"].unique()):
        if strategy not in counts.index:
            values[strategy] = np.nan
            continue
        ti = int(counts.loc[strategy].get("ti", 0))
        tv = int(counts.loc[strategy].get("tv", 0))
        values[strategy] = calc_titv_ratio(ti, tv)
    return pd.Series(values, name="Ti/Tv")


def summarize_strategy_variants(long: pd.DataFrame, count_label: str = "Unique Variants") -> pd.DataFrame:
    if long.empty:
        return pd.DataFrame(
            columns=[
                "Strategy",
                count_label,
                "Ti/Tv",
                "Found in ClinVar",
                "ClinVar found %",
                "ClinVar classified",
                "ClinVar classified %",
                "gnomAD Found",
                "gnomAD found %",
                "P/LP",
                "B/LB",
                "VUS",
                "Other ClinVar",
                "Median gnomAD AF",
            ]
        )

    work = long.copy()
    work["gnomad_found"] = work["gnomad_af"].notna()
    grouped = work.groupby("strategy", observed=True)
    summary = grouped.agg(
        **{
            count_label: ("variant_id", "count"),
            "Genes": ("gene_id", "nunique"),
            "Found in ClinVar": ("clinvar_found", "sum"),
            "ClinVar classified": ("clinvar_classified", "sum"),
            "gnomAD Found": ("gnomad_found", "sum"),
            "Median gnomAD AF": ("gnomad_af", "median"),
        }
    )
    clinvar_counts = pd.crosstab(work["strategy"], work["clinvar_category"])
    for category in CLINVAR_ORDER:
        if category not in clinvar_counts.columns:
            clinvar_counts[category] = 0
    summary["Ti/Tv"] = titv_by_strategy(work)
    summary["P/LP"] = clinvar_counts["P/LP"]
    summary["B/LB"] = clinvar_counts["B/LB"]
    summary["VUS"] = clinvar_counts["VUS"]
    summary["Other ClinVar"] = clinvar_counts["Other"]
    summary["ClinVar found %"] = summary["Found in ClinVar"] / summary[count_label].replace(0, np.nan)
    summary["ClinVar classified %"] = summary["ClinVar classified"] / summary[count_label].replace(0, np.nan)
    summary["gnomAD found %"] = summary["gnomAD Found"] / summary[count_label].replace(0, np.nan)
    summary = summary.reset_index().rename(columns={"strategy": "Strategy"})
    summary["Strategy"] = summary["Strategy"].map(strategy_label)
    ordered = [
        "Strategy",
        count_label,
        "Genes",
        "Ti/Tv",
        "Found in ClinVar",
        "ClinVar found %",
        "ClinVar classified",
        "ClinVar classified %",
        "gnomAD Found",
        "gnomAD found %",
        "P/LP",
        "B/LB",
        "VUS",
        "Other ClinVar",
        "Median gnomAD AF",
    ]
    return summary[ordered].sort_values("Strategy")


def quick_summary_for_report(quick: pd.DataFrame) -> pd.DataFrame:
    if quick.empty:
        return pd.DataFrame()
    report = pd.DataFrame({"Strategy": quick["strategy"].map(strategy_label)})
    if {"aligned_rows", "rows"} <= set(quick.columns):
        report["Aligned orthologs %"] = quick["aligned_rows"] / quick["rows"].replace(0, np.nan)
        report["Aligned orthologs"] = quick["aligned_rows"]
    if "total_event_count" in quick.columns:
        report["Raw support events"] = quick["total_event_count"]
    if "aligned_target_bp" in quick.columns:
        report["Aligned target bp"] = quick["aligned_target_bp"]
    return report


def merge_quick_summary(strategy_stats: pd.DataFrame, quick: pd.DataFrame) -> pd.DataFrame:
    quick_report = quick_summary_for_report(quick)
    if quick_report.empty:
        return strategy_stats
    return strategy_stats.merge(quick_report, on="Strategy", how="left")


def unique_contribution_table(long: pd.DataFrame) -> pd.DataFrame:
    if long.empty:
        return summarize_strategy_variants(long, count_label="Unique To Strategy")
    all_strategies = pd.DataFrame({"Strategy": sorted(long["strategy"].astype(str).unique())})
    all_strategies["Strategy"] = all_strategies["Strategy"].map(strategy_label)
    strategy_counts = long.groupby("variant_id", sort=False)["strategy"].transform("nunique")
    unique_long = long[strategy_counts == 1]
    summary = summarize_strategy_variants(unique_long, count_label="Unique To Strategy")
    merged = all_strategies.merge(summary, on="Strategy", how="left")
    for column in merged.columns:
        if column != "Strategy":
            merged[column] = merged[column].fillna(0)
    return merged


def event_type_counts(long: pd.DataFrame) -> pd.DataFrame:
    if long.empty:
        return pd.DataFrame(columns=["strategy", "event_type", "Variant_Count"])
    counts = (
        long.groupby(["strategy", "event_type"], observed=True)
        .size()
        .reset_index(name="Variant_Count")
        .sort_values(["strategy", "event_type"])
    )
    counts["strategy"] = counts["strategy"].astype(str).map(strategy_label)
    return counts


def clinvar_counts(long: pd.DataFrame) -> pd.DataFrame:
    if long.empty:
        return pd.DataFrame(columns=["strategy", "clinvar_category", "Variant_Count"])
    counts = (
        long.groupby(["strategy", "clinvar_category"], observed=False)
        .size()
        .reset_index(name="Variant_Count")
    )
    counts["clinvar_category"] = pd.Categorical(
        counts["clinvar_category"], categories=CLINVAR_ORDER, ordered=True
    )
    counts["strategy"] = counts["strategy"].astype(str).map(strategy_label)
    return counts.sort_values(["strategy", "clinvar_category"])


def gnomad_found_counts(long: pd.DataFrame) -> pd.DataFrame:
    if long.empty:
        return pd.DataFrame(columns=["strategy", "gnomad_found", "Variant_Count"])
    counts = (
        long.assign(gnomad_found=long["gnomad_af"].notna())
        .groupby(["strategy", "gnomad_found"], observed=True)
        .size()
        .reset_index(name="Variant_Count")
    )
    counts["gnomad_found"] = counts["gnomad_found"].map({True: "Found", False: "Not Found"})
    counts["strategy"] = counts["strategy"].astype(str).map(strategy_label)
    return counts


def binned_gnomad_af(long: pd.DataFrame, bin_count: int = 200) -> pd.DataFrame:
    gnomad = long[long["gnomad_af"].notna() & (long["gnomad_af"] > 0)][["strategy", "gnomad_af"]].copy()
    if gnomad.empty:
        return pd.DataFrame(columns=["strategy", "bin_mid", "Variant_Count", "Density"])

    gnomad["log10_gnomad_af"] = np.log10(gnomad["gnomad_af"])
    min_value = float(gnomad["log10_gnomad_af"].min())
    max_value = float(gnomad["log10_gnomad_af"].max())
    if min_value == max_value:
        min_value -= 0.5
        max_value += 0.5
    bins = np.linspace(min_value, max_value, bin_count + 1)
    gnomad["bin"] = pd.cut(gnomad["log10_gnomad_af"], bins=bins, include_lowest=True)
    counts = gnomad.groupby(["strategy", "bin"], observed=True).size().reset_index(name="Variant_Count")
    counts["bin_mid"] = counts["bin"].apply(lambda value: value.mid).astype(float)
    totals = counts.groupby("strategy", observed=True)["Variant_Count"].transform("sum")
    counts["Density"] = counts["Variant_Count"] / totals.replace(0, np.nan)
    counts["strategy"] = counts["strategy"].astype(str).map(strategy_label)
    return counts[["strategy", "bin_mid", "Variant_Count", "Density"]]


def coverage_summary(cov: pd.DataFrame, feature_types: list[str] | None = None) -> pd.DataFrame:
    if cov.empty:
        return pd.DataFrame()
    if feature_types:
        cov = cov[cov["feature_type"].isin(feature_types)].copy()
    summary = (
        cov.groupby(["strategy", "feature_type"], as_index=False)
        .agg(
            Gene_Count=("gene_id", "nunique"),
            Feature_Count=("feature_id", "count"),
            Total_Length_bp=("length_bp", "sum"),
            Covered_Bases=("covered_bases", "sum"),
            Depth_Bases=("depth_bases", "sum"),
            Median_Breadth=("coverage_breadth", "median"),
            Median_Mean_Depth=("mean_depth", "median"),
            Median_Orthologs_Covered=("orthologs_covered", "median"),
        )
    )
    summary["Breadth_Weighted"] = summary["Covered_Bases"] / summary["Total_Length_bp"].replace(0, np.nan)
    summary["Mean_Depth_Weighted"] = summary["Depth_Bases"] / summary["Total_Length_bp"].replace(0, np.nan)
    summary["feature_type"] = pd.Categorical(summary["feature_type"], categories=FEATURE_ORDER, ordered=True)
    summary["strategy"] = summary["strategy"].astype(str).map(strategy_label)
    return summary.sort_values(["strategy", "feature_type"])


def fig_html(fig, include_plotlyjs: bool = False) -> str:
    return fig.to_html(full_html=False, include_plotlyjs="cdn" if include_plotlyjs else False)


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


def build_overview(
    variants: pd.DataFrame,
    long: pd.DataFrame,
    cov: pd.DataFrame,
    strategy_stats: pd.DataFrame,
    annotation_manifest: dict,
    alignment_manifest: dict,
) -> list[str]:
    unique_variant_count = variants["variant_id"].nunique()
    event_row_count = annotation_manifest.get("event_row_count") or alignment_manifest.get("raw_alignment_event_count") or ""
    clinvar_found = int(variants["clinvar_found"].sum())
    clinvar_classified = int(variants["clinvar_classified"].sum())
    gnomad_found = int(variants["gnomad_af"].notna().sum())
    annotation_warnings = int(annotation_manifest.get("failure_count", 0) or 0)
    cards = [
        ("Raw support events", format_int(event_row_count) if event_row_count != "" else "n/a"),
        ("Unique candidate variants", format_int(unique_variant_count)),
        ("Strategies", format_int(long["strategy"].nunique())),
        ("Genes", format_int(variants["gene_id"].nunique())),
        ("Found in ClinVar", f"{format_int(clinvar_found)} ({format_percent(clinvar_found / unique_variant_count)})"),
        (
            "ClinVar with CLNSIG",
            f"{format_int(clinvar_classified)} ({format_percent(clinvar_classified / unique_variant_count)})",
        ),
        ("Found in gnomAD", f"{format_int(gnomad_found)} ({format_percent(gnomad_found / unique_variant_count)})"),
        ("Annotation warnings", format_int(annotation_warnings)),
    ]
    sections = [metric_cards(cards)]
    sections.append("<h2>Strategy Summary</h2>")
    sections.append(table_html(strategy_stats))
    return sections


def build_variant_sections(
    long: pd.DataFrame,
    strategy_stats: pd.DataFrame,
    include_plotly: bool,
) -> list[str]:
    sections = ["<h2>Variant Profile</h2>"]
    variant_volume = sort_by_metric(strategy_stats[["Strategy", "Unique Variants"]], "Unique Variants")
    fig_volume = px.bar(
        variant_volume,
        x="Strategy",
        y="Unique Variants",
        title="Unique candidate variants by strategy",
        category_orders={"Strategy": variant_volume["Strategy"].tolist()},
    )
    compact_figure(fig_volume)
    sections.append(fig_html(fig_volume, include_plotlyjs=include_plotly))

    titv = sort_by_metric(strategy_stats[["Strategy", "Ti/Tv"]], "Ti/Tv")
    fig_titv = px.bar(
        titv,
        x="Strategy",
        y="Ti/Tv",
        title="Ti/Tv by strategy",
        category_orders={"Strategy": titv["Strategy"].tolist()},
    )
    compact_figure(fig_titv)
    sections.append(fig_html(fig_titv))

    unique_contrib = unique_contribution_table(long)
    unique_contrib_plot = sort_by_metric(unique_contrib[["Strategy", "Unique To Strategy"]], "Unique To Strategy")
    fig_unique = px.bar(
        unique_contrib_plot,
        x="Strategy",
        y="Unique To Strategy",
        title="Variants found only by one strategy",
        category_orders={"Strategy": unique_contrib_plot["Strategy"].tolist()},
    )
    compact_figure(fig_unique)
    sections.append(fig_html(fig_unique))

    counts = event_type_counts(long)
    totals = counts.groupby("strategy", observed=True)["Variant_Count"].transform("sum")
    counts["Fraction"] = counts["Variant_Count"] / totals.replace(0, np.nan)
    snv_order = (
        counts[counts["event_type"].astype(str).str.lower() == "snv"]
        .sort_values("Fraction", ascending=False)
        ["strategy"]
        .tolist()
    )
    order = snv_order + [strategy for strategy in variant_volume["Strategy"].tolist() if strategy not in snv_order]
    fig_events = px.bar(
        counts,
        x="strategy",
        y="Fraction",
        color="event_type",
        barmode="stack",
        title="Variant type composition by strategy",
        category_orders={"strategy": order},
        labels={"strategy": "", "Fraction": "Variant fraction", "event_type": "Variant type"},
    )
    fig_events.update_layout(yaxis_tickformat=".0%")
    compact_figure(fig_events, height=360)
    sections.append(fig_html(fig_events))
    return sections


def build_clinvar_gnomad_sections(
    variants: pd.DataFrame,
    long: pd.DataFrame,
    strategy_stats: pd.DataFrame,
    include_plotly: bool,
) -> list[str]:
    sections = ["<h2>External Evidence</h2>"]
    sections.append(
        metric_cards(
            [
                ("Found in ClinVar", format_int(variants["clinvar_found"].sum())),
                ("ClinVar with CLNSIG", format_int(variants["clinvar_classified"].sum())),
                ("Found in gnomAD", format_int(variants["gnomad_af"].notna().sum())),
            ]
        )
    )

    clinvar_rate = sort_by_metric(strategy_stats[["Strategy", "ClinVar found %"]], "ClinVar found %")
    fig_clin_rate = px.bar(
        clinvar_rate,
        x="Strategy",
        y="ClinVar found %",
        title="ClinVar hit rate by strategy",
        category_orders={"Strategy": clinvar_rate["Strategy"].tolist()},
    )
    fig_clin_rate.update_layout(yaxis_tickformat=".2%")
    compact_figure(fig_clin_rate)
    sections.append(fig_html(fig_clin_rate, include_plotlyjs=include_plotly))

    gnomad_rate = sort_by_metric(strategy_stats[["Strategy", "gnomAD found %"]], "gnomAD found %")
    fig_gnomad_rate = px.bar(
        gnomad_rate,
        x="Strategy",
        y="gnomAD found %",
        title="gnomAD hit rate by strategy",
        category_orders={"Strategy": gnomad_rate["Strategy"].tolist()},
    )
    fig_gnomad_rate.update_layout(yaxis_tickformat=".1%")
    compact_figure(fig_gnomad_rate)
    sections.append(fig_html(fig_gnomad_rate))

    clin_counts = clinvar_counts(long)
    clin_plot = clin_counts[clin_counts["clinvar_category"].astype(str) != "Not Found"].copy()
    totals = clin_plot.groupby("strategy", observed=True)["Variant_Count"].transform("sum")
    clin_plot["Fraction"] = clin_plot["Variant_Count"] / totals.replace(0, np.nan)
    clin_order = (
        clin_plot[clin_plot["clinvar_category"].astype(str) == "B/LB"]
        .sort_values("Fraction", ascending=False)
        ["strategy"]
        .tolist()
    )
    clin_order += [strategy for strategy in clinvar_rate["Strategy"].tolist() if strategy not in clin_order]
    fig_clin = px.bar(
        clin_plot,
        x="strategy",
        y="Fraction",
        color="clinvar_category",
        barmode="stack",
        title="ClinVar classification mix among classified variants",
        category_orders={"strategy": clin_order, "clinvar_category": CLINVAR_ORDER},
        color_discrete_map=CLINVAR_COLORS,
        labels={"strategy": "", "Fraction": "ClinVar class fraction", "clinvar_category": "ClinVar class"},
    )
    fig_clin.update_layout(yaxis_tickformat=".0%")
    compact_figure(fig_clin, height=360)
    sections.append(fig_html(fig_clin))

    gnomad_bins = binned_gnomad_af(long)
    if not gnomad_bins.empty:
        fig_af = px.bar(
            gnomad_bins,
            x="bin_mid",
            y="Density",
            color="strategy",
            barmode="overlay",
            opacity=0.65,
            title="gnomAD AF Distribution by Strategy",
            labels={"bin_mid": "log10 gnomAD AF", "Density": "Within-strategy density", "strategy": ""},
        )
        fig_af.update_layout(yaxis_title="Within-strategy density", xaxis_title="log10 gnomAD AF")
        fig_af.update_traces(marker_line_width=0)
        compact_figure(fig_af, height=380, show_x_title=True)
        sections.append("<h3>gnomAD AF Distribution</h3>")
        sections.append(fig_html(fig_af))
    else:
        sections.append("<p>No non-zero gnomAD AF values were found.</p>")

    return sections


def build_feature_sections(cov: pd.DataFrame, include_plotly: bool) -> list[str]:
    sections = ["<h2>Target Feature Coverage</h2>"]
    if cov.empty:
        sections.append("<p>No feature coverage table was found.</p>")
        return sections

    disjoint_summary = coverage_summary(cov, DISJOINT_FEATURE_ORDER)
    if disjoint_summary.empty:
        sections.append("<p>No CDS/UTR/intron coverage rows were found.</p>")
        return sections

    breadth_cards = []
    for feature_type in DISJOINT_FEATURE_ORDER:
        feature = disjoint_summary[disjoint_summary["feature_type"].astype(str) == feature_type]
        if feature.empty:
            continue
        min_breadth = feature["Breadth_Weighted"].min()
        max_breadth = feature["Breadth_Weighted"].max()
        breadth_cards.append((f"{feature_type.upper()} breadth range", f"{format_percent(min_breadth)}-{format_percent(max_breadth)}"))
    sections.append(metric_cards(breadth_cards))

    cds_depth = (
        disjoint_summary[disjoint_summary["feature_type"].astype(str) == "cds"]
        .sort_values("Mean_Depth_Weighted", ascending=False)
    )
    strategy_order = cds_depth["strategy"].tolist() or sorted(disjoint_summary["strategy"].unique())

    fig_depth = px.bar(
        disjoint_summary,
        x="strategy",
        y="Mean_Depth_Weighted",
        color="feature_type",
        barmode="group",
        title="Weighted mean ortholog depth by target feature",
        category_orders={"strategy": strategy_order, "feature_type": DISJOINT_FEATURE_ORDER},
        labels={
            "strategy": "",
            "Mean_Depth_Weighted": "Weighted mean ortholog depth",
            "feature_type": "Feature",
        },
    )
    compact_figure(fig_depth, height=360)
    sections.append(fig_html(fig_depth, include_plotlyjs=include_plotly))
    return sections


def build_methods_sections(
    inputs: RunInputs,
    out_html: Path,
    variants: pd.DataFrame,
    long: pd.DataFrame,
    cov: pd.DataFrame,
    failures: pd.DataFrame,
    annotation_manifest: dict,
    alignment_manifest: dict,
) -> list[str]:
    files = [
        ("Run Dir", inputs.run_dir),
        ("Variant Annotations", inputs.variant_annotations_tsv),
        ("Feature Coverage", inputs.feature_coverage_tsv),
        ("Strategy Quick Summary", inputs.strategy_quick_summary_tsv),
        ("Annotation Manifest", inputs.annotation_manifest_json),
        ("Alignment Manifest", inputs.alignment_manifest_json),
        ("Output HTML", out_html),
    ]
    file_rows = [
        {"Key": label, "Path": str(path), "Exists": path.exists(), "Size": file_size_label(path)}
        for label, path in files
    ]

    ok_events = int(annotation_manifest.get("event_key_status_counts", {}).get("ok", 0))
    missing_left_anchor = int(annotation_manifest.get("event_key_status_counts", {}).get("missing_left_anchor", 0))
    sections = [
        "<h2>QC</h2>",
        metric_cards(
            [
                ("Event keys normalized", format_int(ok_events)),
                ("Missing left anchor", format_int(missing_left_anchor)),
                ("gnomAD regions failed", format_int(annotation_manifest.get("gnomad_region_failure_count", 0))),
                ("ClinVar cached variants", format_int(annotation_manifest.get("clinvar_cached_variant_count", 0))),
                ("gnomAD cached variants", format_int(annotation_manifest.get("gnomad_cached_variant_count", 0))),
                ("Feature coverage rows", format_int(len(cov))),
            ]
        ),
    ]
    if not failures.empty:
        sections.append("<h3>Annotation Failures</h3>")
        sections.append(table_html(failures, classes="table table-sm table-striped", max_rows=50))
    sections.append("<details><summary>Input files and loaded row counts</summary>")
    sections.append(table_html(pd.DataFrame(file_rows), classes="table table-sm table-striped"))
    sections.append(
        table_html(
            pd.DataFrame(
                [
                    {"Metric": "Unique candidate variants loaded", "Value": len(variants)},
                    {"Metric": "Strategy-supported variant records loaded", "Value": len(long)},
                    {"Metric": "Feature coverage rows loaded", "Value": len(cov)},
                    {"Metric": "Annotation failure rows", "Value": len(failures)},
                    {"Metric": "Alignment event mode", "Value": alignment_manifest.get("alignment_event_mode", "")},
                ]
            ),
            classes="table table-sm table-striped",
        )
    )
    sections.append("</details>")
    return sections


def render_tabs(sections: list[tuple[str, str, list[str]]]) -> str:
    buttons = []
    pages = []
    for index, (tab_id, title, html_parts) in enumerate(sections):
        active = " active" if index == 0 else ""
        buttons.append(f'<button class="tab-button{active}" data-tab="{tab_id}">{title}</button>')
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
            window.dispatchEvent(new Event('resize'));
        }});
    }});
    </script>
    """


def render_html(sections: list[tuple[str, str, list[str]]]) -> str:
    return f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>Alignment Strategies Report</title>
        <style>
            body {{
                padding: 20px;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
                color: #1f2933;
            }}
            h1 {{ margin-bottom: 4px; }}
            h2 {{ margin-top: 22px; border-bottom: 1px solid #d5d9df; padding-bottom: 6px; }}
            h3 {{ margin-top: 16px; }}
            .lead {{ margin-top: 0; color: #52606d; }}
            .tab-bar {{
                display: flex;
                flex-wrap: wrap;
                gap: 8px;
                margin: 16px 0;
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
                gap: 10px;
                margin: 12px 0 18px 0;
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
            th, td {{ border: 1px solid #d5d9df; padding: 6px 8px; text-align: center; }}
            th {{ background: #f5f7fa; }}
            td:first-child, th:first-child {{ text-align: left; }}
            details {{
                margin: 16px 0;
                border: 1px solid #d5d9df;
                border-radius: 6px;
                padding: 10px 12px;
                background: #fbfcfd;
            }}
            summary {{ cursor: pointer; font-weight: 600; }}
            .plotly-graph-div {{ min-height: 300px; }}
        </style>
    </head>
    <body>
        <h1>Alignment Strategies Report</h1>
        <p class="lead">A compact comparison of alignment support, candidate variant profiles, external evidence, and target-feature coverage.</p>
        {render_tabs(sections)}
    </body>
    </html>
    """


def main() -> None:
    args = parse_args()
    inputs = resolve_run_inputs(args.run_dir)
    out_html = resolve_out_html(args, inputs.run_dir)

    variants = read_variant_annotations(inputs.variant_annotations_tsv)
    long = explode_strategy_variants(variants)
    cov = read_feature_coverage(inputs.feature_coverage_tsv)
    quick = read_strategy_quick_summary(inputs.strategy_quick_summary_tsv)
    failures = read_failures(inputs.annotation_failures_tsv)
    annotation_manifest = read_json(inputs.annotation_manifest_json)
    alignment_manifest = read_json(inputs.alignment_manifest_json)

    print("Computing strategy metrics...")
    strategy_stats_full = merge_quick_summary(summarize_strategy_variants(long), quick)
    summary_columns = [
        "Strategy",
        "Unique Variants",
        "Ti/Tv",
        "Found in ClinVar",
        "ClinVar found %",
        "gnomAD Found",
        "gnomAD found %",
        "Aligned orthologs %",
        "Raw support events",
    ]
    strategy_stats = strategy_stats_full[[column for column in summary_columns if column in strategy_stats_full.columns]]

    sections = [
        ("overview", "Overview", build_overview(variants, long, cov, strategy_stats, annotation_manifest, alignment_manifest)),
        ("variants", "Variant Profile", build_variant_sections(long, strategy_stats, include_plotly=True)),
        (
            "external-evidence",
            "External Evidence",
            build_clinvar_gnomad_sections(variants, long, strategy_stats_full, include_plotly=True),
        ),
        ("coverage", "Feature Coverage", build_feature_sections(cov, include_plotly=True)),
        (
            "qc",
            "QC",
            build_methods_sections(
                inputs,
                out_html,
                variants,
                long,
                cov,
                failures,
                annotation_manifest,
                alignment_manifest,
            ),
        ),
    ]

    print(f"Writing report to {out_html}...")
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(render_html(sections))
    print("Done!")


if __name__ == "__main__":
    main()
