#!/usr/bin/env python3
"""Build an HTML report for one completed GAPH run."""

from __future__ import annotations

import argparse
import json
import math
import re
import time
import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from analytics.core.clinvar_validation import build_validation
from analytics.core.conservation import DEFAULT_TRACK_NAMES, build_conservation_annotations
from analytics.core.conservation_validation import (
    CONSEQUENCE_OPTIONS,
    CONSEQUENCE_TERMS,
    PHYLOP_BANDS,
    SCORE_COLUMN,
    SPLINE_DF,
    VARIANT_TYPE_OPTIONS,
    ConservationValidation,
    build_conservation_cohort,
    compute_conservation_validation,
)
from analytics.core.negative_controls import TargetSpaceNullAnalysis, build_target_space_null
from analytics.core.variant_summary import StrategyOverlap, VariantSummary, build_variant_summary


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
REVIEW_STAR_ORDER = ["4", "3", "2", "1", "0", "Unmapped"]
REVIEW_STAR_COLORS = {
    "4": "#08519c",
    "3": "#3182bd",
    "2": "#6baed6",
    "1": "#9ecae1",
    "0": "#fdbb84",
    "Unmapped": "#bdbdbd",
}
CONSEQUENCE_GROUP_ORDER = ["LoF/splice", "Missense/inframe", "Synonymous", "Noncoding/UTR/intron", "Other"]
CONSEQUENCE_GROUP_COLORS = {
    "LoF/splice": "#de2d26",
    "Missense/inframe": "#fb6a4a",
    "Synonymous": "#74add1",
    "Noncoding/UTR/intron": "#abd9e9",
    "Other": "#9e9e9e",
}
CONSEQUENCE_GROUP_TERMS = {
    "LoF/splice": [
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "start_lost",
        "stop_gained",
        "stop_lost",
    ],
    "Missense/inframe": [
        "inframe_deletion",
        "inframe_insertion",
        "missense_variant",
        "protein_altering_variant",
    ],
    "Synonymous": [
        "stop_retained_variant",
        "synonymous_variant",
    ],
    "Noncoding/UTR/intron": [
        "3_prime_UTR_variant",
        "5_prime_UTR_variant",
        "intron_variant",
        "non_coding_transcript_exon_variant",
        "splice_region_variant",
    ],
}
STRATEGY_LABELS = {
    "bwa_pseudoreads": "BWA pseudo",
    "minimap2_asm10": "minimap2 asm10",
    "minimap2_asm20": "minimap2 asm20",
    "minimap2_taxonomy_adaptive": "minimap2 adaptive",
    "nucmer": "nucmer",
    "precomputed_ensembl_92_mammals_epo_extended": "Ensembl EPO",
}

@dataclass(frozen=True)
class RunInputs:
    run_dir: Path
    genes_tsv: Path
    target_features_tsv: Path
    target_sequences_dir: Path
    variant_annotations_tsv: Path
    annotation_manifest_json: Path
    annotation_failures_tsv: Path
    feature_coverage_tsv: Path
    alignment_segments_tsv: Path
    alignment_manifest_json: Path
    strategy_summary_tsv: Path


@dataclass(frozen=True)
class ConservationAnalysis:
    annotations_path: Path
    manifest_path: Path
    manifest: dict
    validation: ConservationValidation


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path, help="Completed GAPH run directory.")
    parser.add_argument(
        "--clinvar-vcf",
        type=Path,
        default=project_root() / "assets" / "reference" / "clinvar" / "clinvar.vcf.gz",
        help="Indexed ClinVar VCF used for validation. Default: assets/reference/clinvar/clinvar.vcf.gz",
    )
    parser.add_argument(
        "--out-html",
        type=Path,
        help="Output HTML path. Default: <run-dir>/reports/strategy_compare.html",
    )
    parser.add_argument(
        "--report-name",
        help="Short report file name inside <run-dir>/reports. '.html' is added if omitted.",
    )
    parser.add_argument(
        "--target-space-null",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Build the consequence-matched target-space null. Disabled by default because it uses Ensembl REST VEP.",
    )
    parser.add_argument(
        "--target-space-null-sample-size",
        type=int,
        default=25_000,
        help="Maximum deterministic focal-SNV sample per strategy for the target-space null.",
    )
    parser.add_argument(
        "--target-space-null-resamples",
        type=int,
        default=1_000,
        help="Target-space-null resampling iterations. Default: 1000.",
    )
    parser.add_argument(
        "--target-space-null-seed",
        type=int,
        default=20_260_721,
        help="Deterministic target-space-null seed.",
    )
    return parser.parse_args()


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


@contextmanager
def timed_stage(name: str, timings: list[dict[str, object]]) -> Iterator[dict[str, object]]:
    entry: dict[str, object] = {"Stage": name, "Status": "completed", "Details": ""}
    started = time.perf_counter()
    try:
        yield entry
    except Exception:
        entry["Status"] = "failed"
        raise
    finally:
        entry["Seconds"] = round(time.perf_counter() - started, 3)
        timings.append(entry)
        print(f"{name}: {entry['Status']} in {entry['Seconds']:.3f} s")


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
        genes_tsv=run_dir / "fetch" / "genes.tsv.gz",
        target_features_tsv=run_dir / "fetch" / "target_features.tsv.gz",
        target_sequences_dir=run_dir / "fetch" / "sequences" / "targets",
        variant_annotations_tsv=run_dir / "annotation" / "variant_annotations.tsv.gz",
        annotation_manifest_json=run_dir / "annotation" / "manifest.json",
        annotation_failures_tsv=run_dir / "annotation" / "failures.tsv.gz",
        feature_coverage_tsv=run_dir / "alignment" / "feature_coverage.tsv.gz",
        alignment_segments_tsv=run_dir / "alignment" / "alignment_segments.tsv.gz",
        alignment_manifest_json=run_dir / "alignment" / "manifest.json",
        strategy_summary_tsv=run_dir / "alignment" / "strategy_summary.tsv.gz",
    )
    if not inputs.variant_annotations_tsv.exists():
        raise FileNotFoundError(
            "Missing annotation/variant_annotations.tsv.gz under --run-dir. "
            "Run the annotation stage before building this report."
        )
    if not inputs.genes_tsv.exists():
        raise FileNotFoundError("Missing fetch/genes.tsv.gz under --run-dir.")
    if not inputs.target_features_tsv.exists():
        raise FileNotFoundError("Missing fetch/target_features.tsv.gz under --run-dir.")
    if not inputs.target_sequences_dir.exists():
        raise FileNotFoundError("Missing fetch/sequences/targets under --run-dir.")
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


def format_pvalue(value) -> str:
    if pd.isna(value):
        return ""
    value = float(value)
    if value == 0:
        return "0"
    if value < 0.001:
        return f"{value:.2e}"
    return f"{value:.3g}"


def format_ratio(value) -> str:
    if pd.isna(value):
        return ""
    value = float(value)
    if np.isposinf(value):
        return "inf"
    if np.isneginf(value):
        return "-inf"
    return format_float(value, 3)


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


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


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


def read_strategy_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError("Missing alignment/strategy_summary.tsv.gz under --run-dir.")
    print(f"Reading {path}...")
    summary = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    required = {"strategy", "summary_row_count", "aligned_summary_row_count", "event_count"}
    missing = required - set(summary.columns)
    if missing:
        raise ValueError(f"Strategy summary missing required columns: {', '.join(sorted(missing))}")
    for column in summary.columns:
        if column != "strategy":
            summary[column] = pd.to_numeric(summary[column], errors="coerce")
    return summary


def read_failures(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    failures = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
    return failures


def alignment_summary_for_report(summary: pd.DataFrame) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame()
    report = pd.DataFrame({"Strategy": summary["strategy"].map(strategy_label)})
    report["Aligned support records %"] = (
        summary["aligned_summary_row_count"] / summary["summary_row_count"].replace(0, np.nan)
    )
    report["Aligned support records"] = summary["aligned_summary_row_count"]
    report["Raw support events"] = summary["event_count"]
    if "aligned_target_bp" in summary.columns:
        report["Aligned target bp"] = summary["aligned_target_bp"]
    return report


def merge_alignment_summary(strategy_stats: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    report_summary = alignment_summary_for_report(summary)
    if report_summary.empty:
        return strategy_stats
    return strategy_stats.merge(report_summary, on="Strategy", how="left")


def strategy_overlap_figure(overlap: StrategyOverlap | None):
    if overlap is None:
        return None
    labels = [strategy_label(strategy) for strategy in overlap.strategies]

    fig = go.Figure(
        data=go.Heatmap(
            z=overlap.jaccard,
            x=labels,
            y=labels,
            text=np.vectorize(lambda value: f"{value:.0%}")(overlap.jaccard),
            customdata=np.dstack([overlap.intersections, overlap.unions]),
            colorscale="Blues",
            zmin=0,
            zmax=1,
            colorbar={"title": "Jaccard"},
            hovertemplate=(
                "%{y} vs %{x}<br>"
                "Jaccard: %{z:.1%}<br>"
                "Shared variants: %{customdata[0]:,}<br>"
                "Union variants: %{customdata[1]:,}<extra></extra>"
            ),
        )
    )
    fig.update_traces(texttemplate="%{text}", textfont_size=12)
    fig.update_layout(
        title="Pairwise strategy overlap",
        height=520,
        margin={"l": 135, "r": 30, "t": 95, "b": 35},
        template="plotly_white",
    )
    fig.update_xaxes(side="top", tickangle=-35, title_text=None, automargin=True)
    fig.update_yaxes(title_text=None, automargin=True)
    return fig


def review_star_category(row: pd.Series) -> str:
    stars = str(row.get("clinvar_review_stars", "") or "").strip()
    if stars in {"0", "1", "2", "3", "4"}:
        return stars
    return "Unmapped"


def consequence_group(value: str) -> str:
    consequence = str(value or "")
    for group, terms in CONSEQUENCE_GROUP_TERMS.items():
        if consequence in terms:
            return group
    return "Other"


def consequence_grouping_table() -> pd.DataFrame:
    rows = [
        {
            "Group": group,
            "gnomAD consequence values": ", ".join(CONSEQUENCE_GROUP_TERMS.get(group, []))
            if group != "Other"
            else "Any non-empty gnomAD consequence not listed above.",
        }
        for group in CONSEQUENCE_GROUP_ORDER
    ]
    return pd.DataFrame(rows)


def clinvar_class_mapping_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Class": "P/LP",
                "Rule": "CLNSIG contains pathogenic and does not contain benign, uncertain/VUS, or conflicting.",
                "Used for": "ClinVar class composition and pathogenic-only evidence plots.",
            },
            {
                "Class": "B/LB",
                "Rule": "CLNSIG contains benign and does not contain pathogenic, uncertain/VUS, or conflicting.",
                "Used for": "ClinVar class composition and validation enrichment counts.",
            },
            {
                "Class": "VUS",
                "Rule": "CLNSIG contains uncertain or VUS, unless the record is marked conflicting.",
                "Used for": "Shown as uncertainty in class composition; excluded from validation denominator.",
            },
            {
                "Class": "Other",
                "Rule": "Conflicting, mixed, ambiguous, or non-empty CLNSIG values outside the clean classes above.",
                "Used for": "Shown separately in class composition; excluded from validation denominator.",
            },
            {
                "Class": "Not Found",
                "Rule": "No CLNSIG value in the annotation row.",
                "Used for": "Not included in classified-variant plots or validation denominator.",
            },
        ]
    )


def clinvar_review_star_mapping_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Stars": "4",
                "ClinVar review status values": "practice_guideline",
                "Interpretation": "Practice guideline.",
            },
            {
                "Stars": "3",
                "ClinVar review status values": "reviewed_by_expert_panel",
                "Interpretation": "Reviewed by expert panel.",
            },
            {
                "Stars": "2",
                "ClinVar review status values": (
                    "criteria_provided,_multiple_submitters,_no_conflicts; "
                    "criteria_provided,_multiple_submitters"
                ),
                "Interpretation": "Multiple submitters with criteria; no-conflict status when provided by ClinVar.",
            },
            {
                "Stars": "1",
                "ClinVar review status values": (
                    "criteria_provided,_single_submitter; "
                    "criteria_provided,_conflicting_classifications; "
                    "criteria_provided,_conflicting_interpretations"
                ),
                "Interpretation": "Criteria provided, but lower review confidence or conflicting submissions.",
            },
            {
                "Stars": "0",
                "ClinVar review status values": (
                    "no_assertion_criteria_provided; no_assertion_provided; "
                    "no_classification_provided; no_classification_for_the_individual_variant"
                ),
                "Interpretation": "No assertion criteria or no classification.",
            },
            {
                "Stars": "Unmapped",
                "ClinVar review status values": "Missing, empty, or unrecognized review status.",
                "Interpretation": "Kept visible only when such records are present.",
            },
        ]
    )


def validation_method_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Step": "Universe",
                "Definition": (
                    "ClinVar alleles overlapping fetched target loci, normalized to the same target-context "
                    "variant_key representation as GAPH annotations."
                ),
            },
            {
                "Step": "Variant types",
                "Definition": "SNV and INDEL validation are computed separately; complex/MNV/symbolic alleles are excluded.",
            },
            {
                "Step": "Included labels",
                "Definition": "Only clean B/LB and clean P/LP ClinVar labels enter the validation denominator.",
            },
            {
                "Step": "Excluded labels",
                "Definition": "VUS, missing CLNSIG, conflicting/other, mixed B/P labels, and unnormalizable alleles.",
            },
            {
                "Step": "Observed",
                "Definition": "A ClinVar allele is observed for a strategy when its normalized variant_key is present in that strategy's variant_annotations rows.",
            },
            {
                "Step": "2x2 table",
                "Definition": "[B/LB observed, P/LP observed; B/LB not observed, P/LP not observed].",
            },
            {
                "Step": "Statistics",
                "Definition": (
                    "Raw odds ratio, approximate 95% CI on log(OR) with Haldane 0.5 correction for zero cells, "
                    "two-sided Fisher exact p-value, and Benjamini-Hochberg FDR within each SNV or INDEL family."
                ),
            },
        ]
    )


def validation_consequence_grouping_table() -> pd.DataFrame:
    rows = []
    for key, label in CONSEQUENCE_OPTIONS:
        if key == "all":
            continue
        terms = CONSEQUENCE_TERMS.get(key)
        rows.append(
            {
                "Consequence subset": label,
                "ClinVar MC terms": ", ".join(sorted(terms))
                if terms
                else "Missing MC or any MC term not assigned to a named subset.",
            }
        )
    return pd.DataFrame(rows)


def conservation_validation_method_table(analysis: ConservationAnalysis | None) -> pd.DataFrame:
    versions = analysis.validation.r_versions if analysis is not None else {}
    return pd.DataFrame(
        [
            {
                "Step": "Variant set",
                "Definition": (
                    "Clean B/LB and P/LP ClinVar SNVs and simple INDELs from the normalized target-locus universe. "
                    "ALT_observed=0 means that the strategy did not report that exact normalized ALT; no callability "
                    "filter or covariate is applied."
                ),
            },
            {
                "Step": "Consequence subsets",
                "Definition": (
                    "Subsets use ClinVar MC Sequence Ontology terms. A record with terms from multiple groups enters "
                    "each matching group once and also enters All consequences; consequence-view counts are not additive."
                ),
            },
            {
                "Step": "Conservation annotation",
                "Definition": (
                    "phyloP100way is read from the hg38 UCSC bigWig. SNVs use the substituted base; deletions use "
                    "the mean across deleted reference bases excluding the VCF padding base; insertions use the mean "
                    "of the two flanking reference bases. All required bases must have a score."
                ),
            },
            {
                "Step": "Fixed bands",
                "Definition": (
                    "<= -1.30103 nominal acceleration band; (-1.30103, 1.30103) central band; "
                    ">= 1.30103 nominal conservation band. The cutoffs equal signed -log10(0.05) "
                    "for a single-base phyloP score and are descriptive, not genome-wide significance claims."
                ),
            },
            {
                "Step": "Fixed-band statistics",
                "Definition": (
                    "Each band receives a B/LB-vs-P/LP ALT-observed 2x2 table, OR, approximate 95% CI, and two-sided "
                    "Fisher test. A Mantel-Haenszel common OR and CMH test summarize across bands. This is a "
                    "sensitivity analysis because residual phyloP differences can remain within a band."
                ),
            },
            {
                "Step": "Continuous analysis",
                "Definition": (
                    f"Firth logistic regression: logit P(B/LB) = intercept + beta_ALT*ALT_observed + natural "
                    f"spline(phyloP100way, df={SPLINE_DF}). exp(beta_ALT) is the adjusted OR. The 95% CI and p-value "
                    "use profile penalized likelihood, not a Wald approximation."
                ),
            },
            {
                "Step": "Continuous estimability",
                "Definition": (
                    "Both clinical classes, both ALT-observed groups, at least four distinct scores, and overlap "
                    "between the groups' observed phyloP ranges are required. Complete outcome separation is handled "
                    "by Firth penalization rather than by discarding the view."
                ),
            },
            {
                "Step": "Multiplicity",
                "Definition": (
                    "Benjamini-Hochberg correction is applied as three prespecified families across displayed "
                    "selector combinations: per-band Fisher tests, pooled CMH tests, and continuous-model PLR tests."
                ),
            },
            {
                "Step": "INDEL interpretation",
                "Definition": (
                    "The fixed thresholds have their nominal single-base p-value interpretation only for SNVs. "
                    "INDEL, insertion, deletion, and All-variant views apply the same bands to an aggregate score "
                    "for descriptive comparability."
                ),
            },
            {
                "Step": "Software",
                "Definition": (
                    f"R {versions.get('R', 'not recorded')}; logistf {versions.get('logistf', 'not recorded')}. "
                    "The report does not fall back to ordinary maximum-likelihood logistic regression."
                ),
            },
        ]
    )


def negative_control_method_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Step": "Focal sample",
                "Definition": (
                    "Normalized GAPH SNVs are selected by a stable hash independently per strategy, up to the "
                    "configured engineering cap. All eligible SNVs are used below the cap."
                ),
            },
            {
                "Step": "Target-space matching",
                "Definition": (
                    "Each focal SNV is matched to up to five SNVs from the same gene and target context, with the "
                    "same genomic REF>ALT substitution and the same primary RefSeq VEP consequence. A control is "
                    "excluded when the same GAPH strategy observed it."
                ),
            },
            {
                "Step": "VEP consequence",
                "Definition": (
                    "Ensembl REST VEP uses RefSeq transcripts and pick_allele_gene. The target Entrez Gene ID is "
                    "selected; the most severe Sequence Ontology term on the picked transcript is the matching key."
                ),
            },
            {
                "Step": "Outcome",
                "Definition": (
                    "phyloP100way is compared as a focal-weighted ECDF and a descriptive median summary. "
                    "Conservation and GAPH callability are not matching variables, and no inferential p-value is reported."
                ),
            },
        ]
    )


def feature_coverage_formula_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Metric": "Per-feature breadth",
                "Formula": "coverage_breadth = covered_bases / length_bp",
                "Notes": "Read directly from alignment/feature_coverage.tsv.gz.",
            },
            {
                "Metric": "Weighted breadth",
                "Formula": "sum(covered_bases) / sum(length_bp)",
                "Notes": "Length-weighted aggregate used for feature breadth ranges.",
            },
            {
                "Metric": "Per-feature mean depth",
                "Formula": "mean_depth = depth_bases / length_bp",
                "Notes": "Uses pipeline-provided depth_bases; the report does not recompute depth from raw aligner output.",
            },
            {
                "Metric": "Weighted mean ortholog depth",
                "Formula": "sum(depth_bases) / sum(length_bp)",
                "Notes": "Main Feature Coverage plot metric.",
            },
            {
                "Metric": "Median feature metrics",
                "Formula": "median(coverage_breadth), median(mean_depth), median(orthologs_covered)",
                "Notes": "Computed over feature rows within each strategy and feature type.",
            },
            {
                "Metric": "Main feature classes",
                "Formula": "CDS, UTR, intron",
                "Notes": "The main plot uses disjoint target feature classes to avoid exon/gene aggregates dominating interpretation.",
            },
        ]
    )


def group_consequence_counts(raw_counts: pd.DataFrame) -> pd.DataFrame:
    if raw_counts.empty:
        return pd.DataFrame(columns=["Strategy", "Consequence group", "Variant_Count", "Fraction"])
    counts = raw_counts.rename(columns={"strategy": "Strategy"}).copy()
    counts["Consequence group"] = counts["value"].map(consequence_group)
    counts = (
        counts.groupby(["Strategy", "Consequence group"], observed=True)["Variant_Count"]
        .sum()
        .reset_index()
    )
    totals = counts.groupby("Strategy", observed=True)["Variant_Count"].transform("sum")
    counts["Fraction"] = counts["Variant_Count"] / totals.replace(0, np.nan)
    counts["Consequence group"] = pd.Categorical(
        counts["Consequence group"], categories=CONSEQUENCE_GROUP_ORDER, ordered=True
    )
    return counts.sort_values(["Strategy", "Consequence group"])


def consequence_strategy_order(counts: pd.DataFrame) -> list[str]:
    if counts.empty:
        return []
    pivot = counts.pivot_table(
        index="Strategy",
        columns="Consequence group",
        values="Fraction",
        aggfunc="sum",
        fill_value=0,
        observed=True,
    )
    for column in CONSEQUENCE_GROUP_ORDER:
        if column not in pivot.columns:
            pivot[column] = 0.0
    pivot["impact_fraction"] = pivot["LoF/splice"] + pivot["Missense/inframe"]
    pivot["total_count"] = counts.groupby("Strategy", observed=True)["Variant_Count"].sum()
    return pivot.sort_values(["impact_fraction", "total_count"], ascending=False).index.tolist()


def compact_list_text(value: str, max_items: int = 2, max_chars: int = 90) -> str:
    items = [item for item in re.split(r"[|,]", str(value or "")) if item and item != "."]
    if not items:
        return ""
    shown = "; ".join(items[:max_items])
    if len(items) > max_items:
        shown += f"; +{len(items) - max_items}"
    if len(shown) > max_chars:
        shown = shown[: max_chars - 1].rstrip() + "..."
    return shown


def format_strategy_list(value: str) -> str:
    strategies = [strategy_label(item.strip()) for item in str(value or "").split(",") if item.strip()]
    return ", ".join(strategies)


def pathogenic_variant_table(variants: pd.DataFrame) -> pd.DataFrame:
    pathogenic = variants[variants["clinvar_category"].astype(str) == "P/LP"].copy()
    if pathogenic.empty:
        return pd.DataFrame()

    pathogenic["Stars"] = pathogenic.apply(review_star_category, axis=1)
    pathogenic["Strategies"] = pathogenic["strategies"].map(format_strategy_list)
    pathogenic["Disease"] = pathogenic["clinvar_disease"].map(compact_list_text)
    pathogenic["HGVS"] = pathogenic["clinvar_hgvs"].map(lambda value: compact_list_text(value, max_items=1, max_chars=70))
    pathogenic["gnomAD AF"] = pathogenic["gnomad_af"]
    table = pd.DataFrame(
        {
            "Key": pathogenic["variant_id"],
            "Gene": pathogenic["gene_id"],
            "Event": pathogenic["event_type"],
            "ClinVar sig": pathogenic["clinvar_sig"],
            "Stars": pathogenic["Stars"],
            "Review status": pathogenic["clinvar_revstat"],
            "SCVs": pathogenic["clinvar_scv_count"],
            "ClinVar ID": pathogenic["clinvar_id"],
            "Allele ID": pathogenic["clinvar_allele_id"],
            "Disease": pathogenic["Disease"],
            "HGVS": pathogenic["HGVS"],
            "ClinVar type": pathogenic["clinvar_variant_type"],
            "gnomAD AF": pathogenic["gnomAD AF"],
            "gnomAD consequence": pathogenic["gnomad_csq"],
            "Orthologs": pathogenic["support_ortholog_count"],
            "Support events": pathogenic["support_row_count"],
            "Strategies": pathogenic["Strategies"],
        }
    )
    star_rank = table["Stars"].map({star: index for index, star in enumerate(REVIEW_STAR_ORDER[::-1])}).fillna(-1)
    table["_star_rank"] = star_rank
    table = table.sort_values(["_star_rank", "Orthologs", "Support events"], ascending=False).drop(columns=["_star_rank"])
    return table


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
    variant_summary: VariantSummary,
    cov: pd.DataFrame,
    strategy_stats: pd.DataFrame,
    annotation_manifest: dict,
    alignment_manifest: dict,
) -> list[str]:
    unique_variant_count = variant_summary.unique_variant_count
    event_row_count = annotation_manifest.get("event_row_count") or alignment_manifest.get("raw_alignment_event_count") or ""
    clinvar_found = variant_summary.clinvar_found
    clinvar_classified = variant_summary.clinvar_classified
    gnomad_found = variant_summary.gnomad_found
    annotation_warnings = int(annotation_manifest.get("failure_count", 0) or 0)
    cards = [
        ("Raw support events", format_int(event_row_count) if event_row_count != "" else "n/a"),
        ("Unique candidate variants", format_int(unique_variant_count)),
        ("Strategies", format_int(len(variant_summary.strategies))),
        ("Genes", format_int(variant_summary.gene_count)),
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
    variant_summary: VariantSummary,
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

    unique_contrib = variant_summary.unique_contribution
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

    fig_overlap = strategy_overlap_figure(variant_summary.overlap)
    if fig_overlap is not None:
        sections.append("<h3>Strategy Overlap</h3>")
        sections.append(fig_html(fig_overlap))

    counts = variant_summary.event_counts.copy()
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
    variant_summary: VariantSummary,
    strategy_stats: pd.DataFrame,
    include_plotly: bool,
) -> list[str]:
    sections = ["<h2>External Evidence</h2>"]
    sections.append(
        metric_cards(
            [
                ("Found in ClinVar", format_int(variant_summary.clinvar_found)),
                ("ClinVar with CLNSIG", format_int(variant_summary.clinvar_classified)),
                ("Found in gnomAD", format_int(variant_summary.gnomad_found)),
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

    clin_counts = variant_summary.clinvar_counts.copy()
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

    gnomad_bins = variant_summary.gnomad_bins
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

    star_counts = variant_summary.pathogenic_star_counts.copy()
    if not star_counts.empty:
        present_stars = [star for star in REVIEW_STAR_ORDER if star in set(star_counts["Review stars"].astype(str))]
        totals = star_counts.groupby("Strategy", observed=True)["Variant_Count"].sum()
        high_conf = star_counts[star_counts["Review stars"].astype(str).isin(["4", "3", "2"])]
        high_conf_totals = high_conf.groupby("Strategy", observed=True)["Variant_Count"].sum()
        star_order = (
            pd.DataFrame({"total": totals, "high_conf": high_conf_totals})
            .fillna(0)
            .sort_values(["high_conf", "total"], ascending=False)
            .index.tolist()
        )
        fig_stars = px.bar(
            star_counts,
            x="Strategy",
            y="Variant_Count",
            color="Review stars",
            barmode="stack",
            title="Pathogenic ClinVar hits by review stars",
            category_orders={"Strategy": star_order, "Review stars": present_stars},
            color_discrete_map=REVIEW_STAR_COLORS,
            labels={"Strategy": "", "Variant_Count": "P/LP ClinVar variants", "Review stars": "Review stars"},
        )
        compact_figure(fig_stars, height=340)
        sections.append("<h3>Pathogenic ClinVar Evidence</h3>")
        sections.append(fig_html(fig_stars))
    else:
        sections.append("<h3>Pathogenic ClinVar Evidence</h3>")
        sections.append("<p>No P/LP ClinVar variants were found in the candidate set.</p>")

    consequence_counts = group_consequence_counts(variant_summary.consequence_counts)
    if not consequence_counts.empty:
        order = consequence_strategy_order(consequence_counts)
        fig_conseq = px.bar(
            consequence_counts,
            x="Strategy",
            y="Fraction",
            color="Consequence group",
            barmode="stack",
            title="gnomAD consequence mix among gnomAD hits",
            category_orders={"Strategy": order, "Consequence group": CONSEQUENCE_GROUP_ORDER},
            color_discrete_map=CONSEQUENCE_GROUP_COLORS,
            labels={"Strategy": "", "Fraction": "Within-strategy fraction", "Consequence group": "Consequence group"},
        )
        fig_conseq.update_layout(yaxis_tickformat=".0%")
        compact_figure(fig_conseq, height=360)
        sections.append("<h3>gnomAD Consequence Profile</h3>")
        sections.append(fig_html(fig_conseq))
    else:
        sections.append("<p>No gnomAD consequences were found.</p>")

    pathogenic_consequence_counts = group_consequence_counts(variant_summary.pathogenic_consequence_counts)
    if not pathogenic_consequence_counts.empty:
        pathogenic_order = (
            pathogenic_consequence_counts.groupby("Strategy", observed=True)["Variant_Count"]
            .sum()
            .sort_values(ascending=False)
            .index.tolist()
        )
        fig_path_conseq = px.bar(
            pathogenic_consequence_counts,
            x="Strategy",
            y="Variant_Count",
            color="Consequence group",
            barmode="stack",
            title="gnomAD consequence groups for pathogenic ClinVar hits",
            category_orders={"Strategy": pathogenic_order, "Consequence group": CONSEQUENCE_GROUP_ORDER},
            color_discrete_map=CONSEQUENCE_GROUP_COLORS,
            labels={"Strategy": "", "Variant_Count": "P/LP ClinVar variants", "Consequence group": "Consequence group"},
        )
        compact_figure(fig_path_conseq, height=320)
        sections.append(fig_html(fig_path_conseq))

    pathogenic_table = pathogenic_variant_table(variant_summary.pathogenic_rows)
    if not pathogenic_table.empty:
        sections.append("<h3>Pathogenic ClinVar Variants Found</h3>")
        sections.append(table_html(pathogenic_table, classes="table table-sm table-striped", max_rows=100))

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


def validation_excluded_count(manifest: dict, variant_kind: str) -> int:
    return (
        int(manifest.get(f"excluded_vus_{variant_kind}_count", 0))
        + int(manifest.get(f"excluded_missing_{variant_kind}_count", 0))
        + int(manifest.get(f"excluded_other_{variant_kind}_count", 0))
        + int(manifest.get(f"excluded_normalization_{variant_kind}_count", 0))
        + int(manifest.get(f"ambiguous_mixed_label_{variant_kind}_count", 0))
    )


def validation_kind_label(variant_kind: str) -> str:
    return "INDEL" if variant_kind == "indel" else "SNV"


def build_validation_sections(validation, include_plotly: bool) -> list[str]:
    manifest = validation.manifest
    results = validation.strategy_results.copy()
    sections = ["<h2>ClinVar Enrichment</h2>"]
    sections.append(
        metric_cards(
            [
                ("ClinVar SNV universe", format_int(manifest.get("usable_snv_allele_count", 0))),
                ("B/LB SNVs", format_int(manifest.get("benign_snv_count", 0))),
                ("P/LP SNVs", format_int(manifest.get("pathogenic_snv_count", 0))),
                ("Excluded SNV alleles", format_int(validation_excluded_count(manifest, "snv"))),
                ("ClinVar INDEL universe", format_int(manifest.get("usable_indel_allele_count", 0))),
                ("B/LB INDELs", format_int(manifest.get("benign_indel_count", 0))),
                ("P/LP INDELs", format_int(manifest.get("pathogenic_indel_count", 0))),
                ("Excluded INDEL alleles", format_int(validation_excluded_count(manifest, "indel"))),
            ]
        )
    )
    sections.append(
        "<p class=\"lead\">ClinVar validation asks whether observed alternate alleles are enriched for B/LB over P/LP labels. SNV and INDEL are computed separately.</p>"
    )

    if results.empty:
        sections.append("<p>No usable ClinVar validation rows were found.</p>")
        return sections

    for variant_kind in ["snv", "indel"]:
        sections.extend(build_validation_kind_sections(results, variant_kind, include_plotly))
    return sections


def build_validation_kind_sections(results: pd.DataFrame, variant_kind: str, include_plotly: bool) -> list[str]:
    label = validation_kind_label(variant_kind)
    subset = results[results["variant_type"].astype(str) == variant_kind].copy()
    sections = [f"<h3>{label} Enrichment</h3>"]
    if subset.empty:
        sections.append(f"<p>No usable ClinVar {label} rows were found.</p>")
        return sections

    subset["Strategy"] = subset["strategy"].map(strategy_label)
    plot_df = subset.dropna(subset=["ci_low", "ci_high"]).copy()
    plot_df = plot_df[(plot_df["ci_low"] > 0) & (plot_df["ci_high"] > 0)]
    plot_df["plot_odds_ratio"] = plot_df["odds_ratio"]
    infinite_or = ~np.isfinite(plot_df["plot_odds_ratio"])
    plot_df.loc[infinite_or, "plot_odds_ratio"] = np.sqrt(
        plot_df.loc[infinite_or, "ci_low"] * plot_df.loc[infinite_or, "ci_high"]
    )
    plot_df = plot_df[np.isfinite(plot_df["plot_odds_ratio"]) & (plot_df["plot_odds_ratio"] > 0)]
    if not plot_df.empty:
        plot_df = plot_df.sort_values("plot_odds_ratio", ascending=False)
        fig = go.Figure(
            data=go.Scatter(
                x=plot_df["plot_odds_ratio"],
                y=plot_df["Strategy"],
                mode="markers",
                marker={"size": 10, "color": "#356d8f" if variant_kind == "snv" else "#6f4aa8"},
                error_x={
                    "type": "data",
                    "symmetric": False,
                    "array": plot_df["ci_high"] - plot_df["plot_odds_ratio"],
                    "arrayminus": plot_df["plot_odds_ratio"] - plot_df["ci_low"],
                    "thickness": 1.4,
                },
                hovertemplate=(
                    "%{y}<br>OR: %{x:.3g}<br>"
                    "Raw OR: %{customdata[0]}<br>"
                    "95% CI: %{customdata[1]:.3g}-%{customdata[2]:.3g}<br>"
                    "Fisher p: %{customdata[3]:.3g}<br>"
                    "BH q: %{customdata[4]:.3g}<extra></extra>"
                ),
                customdata=np.stack(
                    [
                        plot_df["odds_ratio"].map(format_ratio),
                        plot_df["ci_low"],
                        plot_df["ci_high"],
                        plot_df["fisher_p"],
                        plot_df["fisher_q"],
                    ],
                    axis=-1,
                ),
            )
        )
        fig.add_vline(x=1.0, line_dash="dash", line_color="#8c8c8c")
        fig.update_layout(
            title=f"ClinVar B/LB enrichment among observed {label} alternate alleles",
            xaxis_title="Odds ratio (log scale)",
            yaxis_title="",
            xaxis_type="log",
            height=360,
            margin={"l": 140, "r": 30, "t": 52, "b": 58},
            template="plotly_white",
        )
        fig.update_yaxes(categoryorder="array", categoryarray=plot_df["Strategy"].tolist()[::-1])
        if infinite_or.any():
            fig.add_annotation(
                text="Infinite raw ORs are plotted at the Haldane-corrected CI center.",
                xref="paper",
                yref="paper",
                x=0,
                y=-0.22,
                showarrow=False,
                font={"size": 12, "color": "#52606d"},
                align="left",
            )
        sections.append(fig_html(fig, include_plotlyjs=include_plotly))
    else:
        sections.append(f"<p>{label} odds ratios were not finite enough to draw a log-scale forest plot.</p>")

    table = subset.sort_values("odds_ratio", ascending=False, na_position="last").copy()
    table["Odds Ratio"] = table["odds_ratio"].map(format_ratio)
    table["95% CI"] = table.apply(lambda row: f"{format_ratio(row['ci_low'])}-{format_ratio(row['ci_high'])}", axis=1)
    table["Fisher p"] = table["fisher_p"].map(format_pvalue)
    table["BH q"] = table["fisher_q"].map(format_pvalue)
    table = table.rename(
        columns={
            "benign_observed": "B/LB observed",
            "pathogenic_observed": "P/LP observed",
            "benign_not_observed": "B/LB not observed",
            "pathogenic_not_observed": "P/LP not observed",
        }
    )
    sections.append(f"<h4>{label} 2x2 Tables by Strategy</h4>")
    sections.append(
        table_html(
            table[
                [
                    "Strategy",
                    "B/LB observed",
                    "P/LP observed",
                    "B/LB not observed",
                    "P/LP not observed",
                    "Odds Ratio",
                    "95% CI",
                    "Fisher p",
                    "BH q",
                ]
            ],
            classes="table table-sm table-striped",
        )
    )
    return sections


def build_conservation_analysis(
    *,
    inputs: RunInputs,
    validation,
    strategies: list[str],
) -> ConservationAnalysis:
    conservation = build_conservation_annotations(
        universe=validation.universe,
        universe_path=validation.universe_path,
        analytics_dir=inputs.run_dir / "analytics",
        track_names=DEFAULT_TRACK_NAMES,
    )
    cohort = build_conservation_cohort(
        universe=validation.universe,
        conservation=conservation.annotations,
    )
    results = compute_conservation_validation(
        cohort=cohort,
        observed_by_strategy_type=validation.observed_by_strategy_type,
        strategies=strategies,
        analytics_dir=inputs.run_dir / "analytics",
    )
    return ConservationAnalysis(
        annotations_path=conservation.annotations_path,
        manifest_path=conservation.manifest_path,
        manifest=conservation.manifest,
        validation=results,
    )


def build_fixed_conservation_sections(
    analysis: ConservationAnalysis,
) -> list[str]:
    summary = analysis.validation.cohort.summary
    sections = ["<h2>ClinVar Enrichment Within Fixed phyloP Bands</h2>"]
    sections.append(
        "<p class=\"lead\">This sensitivity analysis repeats the B/LB-versus-P/LP enrichment test within "
        "three prespecified phyloP100way bands. The Mantel-Haenszel estimate pools the band-specific 2x2 tables; "
        "it reduces, but cannot eliminate, residual conservation differences inside each band.</p>"
    )
    sections.append(conservation_cohort_cards(summary))
    sections.append(
        conservation_selector_view(
            view_id="fixed-conservation",
            strategies=analysis.validation.fixed_adjusted["strategy"].drop_duplicates().tolist(),
            primary=analysis.validation.fixed_adjusted,
            detail=analysis.validation.fixed_bins,
            mode="fixed",
        )
    )
    return sections


def build_continuous_firth_sections(
    analysis: ConservationAnalysis,
) -> list[str]:
    summary = analysis.validation.cohort.summary
    sections = ["<h2>ClinVar Enrichment with Continuous phyloP Adjustment</h2>"]
    sections.append(
        "<p class=\"lead\">The primary model uses Firth logistic regression with a three-degree-of-freedom "
        "natural spline for phyloP100way. The adjusted OR asks whether exact ALT observation remains associated "
        "with B/LB classification after modeling the continuous, potentially nonlinear conservation relationship.</p>"
    )
    sections.append(conservation_cohort_cards(summary))
    sections.append(
        conservation_selector_view(
            view_id="continuous-conservation",
            strategies=analysis.validation.continuous["strategy"].drop_duplicates().tolist(),
            primary=analysis.validation.continuous,
            detail=analysis.validation.distributions,
            mode="continuous",
        )
    )
    return sections


def conservation_cohort_cards(summary: dict[str, int]) -> str:
    return metric_cards(
        [
            ("ClinVar B/LB + P/LP alleles", format_int(summary.get("allele_count", 0))),
            ("With phyloP100way", format_int(summary.get("scored_allele_count", 0))),
            ("SNVs", format_int(summary.get("snv_count", 0))),
            ("Insertions", format_int(summary.get("insertion_count", 0))),
            ("Deletions", format_int(summary.get("deletion_count", 0))),
        ]
    )


def conservation_selector_view(
    *,
    view_id: str,
    strategies: list[str],
    primary: pd.DataFrame,
    detail: pd.DataFrame,
    mode: str,
) -> str:
    indel_note = (
        "phyloP thresholds have a nominal single-base p-value interpretation only for SNVs. This view applies "
        "them to the prespecified INDEL aggregate score for descriptive comparability."
        if mode == "fixed"
        else "INDEL phyloP is an aggregate allele score (deleted-base mean or insertion-flank mean), not a "
        "single-base phyloP value. Interpret the adjusted INDEL association separately from the primary SNV view."
    )
    payload = {
        "viewId": view_id,
        "mode": mode,
        "strategies": [{"key": value, "label": strategy_label(value)} for value in strategies],
        "variantTypes": [{"key": key, "label": label} for key, label in VARIANT_TYPE_OPTIONS],
        "consequences": [{"key": key, "label": label} for key, label in CONSEQUENCE_OPTIONS],
        "primary": dataframe_records(primary),
        "detail": dataframe_records(detail),
    }
    payload_json = json.dumps(payload, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    return f"""
    <div class="analysis-controls" id="{view_id}-controls">
      <label>Strategy<select data-role="strategy"></select></label>
      <label>Variant type<select data-role="variant-type"></select></label>
      <label>Consequence subset<select data-role="consequence"></select></label>
    </div>
    <div id="{view_id}-indel-note" class="analysis-note" hidden>
      {indel_note}
    </div>
    <div id="{view_id}-status" class="analysis-note" hidden></div>
    <div id="{view_id}-metrics" class="metric-grid"></div>
    <div id="{view_id}-plot" class="analysis-plot"></div>
    <div id="{view_id}-table"></div>
    <script>
    (() => {{
      const config = {payload_json};
      const root = document.getElementById(config.viewId + '-controls');
      const strategySelect = root.querySelector('[data-role="strategy"]');
      const variantSelect = root.querySelector('[data-role="variant-type"]');
      const consequenceSelect = root.querySelector('[data-role="consequence"]');
      const addOptions = (select, values) => values.forEach(value => {{
        const option = document.createElement('option');
        option.value = value.key; option.textContent = value.label; select.appendChild(option);
      }});
      addOptions(strategySelect, config.strategies);
      addOptions(variantSelect, config.variantTypes);
      addOptions(consequenceSelect, config.consequences);
      variantSelect.value = 'snv';
      consequenceSelect.value = 'missense';

      const finite = value => value !== null && Number.isFinite(Number(value));
      const number = value => finite(value) ? Number(value) : null;
      const fmt = value => {{
        if (value === 'inf') return '∞';
        if (value === '-inf') return '-∞';
        const item = number(value);
        if (item === null) return 'NA';
        if (item === 0) return '0';
        if (Math.abs(item) < 0.001 || Math.abs(item) >= 1000) return item.toExponential(2);
        return item.toPrecision(3);
      }};
      const count = value => number(value) === null ? '0' : Math.round(Number(value)).toLocaleString('en-US').replaceAll(',', ' ');
      const ci = row => finite(row.ci_low) && finite(row.ci_high) ? `${{fmt(row.ci_low)}}–${{fmt(row.ci_high)}}` : 'NA';
      const metric = (label, value) => `<div class="metric-card"><div class="metric-label">${{label}}</div><div class="metric-value">${{value}}</div></div>`;
      const matches = row => row.strategy === strategySelect.value && row.variant_type === variantSelect.value && row.consequence === consequenceSelect.value;
      const cell = value => `<td>${{value}}</td>`;

      function renderFixed(row, detailRows) {{
        document.getElementById(config.viewId + '-metrics').innerHTML = [
          metric('MH adjusted OR', fmt(row?.odds_ratio_mh)),
          metric('95% CI', row ? ci(row) : 'NA'),
          metric('CMH p', fmt(row?.cmh_p)),
          metric('BH q', fmt(row?.cmh_q)),
          metric('Scored alleles', count(row?.usable_rows)),
        ].join('');
        const plotRows = detailRows.filter(item => finite(item.odds_ratio) && number(item.odds_ratio) > 0 && finite(item.ci_low) && finite(item.ci_high));
        Plotly.react(config.viewId + '-plot', [{{
          type: 'scatter', mode: 'markers',
          x: plotRows.map(item => number(item.odds_ratio)),
          y: plotRows.map(item => item.band_label),
          error_x: {{type: 'data', symmetric: false,
            array: plotRows.map(item => number(item.ci_high) - number(item.odds_ratio)),
            arrayminus: plotRows.map(item => number(item.odds_ratio) - number(item.ci_low))}},
          marker: {{size: 10, color: '#2f6f62'}},
          customdata: plotRows.map(item => [ci(item), fmt(item.fisher_p), fmt(item.fisher_q), count(item.row_count)]),
          hovertemplate: '%{{y}}<br>OR: %{{x:.3g}}<br>95% CI: %{{customdata[0]}}<br>Fisher p: %{{customdata[1]}}<br>BH q: %{{customdata[2]}}<br>N: %{{customdata[3]}}<extra></extra>'
        }}], {{title: 'Band-specific enrichment', template: 'plotly_white', height: 330, margin: {{l: 155, r: 25, t: 50, b: 55}}, xaxis: {{title: 'Odds ratio (log scale)', type: 'log'}}, yaxis: {{title: ''}}, shapes: [{{type: 'line', x0: 1, x1: 1, y0: 0, y1: 1, yref: 'paper', line: {{dash: 'dash', color: '#8c8c8c'}}}}]}}, {{responsive: true}});
        const rows = detailRows.map(item => `<tr>${{cell(item.band_label)}}${{cell(item.band_range)}}${{cell(count(item.row_count))}}${{cell(count(item.benign_observed))}}${{cell(count(item.pathogenic_observed))}}${{cell(count(item.benign_not_observed))}}${{cell(count(item.pathogenic_not_observed))}}${{cell(fmt(item.odds_ratio))}}${{cell(ci(item))}}${{cell(fmt(item.fisher_p))}}${{cell(fmt(item.fisher_q))}}${{cell(item.status === 'estimated' ? 'Estimated' : item.reason)}}</tr>`).join('');
        document.getElementById(config.viewId + '-table').innerHTML = `<h3>Band-specific 2x2 Tables</h3><table><thead><tr><th>Band</th><th>Range</th><th>N</th><th>B/LB observed</th><th>P/LP observed</th><th>B/LB not observed</th><th>P/LP not observed</th><th>OR</th><th>95% CI</th><th>Fisher p</th><th>BH q</th><th>Status</th></tr></thead><tbody>${{rows}}</tbody></table>`;
      }}

      function renderContinuous(row, detailRows) {{
        document.getElementById(config.viewId + '-metrics').innerHTML = [
          metric('Firth adjusted OR', fmt(row?.odds_ratio)),
          metric('95% profile CI', row ? ci(row) : 'NA'),
          metric('PLR p', fmt(row?.plr_p)),
          metric('BH q', fmt(row?.plr_q)),
          metric('Scored alleles', count(row?.usable_rows)),
        ].join('');
        const groups = ['ALT observed', 'ALT not observed'];
        const colors = {{'ALT observed': '#2166ac', 'ALT not observed': '#8c8c8c'}};
        const traces = groups.map(group => {{
          const values = detailRows.filter(item => item.group === group);
          return {{type: 'scatter', mode: 'lines', name: group, x: values.map(item => item.score), y: values.map(item => item.ecdf), line: {{color: colors[group], width: 2}}, hovertemplate: 'phyloP: %{{x:.3g}}<br>Cumulative fraction: %{{y:.1%}}<extra>' + group + '</extra>'}};
        }}).filter(trace => trace.x.length);
        Plotly.react(config.viewId + '-plot', traces, {{title: 'phyloP100way distributions', template: 'plotly_white', height: 350, margin: {{l: 65, r: 25, t: 50, b: 55}}, xaxis: {{title: 'phyloP100way'}}, yaxis: {{title: 'Cumulative fraction', tickformat: '.0%'}}, legend: {{orientation: 'h', y: 1.12}}}}, {{responsive: true}});
        const overlap = row && finite(row.overlap_low) && finite(row.overlap_high) ? `${{fmt(row.overlap_low)}}–${{fmt(row.overlap_high)}}` : 'None';
        const status = row?.status === 'estimated' ? 'Estimated' : (row?.reason || 'Not estimable');
        document.getElementById(config.viewId + '-table').innerHTML = `<h3>Model Data and Estimability</h3><table><thead><tr><th>N</th><th>B/LB</th><th>P/LP</th><th>ALT observed</th><th>ALT not observed</th><th>Score range</th><th>Range overlap</th><th>Status</th></tr></thead><tbody><tr>${{cell(count(row?.usable_rows))}}${{cell(count(row?.benign_rows))}}${{cell(count(row?.pathogenic_rows))}}${{cell(count(row?.observed_rows))}}${{cell(count(row?.not_observed_rows))}}${{cell(row ? `${{fmt(row.score_min)}}–${{fmt(row.score_max)}}` : 'NA')}}${{cell(overlap)}}${{cell(status)}}</tr></tbody></table>`;
      }}

      function render() {{
        const row = config.primary.find(matches);
        const detailRows = config.detail.filter(matches);
        const indelNote = document.getElementById(config.viewId + '-indel-note');
        indelNote.hidden = variantSelect.value === 'snv';
        const status = document.getElementById(config.viewId + '-status');
        const messages = [];
        if (row && row.status !== 'estimated') messages.push(`<strong>Not estimable:</strong> ${{row.reason || 'insufficient data'}}`);
        if (config.mode === 'continuous' && row?.overlap_warning) messages.push('The ALT groups overlap across less than 10% of their combined phyloP range; the adjusted estimate relies on limited common support.');
        status.innerHTML = messages.join('<br>'); status.hidden = messages.length === 0;
        if (config.mode === 'fixed') renderFixed(row, detailRows); else renderContinuous(row, detailRows);
      }}
      [strategySelect, variantSelect, consequenceSelect].forEach(select => select.addEventListener('change', render));
      render();
    }})();
    </script>
    """


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


def format_ci(low, high) -> str:
    if pd.isna(low) or pd.isna(high):
        return ""
    return f"{format_ratio(low)}-{format_ratio(high)}"


def build_target_space_null_sections(
    analysis: TargetSpaceNullAnalysis | None,
    include_plotly: bool,
    *,
    enabled: bool = True,
) -> list[str]:
    sections = ["<h2>Target-Space Null</h2>"]
    sections.append(
        "<p class=\"lead\">Each sampled GAPH SNV is compared with SNVs from the same gene and target context, "
        "with the same genomic REF&gt;ALT substitution and the same primary RefSeq VEP consequence. Controls "
        "observed by that GAPH strategy are excluded. phyloP is the outcome and is not used for matching.</p>"
    )
    if not enabled:
        sections.append(
            "<p>Target-Space Null was disabled for this report run. Enable it with "
            "<code>--target-space-null</code>; this analysis uses Ensembl REST VEP and may take hours.</p>"
        )
        return sections
    if analysis is None:
        sections.append("<p>No target-space-null analysis is available for this run.</p>")
        return sections
    summary = analysis.summary.copy()
    if summary.empty:
        sections.append("<p>No consequence-matched target-space controls could be constructed.</p>")
        return sections

    focal_vep = analysis.manifest.get("focal_vep", {})
    sections.append(
        metric_cards(
            [
                (
                    "Sample cap per strategy",
                    format_int(analysis.manifest.get("inputs", {}).get("sample_size_per_strategy", 0)),
                ),
                ("Sampled GAPH SNVs", format_int(analysis.manifest.get("sampled_focal_count", 0))),
                ("VEP-annotated focal SNVs", format_int(analysis.manifest.get("vep_annotated_focal_count", 0))),
                ("Matched focal SNVs", format_int(analysis.manifest.get("matched_focal_count", 0))),
                ("VEP release", str(focal_vep.get("release", ""))),
                ("Control resamples", format_int(analysis.resamples)),
            ]
        )
    )
    conservation_status = analysis.manifest.get("conservation", {}).get("status", "")
    if conservation_status != "complete":
        error = analysis.manifest.get("conservation", {}).get("error", "")
        sections.append(f"<p>phyloP annotation was incomplete: {error or conservation_status}</p>")

    sections.append(
        "<p>The ECDF is the primary distribution view. At each phyloP value it shows the fraction of SNVs at or "
        "below that value. A GAPH curve shifted upward and left indicates lower conservation than its matched "
        "target-space null.</p>"
    )
    ecdf = analysis.ecdf.copy()
    if not ecdf.empty:
        ecdf["Strategy"] = ecdf["strategy"].map(strategy_label)
        fig_ecdf = px.line(
            ecdf,
            x="phyloP100way",
            y="fraction_leq",
            color="set",
            facet_col="Strategy",
            facet_col_wrap=2,
            title="phyloP distributions: GAPH and consequence-matched target-space null",
            labels={"fraction_leq": "Cumulative fraction", "set": ""},
            color_discrete_map={"GAPH": "#2166ac", "Matched target-space null": "#8c8c8c"},
        )
        fig_ecdf.for_each_annotation(lambda item: item.update(text=item.text.split("=")[-1]))
        fig_ecdf.update_yaxes(tickformat=".0%")
        compact_figure(fig_ecdf, height=max(420, 260 * math.ceil(ecdf["Strategy"].nunique() / 2)))
        sections.append(fig_html(fig_ecdf, include_plotlyjs=include_plotly))

    plot = summary.copy()
    plot["Strategy"] = plot["strategy"].map(strategy_label)
    plot = plot.sort_values("median_difference", ascending=False, kind="mergesort")
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=plot["null_median"],
            y=plot["Strategy"],
            mode="markers",
            marker={"color": "#8c8c8c", "size": 8},
            error_x={
                "type": "data",
                "symmetric": False,
                "array": plot["null_ci_high"] - plot["null_median"],
                "arrayminus": plot["null_median"] - plot["null_ci_low"],
                "color": "#8c8c8c",
            },
            name="Target-space-null median (95% resampling interval)",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=plot["observed_median"],
            y=plot["Strategy"],
            mode="markers",
            marker={"color": "#2166ac", "size": 10, "symbol": "diamond"},
            name="GAPH median",
        )
    )
    fig.update_layout(
        title="Descriptive phyloP median summary",
        xaxis_title="phyloP100way",
        yaxis={"categoryorder": "array", "categoryarray": plot["Strategy"].tolist()[::-1]},
    )
    compact_figure(fig, height=max(360, 52 * len(plot) + 120), show_x_title=True)
    sections.append(fig_html(fig, include_plotlyjs=False if not ecdf.empty else include_plotly))

    table = summary.rename(
        columns={
            "strategy": "Strategy",
            "matched_focals": "Matched SNVs",
            "observed_median": "GAPH median",
            "null_median": "Target-space median",
            "null_ci_low": "Target-space median Q2.5",
            "null_ci_high": "Target-space median Q97.5",
            "median_difference": "Median difference",
        }
    )
    table["Strategy"] = table["Strategy"].map(strategy_label)
    sections.append("<h3>Strategy Summary</h3>")
    sections.append(table_html(table, classes="table table-sm table-striped"))

    matching = pd.DataFrame(analysis.manifest.get("matching_by_consequence", []))
    if not matching.empty:
        matching = matching.rename(
            columns={
                "strategy": "Strategy",
                "primary_consequence": "Primary VEP consequence",
                "eligible_focals": "VEP-annotated focal SNVs",
                "matched_focals": "Matched focal SNVs",
                "match_rate": "Matched %",
            }
        )
        matching["Strategy"] = matching["Strategy"].map(strategy_label)
        sections.append("<details><summary>Matching yield by consequence</summary>")
        sections.append(
            "<p>Low match yield identifies consequence classes for which the target-space comparison is poorly supported.</p>"
        )
        sections.append(table_html(matching, classes="table table-sm table-striped"))
        sections.append("</details>")

    consequence = analysis.consequence_summary.copy()
    if not consequence.empty:
        consequence = consequence.rename(
            columns={
                "strategy": "Strategy",
                "primary_consequence": "Primary VEP consequence",
                "matched_focals": "Matched SNVs",
                "observed_median": "GAPH median",
                "null_median": "Target-space median",
                "median_difference": "Median difference",
            }
        )
        consequence["Strategy"] = consequence["Strategy"].map(strategy_label)
        consequence = consequence[
            [
                "Strategy",
                "Primary VEP consequence",
                "Matched SNVs",
                "GAPH median",
                "Target-space median",
                "Median difference",
            ]
        ]
        sections.append("<details><summary>Results by primary VEP consequence</summary>")
        sections.append(table_html(consequence, classes="table table-sm table-striped"))
        sections.append("</details>")
    return sections


def build_methods_sections(
    inputs: RunInputs,
    out_html: Path,
    variant_summary: VariantSummary,
    cov: pd.DataFrame,
    failures: pd.DataFrame,
    annotation_manifest: dict,
    alignment_manifest: dict,
    validation=None,
    conservation_analysis: ConservationAnalysis | None = None,
    negative_controls: TargetSpaceNullAnalysis | None = None,
    report_timings: list[dict[str, object]] | None = None,
) -> list[str]:
    files = [
        ("Run Dir", inputs.run_dir),
        ("Variant Annotations", inputs.variant_annotations_tsv),
        ("Target Features", inputs.target_features_tsv),
        ("Target Sequences", inputs.target_sequences_dir),
        ("Feature Coverage", inputs.feature_coverage_tsv),
        ("Alignment Segments", inputs.alignment_segments_tsv),
        ("Strategy Summary", inputs.strategy_summary_tsv),
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
    if report_timings:
        sections.append("<details><summary>Report computation timing</summary>")
        sections.append(
            "<p>Durations describe this report invocation; cache status is shown where the stage exposes it directly.</p>"
        )
        sections.append(
            table_html(
                pd.DataFrame(report_timings),
                classes="table table-sm table-striped",
            )
        )
        sections.append("</details>")
    sections.append("<details><summary>Input files and loaded row counts</summary>")
    sections.append(table_html(pd.DataFrame(file_rows), classes="table table-sm table-striped"))
    sections.append(
        table_html(
            pd.DataFrame(
                [
                    {"Metric": "Variant-context rows read", "Value": variant_summary.input_row_count},
                    {"Metric": "Unique candidate variants", "Value": variant_summary.unique_variant_count},
                    {"Metric": "Strategy-supported variant records", "Value": variant_summary.strategy_record_count},
                    {"Metric": "Feature coverage rows loaded", "Value": len(cov)},
                    {"Metric": "Annotation failure rows", "Value": len(failures)},
                    {"Metric": "Alignment event mode", "Value": alignment_manifest.get("alignment_event_mode", "")},
                ]
            ),
            classes="table table-sm table-striped",
        )
    )
    sections.append("</details>")
    sections.append("<details><summary>ClinVar class mapping</summary>")
    sections.append(
        "<p class=\"lead\">The report collapses raw <code>clinvar_sig</code> values into conservative plotting classes.</p>"
    )
    sections.append(table_html(clinvar_class_mapping_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    sections.append("<details><summary>ClinVar review stars</summary>")
    sections.append(
        "<p class=\"lead\">Review-star plots use the normalized star value written during annotation from ClinVar review status.</p>"
    )
    sections.append(table_html(clinvar_review_star_mapping_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    sections.append("<details><summary>gnomAD consequence grouping</summary>")
    sections.append(
        "<p class=\"lead\">The External Evidence consequence plots group raw values from the "
        "<code>gnomad_csq</code> annotation column as follows.</p>"
    )
    sections.append(table_html(consequence_grouping_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    sections.append("<details><summary>ClinVar validation denominator and statistics</summary>")
    sections.append(
        "<p class=\"lead\">The ClinVar Enrichment tab intentionally uses a stricter ClinVar subset than External Evidence hit-rate plots.</p>"
    )
    sections.append(table_html(validation_method_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    sections.append("<details><summary>Conservation-adjusted ClinVar validation method</summary>")
    sections.append(
        "<p class=\"lead\">The two conservation tabs use the same normalized ClinVar allele universe and "
        "differ only in whether phyloP100way is represented by fixed bands or a continuous spline.</p>"
    )
    sections.append(table_html(conservation_validation_method_table(conservation_analysis), classes="table table-sm table-striped"))
    sections.append("<h4>ClinVar MC consequence subsets</h4>")
    sections.append(table_html(validation_consequence_grouping_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    sections.append("<details><summary>Negative-control construction</summary>")
    sections.append(
        "<p class=\"lead\">The target-space null uses normalized SNVs and preserves gene, target context, "
        "substitution, and allele-specific functional consequence.</p>"
    )
    sections.append(table_html(negative_control_method_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    sections.append("<details><summary>Feature coverage formulas</summary>")
    sections.append(
        "<p class=\"lead\">Feature Coverage uses the normalized feature-level table emitted by the alignment stage.</p>"
    )
    sections.append(table_html(feature_coverage_formula_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    if validation is not None:
        validation_files = [
            ("ClinVar universe", validation.universe_path),
            ("ClinVar universe manifest", validation.manifest_path),
        ]
        regions_bed = validation.manifest.get("regions_bed", "")
        if regions_bed:
            validation_files.append(("ClinVar target regions", Path(regions_bed)))
        sections.append("<details><summary>Validation cache files</summary>")
        sections.append(
            table_html(
                pd.DataFrame(
                    [
                        {"Key": label, "Path": str(path), "Exists": path.exists(), "Size": file_size_label(path)}
                        for label, path in validation_files
                    ]
                ),
                classes="table table-sm table-striped",
            )
        )
        sections.append("</details>")
    if conservation_analysis is not None:
        conservation_files = [
            ("Conservation allele annotations", conservation_analysis.annotations_path),
            ("Conservation annotation manifest", conservation_analysis.manifest_path),
        ]
        sections.append("<details><summary>Conservation cache files</summary>")
        sections.append(
            table_html(
                pd.DataFrame(
                    [
                        {"Key": label, "Path": str(path), "Exists": path.exists(), "Size": file_size_label(path)}
                        for label, path in conservation_files
                    ]
                ),
                classes="table table-sm table-striped",
            )
        )
        track_rows = [
            {
                "Track": item.get("track", ""),
                "Status": item.get("status", ""),
                "Annotated positions": item.get("annotated_positions", ""),
                "Unique positions": item.get("unique_positions", ""),
                "Annotated alleles": item.get("annotated_variants", ""),
                "Missing alleles": item.get("missing_variants", ""),
                "Blocks": item.get("block_count", ""),
                "Failed blocks": item.get("failed_block_count", ""),
                "Open seconds": item.get("open_seconds", ""),
                "Read seconds": item.get("read_seconds", ""),
                "Error": item.get("error", ""),
                "URL": item.get("url", ""),
            }
            for item in conservation_analysis.manifest.get("tracks", [])
        ]
        if track_rows:
            sections.append(table_html(pd.DataFrame(track_rows), classes="table table-sm table-striped"))
        sections.append("</details>")
    if negative_controls is not None:
        control_files = [
            ("Negative-control manifest", negative_controls.manifest_path),
            ("Target-space-null rows", negative_controls.matched_path),
            ("Target-space-null phyloP annotations", negative_controls.conservation_path),
            ("VEP consequence cache", negative_controls.vep_cache_path),
        ]
        sections.append("<details><summary>Negative-control cache files</summary>")
        sections.append(
            table_html(
                pd.DataFrame(
                    [
                        {"Key": label, "Path": str(path), "Exists": path.exists(), "Size": file_size_label(path)}
                        for label, path in control_files
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
        <title>GAPH Variant Analytics Report</title>
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
            .analysis-controls {{
                display: grid;
                grid-template-columns: repeat(3, minmax(180px, 260px));
                gap: 12px;
                margin: 12px 0;
            }}
            .analysis-controls label {{ color: #52606d; font-size: 13px; }}
            .analysis-controls select {{
                display: block;
                width: 100%;
                margin-top: 4px;
                padding: 7px 8px;
                border: 1px solid #cbd2d9;
                border-radius: 4px;
                background: white;
                color: #1f2933;
            }}
            .analysis-note {{
                margin: 10px 0;
                padding: 9px 11px;
                border-left: 3px solid #d99b2b;
                background: #fff8e8;
                color: #594a2a;
                font-size: 13px;
            }}
            .analysis-plot {{ min-height: 330px; max-width: 980px; }}
            @media (max-width: 760px) {{
                .analysis-controls {{ grid-template-columns: 1fr; }}
            }}
        </style>
    </head>
    <body>
        <h1>GAPH Variant Analytics Report</h1>
        <p class="lead">Run-level analytics for candidate variant support, strategy overlap, external evidence, and target-feature coverage.</p>
        {render_tabs(sections)}
    </body>
    </html>
    """


def main() -> None:
    args = parse_args()
    if args.target_space_null and args.target_space_null_sample_size < 1:
        raise ValueError("--target-space-null-sample-size must be >= 1")
    if args.target_space_null and args.target_space_null_resamples < 100:
        raise ValueError("--target-space-null-resamples must be >= 100")
    inputs = resolve_run_inputs(args.run_dir)
    out_html = resolve_out_html(args, inputs.run_dir)
    report_started = time.perf_counter()
    timings: list[dict[str, object]] = []

    print(f"Streaming {inputs.variant_annotations_tsv}...")
    with timed_stage("Variant summary", timings) as timing:
        variant_summary = build_variant_summary(
            inputs.variant_annotations_tsv,
            inputs.run_dir / "analytics",
            strategy_label,
        )
        timing["Details"] = "cache hit" if variant_summary.cache_hit else "cache miss"

    with timed_stage("Run summary inputs", timings):
        cov = read_feature_coverage(inputs.feature_coverage_tsv)
        alignment_summary = read_strategy_summary(inputs.strategy_summary_tsv)
        failures = read_failures(inputs.annotation_failures_tsv)
        annotation_manifest = read_json(inputs.annotation_manifest_json)
        alignment_manifest = read_json(inputs.alignment_manifest_json)

    print("Computing strategy metrics...")
    with timed_stage("Strategy metrics", timings):
        strategy_stats_full = merge_alignment_summary(variant_summary.strategy_stats, alignment_summary)
        summary_columns = [
            "Strategy",
            "Unique Variants",
            "Ti/Tv",
            "Found in ClinVar",
            "ClinVar found %",
            "gnomAD Found",
            "gnomAD found %",
            "Aligned support records %",
            "Raw support events",
        ]
        strategy_stats = strategy_stats_full[
            [column for column in summary_columns if column in strategy_stats_full.columns]
        ]
        strategies = variant_summary.strategies

    print("Computing ClinVar enrichment...")
    with timed_stage("ClinVar enrichment", timings):
        validation = build_validation(
            run_dir=inputs.run_dir,
            variant_annotations_tsv=inputs.variant_annotations_tsv,
            genes_tsv=inputs.genes_tsv,
            target_sequences_dir=inputs.target_sequences_dir,
            clinvar_vcf=args.clinvar_vcf.expanduser().resolve(),
            strategies=strategies,
        )

    print("Computing conservation-adjusted ClinVar validation...")
    with timed_stage("Conservation-adjusted validation", timings):
        conservation_analysis = build_conservation_analysis(
            inputs=inputs,
            validation=validation,
            strategies=strategies,
        )

    negative_controls = None
    if args.target_space_null:
        print("Computing consequence-matched target-space null...")
        with timed_stage("Target-space null", timings):
            negative_controls = build_target_space_null(
                run_dir=inputs.run_dir,
                variant_annotations_tsv=inputs.variant_annotations_tsv,
                target_features_tsv=inputs.target_features_tsv,
                genes_tsv=inputs.genes_tsv,
                target_sequences_dir=inputs.target_sequences_dir,
                strategies=strategies,
                sample_size_per_strategy=args.target_space_null_sample_size,
                resamples=args.target_space_null_resamples,
                seed=args.target_space_null_seed,
            )
    else:
        timings.append(
            {
                "Stage": "Target-space null",
                "Status": "disabled",
                "Details": "Enable with --target-space-null",
                "Seconds": 0.0,
            }
        )

    sections = [
        ("overview", "Overview", build_overview(variant_summary, cov, strategy_stats, annotation_manifest, alignment_manifest)),
        ("variants", "Variant Profile", build_variant_sections(variant_summary, strategy_stats, include_plotly=True)),
        (
            "external-evidence",
            "External Evidence",
            build_clinvar_gnomad_sections(variant_summary, strategy_stats_full, include_plotly=True),
        ),
        ("coverage", "Feature Coverage", build_feature_sections(cov, include_plotly=True)),
        (
            "target-space-null",
            "Target-Space Null",
            build_target_space_null_sections(
                negative_controls,
                include_plotly=True,
                enabled=args.target_space_null,
            ),
        ),
        ("clinvar-enrichment", "ClinVar Enrichment", build_validation_sections(validation, include_plotly=True)),
        (
            "conservation-fixed",
            "Conservation: Fixed Bands",
            build_fixed_conservation_sections(conservation_analysis),
        ),
        (
            "conservation-continuous",
            "Conservation: Continuous",
            build_continuous_firth_sections(conservation_analysis),
        ),
        (
            "qc",
            "QC",
            build_methods_sections(
                inputs,
                out_html,
                variant_summary,
                cov,
                failures,
                annotation_manifest,
                alignment_manifest,
                validation,
                conservation_analysis,
                negative_controls,
                timings,
            ),
        ),
    ]

    print(f"Writing report to {out_html}...")
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(render_html(sections))
    print(f"Done in {time.perf_counter() - report_started:.3f} s")


if __name__ == "__main__":
    main()
