#!/usr/bin/env python3
"""Build an HTML report for one completed GAPH run."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import time
import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from analytics.core.clinvar_validation import build_validation
from analytics.core.candidate_conservation import CandidateConservation, build_candidate_conservation
from analytics.core.conservation import DEFAULT_TRACK_NAMES, build_conservation_annotations, universe_rows
from analytics.core.conservation_validation import (
    CONSEQUENCE_OPTIONS,
    CONSEQUENCE_TERMS,
    PHYLOP_BANDS,
    SCORE_COLUMN,
    SPLINE_DF,
    TARGET_CONTEXT_OPTIONS,
    VARIANT_TYPE_OPTIONS,
    ConservationValidation,
    build_conservation_cohort,
    compute_conservation_validation,
)
from analytics.core.negative_controls import TargetSpaceNullAnalysis, build_target_space_null
from analytics.core.variant_summary import (
    StrategyOverlap,
    VariantSummary,
    build_variant_summary,
    read_taxonomic_ortholog_evidence,
)


warnings.filterwarnings("ignore", r"All-NaN (slice|axis) encountered")
warnings.filterwarnings("ignore", r"Mean of empty slice")


FEATURE_ORDER = ["gene", "exon", "cds", "utr", "intron"]
PROFILE_FEATURE_ORDER = ["cds", "utr", "intron"]
TARGET_CONTEXT_ORDER = ["cds", "utr", "other_exon", "intron", "other"]
TARGET_CONTEXT_LABELS = {
    "cds": "CDS",
    "utr": "UTR",
    "other_exon": "Other exon",
    "intron": "Intron",
    "other": "Other",
}
TARGET_CONTEXT_COLORS = {
    "CDS": "#2166ac",
    "UTR": "#67a9cf",
    "Other exon": "#d1e5f0",
    "Intron": "#ef8a62",
    "Other": "#bdbdbd",
}
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
TAXONOMIC_SCOPE_ORDER = [
    "all",
    "eukaryota",
    "metazoa",
    "vertebrata",
    "tetrapoda",
    "amniota",
    "mammalia",
    "primates",
]
TAXONOMIC_SCOPE_LABELS = {
    "all": "All selected",
    "eukaryota": "Eukaryota",
    "metazoa": "Metazoa",
    "vertebrata": "Vertebrata",
    "tetrapoda": "Tetrapoda",
    "amniota": "Amniota",
    "mammalia": "Mammalia",
    "primates": "Primates",
}
EVIDENCE_UNIT_ORDER = ["ortholog", "species", "genus", "family", "order"]
EVIDENCE_UNIT_LABELS = {
    "ortholog": "Ortholog",
    "species": "Species",
    "genus": "Genus",
    "family": "Family",
    "order": "Order",
}

@dataclass(frozen=True)
class RunInputs:
    run_dir: Path
    fetch_manifest_json: Path
    genes_tsv: Path
    target_features_tsv: Path
    target_sequences_dir: Path
    variant_annotations_tsv: Path
    variant_strategy_support_tsv: Path
    ortholog_evidence_summary_tsv: Path
    annotation_manifest_json: Path
    annotation_failures_tsv: Path
    feature_coverage_tsv: Path
    alignment_segments_tsv: Path
    alignment_manifest_json: Path
    strategy_summary_tsv: Path
    taxonomy_summary_tsv: Path


@dataclass(frozen=True)
class ConservationAnalysis:
    annotations_path: Path
    manifest_path: Path
    manifest: dict
    validation: ConservationValidation
    candidate: CandidateConservation


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path, help="Completed GAPH run directory.")
    parser.add_argument(
        "--annotation-dir",
        type=Path,
        help="Annotation directory override. Default: <run-dir>/annotation.",
    )
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
        help=(
            "Build the consequence-matched target-space null. Disabled by default because it uses Ensembl VEP "
            "and the gnomAD GraphQL API."
        ),
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
    parser.add_argument(
        "--gnomad-cache-dir",
        type=Path,
        default=os.environ.get("GAPH_GNOMAD_CACHE_DIR") or None,
        help="Optional shared directory for resumable gnomAD regional responses.",
    )
    parser.add_argument(
        "--vep-backend",
        choices=("rest", "local"),
        default=os.environ.get("GAPH_VEP_BACKEND", "rest"),
        help="VEP backend for unified consequences and the target-space null. Default: rest.",
    )
    parser.add_argument(
        "--vep-release",
        default=os.environ.get("GAPH_VEP_RELEASE") or None,
        help="Pinned Ensembl VEP release. Required for local VEP; REST detects the current release.",
    )
    parser.add_argument(
        "--vep-executable",
        default=os.environ.get("GAPH_VEP_EXECUTABLE", "vep"),
        help="Local VEP executable or wrapper. Used with --vep-backend local.",
    )
    parser.add_argument(
        "--vep-cache-dir",
        type=Path,
        default=os.environ.get("GAPH_VEP_CACHE_DIR") or None,
        help="Local indexed VEP cache root. Required with --vep-backend local.",
    )
    parser.add_argument(
        "--vep-forks",
        type=int,
        default=int(os.environ.get("GAPH_VEP_FORKS", "4")),
        help="Worker processes for local VEP. Default: 4.",
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


def resolve_run_inputs(run_dir: Path, annotation_dir: Path | None = None) -> RunInputs:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory does not exist: {run_dir}")
    if not run_dir.is_dir():
        raise NotADirectoryError(f"--run-dir is not a directory: {run_dir}")

    annotation_override = annotation_dir is not None
    annotation_dir = (
        annotation_dir.expanduser().resolve()
        if annotation_dir is not None
        else run_dir / "annotation"
    )
    variant_annotations_tsv = annotation_dir / "variant_annotations.tsv.gz"
    if not annotation_override:
        variant_annotations_tsv = resolve_vep_variant_annotations(run_dir, variant_annotations_tsv)
    base_annotation_dir = run_dir / "annotation"
    reuse_base_support = False
    recovery_manifest_path = annotation_dir / "manifest.json"
    if annotation_override and recovery_manifest_path.exists():
        recovery_manifest = json.loads(recovery_manifest_path.read_text())
        recovery_source = (
            recovery_manifest.get("gnomad_completion", {}).get("source_annotation_dir", "")
        )
        if recovery_source:
            reuse_base_support = Path(recovery_source).expanduser().resolve() == base_annotation_dir.resolve()

    def annotation_support_path(filename: str) -> Path:
        candidate = annotation_dir / filename
        if candidate.exists() or not reuse_base_support:
            return candidate
        return base_annotation_dir / filename

    inputs = RunInputs(
        run_dir=run_dir,
        fetch_manifest_json=run_dir / "fetch" / "manifest.json",
        genes_tsv=run_dir / "fetch" / "genes.tsv.gz",
        target_features_tsv=run_dir / "fetch" / "target_features.tsv.gz",
        target_sequences_dir=run_dir / "fetch" / "sequences" / "targets",
        variant_annotations_tsv=variant_annotations_tsv,
        variant_strategy_support_tsv=annotation_support_path("variant_strategy_support.tsv.gz"),
        ortholog_evidence_summary_tsv=annotation_support_path("ortholog_evidence_summary.tsv.gz"),
        annotation_manifest_json=annotation_dir / "manifest.json",
        annotation_failures_tsv=annotation_dir / "failures.tsv.gz",
        feature_coverage_tsv=run_dir / "alignment" / "feature_coverage.tsv.gz",
        alignment_segments_tsv=run_dir / "alignment" / "alignment_segments.tsv.gz",
        alignment_manifest_json=run_dir / "alignment" / "manifest.json",
        strategy_summary_tsv=run_dir / "alignment" / "strategy_summary.tsv.gz",
        taxonomy_summary_tsv=run_dir / "alignment" / "taxonomy_summary.tsv.gz",
    )
    if not inputs.variant_annotations_tsv.exists():
        raise FileNotFoundError(
            f"Missing variant_annotations.tsv.gz under annotation directory: {annotation_dir}. "
            "Run the annotation stage before building this report."
        )
    if not inputs.genes_tsv.exists():
        raise FileNotFoundError("Missing fetch/genes.tsv.gz under --run-dir.")
    if not inputs.target_features_tsv.exists():
        raise FileNotFoundError("Missing fetch/target_features.tsv.gz under --run-dir.")
    if not inputs.target_sequences_dir.exists():
        raise FileNotFoundError("Missing fetch/sequences/targets under --run-dir.")
    return inputs


def resolve_vep_variant_annotations(run_dir: Path, source: Path) -> Path:
    """Use the completed bulk-VEP artifact only when it matches the source table."""

    artifact_dir = run_dir / "analytics" / "vep_consequences"
    manifest_path = artifact_dir / "manifest.json"
    output_path = artifact_dir / "variant_annotations.vep.tsv.gz"
    if not manifest_path.exists() and not output_path.exists():
        return source
    if not manifest_path.exists() or not output_path.exists():
        raise ValueError(f"Incomplete bulk VEP artifact under {artifact_dir}")

    manifest = read_json(manifest_path)
    source_stat = source.stat()
    source_identity = {
        "path": str(source.resolve()),
        "size_bytes": source_stat.st_size,
        "mtime": int(source_stat.st_mtime),
    }
    if manifest.get("status") != "complete" or manifest.get("source") != source_identity:
        raise ValueError(f"Bulk VEP artifact does not match {source}")
    output_stat = output_path.stat()
    output_identity = {
        "size_bytes": output_stat.st_size,
        "mtime_ns": output_stat.st_mtime_ns,
    }
    if manifest.get("output") != output_identity:
        raise ValueError(f"Bulk VEP output metadata changed: {output_path}")
    return output_path


def bulk_vep_release(inputs: RunInputs) -> str | None:
    expected = inputs.run_dir / "analytics" / "vep_consequences" / "variant_annotations.vep.tsv.gz"
    if inputs.variant_annotations_tsv != expected:
        return None
    manifest = read_json(expected.parent / "manifest.json")
    release = manifest.get("config", {}).get("release")
    return str(release) if release else None


def validate_report_inputs(inputs: RunInputs) -> None:
    """Fail before expensive work when a production table contract is incompatible."""
    contracts = {
        inputs.variant_annotations_tsv: {
            "variant_key",
            "gene_id",
            "event_type",
            "ref",
            "alt",
            "lookup_status",
            "strategies",
            "support_row_count",
            "support_ortholog_count",
            "clinvar_id",
            "clinvar_sig",
            "clinvar_review_stars",
            "clinvar_scv_count",
            "gnomad_af",
            "gnomad_csq",
        },
        inputs.genes_tsv: {"gene_id", "chromosome", "begin", "end", "sequence_length"},
        inputs.target_features_tsv: {
            "gene_id",
            "feature_type",
            "target_start0",
            "target_end0",
        },
        inputs.feature_coverage_tsv: {"gene_id", "strategy", "feature_type"},
        inputs.strategy_summary_tsv: {
            "strategy",
            "gene_count",
            "summary_row_count",
            "aligned_summary_row_count",
            "event_count",
        },
        inputs.variant_strategy_support_tsv: {
            "variant_key",
            "gene_id",
            "strategy",
            "alt_support_row_count",
            "alt_support_ortholog_count",
        },
    }
    for path, required in contracts.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing report input: {path}")
        compression = "gzip" if path.suffix == ".gz" else None
        header = set(pd.read_csv(path, sep="\t", compression=compression, nrows=0).columns)
        missing = required - header
        if missing:
            raise ValueError(
                f"Report input {path} is missing required columns: {', '.join(sorted(missing))}"
            )
    optional_contracts = {
        inputs.ortholog_evidence_summary_tsv: {
            "strategy",
            "target_context",
            "taxonomic_scope",
            "evidence_unit",
            "site_aligned_count",
            "alt_support_count",
            "gnomad_found_count",
            "gnomad_not_found_count",
            "gnomad_lookup_failed_count",
        },
        inputs.taxonomy_summary_tsv: {
            "taxonomic_scope",
            "evidence_unit",
            "gene_count",
            "ortholog_count",
            "taxon_count",
            "unit_count",
            "orthologs_per_gene_median",
            "units_per_gene_median",
        },
    }
    for path, required in optional_contracts.items():
        if not path.exists():
            continue
        header = set(pd.read_csv(path, sep="\t", compression="gzip", nrows=0).columns)
        missing = required - header
        if missing:
            raise ValueError(
                f"Report input {path} is missing required columns: {', '.join(sorted(missing))}"
            )


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


def read_input_gene_count(path: Path) -> int:
    genes = pd.read_csv(path, sep="\t", compression="gzip", usecols=["gene_id"])
    return int(genes["gene_id"].astype(str).nunique())


def read_failures(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    failures = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
    return failures


def read_taxonomy_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    summary = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
    numeric_columns = [
        "gene_count",
        "ortholog_count",
        "taxon_count",
        "unit_count",
        "orthologs_per_gene_min",
        "orthologs_per_gene_median",
        "orthologs_per_gene_mean",
        "orthologs_per_gene_max",
        "units_per_gene_min",
        "units_per_gene_median",
        "units_per_gene_mean",
        "units_per_gene_max",
    ]
    for column in numeric_columns:
        if column in summary:
            summary[column] = pd.to_numeric(summary[column], errors="raise")
    return summary


def alignment_summary_for_report(summary: pd.DataFrame) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame()
    report = pd.DataFrame({"Strategy": summary["strategy"].map(strategy_label)})
    if "gene_count" in summary.columns:
        report["Genes with result"] = summary["gene_count"]
    report["Orthologs aligned %"] = (
        summary["aligned_summary_row_count"] / summary["summary_row_count"].replace(0, np.nan)
    )
    report["Orthologs evaluated"] = summary["summary_row_count"]
    report["Orthologs aligned"] = summary["aligned_summary_row_count"]
    report["Raw support events"] = summary["event_count"]
    if "aligned_target_bp" in summary.columns:
        report["Aligned target bp"] = summary["aligned_target_bp"]
    return report


def merge_alignment_summary(strategy_stats: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    report_summary = alignment_summary_for_report(summary)
    if report_summary.empty:
        return strategy_stats
    return strategy_stats.merge(report_summary, on="Strategy", how="left")


def alignment_gene_ids_by_strategy(cov: pd.DataFrame) -> dict[str, set[str]]:
    if cov.empty or not {"strategy", "gene_id"}.issubset(cov.columns):
        return {}
    return {
        str(strategy): set(group["gene_id"].astype(str))
        for strategy, group in cov.groupby("strategy", sort=False)
    }


def strategy_overlap_figure(overlap: StrategyOverlap | None):
    if overlap is None:
        return None
    order = np.arange(len(overlap.strategies))
    if len(order) > 2:
        distance = np.clip(1.0 - overlap.jaccard, 0.0, 1.0)
        np.fill_diagonal(distance, 0.0)
        order = leaves_list(linkage(squareform(distance, checks=False), method="average", optimal_ordering=True))
        left, right = int(order[0]), int(order[-1])
        mean_distance = distance.mean(axis=1)
        if (mean_distance[left], overlap.strategies[left]) < (mean_distance[right], overlap.strategies[right]):
            order = order[::-1]
    labels = [strategy_label(overlap.strategies[index]) for index in order]
    jaccard = overlap.jaccard[np.ix_(order, order)]
    intersections = overlap.intersections[np.ix_(order, order)]
    unions = overlap.unions[np.ix_(order, order)]

    fig = go.Figure(
        data=go.Heatmap(
            z=jaccard,
            x=labels,
            y=labels,
            text=np.vectorize(lambda value: f"{value:.0%}")(jaccard),
            customdata=np.dstack([intersections, unions]),
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


def consequence_grouping_table(source: str) -> pd.DataFrame:
    rows = [
        {
            "Group": group,
            f"{source} consequence values": ", ".join(CONSEQUENCE_GROUP_TERMS.get(group, []))
            if group != "Other"
            else f"Any non-empty {source} consequence not listed above.",
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
                    "and two-sided Fisher exact p-value. Benjamini-Hochberg FDR is computed across strategies "
                    "within each variant-type, target-context, and consequence selection."
                ),
            },
        ]
    )


def validation_consequence_grouping_table(source: str = "ClinVar MC") -> pd.DataFrame:
    rows = []
    for key, label in CONSEQUENCE_OPTIONS:
        if key == "all":
            continue
        terms = CONSEQUENCE_TERMS.get(key)
        rows.append(
            {
                "Consequence subset": label,
                f"{source} terms": ", ".join(sorted(terms))
                if terms
                else f"Missing {source} consequence or any term not assigned to a named subset.",
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
                    "For each strategy, the denominator is restricted to genes with an alignment result for that "
                    "strategy. Within those genes, ALT_observed=0 means that the strategy did not report that exact "
                    "normalized ALT; no per-base callability filter is applied."
                ),
            },
            {
                "Step": "Target contexts",
                "Definition": (
                    "Each normalized ClinVar allele is assigned to one target-locus context using the affected "
                    "target position and the precedence CDS > UTR > exon > intron > other. If the same allele "
                    "overlaps multiple target genes, the highest-priority context is used once, without duplicating "
                    "the allele."
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
                    "For each analysis, variant-type, target-context, and consequence selection, "
                    "Benjamini-Hochberg correction is applied across strategies. Band-specific Fisher tests are "
                    "corrected across strategies within the same band."
                ),
            },
            {
                "Step": "INDEL interpretation",
                "Definition": (
                    "The fixed thresholds have their nominal single-base p-value interpretation only for SNVs. "
                    "INDEL views apply the same bands to an aggregate score for descriptive comparability."
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
                    "Ensembl VEP uses RefSeq transcripts and pick_allele_gene. The target Entrez Gene ID is "
                    "selected; the most severe Sequence Ontology term on the picked transcript is the matching key."
                ),
            },
            {
                "Step": "Outcomes",
                "Definition": (
                    "The same matched sets compare phyloP100way, exact-allele gnomAD overlap and AF, and exact-allele "
                    "ClinVar overlap and class composition. Conservation and external evidence are outcomes, not "
                    "matching variables."
                ),
            },
            {
                "Step": "Resampling",
                "Definition": (
                    "Each iteration resamples matched sets with replacement and selects one available control from "
                    "each selected set. GAPH, matched-control, and paired-difference statistics use the same draws. "
                    "The report shows descriptive 95% paired matched-set bootstrap intervals and no inferential "
                    "p-value. ClinVar class proportions exclude records with missing CLNSIG; failed gnomAD regions "
                    "remain missing rather than absent."
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
                "Notes": "Length-weighted aggregate used in the target-bases-covered plot.",
            },
            {
                "Metric": "Per-feature mean depth",
                "Formula": "mean_depth = depth_bases / length_bp",
                "Notes": "Uses pipeline-provided depth_bases; the report does not recompute depth from raw aligner output.",
            },
            {
                "Metric": "Weighted mean ortholog depth",
                "Formula": "sum(depth_bases) / sum(length_bp)",
                "Notes": "Main ortholog-depth plot metric.",
            },
            {
                "Metric": "Median feature metrics",
                "Formula": "median(coverage_breadth), median(mean_depth), median(orthologs_covered)",
                "Notes": "Computed over feature rows within each strategy and feature type.",
            },
            {
                "Metric": "Main feature classes",
                "Formula": "CDS, UTR, intron",
                "Notes": "Exon and gene aggregates are omitted because they overlap CDS/UTR/intron and would be redundant.",
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
    pathogenic["Ortholog support / strategy"] = pathogenic.apply(
        lambda row: (
            ""
            if pd.isna(row.get("support_ortholog_mean"))
            else (
                f"{float(row['support_ortholog_mean']):.1f} "
                f"({int(row['support_ortholog_min'])}-{int(row['support_ortholog_max'])})"
            )
        ),
        axis=1,
    )
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
            **(
                {"VEP consequence": pathogenic["vep_primary_consequence"]}
                if "vep_primary_consequence" in pathogenic.columns
                else {}
            ),
            "Ortholog support / strategy": pathogenic["Ortholog support / strategy"],
            "Strategies": pathogenic["Strategies"],
        }
    )
    star_rank = table["Stars"].map({star: index for index, star in enumerate(REVIEW_STAR_ORDER[::-1])}).fillna(-1)
    table["_star_rank"] = star_rank
    table = table.sort_values(
        ["_star_rank", "SCVs", "Key"],
        ascending=[False, False, True],
        kind="mergesort",
    ).drop(columns=["_star_rank"])
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


def ortholog_evidence_figure(
    cells: pd.DataFrame,
    strategy: str,
    quantile_count: int,
    taxonomic_scope: str = "all",
    evidence_unit: str = "ortholog",
):
    contexts = [("cds", "CDS"), ("utr", "UTR"), ("intron", "Intron")]
    figure = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=[label for _context, label in contexts],
        horizontal_spacing=0.08,
    )
    selected = cells[
        cells["strategy"].astype(str).eq(strategy)
        & cells["quantile_count"].astype(int).eq(quantile_count)
        & cells["taxonomic_scope"].astype(str).eq(taxonomic_scope)
        & cells["evidence_unit"].astype(str).eq(evidence_unit)
    ]
    for column, (context, _label) in enumerate(contexts, start=1):
        subset = selected[selected["target_context"].astype(str).eq(context)]
        depth_labels = {int(row.depth_bin): str(row.depth_label) for row in subset.itertuples()}
        alt_labels = {
            int(row.alt_bin): str(row.alt_label)
            for row in subset.itertuples()
        }
        x = [depth_labels.get(index, f"Q{index + 1} (empty)") for index in range(quantile_count)]
        y = [
            alt_labels.get(index, f"Q{index + 1} (empty)")
            for index in range(quantile_count)
        ]
        values = {
            (int(row.depth_bin), int(row.alt_bin)): row
            for row in subset.itertuples()
        }
        z = []
        customdata = []
        for y_index in range(quantile_count):
            z_row = []
            custom_row = []
            for x_index in range(quantile_count):
                row = values.get((x_index, y_index))
                z_row.append(None if row is None else float(row.gnomad_found_fraction))
                custom_row.append(
                    [0, 0]
                    if row is None
                    else [int(row.gnomad_found_count), int(row.gnomad_eligible_count)]
                )
            z.append(z_row)
            customdata.append(custom_row)
        figure.add_trace(
            go.Heatmap(
                x=x,
                y=y,
                z=z,
                customdata=customdata,
                zmin=0,
                zmax=1,
                colorscale="Viridis",
                colorbar={"title": "gnomAD found", "tickformat": ".0%"},
                showscale=column == len(contexts),
                hoverongaps=False,
                hovertemplate=(
                    "Site-aligned orthologs: %{x}<br>"
                    "Exact-ALT support: %{y}<br>"
                    "gnomAD found: %{customdata[0]:,} / %{customdata[1]:,} "
                    "(%{z:.1%})<extra></extra>"
                ),
            ),
            row=1,
            col=column,
        )
        figure.update_xaxes(title_text="Site-aligned evidence units", row=1, col=column)
        if column == 1:
            figure.update_yaxes(title_text="Exact-ALT evidence units", row=1, col=column)
    figure.update_layout(
        height=500,
        margin={"l": 70, "r": 90, "t": 55, "b": 100},
        template="plotly_white",
    )
    figure.update_xaxes(automargin=True)
    figure.update_yaxes(automargin=True)
    return figure


def build_ortholog_evidence_sections(
    variant_summary: VariantSummary,
    include_plotly: bool,
    taxonomy_summary: pd.DataFrame | None = None,
) -> list[str]:
    sections = [
        "<h2>Ortholog Evidence</h2>",
        "<p class=\"lead\">SNV evidence strength by the number of selected evidence units aligned "
        "at the variant site and carrying the exact ALT allele. Cell color is the "
        "gnomAD found fraction; failed lookups are excluded.</p>",
    ]
    if not variant_summary.ortholog_evidence_available:
        sections.append(
            "<p class=\"analysis-note\">Unavailable for this run: "
            "variant_strategy_support.tsv.gz predates site-aligned ortholog depth.</p>"
        )
        return sections
    cells = variant_summary.ortholog_evidence_cells
    if cells.empty:
        sections.append("<p>No eligible SNVs with ortholog evidence and successful gnomAD lookup.</p>")
        return sections

    supported_strategies = [
        strategy
        for strategy in variant_summary.strategies
        if strategy in set(cells["strategy"].astype(str))
    ]
    if not supported_strategies:
        sections.append("<p>No strategies expose taxonomically identified ortholog evidence.</p>")
        return sections
    default_strategy = supported_strategies[0]
    quantile_options = {2: "Median", 4: "Quartiles", 10: "Deciles"}
    available_scopes = set(cells["taxonomic_scope"].astype(str))
    taxonomy_summary = taxonomy_summary if taxonomy_summary is not None else pd.DataFrame()
    visible_scopes = []
    seen_scope_signatures = set()
    for scope in TAXONOMIC_SCOPE_ORDER:
        if scope not in available_scopes:
            continue
        signature = None
        if not taxonomy_summary.empty:
            row = taxonomy_summary[
                taxonomy_summary["taxonomic_scope"].astype(str).eq(scope)
                & taxonomy_summary["evidence_unit"].astype(str).eq("ortholog")
            ]
            if not row.empty:
                signature = (
                    int(row.iloc[0]["ortholog_count"]),
                    int(row.iloc[0]["taxon_count"]),
                    float(row.iloc[0]["orthologs_per_gene_median"]),
                )
        if signature is not None and signature in seen_scope_signatures:
            continue
        visible_scopes.append(scope)
        if signature is not None:
            seen_scope_signatures.add(signature)
    visible_scopes.extend(sorted(available_scopes - set(visible_scopes) - set(TAXONOMIC_SCOPE_ORDER)))
    available_units = [
        unit
        for unit in EVIDENCE_UNIT_ORDER
        if unit in set(cells["evidence_unit"].astype(str))
    ]
    default_scope = "all" if "all" in visible_scopes else visible_scopes[0]
    default_unit = "ortholog" if "ortholog" in available_units else available_units[0]

    figures = {}
    for strategy in supported_strategies:
        figures[strategy] = {}
        for scope in visible_scopes:
            scoped = cells[
                cells["strategy"].astype(str).eq(strategy)
                & cells["taxonomic_scope"].astype(str).eq(scope)
            ]
            if scoped.empty:
                continue
            figures[strategy][scope] = {}
            for unit in available_units:
                if scoped["evidence_unit"].astype(str).eq(unit).sum() == 0:
                    continue
                figures[strategy][scope][unit] = {}
                for quantile_count in quantile_options:
                    figure = ortholog_evidence_figure(
                        cells,
                        strategy,
                        quantile_count,
                        scope,
                        unit,
                    )
                    figures[strategy][scope][unit][str(quantile_count)] = json.loads(
                        figure.to_json()
                    )

    initial = ortholog_evidence_figure(
        cells,
        default_strategy,
        4,
        default_scope,
        default_unit,
    )
    initial_html = initial.to_html(
        full_html=False,
        include_plotlyjs="cdn" if include_plotly else False,
        div_id="ortholog-evidence-plot",
    )
    strategy_options = "".join(
        f'<option value="{strategy}">{strategy_label(strategy)}</option>'
        for strategy in supported_strategies
    )
    unsupported_options = "".join(
        f'<option disabled>{strategy_label(strategy)} (taxonomy unavailable)</option>'
        for strategy in variant_summary.strategies
        if strategy not in supported_strategies
    )
    scope_options = "".join(
        f'<option value="{scope}"{" selected" if scope == default_scope else ""}>'
        f'{TAXONOMIC_SCOPE_LABELS.get(scope, scope)}</option>'
        for scope in visible_scopes
    )
    unit_options = "".join(
        f'<option value="{unit}"{" selected" if unit == default_unit else ""}>'
        f'{EVIDENCE_UNIT_LABELS.get(unit, unit)}</option>'
        for unit in available_units
    )
    quantile_html = "".join(
        f'<option value="{count}"{" selected" if count == 4 else ""}>{label}</option>'
        for count, label in quantile_options.items()
    )
    payload = json.dumps(figures, separators=(",", ":"))
    stats = {}
    if not taxonomy_summary.empty:
        for scope in visible_scopes:
            ortholog_row = taxonomy_summary[
                taxonomy_summary["taxonomic_scope"].astype(str).eq(scope)
                & taxonomy_summary["evidence_unit"].astype(str).eq("ortholog")
            ]
            if ortholog_row.empty:
                continue
            stats[scope] = {}
            ortholog_median = float(ortholog_row.iloc[0]["orthologs_per_gene_median"])
            for unit in available_units:
                unit_row = taxonomy_summary[
                    taxonomy_summary["taxonomic_scope"].astype(str).eq(scope)
                    & taxonomy_summary["evidence_unit"].astype(str).eq(unit)
                ]
                if unit_row.empty:
                    continue
                row = unit_row.iloc[0]
                stats[scope][unit] = (
                    f"Median selected orthologs/gene: {ortholog_median:,.1f}; "
                    f"median {EVIDENCE_UNIT_LABELS.get(unit, unit).lower()} units/gene: "
                    f"{float(row['units_per_gene_median']):,.1f}; "
                    f"distinct units in run: {int(row['unit_count']):,}."
                )
    stats_payload = json.dumps(stats, separators=(",", ":"))
    sections.append(
        f"""
        <div class="analysis-controls" id="ortholog-evidence-controls">
            <label>Strategy<select id="ortholog-evidence-strategy">{strategy_options}{unsupported_options}</select></label>
            <label>Taxonomic scope<select id="ortholog-evidence-scope">{scope_options}</select></label>
            <label>Evidence unit<select id="ortholog-evidence-unit">{unit_options}</select></label>
            <label>Groups<select id="ortholog-evidence-quantiles">{quantile_html}</select></label>
        </div>
        <p class="analysis-note" id="ortholog-evidence-stats"></p>
        {initial_html}
        <script>
        (() => {{
            const figures = {payload};
            const stats = {stats_payload};
            const strategy = document.getElementById('ortholog-evidence-strategy');
            const scope = document.getElementById('ortholog-evidence-scope');
            const unit = document.getElementById('ortholog-evidence-unit');
            const quantiles = document.getElementById('ortholog-evidence-quantiles');
            const summary = document.getElementById('ortholog-evidence-stats');
            const firstKey = value => Object.keys(value)[0];
            const render = () => {{
                const strategyFigures = figures[strategy.value];
                if (!strategyFigures[scope.value]) scope.value = firstKey(strategyFigures);
                const scopeFigures = strategyFigures[scope.value];
                if (!scopeFigures[unit.value]) unit.value = firstKey(scopeFigures);
                const figure = scopeFigures[unit.value][quantiles.value];
                summary.textContent = stats[scope.value]?.[unit.value] || '';
                Plotly.react('ortholog-evidence-plot', figure.data, figure.layout, {{responsive: true}});
            }};
            strategy.addEventListener('change', render);
            scope.addEventListener('change', render);
            unit.addEventListener('change', render);
            quantiles.addEventListener('change', render);
            render();
        }})();
        </script>
        """
    )
    return sections


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


def target_gene_coverage_for_report(cov: pd.DataFrame) -> pd.DataFrame:
    required = {"strategy", "feature_type", "length_bp", "covered_bases"}
    if cov.empty or not required.issubset(cov.columns):
        return pd.DataFrame(columns=["Strategy", "Target bases covered %"])
    genes = cov[cov["feature_type"].astype(str).str.lower().eq("gene")]
    if genes.empty:
        return pd.DataFrame(columns=["Strategy", "Target bases covered %"])
    coverage = (
        genes.groupby("strategy", as_index=False)
        .agg(Target_Length_bp=("length_bp", "sum"), Covered_Bases=("covered_bases", "sum"))
    )
    coverage["Target bases covered %"] = (
        coverage["Covered_Bases"] / coverage["Target_Length_bp"].replace(0, np.nan)
    )
    coverage["Strategy"] = coverage["strategy"].map(strategy_label)
    return coverage[["Strategy", "Target bases covered %"]]


def overview_strategy_table(
    variant_summary: VariantSummary,
    cov: pd.DataFrame,
    strategy_stats: pd.DataFrame,
    input_gene_count: int,
) -> pd.DataFrame:
    stats = strategy_stats.merge(variant_summary.unique_contribution, on="Strategy", how="left")
    stats = stats.merge(target_gene_coverage_for_report(cov), on="Strategy", how="left")
    stats["Unique To Strategy"] = stats["Unique To Strategy"].fillna(0)
    stats = stats.sort_values("Unique Variants", ascending=False, kind="mergesort")

    table = pd.DataFrame(
        {
            "Strategy": stats["Strategy"],
            "Genes with result": [
                format_count_ratio(count, input_gene_count)
                for count in stats["Genes with result"]
            ],
            "Candidate variants": stats["Unique Variants"],
            "Only this strategy": [
                format_count_share(count, total)
                for count, total in zip(stats["Unique To Strategy"], stats["Unique Variants"])
            ],
            "gnomAD matches": [
                format_count_share(count, total)
                for count, total in zip(stats["gnomAD Found"], stats["gnomAD Eligible"])
            ],
            "ClinVar matches": [
                format_count_share(count, total)
                for count, total in zip(stats["Found in ClinVar"], stats["Unique Variants"])
            ],
            "Orthologs aligned": [
                format_count_ratio(count, total)
                for count, total in zip(stats["Orthologs aligned"], stats["Orthologs evaluated"])
            ],
            "Target bases covered %": stats["Target bases covered %"],
        }
    )
    return table.reset_index(drop=True)


def build_overview(
    variant_summary: VariantSummary,
    cov: pd.DataFrame,
    strategy_stats: pd.DataFrame,
    annotation_manifest: dict,
    input_gene_count: int,
) -> list[str]:
    unique_variant_count = variant_summary.unique_variant_count
    all_strategy_count = variant_summary.all_strategy_variant_count
    annotation_warnings = int(annotation_manifest.get("failure_count", 0) or 0)
    cards = [
        ("Unique candidate variants", format_int(unique_variant_count)),
        ("Strategies", format_int(len(variant_summary.strategies))),
        ("Input genes", format_int(input_gene_count)),
        ("Candidates found by all strategies", format_count_share(all_strategy_count, unique_variant_count)),
        ("Annotation warnings", format_int(annotation_warnings)),
    ]
    sections = [metric_cards(cards)]
    sections.append("<h2>Strategies</h2>")
    sections.append(
        table_html(
            overview_strategy_table(variant_summary, cov, strategy_stats, input_gene_count),
            classes="table overview-table",
        )
    )
    return sections


def build_variant_sections(
    variant_summary: VariantSummary,
    strategy_stats: pd.DataFrame,
    include_plotly: bool,
) -> list[str]:
    sections = ["<h2>Strategy Concordance</h2>"]
    plotly_pending = include_plotly
    fig_overlap = strategy_overlap_figure(variant_summary.overlap)
    if fig_overlap is not None:
        sections.append(fig_html(fig_overlap, include_plotlyjs=plotly_pending))
        plotly_pending = False

    sections.append("<h2>Variant Composition</h2>")
    counts = variant_summary.event_counts.copy()
    totals = counts.groupby("strategy", observed=True)["Variant_Count"].transform("sum")
    counts["Fraction"] = counts["Variant_Count"] / totals.replace(0, np.nan)
    snv_order = (
        counts[counts["event_type"].astype(str).str.lower() == "snv"]
        .sort_values("Fraction", ascending=False)
        ["strategy"]
        .tolist()
    )
    all_strategies = strategy_stats["Strategy"].tolist()
    order = snv_order + [strategy for strategy in all_strategies if strategy not in snv_order]
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
    sections.append(fig_html(fig_events, include_plotlyjs=plotly_pending))
    plotly_pending = False

    contexts = variant_summary.target_context_counts.copy()
    if not contexts.empty:
        totals = contexts.groupby("strategy", observed=True)["Variant_Count"].transform("sum")
        contexts["Fraction"] = contexts["Variant_Count"] / totals.replace(0, np.nan)
        contexts["Target context"] = contexts["target_context"].map(TARGET_CONTEXT_LABELS).fillna("Other")
        cds_order = (
            contexts[contexts["target_context"].astype(str).eq("cds")]
            .sort_values("Fraction", ascending=False)["strategy"]
            .tolist()
        )
        context_order = cds_order + [strategy for strategy in all_strategies if strategy not in cds_order]
        fig_context = px.bar(
            contexts,
            x="strategy",
            y="Fraction",
            color="Target context",
            barmode="stack",
            title="Target context composition by strategy",
            category_orders={
                "strategy": context_order,
                "Target context": [TARGET_CONTEXT_LABELS[item] for item in TARGET_CONTEXT_ORDER],
            },
            color_discrete_map=TARGET_CONTEXT_COLORS,
            labels={"strategy": "", "Fraction": "Variant fraction"},
            custom_data=["Variant_Count"],
        )
        fig_context.update_layout(yaxis_tickformat=".0%")
        fig_context.update_traces(
            hovertemplate="%{x}<br>%{fullData.name}: %{customdata[0]:,} (%{y:.1%})<extra></extra>"
        )
        compact_figure(fig_context, height=360)
        sections.append(fig_html(fig_context))
    return sections


def build_clinvar_gnomad_sections(
    variant_summary: VariantSummary,
    strategy_stats: pd.DataFrame,
    include_plotly: bool,
) -> list[str]:
    sections = ["<h2>Population Evidence</h2>"]
    gnomad_rate = sort_by_metric(strategy_stats[["Strategy", "gnomAD found %"]], "gnomAD found %")
    gnomad_rate = gnomad_rate.merge(
        strategy_stats[["Strategy", "gnomAD Found", "gnomAD Eligible"]],
        on="Strategy",
        how="left",
    )
    fig_gnomad_rate = px.bar(
        gnomad_rate,
        x="Strategy",
        y="gnomAD found %",
        title="gnomAD hit rate by strategy",
        category_orders={"Strategy": gnomad_rate["Strategy"].tolist()},
        custom_data=["gnomAD Found", "gnomAD Eligible"],
    )
    fig_gnomad_rate.update_layout(yaxis_tickformat=".1%")
    fig_gnomad_rate.update_traces(
        hovertemplate=(
            "%{x}<br>Found in gnomAD: %{customdata[0]:,}<br>"
            "Completed lookups: %{customdata[1]:,}<br>Hit rate: %{y:.2%}<extra></extra>"
        )
    )
    compact_figure(fig_gnomad_rate)
    sections.append(fig_html(fig_gnomad_rate, include_plotlyjs=include_plotly))

    af_summary = variant_summary.gnomad_af_summary.sort_values("Median", ascending=False)
    if not af_summary.empty:
        fig_af = go.Figure()
        for width, low, high, color, name in [
            (3, "Q05", "Q95", "#9ecae1", "5-95% interval"),
            (10, "Q25", "Q75", "#3182bd", "Interquartile interval"),
        ]:
            x_values, y_values = [], []
            for row in af_summary.itertuples(index=False):
                x_values.extend([getattr(row, low), getattr(row, high), None])
                y_values.extend([row.Strategy, row.Strategy, None])
            fig_af.add_trace(go.Scatter(
                x=x_values, y=y_values, mode="lines", line={"width": width, "color": color},
                name=name, hoverinfo="skip",
            ))
        fig_af.add_trace(go.Scatter(
            x=af_summary["Median"],
            y=af_summary["Strategy"],
            mode="markers",
            marker={"size": 9, "color": "#08306b"},
            name="Median",
            customdata=af_summary[["Count", "Q05", "Q25", "Q75", "Q95"]],
            hovertemplate=(
                "%{y}<br>Variants with AF &gt; 0: %{customdata[0]:,}<br>"
                "5th percentile: %{customdata[1]:.3f}<br>Q1: %{customdata[2]:.3f}<br>"
                "Median: %{x:.3f}<br>Q3: %{customdata[3]:.3f}<br>"
                "95th percentile: %{customdata[4]:.3f}<extra></extra>"
            ),
        ))
        fig_af.update_layout(title="gnomAD allele frequency among exact hits", xaxis_title="log10 gnomAD AF")
        compact_figure(fig_af, height=380, show_x_title=True)
        sections.append(fig_html(fig_af))
    else:
        sections.append("<p>No non-zero gnomAD AF values were found.</p>")

    consequence_counts = group_consequence_counts(variant_summary.consequence_counts)
    if not consequence_counts.empty:
        order = consequence_strategy_order(consequence_counts)
        fig_conseq = px.bar(
            consequence_counts,
            x="Strategy",
            y="Fraction",
            color="Consequence group",
            barmode="stack",
            title=f"{variant_summary.consequence_source} consequence mix among candidates",
            category_orders={"Strategy": order, "Consequence group": CONSEQUENCE_GROUP_ORDER},
            color_discrete_map=CONSEQUENCE_GROUP_COLORS,
            labels={"Strategy": "", "Fraction": "Within-strategy fraction", "Consequence group": "Consequence group"},
        )
        fig_conseq.update_layout(yaxis_tickformat=".0%")
        compact_figure(fig_conseq, height=360)
        sections.append(fig_html(fig_conseq))

    sections.append("<h2>Clinical Evidence</h2>")
    clinvar_rate = sort_by_metric(strategy_stats[["Strategy", "ClinVar found %"]], "ClinVar found %")
    clinvar_rate = clinvar_rate.merge(
        strategy_stats[["Strategy", "Found in ClinVar", "Unique Variants"]], on="Strategy", how="left"
    )
    fig_clin_rate = px.bar(
        clinvar_rate,
        x="Strategy",
        y="ClinVar found %",
        title="ClinVar hit rate by strategy",
        category_orders={"Strategy": clinvar_rate["Strategy"].tolist()},
        custom_data=["Found in ClinVar", "Unique Variants"],
    )
    fig_clin_rate.update_layout(yaxis_tickformat=".2%")
    fig_clin_rate.update_traces(
        hovertemplate=(
            "%{x}<br>Found in ClinVar: %{customdata[0]:,}<br>"
            "Candidate variants: %{customdata[1]:,}<br>Hit rate: %{y:.3%}<extra></extra>"
        )
    )
    compact_figure(fig_clin_rate)
    sections.append(fig_html(fig_clin_rate))

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

    star_counts = variant_summary.pathogenic_star_counts.copy()
    if not star_counts.empty:
        present_stars = [star for star in REVIEW_STAR_ORDER if star != "Unmapped"]
        if "Unmapped" in set(star_counts["Review stars"].astype(str)):
            present_stars.append("Unmapped")
        complete_index = pd.MultiIndex.from_product(
            [strategy_stats["Strategy"].tolist(), present_stars], names=["Strategy", "Review stars"]
        )
        star_counts = (
            star_counts.set_index(["Strategy", "Review stars"])
            .reindex(complete_index, fill_value=0)
            .reset_index()
        )
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
            title=f"{variant_summary.consequence_source} consequence groups for pathogenic ClinVar hits",
            category_orders={"Strategy": pathogenic_order, "Consequence group": CONSEQUENCE_GROUP_ORDER},
            color_discrete_map=CONSEQUENCE_GROUP_COLORS,
            labels={"Strategy": "", "Variant_Count": "P/LP ClinVar variants", "Consequence group": "Consequence group"},
        )
        compact_figure(fig_path_conseq, height=320)
        sections.append(fig_html(fig_path_conseq))

    return sections


def build_feature_sections(cov: pd.DataFrame, include_plotly: bool) -> list[str]:
    sections = ["<h2>Alignment Coverage</h2>"]
    if cov.empty:
        sections.append("<p>No feature coverage table was found.</p>")
        return sections

    profile_summary = coverage_summary(cov, PROFILE_FEATURE_ORDER)
    if profile_summary.empty:
        sections.append("<p>No CDS/UTR/intron coverage rows were found.</p>")
        return sections

    cds_depth = (
        profile_summary[profile_summary["feature_type"].astype(str) == "cds"]
        .sort_values("Mean_Depth_Weighted", ascending=False)
    )
    strategy_order = cds_depth["strategy"].tolist() or sorted(profile_summary["strategy"].unique())

    fig_breadth = px.bar(
        profile_summary,
        x="strategy",
        y="Breadth_Weighted",
        color="feature_type",
        barmode="group",
        title="Target bases covered by one or more orthologs",
        category_orders={"strategy": strategy_order, "feature_type": PROFILE_FEATURE_ORDER},
        labels={"strategy": "", "Breadth_Weighted": "Target bases covered", "feature_type": "Feature"},
    )
    fig_breadth.update_layout(yaxis_tickformat=".0%")
    compact_figure(fig_breadth, height=360)
    sections.append(fig_html(fig_breadth, include_plotlyjs=include_plotly))

    fig_depth = px.bar(
        profile_summary,
        x="strategy",
        y="Mean_Depth_Weighted",
        color="feature_type",
        barmode="group",
        title="Weighted mean ortholog depth by target feature",
        category_orders={"strategy": strategy_order, "feature_type": PROFILE_FEATURE_ORDER},
        labels={
            "strategy": "",
            "Mean_Depth_Weighted": "Weighted mean ortholog depth",
            "feature_type": "Feature",
        },
    )
    compact_figure(fig_depth, height=360)
    sections.append(fig_html(fig_depth, include_plotlyjs=include_plotly))
    return sections


def gnomad_stratification_figure(
    counts: pd.DataFrame,
    category_column: str,
    category_order: list[str],
    strategy_order: list[str],
    title: str,
    color_map: dict[str, str] | None = None,
):
    if counts.empty:
        return None
    plot = counts.copy()
    plot["strategy"] = plot["strategy"].astype(str)
    plot["gnomad_status"] = plot["gnomad_status"].astype(str)
    plot[category_column] = plot[category_column].astype(str)
    present_categories = set(plot[category_column])
    categories = [category for category in category_order if category in present_categories]
    categories += sorted(present_categories - set(categories))
    statuses = [("found", "Found"), ("not_found", "Not found")]
    combinations = [
        (strategy, status, label)
        for strategy in strategy_order
        for status, label in statuses
    ]
    indexed = plot.set_index(["strategy", "gnomad_status", category_column])["Variant_Count"]
    totals = plot.groupby(["strategy", "gnomad_status"], observed=True)["Variant_Count"].sum()
    x_strategy = [strategy for strategy, _status, _label in combinations]
    x_status = [label for _strategy, _status, label in combinations]
    fig = go.Figure()
    for category in categories:
        counts_for_category = np.asarray(
            [int(indexed.get((strategy, status, category), 0)) for strategy, status, _label in combinations]
        )
        totals_for_group = np.asarray(
            [int(totals.get((strategy, status), 0)) for strategy, status, _label in combinations]
        )
        fractions = np.divide(
            counts_for_category,
            totals_for_group,
            out=np.zeros(len(combinations), dtype=float),
            where=totals_for_group > 0,
        )
        fig.add_trace(
            go.Bar(
                x=[x_strategy, x_status],
                y=fractions,
                name=category,
                marker_color=(color_map or {}).get(category),
                customdata=np.column_stack(
                    [x_strategy, x_status, counts_for_category, totals_for_group]
                ),
                hovertemplate=(
                    "%{customdata[0]}<br>%{customdata[1]} in gnomAD<br>"
                    + category + ": %{customdata[2]:,} / %{customdata[3]:,} "
                    "(%{y:.1%})<extra></extra>"
                ),
            )
        )
    fig.update_layout(title=title, barmode="stack", height=440, bargap=0.18)
    fig.update_yaxes(tickformat=".0%", title_text="Within-stratum fraction", range=[0, 1])
    compact_figure(fig, height=440)
    return fig


def build_gnomad_stratification_sections(
    variant_summary: VariantSummary,
    strategy_stats: pd.DataFrame,
    candidate_conservation: CandidateConservation,
    include_plotly: bool,
) -> list[str]:
    sections = [
        "<h2>gnomAD Stratification</h2>",
        "<p class=\"lead\">Descriptive comparison of candidate alleles found and not found in gnomAD. "
        "This is not a matched-control analysis.</p>",
    ]
    strategy_order = sort_by_metric(
        strategy_stats[["Strategy", "gnomAD found %"]], "gnomAD found %"
    )["Strategy"].tolist()

    event_counts = variant_summary.gnomad_event_counts.copy()
    event_order = [
        item for item in ["snv", "ins", "del", "mnv", "complex"]
        if item in set(event_counts.get("event_type", pd.Series(dtype=str)).astype(str))
    ]
    event_order += sorted(set(event_counts.get("event_type", pd.Series(dtype=str)).astype(str)) - set(event_order))
    plotly_pending = include_plotly
    event_fig = gnomad_stratification_figure(
        event_counts,
        "event_type",
        event_order,
        strategy_order,
        "Variant type: gnomAD hits versus non-hits",
    )
    if event_fig is not None:
        sections.append(fig_html(event_fig, include_plotlyjs=plotly_pending))
        plotly_pending = False

    context_counts = variant_summary.gnomad_context_counts.copy()
    context_counts["Target context"] = context_counts.get("target_context", pd.Series(dtype=str)).map(
        TARGET_CONTEXT_LABELS
    ).fillna("Other")
    context_fig = gnomad_stratification_figure(
        context_counts,
        "Target context",
        [TARGET_CONTEXT_LABELS[item] for item in TARGET_CONTEXT_ORDER],
        strategy_order,
        "Target context: gnomAD hits versus non-hits",
        TARGET_CONTEXT_COLORS,
    )
    if context_fig is not None:
        sections.append(fig_html(context_fig, include_plotlyjs=plotly_pending))
        plotly_pending = False

    sections.append("<h3>Conservation</h3>")
    phylop_fig = candidate_phylop_figure(candidate_conservation, strategy_order)
    if phylop_fig is not None:
        sections.append(
            "<p class=\"lead\">Candidate phyloP100way distributions are shown separately for exact gnomAD hits "
            "and alleles without a gnomAD hit. Select one strategy to compare the two strata.</p>"
        )
        sections.append(fig_html(phylop_fig, include_plotlyjs=plotly_pending))
        plotly_pending = False
        phylop_summary_fig = candidate_phylop_summary_figure(candidate_conservation, strategy_order)
        if phylop_summary_fig is not None:
            sections.append(fig_html(phylop_summary_fig, include_plotlyjs=False))
    else:
        sections.append("<p>No candidate-wide phyloP100way scores were available.</p>")

    sections.extend(
        [
            "<h3>Functional Consequence</h3>",
            "<p class=\"analysis-note\">Not computed. A defensible comparison requires the same VEP release, "
            "transcript set, and consequence-selection rule for both gnomAD strata.</p>",
        ]
    )
    return sections


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


def build_conservation_analysis(
    *,
    inputs: RunInputs,
    validation,
    strategies: list[str],
    eligible_gene_ids_by_strategy: dict[str, set[str]],
) -> ConservationAnalysis:
    candidate = build_candidate_conservation(
        variant_annotations_tsv=inputs.variant_annotations_tsv,
        analytics_dir=inputs.run_dir / "analytics",
        annotation_failures_tsv=inputs.annotation_failures_tsv,
        additional_rows=universe_rows(validation.universe),
        track_names=DEFAULT_TRACK_NAMES,
    )
    conservation = build_conservation_annotations(
        universe=validation.universe,
        universe_path=validation.universe_path,
        analytics_dir=inputs.run_dir / "analytics",
        track_names=DEFAULT_TRACK_NAMES,
        position_scores=candidate.position_scores,
    )
    cohort = build_conservation_cohort(
        universe=validation.universe,
        conservation=conservation.annotations,
        genes_tsv=inputs.genes_tsv,
        target_features_tsv=inputs.target_features_tsv,
        consequence_column=validation.consequence_column,
    )
    results = compute_conservation_validation(
        cohort=cohort,
        observed_by_strategy_type=validation.observed_by_strategy_type,
        strategies=strategies,
        analytics_dir=inputs.run_dir / "analytics",
        eligible_gene_ids_by_strategy=eligible_gene_ids_by_strategy,
    )
    return ConservationAnalysis(
        annotations_path=conservation.annotations_path,
        manifest_path=conservation.manifest_path,
        manifest=conservation.manifest,
        validation=results,
        candidate=candidate.without_position_scores(),
    )


def build_clinvar_association_sections(analysis: ConservationAnalysis) -> list[str]:
    return [
        "<h2>ClinVar Association</h2>",
        clinvar_association_view(analysis.validation),
    ]


def clinvar_association_view(validation: ConservationValidation) -> str:
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
        "viewId": "clinvar-association",
        "modes": [
            {"key": "unadjusted", "label": "Unadjusted"},
            {"key": "fixed", "label": "phyloP fixed bands"},
            {"key": "continuous", "label": "phyloP continuous"},
        ],
        "strategies": [{"key": value, "label": strategy_label(value)} for value in strategies],
        "variantTypes": [{"key": key, "label": label} for key, label in VARIANT_TYPE_OPTIONS],
        "targetContexts": [{"key": key, "label": label} for key, label in TARGET_CONTEXT_OPTIONS],
        "consequences": [{"key": key, "label": label} for key, label in CONSEQUENCE_OPTIONS],
        "primary": dataframe_records(primary),
        "fixedDetail": dataframe_records(validation.fixed_bins),
        "continuousDetail": dataframe_records(validation.distributions),
    }
    payload_json = json.dumps(payload, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    return f"""
    <div class="analysis-controls" id="clinvar-association-controls">
      <label>Analysis<select data-role="mode"></select></label>
      <label>Variant type<select data-role="variant-type"></select></label>
      <label>Target context<select data-role="target-context"></select></label>
      <label>Consequence subset<select data-role="consequence"></select></label>
    </div>
    <div id="clinvar-association-status" class="analysis-note" hidden></div>
    <div id="clinvar-association-forest" class="analysis-plot"></div>
    <div id="clinvar-association-results"></div>
    <div class="analysis-controls analysis-controls-single" id="clinvar-association-strategy-control">
      <label>Inspect strategy<select data-role="strategy"></select></label>
    </div>
    <div id="clinvar-association-detail-plot" class="analysis-plot"></div>
    <div id="clinvar-association-detail-table"></div>
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


def target_null_interval_figure(
    frame: pd.DataFrame,
    *,
    title: str,
    xaxis_title: str,
    observed_name: str,
    null_name: str,
    tickformat: str | None = None,
    log_x: bool = False,
) -> go.Figure:
    plot = frame.copy()
    plot["Strategy"] = plot["strategy"].map(strategy_label)
    plot["difference"] = plot["observed_value"] - plot["null_value"]
    plot = plot.sort_values("difference", ascending=False, kind="mergesort")
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=plot["null_value"],
            y=plot["Strategy"],
            mode="markers",
            marker={"color": "#8c8c8c", "size": 8},
            error_x={
                "type": "data",
                "symmetric": False,
                "array": plot["null_ci_high"] - plot["null_value"],
                "arrayminus": plot["null_value"] - plot["null_ci_low"],
                "color": "#8c8c8c",
            },
            name=null_name,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=plot["observed_value"],
            y=plot["Strategy"],
            mode="markers",
            marker={"color": "#2166ac", "size": 10, "symbol": "diamond"},
            error_x={
                "type": "data",
                "symmetric": False,
                "array": plot["observed_ci_high"] - plot["observed_value"],
                "arrayminus": plot["observed_value"] - plot["observed_ci_low"],
                "color": "#2166ac",
            },
            name=observed_name,
        )
    )
    xaxis: dict[str, object] = {"title": xaxis_title}
    if tickformat:
        xaxis["tickformat"] = tickformat
    if log_x:
        xaxis["type"] = "log"
    fig.update_layout(
        title=title,
        xaxis=xaxis,
        yaxis={"categoryorder": "array", "categoryarray": plot["Strategy"].tolist()[::-1]},
    )
    compact_figure(fig, height=max(360, 52 * len(plot) + 120), show_x_title=True)
    return fig


def clinvar_class_null_figure(frame: pd.DataFrame) -> go.Figure:
    class_order = ["B/LB", "P/LP", "VUS", "Other"]
    colors = [CLINVAR_COLORS[category] for category in class_order]
    benign = frame[frame["clinvar_class"].eq("B/LB")].copy()
    benign["difference"] = benign["observed_value"] - benign["null_value"]
    strategies = benign.sort_values("difference", ascending=False, kind="mergesort")["strategy"].tolist()
    if not strategies:
        strategies = sorted(frame["strategy"].astype(str).unique())
    columns = 2 if len(strategies) > 1 else 1
    rows = math.ceil(len(strategies) / columns)
    fig = make_subplots(
        rows=rows,
        cols=columns,
        subplot_titles=[strategy_label(strategy) for strategy in strategies],
        shared_yaxes=True,
        vertical_spacing=min(0.16, 0.35 / max(rows, 1)),
    )
    for index, strategy in enumerate(strategies):
        row = index // columns + 1
        column = index % columns + 1
        values = frame[frame["strategy"].eq(strategy)].set_index("clinvar_class").reindex(class_order)
        fig.add_trace(
            go.Bar(
                x=class_order,
                y=values["observed_value"],
                marker_color=colors,
                error_y={
                    "type": "data",
                    "symmetric": False,
                    "array": values["observed_ci_high"] - values["observed_value"],
                    "arrayminus": values["observed_value"] - values["observed_ci_low"],
                    "color": "#2166ac",
                },
                name="GAPH (95% paired bootstrap interval)",
                showlegend=index == 0,
                hovertemplate="%{x}<br>GAPH: %{y:.1%}<extra></extra>",
            ),
            row=row,
            col=column,
        )
        fig.add_trace(
            go.Scatter(
                x=class_order,
                y=values["null_value"],
                mode="markers",
                marker={"color": "#4d4d4d", "size": 8},
                error_y={
                    "type": "data",
                    "symmetric": False,
                    "array": values["null_ci_high"] - values["null_value"],
                    "arrayminus": values["null_value"] - values["null_ci_low"],
                    "color": "#4d4d4d",
                },
                name="Matched control (95% paired bootstrap interval)",
                showlegend=index == 0,
                hovertemplate="%{x}<br>Matched null: %{y:.1%}<extra></extra>",
            ),
            row=row,
            col=column,
        )
    fig.update_yaxes(tickformat=".0%", rangemode="tozero", title_text="Classified ClinVar hits")
    fig.update_layout(
        title="ClinVar class composition",
        barmode="group",
        height=max(390, 265 * rows),
        margin={"l": 65, "r": 25, "t": 75, "b": 90},
        legend={"orientation": "h", "x": 0.5, "xanchor": "center", "y": -0.14, "yanchor": "top"},
    )
    return fig


def build_target_space_null_sections(
    analysis: TargetSpaceNullAnalysis | None,
    include_plotly: bool,
    *,
    enabled: bool = True,
) -> list[str]:
    sections = ["<h2>Matched Control</h2>"]
    if not enabled:
        sections.append(
            "<p>Matched Control was disabled for this report run. Enable it with "
            "<code>--target-space-null</code>; this analysis uses Ensembl VEP and may take hours.</p>"
        )
        return sections
    if analysis is None:
        sections.append("<p>No target-space-null analysis is available for this run.</p>")
        return sections
    summary = analysis.summary.copy()
    if summary.empty:
        sections.append("<p>No consequence-matched target-space controls could be constructed.</p>")
        return sections

    sampled = analysis.manifest.get("sampled_focal_count", 0)
    matched = analysis.manifest.get("matched_focal_count", 0)
    sections.append(
        metric_cards(
            [
                (
                    "Sample cap per strategy (input)",
                    format_int(analysis.manifest.get("inputs", {}).get("sample_size_per_strategy", 0)),
                ),
                ("Sampled / matched focal SNVs", f"{format_int(sampled)} / {format_int(matched)}"),
                ("Matched-set bootstrap resamples (input)", format_int(analysis.resamples)),
            ]
        )
    )
    conservation_status = analysis.manifest.get("conservation", {}).get("status", "")
    if conservation_status != "complete":
        error = analysis.manifest.get("conservation", {}).get("error", "")
        sections.append(f"<p>phyloP annotation was incomplete: {error or conservation_status}</p>")

    sections.append("<h3>Conservation</h3>")
    ecdf = analysis.ecdf.copy()
    if not ecdf.empty:
        ecdf["Strategy"] = ecdf["strategy"].map(strategy_label)
        ecdf["set"] = ecdf["set"].replace({"Matched target-space null": "Matched control"})
        fig_ecdf = px.line(
            ecdf,
            x="phyloP100way",
            y="fraction_leq",
            color="set",
            facet_col="Strategy",
            facet_col_wrap=2,
            title="phyloP100way distributions: GAPH and matched control",
            labels={"fraction_leq": "Cumulative fraction", "set": ""},
            color_discrete_map={"GAPH": "#2166ac", "Matched control": "#8c8c8c"},
        )
        fig_ecdf.for_each_annotation(lambda item: item.update(text=item.text.split("=")[-1]))
        fig_ecdf.update_yaxes(tickformat=".0%")
        compact_figure(fig_ecdf, height=max(420, 260 * math.ceil(ecdf["Strategy"].nunique() / 2)))
        sections.append(fig_html(fig_ecdf, include_plotlyjs=include_plotly))

    plot = summary.copy()
    plot["Strategy"] = plot["strategy"].map(strategy_label)
    plot = plot.sort_values("median_difference", ascending=False, kind="mergesort")
    fig = go.Figure()
    for row in plot.itertuples(index=False):
        fig.add_trace(
            go.Scatter(
                x=[row.null_median, row.observed_median],
                y=[row.Strategy, row.Strategy],
                mode="lines",
                line={"color": "#c7c7c7", "width": 2},
                hoverinfo="skip",
                showlegend=False,
            )
        )
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
            name="Matched-control median (95% paired bootstrap interval)",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=plot["observed_median"],
            y=plot["Strategy"],
            mode="markers",
            marker={"color": "#2166ac", "size": 10, "symbol": "diamond"},
            error_x={
                "type": "data",
                "symmetric": False,
                "array": plot["observed_ci_high"] - plot["observed_median"],
                "arrayminus": plot["observed_median"] - plot["observed_ci_low"],
                "color": "#2166ac",
            },
            name="GAPH median (95% paired bootstrap interval)",
        )
    )
    fig.update_layout(
        title="phyloP100way median: GAPH vs matched control",
        xaxis_title="phyloP100way",
        yaxis={"categoryorder": "array", "categoryarray": plot["Strategy"].tolist()[::-1]},
    )
    compact_figure(fig, height=max(360, 52 * len(plot) + 120), show_x_title=True)
    sections.append(fig_html(fig, include_plotlyjs=False if not ecdf.empty else include_plotly))

    gnomad = analysis.gnomad_summary.copy()
    if not gnomad.empty:
        sections.append("<h3>gnomAD</h3>")
        found = gnomad[gnomad["metric"].eq("found_fraction")]
        if not found.empty:
            fig_gnomad_found = target_null_interval_figure(
                found,
                title="Exact alleles found in gnomAD",
                xaxis_title="Fraction found",
                observed_name="GAPH fraction (95% paired bootstrap interval)",
                null_name="Matched-control fraction (95% paired bootstrap interval)",
                tickformat=".0%",
            )
            sections.append(fig_html(fig_gnomad_found, include_plotlyjs=False))
        af = gnomad[gnomad["metric"].eq("median_af")]
        if not af.empty and (af[["observed_value", "null_value"]].gt(0).any(axis=None)):
            fig_gnomad_af = target_null_interval_figure(
                af,
                title="gnomAD allele frequency among exact hits",
                xaxis_title="Median allele frequency (log scale)",
                observed_name="GAPH median AF (95% paired bootstrap interval)",
                null_name="Matched-control median AF (95% paired bootstrap interval)",
                log_x=True,
            )
            sections.append(fig_html(fig_gnomad_af, include_plotlyjs=False))
        gnomad_manifest = analysis.manifest.get("external_evidence", {}).get("gnomad", {})
        if gnomad_manifest.get("failed_region_count", 0):
            sections.append(
                f"<p>gnomAD evidence is incomplete: {format_int(gnomad_manifest['failed_region_count'])} "
                "region request(s) failed. Failed regions were treated as missing, not as absence.</p>"
            )

    clinvar = analysis.clinvar_summary.copy()
    if not clinvar.empty:
        sections.append("<h3>ClinVar</h3>")
        fig_clinvar_found = target_null_interval_figure(
            clinvar,
            title="Exact alleles found in ClinVar",
            xaxis_title="Fraction found",
            observed_name="GAPH fraction (95% paired bootstrap interval)",
            null_name="Matched-control fraction (95% paired bootstrap interval)",
            tickformat=".0%",
        )
        sections.append(fig_html(fig_clinvar_found, include_plotlyjs=False))
        if not analysis.clinvar_class_summary.empty:
            fig_clinvar_class = clinvar_class_null_figure(analysis.clinvar_class_summary)
            sections.append(fig_html(fig_clinvar_class, include_plotlyjs=False))

    return sections


def build_target_space_null_qc_sections(analysis: TargetSpaceNullAnalysis) -> list[str]:
    sections = ["<details><summary>Matched-control QC</summary>"]
    focal_vep = analysis.manifest.get("focal_vep", {})
    sections.append(
        table_html(
            pd.DataFrame(
                [
                    {"Metric": "VEP backend", "Value": focal_vep.get("backend", "")},
                    {"Metric": "VEP release", "Value": focal_vep.get("release", "")},
                    {
                        "Metric": "VEP-annotated focal SNVs",
                        "Value": analysis.manifest.get("vep_annotated_focal_count", 0),
                    },
                ]
            ),
            classes="table table-sm table-striped",
        )
    )

    table = analysis.summary.rename(
        columns={
            "strategy": "Strategy",
            "matched_focals": "Matched SNVs",
            "observed_median": "GAPH median",
            "observed_ci_low": "GAPH median Q2.5",
            "observed_ci_high": "GAPH median Q97.5",
            "null_median": "Matched-control median",
            "null_ci_low": "Matched-control median Q2.5",
            "null_ci_high": "Matched-control median Q97.5",
            "median_difference": "Median difference",
            "difference_ci_low": "Difference Q2.5",
            "difference_ci_high": "Difference Q97.5",
            "valid_resamples": "Valid bootstrap resamples",
        }
    )
    if not table.empty:
        table["Strategy"] = table["Strategy"].map(strategy_label)
        sections.append("<h4>Strategy summary</h4>")
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
        sections.append("<h4>Matching yield by consequence</h4>")
        sections.append(table_html(matching, classes="table table-sm table-striped"))

    consequence = analysis.consequence_summary.copy()
    if not consequence.empty:
        consequence = consequence.rename(
            columns={
                "strategy": "Strategy",
                "primary_consequence": "Primary VEP consequence",
                "matched_focals": "Matched SNVs",
                "observed_median": "GAPH median",
                "observed_ci_low": "GAPH median Q2.5",
                "observed_ci_high": "GAPH median Q97.5",
                "null_median": "Matched-control median",
                "null_ci_low": "Matched-control median Q2.5",
                "null_ci_high": "Matched-control median Q97.5",
                "median_difference": "Median difference",
                "difference_ci_low": "Difference Q2.5",
                "difference_ci_high": "Difference Q97.5",
                "valid_resamples": "Valid bootstrap resamples",
            }
        )
        consequence["Strategy"] = consequence["Strategy"].map(strategy_label)
        consequence = consequence[
            [
                "Strategy",
                "Primary VEP consequence",
                "Matched SNVs",
                "GAPH median",
                "GAPH median Q2.5",
                "GAPH median Q97.5",
                "Matched-control median",
                "Matched-control median Q2.5",
                "Matched-control median Q97.5",
                "Median difference",
                "Difference Q2.5",
                "Difference Q97.5",
                "Valid bootstrap resamples",
            ]
        ]
        sections.append("<h4>Results by primary VEP consequence</h4>")
        sections.append(table_html(consequence, classes="table table-sm table-striped"))

    outcome_frames = []
    if not analysis.gnomad_summary.empty:
        gnomad = analysis.gnomad_summary.copy()
        gnomad["Outcome"] = gnomad["metric"].map(
            {
                "found_fraction": "gnomAD found fraction",
                "median_af": "gnomAD median AF among exact hits",
            }
        )
        outcome_frames.append(gnomad)
    if not analysis.clinvar_summary.empty:
        clinvar = analysis.clinvar_summary.copy()
        clinvar["Outcome"] = "ClinVar found fraction"
        outcome_frames.append(clinvar)
    if not analysis.clinvar_class_summary.empty:
        clinvar_class = analysis.clinvar_class_summary.copy()
        clinvar_class["Outcome"] = "ClinVar class: " + clinvar_class["clinvar_class"].astype(str)
        outcome_frames.append(clinvar_class)
    if outcome_frames:
        outcomes = pd.concat(outcome_frames, ignore_index=True)
        outcomes = outcomes.rename(
            columns={
                "strategy": "Strategy",
                "matched_focals": "Matched SNVs",
                "observed_value": "GAPH statistic",
                "observed_ci_low": "GAPH Q2.5",
                "observed_ci_high": "GAPH Q97.5",
                "null_value": "Matched-control statistic",
                "null_ci_low": "Matched-control Q2.5",
                "null_ci_high": "Matched-control Q97.5",
                "difference": "Paired difference",
                "difference_ci_low": "Difference Q2.5",
                "difference_ci_high": "Difference Q97.5",
                "valid_resamples": "Valid bootstrap resamples",
            }
        )
        outcomes["Strategy"] = outcomes["Strategy"].map(strategy_label)
        outcomes = outcomes[
            [
                "Strategy",
                "Outcome",
                "Matched SNVs",
                "GAPH statistic",
                "GAPH Q2.5",
                "GAPH Q97.5",
                "Matched-control statistic",
                "Matched-control Q2.5",
                "Matched-control Q97.5",
                "Paired difference",
                "Difference Q2.5",
                "Difference Q97.5",
                "Valid bootstrap resamples",
            ]
        ]
        sections.append("<h4>External-evidence bootstrap summary</h4>")
        sections.append(table_html(outcomes, classes="table table-sm table-striped"))
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
    taxonomy_summary: pd.DataFrame | None = None,
) -> list[str]:
    files = [
        ("Run Dir", inputs.run_dir),
        ("Fetch Manifest", inputs.fetch_manifest_json),
        ("Variant Annotations", inputs.variant_annotations_tsv),
        ("Variant Strategy Support", inputs.variant_strategy_support_tsv),
        ("Ortholog Evidence Summary", inputs.ortholog_evidence_summary_tsv),
        ("Target Features", inputs.target_features_tsv),
        ("Target Sequences", inputs.target_sequences_dir),
        ("Feature Coverage", inputs.feature_coverage_tsv),
        ("Alignment Segments", inputs.alignment_segments_tsv),
        ("Strategy Summary", inputs.strategy_summary_tsv),
        ("Taxonomy Summary", inputs.taxonomy_summary_tsv),
        ("Annotation Manifest", inputs.annotation_manifest_json),
        ("Alignment Manifest", inputs.alignment_manifest_json),
        ("Output HTML", out_html),
    ]
    if conservation_analysis is not None:
        files.extend(
            [
                ("Candidate phyloP distributions", conservation_analysis.candidate.distributions_path),
                ("Candidate phyloP histograms", conservation_analysis.candidate.histograms_path),
                ("Candidate phyloP manifest", conservation_analysis.candidate.manifest_path),
            ]
        )
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
                ("Candidate contexts excluded from gnomAD", format_int(variant_summary.gnomad_lookup_failed)),
                ("ClinVar cached variants", format_int(annotation_manifest.get("clinvar_cached_variant_count", 0))),
                ("gnomAD cached variants", format_int(annotation_manifest.get("gnomad_cached_variant_count", 0))),
                ("Feature coverage rows", format_int(len(cov))),
            ]
        ),
    ]
    if not failures.empty:
        sections.append("<h3>Annotation Failures</h3>")
        sections.append(table_html(failures, classes="table table-sm table-striped", max_rows=50))
    if taxonomy_summary is not None and not taxonomy_summary.empty:
        shown = taxonomy_summary.copy()
        shown["Taxonomic scope"] = shown["taxonomic_scope"].map(
            lambda value: TAXONOMIC_SCOPE_LABELS.get(str(value), str(value))
        )
        shown["Evidence unit"] = shown["evidence_unit"].map(
            lambda value: EVIDENCE_UNIT_LABELS.get(str(value), str(value))
        )
        shown = shown.rename(
            columns={
                "gene_count": "Genes",
                "ortholog_count": "Selected ortholog rows",
                "taxon_count": "Distinct taxa",
                "unit_count": "Distinct units",
                "orthologs_per_gene_median": "Median orthologs/gene",
                "units_per_gene_median": "Median units/gene",
            }
        )[
            [
                "Taxonomic scope",
                "Evidence unit",
                "Genes",
                "Selected ortholog rows",
                "Distinct taxa",
                "Distinct units",
                "Median orthologs/gene",
                "Median units/gene",
            ]
        ]
        sections.append("<details><summary>Taxonomic evidence scope and grouping</summary>")
        sections.append(table_html(shown, classes="table table-sm table-striped"))
        sections.append("</details>")
    pathogenic_table = pathogenic_variant_table(variant_summary.pathogenic_rows)
    if not pathogenic_table.empty:
        shown = min(len(pathogenic_table), 100)
        sections.append(
            "<details><summary>Top "
            f"{format_int(shown)} of {format_int(variant_summary.pathogenic_variant_count)} "
            "unique P/LP variants</summary>"
        )
        sections.append(
            "<p>Sorted by ClinVar review stars, then supporting SCV count.</p>"
        )
        sections.append(
            table_html(pathogenic_table, classes="table table-sm table-striped", max_rows=100)
        )
        sections.append("</details>")
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
    sections.append(f"<details><summary>{variant_summary.consequence_source} consequence grouping</summary>")
    sections.append(
        f"<p class=\"lead\">The Candidate Profile consequence plots use "
        f"{variant_summary.consequence_source} annotations and group them as follows. "
        "Raw <code>gnomad_csq</code> remains available as provenance.</p>"
    )
    sections.append(
        table_html(
            consequence_grouping_table(variant_summary.consequence_source),
            classes="table table-sm table-striped",
        )
    )
    sections.append("</details>")
    sections.append("<details><summary>Target-context assignment</summary>")
    sections.append(
        "<p class=\"lead\">Candidate and ClinVar validation positions are assigned to one exclusive target context. "
        "Overlapping transcript features use the precedence CDS &gt; UTR &gt; exon &gt; intron; exon sequence "
        "outside CDS/UTR is labelled Other exon, and remaining target sequence is labelled Other. "
        "The ClinVar Association selector exposes All, CDS, UTR, and Intron; Other exon and Other remain included "
        "in the All denominator.</p>"
    )
    sections.append("</details>")
    sections.append("<details><summary>ClinVar validation denominator and statistics</summary>")
    if validation is not None:
        manifest = validation.manifest
        cohort_flow = pd.DataFrame(
            [
                {"Cohort step": "Raw SNV/INDEL alleles", "Alleles": manifest.get("raw_allele_count", 0)},
                {"Cohort step": "Excluded VUS", "Alleles": manifest.get("excluded_vus_count", 0)},
                {"Cohort step": "Excluded missing CLNSIG", "Alleles": manifest.get("excluded_missing_count", 0)},
                {"Cohort step": "Excluded other/conflicting", "Alleles": manifest.get("excluded_other_count", 0)},
                {"Cohort step": "Included B/LB", "Alleles": manifest.get("benign_count", 0)},
                {"Cohort step": "Included P/LP", "Alleles": manifest.get("pathogenic_count", 0)},
                {"Cohort step": "Final validation cohort", "Alleles": manifest.get("usable_allele_count", 0)},
            ]
        )
        sections.append(table_html(cohort_flow, classes="table table-sm table-striped"))
    sections.append(table_html(validation_method_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    sections.append("<details><summary>ClinVar association modes and consequence subsets</summary>")
    sections.append(
        "<p class=\"lead\">All three modes use the same normalized ClinVar allele cohort. They differ only "
        "in whether and how phyloP100way is included.</p>"
    )
    sections.append(table_html(conservation_validation_method_table(conservation_analysis), classes="table table-sm table-striped"))
    consequence_source = validation.consequence_source if validation is not None else "ClinVar MC"
    sections.append(f"<h4>{consequence_source} consequence subsets</h4>")
    sections.append(
        table_html(
            validation_consequence_grouping_table(consequence_source),
            classes="table table-sm table-striped",
        )
    )
    if conservation_analysis is not None:
        visibility_summary, hidden_views = hidden_clinvar_association_views(
            conservation_analysis.validation
        )
        sections.append("<h4>Adaptive selector visibility</h4>")
        sections.append(
            "<p>Consequence options are hidden from the interactive view only when no strategy has an "
            "estimable result for the selected analysis, variant type, and target context. Hidden combinations "
            "remain listed here; no minimum sample-size threshold is applied.</p>"
        )
        sections.append(
            table_html(visibility_summary, classes="table table-sm table-striped")
        )
        if not hidden_views.empty:
            sections.append(
                table_html(hidden_views, classes="table table-sm table-striped")
            )
    sections.append("</details>")
    sections.append("<details><summary>Candidate-wide phyloP stratification</summary>")
    sections.append(
        table_html(
            pd.DataFrame(
                [
                    {
                        "Step": "Eligible candidate alleles",
                        "Definition": "Normalized lookup_status=ok SNVs and indels with a defined phyloP score basis.",
                    },
                    {
                        "Step": "Allele score",
                        "Definition": (
                            "SNV substituted base; deletion mean across deleted reference bases without VCF padding; "
                            "insertion mean across the two flanking bases. All required bases must have a score."
                        ),
                    },
                    {
                        "Step": "Unit and strata",
                        "Definition": (
                            "Unique variant_key x strategy records with a completed gnomAD lookup, split by presence "
                            "of an exact gnomAD AF annotation. Failed lookups are excluded from both strata."
                        ),
                    },
                    {
                        "Step": "Distribution",
                        "Definition": (
                            "Exact percentiles from 0 through 100, plus relative-frequency histograms using shared "
                            "Found/Not-found bins selected by the Freedman-Diaconis rule and capped at 80 bins. "
                            "Box plots use Tukey 1.5-IQR whiskers."
                        ),
                    },
                    {
                        "Step": "Shared read",
                        "Definition": (
                            "A cold report reads the union of candidate and ClinVar-required positions from bigWig once; "
                            "both analyses reuse that positional score map."
                        ),
                    },
                ]
            ),
            classes="table table-sm table-striped",
        )
    )
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
        "<p class=\"lead\">Candidate Profile coverage plots use the normalized feature-level table emitted by the alignment stage.</p>"
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
        sections.extend(build_target_space_null_qc_sections(negative_controls))
        control_files = [
            ("Negative-control manifest", negative_controls.manifest_path),
            ("Target-space-null rows", negative_controls.matched_path),
            ("Target-space-null phyloP annotations", negative_controls.conservation_path),
            ("VEP consequence cache", negative_controls.vep_cache_path),
            ("Target-space-null external evidence", negative_controls.external_evidence_path),
            ("External-evidence manifest", negative_controls.external_evidence_manifest_path),
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
        <title>GAPH Variant Analysis</title>
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
            .overview-table {{ width: 100%; font-size: 14px; }}
            .overview-table th, .overview-table td {{ padding: 9px 10px; }}
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
            .analysis-controls-single {{ grid-template-columns: minmax(180px, 260px); }}
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
        <h1>GAPH Variant Analysis</h1>
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
    if args.target_space_null and args.vep_forks < 1:
        raise ValueError("--vep-forks must be >= 1")
    if args.target_space_null and args.vep_backend == "local" and args.vep_cache_dir is None:
        raise ValueError("--vep-cache-dir is required with --vep-backend local")
    if args.target_space_null and args.vep_backend == "local" and not args.vep_release:
        raise ValueError("--vep-release is required with --vep-backend local")
    inputs = resolve_run_inputs(args.run_dir, args.annotation_dir)
    validate_report_inputs(inputs)
    annotation_columns = set(
        pd.read_csv(inputs.variant_annotations_tsv, sep="\t", compression="gzip", nrows=0).columns
    )
    use_vep_consequences = {
        "vep_status",
        "vep_primary_consequence",
        "vep_consequence_terms",
    }.issubset(annotation_columns)
    artifact_release = bulk_vep_release(inputs) if use_vep_consequences else None
    if artifact_release and args.vep_release and str(args.vep_release) != artifact_release:
        raise ValueError(
            f"Bulk VEP artifact uses release {artifact_release}, not {args.vep_release}"
        )
    if not args.vep_release:
        args.vep_release = artifact_release
    if use_vep_consequences and not args.vep_release:
        raise ValueError("--vep-release is required when the report uses bulk VEP consequences")
    if use_vep_consequences and args.vep_backend == "local" and args.vep_cache_dir is None:
        raise ValueError("--vep-cache-dir is required with --vep-backend local")
    out_html = resolve_out_html(args, inputs.run_dir)
    report_started = time.perf_counter()
    timings: list[dict[str, object]] = []

    print(f"Streaming {inputs.variant_annotations_tsv}...")
    with timed_stage("Variant summary", timings) as timing:
        variant_summary = build_variant_summary(
            inputs.variant_annotations_tsv,
            inputs.run_dir / "analytics",
            strategy_label,
            target_features_path=inputs.target_features_tsv,
            genes_path=inputs.genes_tsv,
            annotation_failures_path=inputs.annotation_failures_tsv,
            variant_strategy_support_path=inputs.variant_strategy_support_tsv,
        )
        if inputs.ortholog_evidence_summary_tsv.exists():
            available, cells = read_taxonomic_ortholog_evidence(
                inputs.ortholog_evidence_summary_tsv
            )
            variant_summary = replace(
                variant_summary,
                ortholog_evidence_available=available,
                ortholog_evidence_cells=cells,
            )
        timing["Details"] = "cache hit" if variant_summary.cache_hit else "cache miss"

    with timed_stage("Run summary inputs", timings):
        cov = read_feature_coverage(inputs.feature_coverage_tsv)
        alignment_summary = read_strategy_summary(inputs.strategy_summary_tsv)
        fetch_manifest = read_json(inputs.fetch_manifest_json)
        input_gene_count = int(
            fetch_manifest.get("unique_gene_count") or read_input_gene_count(inputs.genes_tsv)
        )
        failures = read_failures(inputs.annotation_failures_tsv)
        annotation_manifest = read_json(inputs.annotation_manifest_json)
        alignment_manifest = read_json(inputs.alignment_manifest_json)
        taxonomy_summary = read_taxonomy_summary(inputs.taxonomy_summary_tsv)

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
            "gnomAD Eligible",
            "gnomAD lookup failed",
            "gnomAD found %",
            "Genes with result",
            "Orthologs evaluated",
            "Orthologs aligned",
            "Orthologs aligned %",
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
            use_vep_consequences=use_vep_consequences,
            vep_backend=args.vep_backend,
            vep_release=args.vep_release,
            vep_executable=args.vep_executable,
            vep_cache_dir=args.vep_cache_dir,
            vep_forks=args.vep_forks,
        )

    print("Computing conservation-adjusted ClinVar validation...")
    with timed_stage("Conservation-adjusted validation", timings):
        conservation_analysis = build_conservation_analysis(
            inputs=inputs,
            validation=validation,
            strategies=strategies,
            eligible_gene_ids_by_strategy=alignment_gene_ids_by_strategy(cov),
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
                clinvar_vcf=args.clinvar_vcf.expanduser().resolve(),
                strategies=strategies,
                sample_size_per_strategy=args.target_space_null_sample_size,
                resamples=args.target_space_null_resamples,
                seed=args.target_space_null_seed,
                gnomad_cache_dir=args.gnomad_cache_dir,
                vep_backend=args.vep_backend,
                vep_release=args.vep_release,
                vep_executable=args.vep_executable,
                vep_cache_dir=args.vep_cache_dir,
                vep_forks=args.vep_forks,
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

    candidate_sections = build_variant_sections(variant_summary, strategy_stats, include_plotly=True)
    candidate_sections.extend(
        build_clinvar_gnomad_sections(variant_summary, strategy_stats_full, include_plotly=False)
    )
    candidate_sections.extend(build_feature_sections(cov, include_plotly=False))
    sections = [
        (
            "overview",
            "Overview",
            build_overview(
                variant_summary,
                cov,
                strategy_stats,
                annotation_manifest,
                input_gene_count,
            ),
        ),
        ("candidates", "Candidate Profile", candidate_sections),
        (
            "ortholog-evidence",
            "Ortholog Evidence",
            build_ortholog_evidence_sections(
                variant_summary,
                include_plotly=False,
                taxonomy_summary=taxonomy_summary,
            ),
        ),
        (
            "gnomad-stratification",
            "gnomAD Stratification",
            build_gnomad_stratification_sections(
                variant_summary,
                strategy_stats_full,
                conservation_analysis.candidate,
                include_plotly=True,
            ),
        ),
        (
            "target-space-null",
            "Matched Control",
            build_target_space_null_sections(
                negative_controls,
                include_plotly=True,
                enabled=args.target_space_null,
            ),
        ),
        ("clinvar-association", "ClinVar Association", build_clinvar_association_sections(conservation_analysis)),
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
                taxonomy_summary,
            ),
        ),
    ]

    print(f"Writing report to {out_html}...")
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(render_html(sections))
    print(f"Done in {time.perf_counter() - report_started:.3f} s")


if __name__ == "__main__":
    main()
