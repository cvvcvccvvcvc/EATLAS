"""Discovery and validation of completed-run inputs."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from analytics.io.alignment_aggregates import resolve_alignment_aggregate_paths
from analytics.io.annotation_support import resolve_annotation_support_paths
from analytics.io.artifacts import file_identity, path_metadata
from analytics.io.taxonomy_summary import resolve_taxonomy_summary_path


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
    alignment_manifest_json: Path
    strategy_summary_tsv: Path
    taxonomy_summary_tsv: Path
    cohort_manifest_json: Path | None = None
    source_run_dirs: tuple[Path, ...] = ()
    cohort_id: str | None = None

    @property
    def is_cohort(self) -> bool:
        return self.cohort_id is not None


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

    annotation_dir = run_dir / "annotation"
    source_annotations_tsv = annotation_dir / "variant_annotations.tsv.gz"
    if not source_annotations_tsv.exists():
        raise FileNotFoundError(
            f"Missing variant_annotations.tsv.gz under {annotation_dir}. "
            "Run the annotation stage before building this report."
        )
    variant_annotations_tsv = resolve_vep_variant_annotations(
        run_dir,
        source_annotations_tsv,
    )
    alignment_aggregates = resolve_alignment_aggregate_paths(run_dir)
    annotation_support = resolve_annotation_support_paths(run_dir)

    inputs = RunInputs(
        run_dir=run_dir,
        fetch_manifest_json=run_dir / "fetch" / "manifest.json",
        genes_tsv=run_dir / "fetch" / "genes.tsv.gz",
        target_features_tsv=run_dir / "fetch" / "target_features.tsv.gz",
        target_sequences_dir=run_dir / "fetch" / "sequences" / "targets",
        variant_annotations_tsv=variant_annotations_tsv,
        variant_strategy_support_tsv=annotation_support.variant_strategy_support_tsv,
        ortholog_evidence_summary_tsv=annotation_support.ortholog_evidence_summary_tsv,
        annotation_manifest_json=annotation_dir / "manifest.json",
        annotation_failures_tsv=annotation_dir / "failures.tsv.gz",
        feature_coverage_tsv=alignment_aggregates.feature_coverage_tsv,
        alignment_manifest_json=run_dir / "alignment" / "manifest.json",
        strategy_summary_tsv=alignment_aggregates.strategy_summary_tsv,
        taxonomy_summary_tsv=resolve_taxonomy_summary_path(run_dir),
    )
    if not inputs.genes_tsv.exists():
        raise FileNotFoundError("Missing fetch/genes.tsv.gz under --run-dir.")
    if not inputs.target_features_tsv.exists():
        raise FileNotFoundError("Missing fetch/target_features.tsv.gz under --run-dir.")
    if not inputs.target_sequences_dir.exists():
        raise FileNotFoundError("Missing fetch/sequences/targets under --run-dir.")
    return inputs


def resolve_vep_variant_annotations(run_dir: Path, source: Path) -> Path:
    """Require the finalized bulk-VEP artifact matching the source table."""

    artifact_dir = run_dir / "analytics" / "vep_consequences"
    manifest_path = artifact_dir / "manifest.json"
    output_path = artifact_dir / "variant_annotations.vep.tsv.gz"
    if not manifest_path.exists() and not output_path.exists():
        raise FileNotFoundError(
            f"Missing finalized bulk VEP artifact under {artifact_dir}. "
            "Build it with `python -m analytics.vep_annotation prepare --run-dir <run>`, "
            "annotate all declared partitions, and run "
            "`python -m analytics.vep_annotation finalize --run-dir <run>` before "
            "generating the report. Individual non-ok vep_status values are allowed in "
            "a finalized artifact."
        )
    if not manifest_path.exists() or not output_path.exists():
        raise ValueError(f"Incomplete bulk VEP artifact under {artifact_dir}")

    manifest = read_json(manifest_path)
    source_identity = path_metadata(source)
    if manifest.get("status") != "complete" or manifest.get("source") != source_identity:
        raise ValueError(f"Bulk VEP artifact does not match {source}")
    output_identity = file_identity(output_path)
    if manifest.get("output") != output_identity:
        raise ValueError(f"Bulk VEP output metadata changed: {output_path}")
    return output_path


def bulk_vep_manifest(inputs: RunInputs) -> dict:
    if inputs.is_cohort:
        manifest = read_json(
            inputs.run_dir / "analytics" / "vep_consequences" / "manifest.json"
        )
        if manifest.get("status") != "complete" or manifest.get("cohort") is not True:
            raise ValueError("Cohort bulk VEP manifest is incomplete")
        return manifest
    expected = inputs.run_dir / "analytics" / "vep_consequences" / "variant_annotations.vep.tsv.gz"
    if inputs.variant_annotations_tsv != expected:
        raise ValueError("Report inputs do not use the finalized bulk VEP artifact")
    return read_json(expected.parent / "manifest.json")


def bulk_vep_release(inputs: RunInputs) -> str:
    manifest = bulk_vep_manifest(inputs)
    release = manifest.get("config", {}).get("release")
    if not release:
        raise ValueError("Bulk VEP manifest is missing config.release")
    return str(release)


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
            "vep_status",
            "vep_primary_consequence",
            "vep_consequence_terms",
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
            "alt_support_genus_count",
        },
    }
    for path, required in contracts.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing report input: {path}")
        if path == inputs.variant_annotations_tsv and path.suffix == ".json":
            from analytics.analyses.variant_summary_aggregation import (
                resolve_variant_aggregation_source,
            )

            header = set(resolve_variant_aggregation_source(path).columns)
        else:
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
        raise FileNotFoundError(f"Missing resolved strategy-summary input: {path}")
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
