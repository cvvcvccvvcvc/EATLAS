"""Discovery and validation of completed-run inputs."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from analytics.analyses.variant_summary_aggregation import (
    resolve_variant_aggregation_source,
)
from analytics.io.alignment_aggregates import resolve_alignment_aggregate_paths
from analytics.io.annotation_support import resolve_annotation_support_paths
from analytics.io.taxonomy_summary import resolve_taxonomy_summary_path


@dataclass(frozen=True)
class RunInputs:
    run_dir: Path
    fetch_manifest_json: Path
    genes_tsv: Path
    target_features_tsv: Path
    target_sequences_dir: Path
    variant_annotations_source: Path
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
    annotation_manifest_json = annotation_dir / "manifest.json"
    variant_annotations_source = resolve_variant_annotations_source(
        annotation_manifest_json
    )
    alignment_aggregates = resolve_alignment_aggregate_paths(run_dir)
    annotation_support = resolve_annotation_support_paths(run_dir)

    inputs = RunInputs(
        run_dir=run_dir,
        fetch_manifest_json=run_dir / "fetch" / "manifest.json",
        genes_tsv=run_dir / "fetch" / "genes.tsv.gz",
        target_features_tsv=run_dir / "fetch" / "target_features.tsv.gz",
        target_sequences_dir=run_dir / "fetch" / "sequences" / "targets",
        variant_annotations_source=variant_annotations_source,
        variant_strategy_support_tsv=annotation_support.variant_strategy_support_tsv,
        ortholog_evidence_summary_tsv=annotation_support.ortholog_evidence_summary_tsv,
        annotation_manifest_json=annotation_manifest_json,
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


def resolve_variant_annotations_source(annotation_manifest_path: Path) -> Path:
    """Resolve the pipeline-owned partitioned variant-annotation dataset."""

    descriptor = _variant_annotation_descriptor(annotation_manifest_path)
    declared = Path(str(descriptor.get("path") or ""))
    if declared.as_posix() != "variant_annotations/manifest.json":
        raise ValueError(
            "Annotation manifest has an unsupported variant_annotations path: "
            f"{annotation_manifest_path}"
        )
    source = (annotation_manifest_path.parent / declared).resolve()
    if not source.is_file():
        raise FileNotFoundError(
            f"Missing pipeline variant-annotation dataset manifest: {source}"
        )
    dataset_manifest = read_json(source)
    if dataset_manifest != descriptor:
        raise ValueError(
            "Annotation variant_annotations descriptor does not match its dataset "
            f"manifest: {source}"
        )
    return source


def variant_annotation_descriptor(inputs: RunInputs) -> dict[str, object]:
    descriptor = _variant_annotation_descriptor(inputs.annotation_manifest_json)
    declared_source = resolve_variant_annotations_source(inputs.annotation_manifest_json)
    if inputs.variant_annotations_source.resolve() != declared_source:
        raise ValueError(
            "Report inputs do not use the variant_annotations dataset declared by "
            f"{inputs.annotation_manifest_json}"
        )
    return descriptor


def variant_annotation_release(descriptor: dict[str, object]) -> str:
    config = descriptor["vep_config"]
    release = config.get("release") if isinstance(config, dict) else None
    if not release:
        raise ValueError("Variant annotation dataset is missing vep_config.release")
    return str(release)


def _variant_annotation_descriptor(annotation_manifest_path: Path) -> dict[str, object]:
    if not annotation_manifest_path.is_file():
        raise FileNotFoundError(
            f"Missing pipeline annotation manifest: {annotation_manifest_path}"
        )
    manifest = read_json(annotation_manifest_path)
    if (
        manifest.get("stage") != "annotation"
        or manifest.get("schema") != "normalized_annotation_evidence_v3"
    ):
        raise ValueError(
            f"Unsupported pipeline annotation contract: {annotation_manifest_path}"
        )
    descriptor = manifest.get("variant_annotations")
    if not isinstance(descriptor, dict):
        raise ValueError(
            "Annotation manifest does not declare variant_annotations: "
            f"{annotation_manifest_path}"
        )
    config = descriptor.get("vep_config")
    status_counts = descriptor.get("vep_status_counts")
    if not isinstance(config, dict) or not config:
        raise ValueError("Variant annotation dataset has no VEP configuration")
    if config.get("backend") not in {"rest", "local"}:
        raise ValueError("Variant annotation dataset has invalid vep_config.backend")
    if not isinstance(status_counts, dict):
        raise ValueError("Variant annotation dataset has no VEP status counts")
    try:
        raw_row_count = descriptor["row_count"]
        if isinstance(raw_row_count, bool):
            raise TypeError
        row_count = int(raw_row_count)
        normalized_counts = {}
        for status, raw_count in status_counts.items():
            if not str(status) or isinstance(raw_count, bool):
                raise TypeError
            normalized_counts[str(status)] = int(raw_count)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Variant annotation dataset has invalid VEP status counts") from exc
    if (
        row_count < 0
        or any(count < 0 for count in normalized_counts.values())
        or sum(normalized_counts.values()) != row_count
    ):
        raise ValueError("Variant annotation VEP status counts do not match row_count")
    variant_annotation_release(descriptor)
    return descriptor


def validate_report_inputs(inputs: RunInputs) -> None:
    """Fail before expensive work when a production table contract is incompatible."""
    descriptor = variant_annotation_descriptor(inputs)
    resolve_variant_aggregation_source(inputs.variant_annotations_source)
    contracts = {
        inputs.variant_annotations_source: {
            "variant_key",
            "gene_id",
            "event_type",
            "ref",
            "alt",
            "lookup_status",
            "strategies",
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
    for path, required in contracts.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing report input: {path}")
        if path == inputs.variant_annotations_source:
            header = set(str(field) for field in descriptor.get("fields", []))
        else:
            compression = "gzip" if path.suffix == ".gz" else None
            header = set(pd.read_csv(path, sep="\t", compression=compression, nrows=0).columns)
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
        raise FileNotFoundError(f"Missing taxonomy summary: {path}")
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
