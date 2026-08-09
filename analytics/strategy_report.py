#!/usr/bin/env python3
"""Build an HTML report for one completed GAPH run."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd

from analytics.analyses.clinvar_validation import build_validation
from analytics.analyses.conservation_analysis import (
    alignment_gene_ids_by_strategy,
    build_conservation_analysis,
)
from analytics.analyses.conservation_validation import validate_firth_runtime
from analytics.analyses.matched_control import build_target_space_null
from analytics.analyses.observed_variant_store import (
    build_or_load_observed_variant_store,
)
from analytics.analyses.variant_summary import build_variant_summary
from analytics.io.run_inputs import (
    bulk_vep_release,
    read_failures,
    read_feature_coverage,
    read_input_gene_count,
    read_json,
    read_strategy_summary,
    read_taxonomy_summary,
    resolve_out_html,
    resolve_run_inputs,
    validate_report_inputs,
)
from analytics.io.artifacts import write_text_atomic
from analytics.io.performance import PerformanceProfile
from analytics.reporting.components import strategy_label
from analytics.reporting.conservation import build_clinvar_association_sections
from analytics.reporting.document import render_html
from analytics.reporting.matched_control import build_target_space_null_sections
from analytics.reporting.ortholog_evidence import build_ortholog_evidence_sections
from analytics.reporting.overview import build_overview, merge_alignment_summary
from analytics.reporting.qc import build_methods_sections
from analytics.reporting.variant_profile import (
    build_clinvar_gnomad_sections,
    build_feature_sections,
    build_gnomad_stratification_sections,
    build_variant_sections,
)


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
        "--vep-result-cache-dir",
        type=Path,
        default=_default_vep_result_cache_dir(),
        help="Shared cross-run cache for completed VEP results.",
    )
    parser.add_argument(
        "--vep-result-cache-tile-size-bp",
        type=int,
        default=int(os.environ.get("GAPH_VEP_RESULT_CACHE_TILE_SIZE_BP", "1000000")),
        help="Genomic tile size for the shared VEP result cache. Default: 1000000.",
    )
    parser.add_argument(
        "--vep-forks",
        type=int,
        default=int(os.environ.get("GAPH_VEP_FORKS", "4")),
        help="Worker processes for local VEP. Default: 4.",
    )
    parser.add_argument(
        "--firth-workers",
        type=int,
        default=_default_firth_workers(),
        help=(
            "Parallel workers for independent Firth models. Defaults to available "
            "Slurm/host CPUs capped at 8."
        ),
    )
    return parser.parse_args()


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _default_vep_result_cache_dir() -> Path | None:
    configured = os.environ.get("GAPH_VEP_RESULT_CACHE_DIR")
    if configured:
        return Path(configured)
    gaph_root = os.environ.get("GAPH_ROOT")
    return Path(gaph_root) / "cache" / "vep_results" if gaph_root else None


def _default_firth_workers() -> int:
    configured = os.environ.get("GAPH_FIRTH_WORKERS")
    if configured:
        return int(configured)
    allocated = os.environ.get("SLURM_CPUS_PER_TASK")
    available = int(allocated) if allocated and allocated.isdigit() else (os.cpu_count() or 1)
    return max(1, min(available, 8))


def main() -> None:
    args = parse_args()
    if args.target_space_null and args.target_space_null_sample_size < 1:
        raise ValueError("--target-space-null-sample-size must be >= 1")
    if args.target_space_null and args.target_space_null_resamples < 100:
        raise ValueError("--target-space-null-resamples must be >= 100")
    if args.target_space_null and args.vep_forks < 1:
        raise ValueError("--vep-forks must be >= 1")
    if args.firth_workers < 1:
        raise ValueError("--firth-workers must be >= 1")
    if args.vep_result_cache_tile_size_bp < 1:
        raise ValueError("--vep-result-cache-tile-size-bp must be >= 1")
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
    analytics_dir = inputs.run_dir / "analytics"
    performance_path = analytics_dir / "performance" / f"{out_html.stem}.json"
    performance = PerformanceProfile(
        performance_path,
        run_dir=inputs.run_dir,
        report_path=out_html,
        tracked_directory=analytics_dir,
    )

    with performance.stage("Firth runtime preflight") as timing:
        timing["metrics"] = validate_firth_runtime()

    if args.vep_result_cache_dir is None:
        print("Shared VEP result cache: disabled")
    else:
        print(
            "Shared VEP result cache: "
            f"{args.vep_result_cache_dir.expanduser()} "
            f"(tile size {args.vep_result_cache_tile_size_bp} bp)"
        )
    print(f"Streaming {inputs.variant_annotations_tsv}...")
    with performance.stage("Variant summary") as timing:
        variant_summary = build_variant_summary(
            inputs.variant_annotations_tsv,
            analytics_dir,
            strategy_label,
            target_features_path=inputs.target_features_tsv,
            genes_path=inputs.genes_tsv,
            annotation_failures_path=inputs.annotation_failures_tsv,
            variant_strategy_support_path=inputs.variant_strategy_support_tsv,
            ortholog_evidence_summary_path=(
                inputs.ortholog_evidence_summary_tsv
                if inputs.ortholog_evidence_summary_tsv.exists()
                else None
            ),
            performance_profile=performance,
        )
        timing["details"] = "cache hit" if variant_summary.cache_hit else "cache miss"

    with performance.stage("Run summary inputs"):
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
    with performance.stage("Strategy metrics"):
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

    print("Building observed variant store...")
    with performance.stage("Observed variant store") as timing:
        observed_store = build_or_load_observed_variant_store(
            variant_annotations_tsv=inputs.variant_annotations_tsv,
            analytics_dir=analytics_dir,
            strategies=strategies,
        )
        timing["details"] = "cache hit" if observed_store.cache_hit else "cache miss"
        timing["metrics"] = {
            "source_rows": int(observed_store.manifest["source_row_count"]),
            "allele_gene_rows": int(observed_store.manifest["allele_gene_count"]),
            "alleles": int(observed_store.manifest["allele_count"]),
            "store_bytes": int(observed_store.allele_gene_path.stat().st_size)
            + int(observed_store.allele_path.stat().st_size),
        }

    print("Computing ClinVar enrichment...")
    with performance.stage("ClinVar enrichment"):
        validation = build_validation(
            run_dir=inputs.run_dir,
            genes_tsv=inputs.genes_tsv,
            target_sequences_dir=inputs.target_sequences_dir,
            clinvar_vcf=args.clinvar_vcf.expanduser().resolve(),
            strategies=strategies,
            observed_store=observed_store,
            use_vep_consequences=use_vep_consequences,
            vep_backend=args.vep_backend,
            vep_release=args.vep_release,
            vep_executable=args.vep_executable,
            vep_cache_dir=args.vep_cache_dir,
            vep_forks=args.vep_forks,
            vep_result_cache_dir=args.vep_result_cache_dir,
            vep_result_cache_tile_size_bp=args.vep_result_cache_tile_size_bp,
            performance_profile=performance,
        )

    print("Computing conservation-adjusted ClinVar validation...")
    with performance.stage("Conservation-adjusted validation"):
        conservation_analysis = build_conservation_analysis(
            inputs=inputs,
            validation=validation,
            strategies=strategies,
            eligible_gene_ids_by_strategy=alignment_gene_ids_by_strategy(cov),
            firth_workers=args.firth_workers,
            performance_profile=performance,
        )

    negative_controls = None
    if args.target_space_null:
        print("Computing consequence-matched target-space null...")
        with performance.stage("Target-space null"):
            negative_controls = build_target_space_null(
                run_dir=inputs.run_dir,
                variant_annotations_tsv=inputs.variant_annotations_tsv,
                target_features_tsv=inputs.target_features_tsv,
                genes_tsv=inputs.genes_tsv,
                target_sequences_dir=inputs.target_sequences_dir,
                clinvar_vcf=args.clinvar_vcf.expanduser().resolve(),
                strategies=strategies,
                observed_store=observed_store,
                sample_size_per_strategy=args.target_space_null_sample_size,
                resamples=args.target_space_null_resamples,
                seed=args.target_space_null_seed,
                gnomad_cache_dir=args.gnomad_cache_dir,
                vep_backend=args.vep_backend,
                vep_release=args.vep_release,
                vep_executable=args.vep_executable,
                vep_cache_dir=args.vep_cache_dir,
                vep_forks=args.vep_forks,
                vep_result_cache_dir=args.vep_result_cache_dir,
                vep_result_cache_tile_size_bp=args.vep_result_cache_tile_size_bp,
                performance_profile=performance,
            )
    else:
        performance.disabled_stage(
            "Target-space null",
            "Enable with --target-space-null",
        )

    with performance.stage("Report sections"):
        candidate_sections = build_variant_sections(variant_summary, strategy_stats)
        candidate_sections.extend(
            build_clinvar_gnomad_sections(variant_summary, strategy_stats_full)
        )
        candidate_sections.extend(build_feature_sections(cov))
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
                ),
            ),
            (
                "target-space-null",
                "Matched Control",
                build_target_space_null_sections(
                    negative_controls,
                    enabled=args.target_space_null,
                ),
            ),
            (
                "clinvar-association",
                "ClinVar Association",
                build_clinvar_association_sections(conservation_analysis),
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
                    performance.table_rows(),
                    taxonomy_summary,
                    performance_path,
                ),
            ),
        ]

    print(f"Writing report to {out_html}...")
    with performance.stage("HTML rendering"):
        html = render_html(sections)
    with performance.stage("HTML write"):
        write_text_atomic(out_html, html)
    performance.finish(artifacts=[out_html])
    print(f"Performance profile: {performance_path}")
    print(f"Done in {performance.total_wall_seconds:.3f} s")


if __name__ == "__main__":
    main()
