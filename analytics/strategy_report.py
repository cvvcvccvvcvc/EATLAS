#!/usr/bin/env python3
"""Build one HTML report from one or more compatible completed GAPH runs."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from analytics.analyses.clinvar_validation import build_validation
from analytics.analyses.basic_filtering import build_basic_filtering_analysis
from analytics.analyses.conservation_analysis import (
    alignment_gene_ids_by_strategy,
    build_conservation_analysis,
)
from analytics.analyses.conservation_validation import validate_firth_runtime
from analytics.analyses.matched_control import build_target_space_null
from analytics.analyses.observed_variant_store import (
    build_or_load_observed_variant_store,
)
from analytics.analyses.pathogenic_variants import build_pathogenic_variant_analysis
from analytics.analyses.variant_summary import build_variant_summary
from analytics.io.run_inputs import (
    build_analysis_inputs,
    read_failures,
    read_feature_coverage,
    read_input_gene_count,
    read_strategy_summary,
    read_taxonomy_summary,
    resolve_analysis_workspace,
    resolve_report_html,
    resolve_source_runs,
    variant_annotation_descriptor,
    variant_annotation_release,
)
from analytics.io.artifacts import path_metadata, write_text_atomic
from analytics.io.performance import PerformanceProfile
from analytics.reporting.components import strategy_label
from analytics.reporting.basic_filtering import build_basic_filtering_sections
from analytics.reporting.conservation import build_clinvar_association_sections
from analytics.reporting.document import render_html
from analytics.reporting.matched_control import build_target_space_null_sections
from analytics.reporting.ortholog_evidence import build_ortholog_evidence_sections
from analytics.reporting.overview import build_overview, merge_alignment_summary
from analytics.reporting.pathogenic_variants import build_pathogenic_variant_sections
from analytics.reporting.qc import build_methods_sections
from analytics.reporting.variant_profile import (
    build_clinvar_gnomad_sections,
    build_feature_sections,
    build_gnomad_stratification_sections,
    build_variant_sections,
)
from genomics.vep.local_runtime import probe_local_vep


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analytics-root",
        type=Path,
        required=True,
        help="External workspace for analytics caches, derived data, and reports.",
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        action="append",
        required=True,
        help="Completed GAPH source run. Repeat for a multi-run analysis.",
    )
    parser.add_argument(
        "--clinvar-vcf",
        type=Path,
        default=_default_clinvar_vcf(),
        help=(
            "Indexed ClinVar VCF used for validation. Default: "
            "$CLINVAR_VCF, then assets/reference/clinvar/clinvar.vcf.gz"
        ),
    )
    parser.add_argument(
        "--report-name",
        required=True,
        help="Report file name inside the analysis reports directory.",
    )
    parser.add_argument(
        "--target-space-null",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Build the consequence-matched target-space null. Disabled by default "
            "because it uses Ensembl VEP and the gnomAD GraphQL API."
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
        "--phylop-bigwig",
        type=Path,
        default=_default_phylop_bigwig(),
        help=(
            "Optional local hg38 phyloP100way BigWig. Defaults to "
            "$GAPH_PHYLOP_BIGWIG or an existing file under $GAPH_ROOT/reference/ucsc."
        ),
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
        help=(
            "Pinned Ensembl VEP release. Local VEP requires it; REST verifies "
            "that the server reports the same release."
        ),
    )
    parser.add_argument(
        "--vep-executable",
        default=os.environ.get("GAPH_VEP_EXECUTABLE", "vep"),
        help="Local VEP executable or wrapper. Used with --vep-backend local.",
    )
    parser.add_argument(
        "--vep-cache-dir",
        type=_non_empty_path,
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
        "--annotation-support-workers",
        type=int,
        default=_default_annotation_support_workers(),
        help=(
            "Worker processes for partition-local annotation support. "
            "Defaults to 4 within a Slurm allocation and 1 locally."
        ),
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


def _non_empty_path(raw: str) -> Path:
    if not raw.strip():
        raise argparse.ArgumentTypeError("path must not be empty")
    return Path(raw)


def _default_clinvar_vcf() -> Path:
    configured = os.environ.get("CLINVAR_VCF")
    if configured:
        return Path(configured)
    return project_root() / "assets" / "reference" / "clinvar" / "clinvar.vcf.gz"


def _default_vep_result_cache_dir() -> Path | None:
    configured = os.environ.get("GAPH_VEP_RESULT_CACHE_DIR")
    if configured:
        return Path(configured)
    gaph_root = os.environ.get("GAPH_ROOT")
    return Path(gaph_root) / "cache" / "vep_results" if gaph_root else None


def _default_phylop_bigwig() -> Path | None:
    configured = os.environ.get("GAPH_PHYLOP_BIGWIG")
    if configured:
        return Path(configured)
    gaph_root = os.environ.get("GAPH_ROOT")
    if not gaph_root:
        return None
    candidate = Path(gaph_root) / "reference" / "ucsc" / "hg38.phyloP100way.bw"
    return candidate if candidate.is_file() else None


def _default_firth_workers() -> int:
    configured = os.environ.get("GAPH_FIRTH_WORKERS")
    if configured:
        return int(configured)
    allocated = os.environ.get("SLURM_CPUS_PER_TASK")
    available = int(allocated) if allocated and allocated.isdigit() else (os.cpu_count() or 1)
    return max(1, min(available, 8))


def _default_annotation_support_workers() -> int:
    configured = os.environ.get("GAPH_ANNOTATION_SUPPORT_WORKERS")
    if configured:
        return int(configured)
    allocated = os.environ.get("SLURM_CPUS_PER_TASK")
    if allocated and allocated.isdigit():
        return max(1, min(int(allocated), 4))
    return 1


def _require_local_vep_runtime(
    *,
    backend: str,
    release: str,
    executable: str,
    cache_dir: Path | None,
) -> None:
    if backend != "local":
        return
    errors: list[str] = []
    probe_local_vep(
        release=release,
        executable=executable,
        cache_dir=cache_dir,
        errors=errors,
    )
    if errors:
        raise RuntimeError("Local VEP runtime check failed:\n- " + "\n- ".join(errors))


def main() -> None:
    args = parse_args()
    args.clinvar_vcf = args.clinvar_vcf.expanduser().resolve()
    if args.phylop_bigwig is not None:
        args.phylop_bigwig = args.phylop_bigwig.expanduser().resolve()
        if not args.phylop_bigwig.is_file():
            raise FileNotFoundError(
                f"Local phyloP BigWig does not exist: {args.phylop_bigwig}"
            )
    if args.target_space_null and args.target_space_null_sample_size < 1:
        raise ValueError("--target-space-null-sample-size must be >= 1")
    if args.target_space_null and args.target_space_null_resamples < 100:
        raise ValueError("--target-space-null-resamples must be >= 100")
    if args.vep_forks < 1:
        raise ValueError("--vep-forks must be >= 1")
    if args.annotation_support_workers < 1:
        raise ValueError("--annotation-support-workers must be >= 1")
    if args.firth_workers < 1:
        raise ValueError("--firth-workers must be >= 1")
    if args.vep_result_cache_tile_size_bp < 1:
        raise ValueError("--vep-result-cache-tile-size-bp must be >= 1")
    source_runs = resolve_source_runs(args.run_dir, clinvar_vcf=args.clinvar_vcf)
    artifact_release = variant_annotation_release(
        source_runs[0].variant_annotation_descriptor
    )
    if args.vep_release and str(args.vep_release) != artifact_release:
        raise ValueError(
            "Pipeline variant annotations use VEP release "
            f"{artifact_release}, not {args.vep_release}"
        )
    if not args.vep_release:
        args.vep_release = artifact_release
    _require_local_vep_runtime(
        backend=args.vep_backend,
        release=str(args.vep_release),
        executable=args.vep_executable,
        cache_dir=args.vep_cache_dir,
    )
    scientific_config = {
        "clinvar": {
            "vcf": path_metadata(args.clinvar_vcf),
            "tbi": path_metadata(Path(f"{args.clinvar_vcf}.tbi")),
        },
        "phylop": path_metadata(args.phylop_bigwig)
        if args.phylop_bigwig is not None
        else None,
        "target_space_null": (
            {
                "enabled": True,
                "sample_size": args.target_space_null_sample_size,
                "resamples": args.target_space_null_resamples,
                "seed": args.target_space_null_seed,
            }
            if args.target_space_null
            else {"enabled": False}
        ),
        "vep": {"backend": args.vep_backend, "release": args.vep_release},
    }
    workspace = resolve_analysis_workspace(
        source_runs,
        analytics_root=args.analytics_root,
        scientific_config=scientific_config,
    )
    out_html = resolve_report_html(workspace, args.report_name)
    performance_path = workspace.analysis_dir / "performance" / f"{out_html.stem}.json"
    performance = PerformanceProfile(
        performance_path,
        analysis_dir=workspace.analysis_dir,
        analysis_id=workspace.analysis_id,
        report_path=out_html,
        tracked_directory=workspace.analysis_dir,
        source_run_dirs=tuple(source.run_dir for source in source_runs),
    )
    with performance.stage("Prepare analysis inputs") as timing:
        inputs = build_analysis_inputs(
            source_runs,
            analytics_root=args.analytics_root,
            scientific_config=scientific_config,
            annotation_support_workers=args.annotation_support_workers,
            workspace=workspace,
            performance_profile=performance,
        )
        timing["metrics"] = {
            "source_runs": len(inputs.source_runs),
            "annotation_support_workers": args.annotation_support_workers,
        }
    candidate_annotation_descriptor = variant_annotation_descriptor(inputs)
    analytics_dir = inputs.derived_dir

    print(
        f"Analysis {inputs.analysis_id}: {len(inputs.source_runs)} source run(s), "
        f"workspace {inputs.analysis_dir}"
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
    print(
        "Streaming variant annotations from "
        f"{len(inputs.variant_annotation_sources)} source run(s)..."
    )
    with performance.stage("Variant summary") as timing:
        variant_summary = build_variant_summary(
            inputs.variant_annotation_sources,
            analytics_dir,
            strategy_label,
            target_features_path=inputs.target_features_tsvs,
            genes_path=inputs.genes_tsvs,
            annotation_failures_path=inputs.annotation_failure_tsvs,
            variant_strategy_support_path=inputs.variant_strategy_support_tsvs,
            ortholog_evidence_summary_path=inputs.ortholog_evidence_summary_tsvs,
            performance_profile=performance,
        )
        timing["details"] = "cache hit" if variant_summary.cache_hit else "cache miss"

    with performance.stage("Run summary inputs"):
        cov = read_feature_coverage(inputs.feature_coverage_tsvs)
        alignment_summary = read_strategy_summary(inputs.strategy_summary_tsvs)
        fetch_manifest = inputs.fetch_manifest
        input_gene_count = int(
            fetch_manifest.get("unique_gene_count")
            or read_input_gene_count(inputs.genes_tsvs)
        )
        failures = read_failures(inputs.annotation_failure_tsvs)
        annotation_manifest = inputs.annotation_manifest
        alignment_manifest = inputs.alignment_manifest
        taxonomy_summary = read_taxonomy_summary(inputs.taxonomy_summary_tsv)

    print("Computing strategy metrics...")
    with performance.stage("Strategy metrics"):
        strategy_stats_full = merge_alignment_summary(
            variant_summary.strategy_stats,
            alignment_summary,
        )
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
            variant_annotations_source=inputs.variant_annotation_sources,
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
            analytics_dir=analytics_dir,
            genes_tsv=inputs.genes_tsvs,
            target_sequences_dir=inputs.target_sequence_dirs,
            clinvar_vcf=args.clinvar_vcf,
            strategies=strategies,
            observed_store=observed_store,
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
    eligible_gene_ids_by_strategy = alignment_gene_ids_by_strategy(cov)
    with performance.stage("Conservation-adjusted validation"):
        conservation_analysis = build_conservation_analysis(
            inputs=inputs,
            validation=validation,
            strategies=strategies,
            eligible_gene_ids_by_strategy=eligible_gene_ids_by_strategy,
            firth_workers=args.firth_workers,
            performance_profile=performance,
            phylop_bigwig=args.phylop_bigwig,
        )

    print("Characterizing pathogenic ClinVar hits...")
    with performance.stage("Pathogenic ClinVar hits") as timing:
        pathogenic_analysis = build_pathogenic_variant_analysis(
            summary=variant_summary,
            clinvar_universe=validation.universe,
            clinvar_vcf=args.clinvar_vcf,
            condition_cache_dir=args.analytics_root / "cache" / "clinvar_conditions",
            eligible_gene_ids_by_strategy=eligible_gene_ids_by_strategy,
            analytics_dir=analytics_dir,
        )
        timing["metrics"] = {
            "variants": int(len(pathogenic_analysis.variants)),
            "conditions": int(
                pathogenic_analysis.condition_counts["condition"].nunique()
                if not pathogenic_analysis.condition_counts.empty
                else 0
            ),
            "support_rows": int(len(pathogenic_analysis.support_rows)),
        }

    print("Computing basic support-filter curves...")
    with performance.stage("Basic filtering") as timing:
        basic_filtering = build_basic_filtering_analysis(
            variant_annotations_source=inputs.variant_annotation_sources,
            variant_strategy_support_tsv=inputs.variant_strategy_support_tsvs,
            annotation_failures_tsv=inputs.annotation_failure_tsvs,
            analytics_dir=analytics_dir,
            cohort=conservation_analysis.validation.cohort,
            strategies=strategies,
            eligible_gene_ids_by_strategy=eligible_gene_ids_by_strategy,
        )
        timing["details"] = "cache hit" if basic_filtering.cache_hit else "cache miss"
        timing["metrics"] = {
            "candidate_curve_rows": int(len(basic_filtering.candidate_curves)),
            "clinvar_curve_rows": int(len(basic_filtering.clinvar_curves)),
        }

    negative_controls = None
    if args.target_space_null:
        print("Computing consequence-matched target-space null...")
        with performance.stage("Target-space null"):
            negative_controls = build_target_space_null(
                analytics_dir=analytics_dir,
                variant_annotations_source=inputs.variant_annotation_sources,
                target_features_tsv=inputs.target_features_tsvs,
                genes_tsv=inputs.genes_tsvs,
                target_sequences_dir=inputs.target_sequence_dirs,
                clinvar_vcf=args.clinvar_vcf,
                strategies=strategies,
                observed_store=observed_store,
                sample_size_per_strategy=args.target_space_null_sample_size,
                resamples=args.target_space_null_resamples,
                seed=args.target_space_null_seed,
                gnomad_cache_dir=args.gnomad_cache_dir,
                phylop_bigwig=args.phylop_bigwig,
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
            (
                "pathogenic-clinvar-hits",
                "Pathogenic ClinVar Hits",
                build_pathogenic_variant_sections(pathogenic_analysis),
            ),
            ("candidates", "Candidate Profile", candidate_sections),
            (
                "basic-filtering",
                "Basic Filtering",
                build_basic_filtering_sections(basic_filtering),
            ),
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
                    taxonomy_summary,
                    candidate_annotation_descriptor,
                    validation=validation,
                    conservation_analysis=conservation_analysis,
                    negative_controls=negative_controls,
                    report_timings=performance.table_rows(),
                    report_profile_path=performance_path,
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
