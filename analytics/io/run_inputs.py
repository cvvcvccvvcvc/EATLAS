"""Read-only pipeline inputs and explicit analytics workspace materialization."""

from __future__ import annotations

import hashlib
import json
import re
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import duckdb
import pandas as pd

from analytics.analyses.variant_summary_aggregation import (
    allele_evidence_comparison_sql,
    resolve_variant_aggregation_source,
)
from analytics.io.alignment_aggregates import resolve_alignment_aggregate_paths
from analytics.io.annotation_support import resolve_annotation_support_paths
from analytics.io.artifacts import content_identity, write_json_atomic
from analytics.io.evidence_verification import verify_source_evidence
from analytics.io.performance import PerformanceProfile, profile_stage
from analytics.io.taxonomy_summary import (
    build_or_load_taxonomy_summary_many,
    resolve_taxonomy_summary_path,
)
from analytics.io.variant_source import (
    resolve_variant_table_source,
    sql_string,
    variant_source_sql,
)
from genomics.gnomad import (
    GNOMAD_DATASET,
    merge_observation_windows,
    validate_observation_window,
)
from genomics.variants import ALLELE_ANNOTATION_FIELDS
from provenance.evidence_inventory import (
    INVENTORY_FILENAME,
    load_bound_evidence_inventory,
)


ANALYSIS_CONTRACT_VERSION = 2


@dataclass(frozen=True)
class SourceRun:
    """One completed pipeline run. Resolution never writes below ``run_dir``."""

    run_dir: Path
    source_id: str
    root_manifest_json: Path
    fetch_manifest_json: Path
    genes_tsv: Path
    target_features_tsv: Path
    target_sequences_dir: Path
    taxonomy_tsv: Path
    orthologs_selected_tsv: Path
    variant_annotations_source: Path
    annotation_manifest_json: Path
    annotation_failures_tsv: Path
    alignment_manifest_json: Path
    requested_gene_ids: frozenset[str]
    target_gene_ids: frozenset[str]
    root_manifest: dict[str, object]
    fetch_manifest: dict[str, object]
    alignment_manifest: dict[str, object]
    annotation_manifest: dict[str, object]
    variant_annotation_descriptor: dict[str, object]
    evidence_inventory: dict[str, object]
    evidence_inventory_descriptor: dict[str, object]


@dataclass(frozen=True)
class AnalysisInputs:
    """Resolved sources and derived paths for one analysis of one or more runs."""

    analytics_root: Path
    analysis_id: str
    analysis_dir: Path
    derived_dir: Path
    source_runs: tuple[SourceRun, ...]
    genes_tsvs: tuple[Path, ...]
    target_features_tsvs: tuple[Path, ...]
    target_sequence_dirs: tuple[Path, ...]
    variant_annotation_sources: tuple[Path, ...]
    variant_strategy_support_tsvs: tuple[Path, ...]
    ortholog_evidence_summary_tsvs: tuple[Path, ...]
    annotation_failure_tsvs: tuple[Path, ...]
    feature_coverage_tsvs: tuple[Path, ...]
    strategy_summary_tsvs: tuple[Path, ...]
    taxonomy_summary_tsv: Path
    fetch_manifest: dict[str, object]
    alignment_manifest: dict[str, object]
    annotation_manifest: dict[str, object]
    variant_annotation_descriptor: dict[str, object]
    analysis_manifest_json: Path

    @property
    def source_run_dirs(self) -> tuple[Path, ...]:
        return tuple(source.run_dir for source in self.source_runs)


@dataclass(frozen=True)
class AnalysisWorkspace:
    analytics_root: Path
    analysis_id: str
    analysis_dir: Path
    derived_dir: Path
    contract: dict[str, object]


def safe_report_name(name: str) -> str:
    value = name.strip()
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", value) is None:
        raise ValueError(
            "Report name may contain only letters, digits, dot, underscore, and hyphen"
        )
    return value


def resolve_analytics_root(
    analytics_root: Path,
    run_dirs: Sequence[Path],
) -> Path:
    """Resolve an analytics root and reject writes below immutable source runs."""

    root = analytics_root.expanduser().resolve()
    for run_dir in run_dirs:
        source = run_dir.expanduser().resolve()
        if root == source or source in root.parents:
            raise ValueError(
                "--analytics-root must be outside every immutable source run: "
                f"{source}"
            )
    return root


def resolve_source_runs(
    run_dirs: Sequence[Path],
    *,
    clinvar_vcf: Path,
) -> tuple[SourceRun, ...]:
    """Validate source runs without creating analytics artifacts."""

    if not run_dirs:
        raise ValueError("At least one --run-dir is required")
    sources = tuple(_resolve_source_run(path) for path in run_dirs)
    if len({source.run_dir for source in sources}) != len(sources):
        raise ValueError("The same --run-dir was supplied more than once")
    if len({source.source_id for source in sources}) != len(sources):
        raise ValueError("Source runs repeat the same completed pipeline identity")
    _require_disjoint_requested_genes(sources)
    _validate_compatibility(sources, clinvar_vcf.expanduser().resolve())
    return tuple(sorted(sources, key=lambda source: source.source_id))


def build_analysis_inputs(
    source_runs: Sequence[SourceRun],
    *,
    analytics_root: Path,
    scientific_config: dict[str, object],
    annotation_support_workers: int = 1,
    workspace: AnalysisWorkspace | None = None,
    performance_profile: PerformanceProfile | None = None,
) -> AnalysisInputs:
    """Materialize reusable analytics data outside immutable source runs."""

    if not source_runs:
        raise ValueError("Analysis requires at least one source run")
    source_runs = tuple(sorted(source_runs, key=lambda source: source.source_id))
    expected_workspace = resolve_analysis_workspace(
        source_runs,
        analytics_root=analytics_root,
        scientific_config=scientific_config,
    )
    if workspace is not None and workspace != expected_workspace:
        raise ValueError("Analysis workspace does not match the requested inputs")
    workspace = workspace or expected_workspace
    analytics_root = workspace.analytics_root
    analysis_id = workspace.analysis_id
    analysis_dir = workspace.analysis_dir
    derived_dir = workspace.derived_dir
    contract = workspace.contract
    for source in source_runs:
        with profile_stage(
            performance_profile, f"Verify source evidence [{source.run_dir.name}]"
        ) as timing:
            cached = verify_source_evidence(
                run_dir=source.run_dir,
                source_id=source.source_id,
                inventory=source.evidence_inventory,
                inventory_descriptor=source.evidence_inventory_descriptor,
                cache_dir=analytics_root / "cache" / source.source_id,
            )
            timing["details"] = "unchanged verified source" if cached else "full SHA-256 verification"
    derived_dir.mkdir(parents=True, exist_ok=True)
    _validate_shared_allele_evidence(source_runs, derived_dir)

    alignment_paths = []
    annotation_paths = []
    taxonomy_paths = []
    for source in source_runs:
        source_cache = analytics_root / "cache" / source.source_id
        with profile_stage(
            performance_profile,
            f"Alignment aggregates [{source.run_dir.name}]",
        ):
            alignment_path = resolve_alignment_aggregate_paths(
                source.run_dir,
                analytics_dir=source_cache,
            )
        alignment_paths.append(alignment_path)
        with profile_stage(
            performance_profile,
            f"Annotation support [{source.run_dir.name}]",
        ):
            annotation_path = resolve_annotation_support_paths(
                source.run_dir,
                analytics_dir=source_cache,
                workers=annotation_support_workers,
                progress_callback=(
                    performance_profile.checkpoint
                    if performance_profile is not None
                    else None
                ),
            )
        annotation_paths.append(annotation_path)
        with profile_stage(
            performance_profile,
            f"Taxonomy summary [{source.run_dir.name}]",
        ):
            taxonomy_path = resolve_taxonomy_summary_path(
                source.run_dir,
                analytics_dir=source_cache,
            )
        taxonomy_paths.append(taxonomy_path)

    taxonomy_summary = (
        taxonomy_paths[0]
        if len(source_runs) == 1
        else build_or_load_taxonomy_summary_many(
            taxonomy_tsvs=tuple(source.taxonomy_tsv for source in source_runs),
            orthologs_tsvs=tuple(
                source.orthologs_selected_tsv for source in source_runs
            ),
            analytics_dir=derived_dir,
        )
    )
    variant_descriptor = _combined_variant_descriptor(source_runs)
    inputs = AnalysisInputs(
        analytics_root=analytics_root,
        analysis_id=analysis_id,
        analysis_dir=analysis_dir,
        derived_dir=derived_dir,
        source_runs=source_runs,
        genes_tsvs=tuple(source.genes_tsv for source in source_runs),
        target_features_tsvs=tuple(source.target_features_tsv for source in source_runs),
        target_sequence_dirs=tuple(source.target_sequences_dir for source in source_runs),
        variant_annotation_sources=tuple(
            source.variant_annotations_source for source in source_runs
        ),
        variant_strategy_support_tsvs=tuple(
            value.variant_strategy_support_tsv for value in annotation_paths
        ),
        ortholog_evidence_summary_tsvs=tuple(
            value.ortholog_evidence_summary_tsv for value in annotation_paths
        ),
        annotation_failure_tsvs=tuple(
            source.annotation_failures_tsv for source in source_runs
        ),
        feature_coverage_tsvs=tuple(
            value.feature_coverage_tsv for value in alignment_paths
        ),
        strategy_summary_tsvs=tuple(
            value.strategy_summary_tsv for value in alignment_paths
        ),
        taxonomy_summary_tsv=taxonomy_summary,
        fetch_manifest=_combined_fetch_manifest(source_runs),
        alignment_manifest=_combined_alignment_manifest(source_runs),
        annotation_manifest=_combined_annotation_manifest(
            source_runs,
            variant_descriptor,
        ),
        variant_annotation_descriptor=variant_descriptor,
        analysis_manifest_json=analysis_dir / "manifest.json",
    )
    validate_report_inputs(inputs)
    write_json_atomic(
        inputs.analysis_manifest_json,
        {
            "schema_version": ANALYSIS_CONTRACT_VERSION,
            "status": "ready",
            "analysis_id": analysis_id,
            "contract": contract,
            "sources": [
                {
                    "source_id": source.source_id,
                    "run_dir": str(source.run_dir),
                    "requested_gene_count": len(source.requested_gene_ids),
                    "target_gene_count": len(source.target_gene_ids),
                    "evidence_inventory": source.evidence_inventory_descriptor,
                    "evidence_tree_sha256": source.evidence_inventory["tree_sha256"],
                }
                for source in source_runs
            ],
        },
    )
    return inputs


def resolve_analysis_workspace(
    source_runs: Sequence[SourceRun],
    *,
    analytics_root: Path,
    scientific_config: dict[str, object],
) -> AnalysisWorkspace:
    """Resolve the stable report location before expensive cache preparation."""

    if not source_runs:
        raise ValueError("Analysis requires at least one source run")
    source_runs = tuple(sorted(source_runs, key=lambda source: source.source_id))
    analytics_root = resolve_analytics_root(
        analytics_root,
        [source.run_dir for source in source_runs],
    )
    contract = {
        "contract_version": ANALYSIS_CONTRACT_VERSION,
        "source_ids": [source.source_id for source in source_runs],
        "scientific_config": scientific_config,
    }
    analysis_id = hashlib.sha256(_canonical_json(contract)).hexdigest()[:24]
    analysis_dir = analytics_root / "analyses" / analysis_id
    derived_dir = analysis_dir / "derived"
    return AnalysisWorkspace(
        analytics_root=analytics_root,
        analysis_id=analysis_id,
        analysis_dir=analysis_dir,
        derived_dir=derived_dir,
        contract=contract,
    )


def resolve_report_html(
    inputs: AnalysisInputs | AnalysisWorkspace,
    report_name: str,
) -> Path:
    name = safe_report_name(report_name)
    if not name.endswith(".html"):
        name += ".html"
    return inputs.analysis_dir / "reports" / name


def variant_annotation_descriptor(inputs: AnalysisInputs) -> dict[str, object]:
    return dict(inputs.variant_annotation_descriptor)


def variant_annotation_release(descriptor: dict[str, object]) -> str:
    config = descriptor["vep_config"]
    release = config.get("release") if isinstance(config, dict) else None
    if not release:
        raise ValueError("Variant annotation dataset is missing vep_config.release")
    return str(release)


def validate_report_inputs(inputs: AnalysisInputs) -> None:
    """Fail before report work when a current source contract is incompatible."""

    resolve_variant_aggregation_source(inputs.variant_annotation_sources)
    contracts: list[tuple[tuple[Path, ...], set[str]]] = [
        (
            inputs.genes_tsvs,
            {"gene_id", "chromosome", "begin", "end", "sequence_length"},
        ),
        (
            inputs.target_features_tsvs,
            {"gene_id", "feature_type", "target_start0", "target_end0"},
        ),
        (inputs.feature_coverage_tsvs, {"gene_id", "strategy", "feature_type"}),
        (
            inputs.strategy_summary_tsvs,
            {
                "strategy",
                "gene_count",
                "summary_row_count",
                "aligned_summary_row_count",
                "event_count",
            },
        ),
        (
            inputs.variant_strategy_support_tsvs,
            {
                "variant_key",
                "gene_id",
                "strategy",
                "alt_support_row_count",
                "alt_support_ortholog_count",
                "alt_support_family_count",
            },
        ),
        (
            inputs.ortholog_evidence_summary_tsvs,
            {
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
        ),
    ]
    for paths, required in contracts:
        for path in paths:
            _require_header(path, required)
    _require_header(
        inputs.taxonomy_summary_tsv,
        {
            "taxonomic_scope",
            "evidence_unit",
            "gene_count",
            "ortholog_count",
            "taxon_count",
            "unit_count",
            "orthologs_per_gene_median",
            "units_per_gene_median",
        },
    )


def resolve_variant_annotations_source(annotation_manifest_path: Path) -> Path:
    """Resolve one pipeline-owned partitioned variant-annotation dataset."""

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
    if _required_json(source) != descriptor:
        raise ValueError(
            "Annotation variant_annotations descriptor does not match its dataset "
            f"manifest: {source}"
        )
    resolve_variant_table_source(source, required_columns=set())
    return source


def _resolve_source_run(run_dir: Path) -> SourceRun:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory does not exist: {run_dir}")
    if not run_dir.is_dir():
        raise NotADirectoryError(f"--run-dir is not a directory: {run_dir}")
    root_path = run_dir / "run_manifest.json"
    fetch_path = run_dir / "fetch" / "manifest.json"
    alignment_path = run_dir / "alignment" / "manifest.json"
    annotation_path = run_dir / "annotation" / "manifest.json"
    root = _required_json(root_path)
    fetch = _required_json(fetch_path)
    alignment = _required_json(alignment_path)
    annotation = _required_json(annotation_path)
    _validate_completed_run(root, run_dir)
    inventory_path = run_dir / INVENTORY_FILENAME
    inventory, inventory_descriptor = load_bound_evidence_inventory(
        inventory_path, root.get("evidence_inventory")
    )
    variant_source = resolve_variant_annotations_source(annotation_path)
    genes = run_dir / "fetch" / "genes.tsv.gz"
    features = run_dir / "fetch" / "target_features.tsv.gz"
    targets = run_dir / "fetch" / "sequences" / "targets"
    taxonomy = run_dir / "fetch" / "taxonomy.tsv.gz"
    selected = run_dir / "fetch" / "orthologs.selected.tsv.gz"
    failures = run_dir / "annotation" / "failures.tsv.gz"
    missing = [
        str(path)
        for path in (genes, features, targets, taxonomy, selected, failures)
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError(
            "Incomplete source run required for analytics; missing: "
            + ", ".join(missing)
        )
    requested_gene_ids = frozenset(
        _read_requested_gene_ids(run_dir / "fetch" / "input.ids.tsv")
    )
    target_gene_ids = frozenset(_read_gene_ids(genes))
    if not target_gene_ids.issubset(requested_gene_ids):
        raise ValueError(f"Target genes fall outside accepted input IDs: {run_dir}")
    manifest_gene_ids = {str(value) for value in alignment.get("gene_ids", [])}
    if manifest_gene_ids != target_gene_ids:
        raise ValueError(
            f"Alignment manifest gene_ids do not match fetch/genes.tsv.gz: {run_dir}"
        )
    identity = {
        "contract_version": ANALYSIS_CONTRACT_VERSION,
        "root_manifest": content_identity(root_path),
        "evidence_tree_sha256": inventory["tree_sha256"],
    }
    source_id = hashlib.sha256(_canonical_json(identity)).hexdigest()[:24]
    return SourceRun(
        run_dir=run_dir,
        source_id=source_id,
        root_manifest_json=root_path,
        fetch_manifest_json=fetch_path,
        genes_tsv=genes,
        target_features_tsv=features,
        target_sequences_dir=targets,
        taxonomy_tsv=taxonomy,
        orthologs_selected_tsv=selected,
        variant_annotations_source=variant_source,
        annotation_manifest_json=annotation_path,
        annotation_failures_tsv=failures,
        alignment_manifest_json=alignment_path,
        requested_gene_ids=requested_gene_ids,
        target_gene_ids=target_gene_ids,
        root_manifest=root,
        fetch_manifest=fetch,
        alignment_manifest=alignment,
        annotation_manifest=annotation,
        variant_annotation_descriptor=_variant_annotation_descriptor(annotation_path),
        evidence_inventory=inventory,
        evidence_inventory_descriptor=inventory_descriptor,
    )


def _validate_completed_run(manifest: dict[str, object], run_dir: Path) -> None:
    if manifest.get("pipeline") != "gaph_v2":
        raise ValueError(f"Not a gaph_v2 run: {run_dir}")
    if manifest.get("schema_version") != 3:
        raise ValueError(f"Unsupported run manifest schema: {run_dir}")
    if manifest.get("status") != "complete" or manifest.get("success") is not True:
        raise ValueError(f"Run is not successfully complete: {run_dir}")
    if manifest.get("git_dirty") is not False:
        raise ValueError(f"Analytics requires a clean pipeline revision: {run_dir}")


def _require_disjoint_requested_genes(sources: Sequence[SourceRun]) -> None:
    owners: dict[str, list[str]] = {}
    for source in sources:
        for gene_id in source.requested_gene_ids:
            owners.setdefault(gene_id, []).append(source.run_dir.name)
    conflicts = {gene: labels for gene, labels in owners.items() if len(labels) > 1}
    if conflicts:
        examples = "; ".join(
            f"{gene}: {', '.join(labels)}"
            for gene, labels in sorted(
                conflicts.items(), key=lambda item: _gene_sort_key(item[0])
            )[:20]
        )
        raise ValueError(
            "Source runs contain overlapping accepted Gene IDs; implicit "
            f"deduplication would change scientific weights. Conflicts: {examples}"
        )


def _validate_compatibility(
    sources: Sequence[SourceRun],
    clinvar_vcf: Path,
) -> None:
    if not clinvar_vcf.is_file():
        raise FileNotFoundError(f"ClinVar VCF does not exist: {clinvar_vcf}")
    clinvar_tbi = Path(f"{clinvar_vcf}.tbi")
    if not clinvar_tbi.is_file():
        raise FileNotFoundError(f"ClinVar index does not exist: {clinvar_tbi}")
    expected_clinvar = {
        "vcf": content_identity(clinvar_vcf),
        "tbi": content_identity(clinvar_tbi),
    }
    contracts = [_source_compatibility(source) for source in sources]
    baseline = contracts[0]
    for source, observed in zip(sources[1:], contracts[1:]):
        differing = [key for key in baseline if observed.get(key) != baseline.get(key)]
        if differing:
            raise ValueError(
                f"Incompatible source run {source.run_dir}; contracts differ: "
                + ", ".join(differing)
            )
    missing = [
        key
        for key, value in baseline.items()
        if value is None or value == "" or value == [] or value == {}
    ]
    if missing:
        raise ValueError(
            "Source runs lack required current compatibility provenance: "
            + ", ".join(missing)
        )
    required = {
        "target_assembly_accession": "GCF_000001405.40",
        "target_assembly_name": "GRCh38.p14",
        "ortholog_scope": "all",
        "alignment_event_mode": "compact_support",
        "event_ortholog_support_format": "event_group_id_v2",
        "annotation_schema": "normalized_annotation_evidence_v5",
        "gnomad_dataset": GNOMAD_DATASET,
    }
    invalid = [
        key for key, value in required.items() if str(baseline.get(key)) != value
    ]
    if str(baseline.get("target_tax_id")) != "9606":
        invalid.append("target_tax_id")
    if invalid:
        raise ValueError(
            "Source runs do not use the current scientific constants: "
            + ", ".join(invalid)
        )
    for source in sources:
        if _declared_clinvar_identity(source) != expected_clinvar:
            raise ValueError(
                f"Run {source.run_dir} was annotated with different ClinVar contents"
            )


def _source_compatibility(source: SourceRun) -> dict[str, object]:
    descriptor = source.variant_annotation_descriptor
    config = descriptor.get("vep_config")
    if not isinstance(config, dict):
        config = {}
    return {
        "pipeline_git_commit": source.root_manifest.get("git_commit"),
        "target_assembly_accession": source.fetch_manifest.get("target_assembly_accession"),
        "target_assembly_name": source.fetch_manifest.get("target_assembly_name"),
        "target_tax_id": source.fetch_manifest.get("target_tax_id"),
        "target_annotation_gff3_sha256": source.fetch_manifest.get(
            "target_annotation_gff3_sha256"
        ),
        "ortholog_scope": source.fetch_manifest.get("ortholog_scope"),
        "datasets_versions": source.fetch_manifest.get("datasets_versions"),
        "strategies": source.alignment_manifest.get("strategies"),
        "strategy_parameters": source.alignment_manifest.get("strategy_parameters"),
        "alignment_event_mode": source.alignment_manifest.get("alignment_event_mode"),
        "event_ortholog_support_format": source.alignment_manifest.get(
            "event_ortholog_support_format"
        ),
        "annotation_schema": source.annotation_manifest.get("schema"),
        "gnomad_api_url": source.annotation_manifest.get("gnomad_api_url"),
        "gnomad_dataset": source.annotation_manifest.get("gnomad_dataset"),
        "vep_backend": config.get("backend"),
        "vep_release": config.get("release"),
        "vep_columns": descriptor.get("fields"),
    }


def _declared_clinvar_identity(source: SourceRun) -> dict[str, object]:
    identities = {}
    for key, manifest_key in (("vcf", "clinvar_vcf"), ("tbi", "clinvar_tbi")):
        metadata = source.annotation_manifest.get(manifest_key)
        if not isinstance(metadata, dict) or set(metadata) != {"size_bytes", "sha256"}:
            raise ValueError(
                f"Run {source.run_dir} annotation manifest has invalid "
                f"{manifest_key} content identity"
            )
        size_bytes = metadata.get("size_bytes")
        sha256 = metadata.get("sha256")
        if (
            isinstance(size_bytes, bool)
            or not isinstance(size_bytes, int)
            or size_bytes < 0
            or not isinstance(sha256, str)
            or re.fullmatch(r"[0-9a-f]{64}", sha256) is None
        ):
            raise ValueError(
                f"Run {source.run_dir} annotation manifest has invalid "
                f"{manifest_key} content identity"
            )
        identities[key] = {"size_bytes": size_bytes, "sha256": sha256}
    return identities


def _combined_variant_descriptor(sources: Sequence[SourceRun]) -> dict[str, object]:
    descriptors = [source.variant_annotation_descriptor for source in sources]
    baseline = descriptors[0]
    exact = ("schema", "status", "layout", "format", "fields", "vep_config")
    for source, descriptor in zip(sources[1:], descriptors[1:]):
        differing = [
            field for field in exact if descriptor.get(field) != baseline.get(field)
        ]
        if differing:
            raise ValueError(
                f"Incompatible variant annotations in {source.run_dir}: "
                + ", ".join(differing)
            )
    statuses: Counter[str] = Counter()
    for descriptor in descriptors:
        statuses.update(
            {
                str(key): int(value)
                for key, value in dict(descriptor.get("vep_status_counts", {})).items()
            }
        )
    return {
        **{field: baseline.get(field) for field in exact},
        "layout": (
            "multi_run_partitioned" if len(sources) > 1 else baseline.get("layout")
        ),
        "run_count": len(sources),
        "partition_count": sum(
            int(value.get("partition_count", 0)) for value in descriptors
        ),
        "shard_count": sum(int(value.get("shard_count", 0)) for value in descriptors),
        "row_count": sum(int(value.get("row_count", 0)) for value in descriptors),
        "vep_status_counts": dict(sorted(statuses.items())),
    }


def _validate_shared_allele_evidence(
    sources: Sequence[SourceRun],
    temporary_root: Path,
) -> None:
    if len(sources) < 2:
        return
    required = {
        "variant_key",
        "gene_id",
        "lookup_status",
        "gnomad_af",
        "clinvar_id",
        "clinvar_sig",
        "clinvar_review_stars",
        "clinvar_scv_count",
    }
    source = resolve_variant_table_source(
        tuple(item.variant_annotations_source for item in sources),
        required_columns=required,
    )
    checks = {
        field: allele_evidence_comparison_sql(field)
        for field in ALLELE_ANNOTATION_FIELDS
        if field in source.columns
    }
    count_sql = ", ".join(
        f"count(DISTINCT {expression}) FILTER "
        f"(WHERE {expression} IS NOT NULL) AS {field}_count"
        for field, expression in checks.items()
    )
    conflict_sql = " OR ".join(
        [
            *(f"{field}_count > 1" for field in checks),
            "gnomad_af_invalid_count > 0",
        ]
    )
    with tempfile.TemporaryDirectory(
        prefix=".shared_allele_preflight.",
        dir=temporary_root,
    ) as temporary:
        with duckdb.connect() as connection:
            connection.execute("SET preserve_insertion_order=false")
            connection.execute(f"SET temp_directory={sql_string(temporary)}")
            connection.execute(
                f"CREATE VIEW source_rows AS SELECT * FROM {variant_source_sql(source)}"
            )
            rows = connection.execute(
                "WITH evidence_counts AS (SELECT variant_key, "
                f"{count_sql}, count(*) FILTER (WHERE nullif(gnomad_af, '') "
                "IS NOT NULL AND try_cast(nullif(gnomad_af, '') AS DOUBLE) IS NULL) "
                "AS gnomad_af_invalid_count FROM source_rows WHERE variant_key <> '' "
                "GROUP BY variant_key) SELECT * FROM evidence_counts WHERE "
                f"{conflict_sql} LIMIT 10"
            ).fetchall()
    conflicts = [
        (field, str(row[0]))
        for row in rows
        for field, count in zip(checks, row[1:])
        if int(count) > 1
    ]
    if conflicts:
        examples = ", ".join(
            f"{variant} ({field})" for field, variant in conflicts[:10]
        )
        raise ValueError(
            "Source runs disagree on successful allele-level external evidence: "
            + examples
        )
    invalid_af = [str(row[0]) for row in rows if int(row[-1]) > 0]
    if invalid_af:
        raise ValueError(
            "Source runs contain non-numeric gnomad_af allele evidence: "
            + ", ".join(invalid_af[:10])
        )


def _combined_fetch_manifest(sources: Sequence[SourceRun]) -> dict[str, object]:
    if len(sources) == 1:
        return dict(sources[0].fetch_manifest)
    baseline = sources[0].fetch_manifest
    return {
        "stage": "fetch",
        "status": "complete",
        "run_count": len(sources),
        "unique_gene_count": sum(len(source.requested_gene_ids) for source in sources),
        "target_gene_count": sum(len(source.target_gene_ids) for source in sources),
        **{
            key: baseline.get(key)
            for key in (
                "target_assembly_accession",
                "target_assembly_name",
                "target_tax_id",
                "target_annotation_gff3_sha256",
                "ortholog_scope",
                "datasets_versions",
            )
        },
    }


def _combined_alignment_manifest(sources: Sequence[SourceRun]) -> dict[str, object]:
    if len(sources) == 1:
        return dict(sources[0].alignment_manifest)
    baseline = sources[0].alignment_manifest
    return {
        "stage": "alignment",
        "schema": baseline.get("schema"),
        "run_count": len(sources),
        "gene_ids": sorted(
            {gene for source in sources for gene in source.target_gene_ids},
            key=_gene_sort_key,
        ),
        "gene_count": sum(len(source.target_gene_ids) for source in sources),
        "strategies": baseline.get("strategies"),
        "strategy_parameters": baseline.get("strategy_parameters"),
        "alignment_event_mode": baseline.get("alignment_event_mode"),
        "event_ortholog_support_format": baseline.get("event_ortholog_support_format"),
    }


def _combined_annotation_manifest(
    sources: Sequence[SourceRun],
    descriptor: dict[str, object],
) -> dict[str, object]:
    if len(sources) == 1:
        return dict(sources[0].annotation_manifest)
    manifests = [source.annotation_manifest for source in sources]
    baseline = manifests[0]
    event_statuses: Counter[str] = Counter()
    for manifest in manifests:
        event_statuses.update(
            {
                str(key): int(value)
                for key, value in dict(
                    manifest.get("event_key_status_counts", {})
                ).items()
            }
        )
    return {
        "stage": "annotation",
        "schema": baseline.get("schema"),
        "run_count": len(sources),
        "variant_annotations": descriptor,
        "gnomad_api_url": baseline.get("gnomad_api_url"),
        "gnomad_dataset": baseline.get("gnomad_dataset"),
        "gnomad_observation_window": merge_observation_windows(
            manifest.get("gnomad_observation_window") for manifest in manifests
        ),
        "failure_count": sum(
            int(value.get("failure_count", 0) or 0) for value in manifests
        ),
        "gnomad_region_failure_count": sum(
            int(value.get("gnomad_region_failure_count", 0) or 0)
            for value in manifests
        ),
        "event_key_status_counts": dict(sorted(event_statuses.items())),
        "clinvar_cached_variant_count": sum(
            int(value.get("clinvar_cached_variant_count", 0) or 0)
            for value in manifests
        ),
        "gnomad_cached_variant_count": sum(
            int(value.get("gnomad_cached_variant_count", 0) or 0)
            for value in manifests
        ),
    }


def _variant_annotation_descriptor(path: Path) -> dict[str, object]:
    manifest = _required_json(path)
    if (
        manifest.get("stage") != "annotation"
        or manifest.get("schema") != "normalized_annotation_evidence_v5"
    ):
        raise ValueError(f"Unsupported pipeline annotation contract: {path}")
    if "gnomad_observation_window" not in manifest:
        raise ValueError(f"Annotation manifest lacks gnomAD observation provenance: {path}")
    observed = manifest["gnomad_observation_window"]
    if observed is not None:
        validate_observation_window(observed)
    descriptor = manifest.get("variant_annotations")
    if not isinstance(descriptor, dict):
        raise ValueError(f"Annotation manifest does not declare variant_annotations: {path}")
    config = descriptor.get("vep_config")
    counts = descriptor.get("vep_status_counts")
    if not isinstance(config, dict) or config.get("backend") not in {"rest", "local"}:
        raise ValueError("Variant annotation dataset has invalid VEP configuration")
    if not isinstance(counts, dict):
        raise ValueError("Variant annotation dataset has no VEP status counts")
    try:
        raw_row_count = descriptor["row_count"]
        if isinstance(raw_row_count, bool):
            raise TypeError
        row_count = int(raw_row_count)
        normalized = {}
        for key, value in counts.items():
            if not str(key) or isinstance(value, bool):
                raise TypeError
            normalized[str(key)] = int(value)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Variant annotation dataset has invalid VEP status counts") from exc
    if (
        row_count < 0
        or any(value < 0 for value in normalized.values())
        or sum(normalized.values()) != row_count
    ):
        raise ValueError("Variant annotation VEP status counts do not match row_count")
    variant_annotation_release(descriptor)
    return descriptor


def _required_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing JSON manifest: {path}")
    try:
        value = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON manifest: {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _read_requested_gene_ids(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing normalized input IDs: {path}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"gene_id", "accepted"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Normalized input IDs have an incompatible header: {path}")
    return frame.loc[frame["accepted"].eq("true"), "gene_id"].astype(str).tolist()


def _read_gene_ids(path: Path) -> list[str]:
    frame = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=["gene_id"],
        dtype=str,
    )
    values = frame["gene_id"].astype(str).tolist()
    if len(values) != len(set(values)):
        raise ValueError(f"Run contains duplicate target Gene IDs: {path}")
    return values


def _require_header(path: Path, required: set[str]) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing report input: {path}")
    compression = "gzip" if path.suffix == ".gz" else None
    header = set(
        pd.read_csv(path, sep="\t", compression=compression, nrows=0).columns
    )
    missing = required - header
    if missing:
        raise ValueError(
            f"Report input {path} is missing required columns: "
            + ", ".join(sorted(missing))
        )


def _canonical_json(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def _gene_sort_key(value: str) -> tuple[int, int | str]:
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)


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


def read_feature_coverage(paths: Sequence[Path]) -> pd.DataFrame:
    frames = [
        pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
        for path in paths
    ]
    cov = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    for col in (
        "length_bp",
        "ortholog_count",
        "orthologs_covered",
        "covered_bases",
        "coverage_breadth",
        "depth_bases",
        "mean_depth",
    ):
        if col in cov.columns:
            cov[col] = pd.to_numeric(cov[col], errors="coerce")
    return cov


def read_strategy_summary(paths: Sequence[Path]) -> pd.DataFrame:
    frames = [
        pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
        for path in paths
    ]
    if not frames:
        raise ValueError("Strategy summary requires at least one source table")
    columns = list(frames[0].columns)
    if any(list(frame.columns) != columns for frame in frames[1:]):
        raise ValueError("Strategy summary columns differ across source runs")
    combined = pd.concat(frames, ignore_index=True)
    numeric = [column for column in columns if column != "strategy"]
    for column in numeric:
        combined[column] = pd.to_numeric(combined[column], errors="raise")
    return combined.groupby("strategy", as_index=False, sort=True)[numeric].sum()[columns]


def read_input_gene_count(paths: Sequence[Path]) -> int:
    values: set[str] = set()
    for path in paths:
        frame = pd.read_csv(
            path,
            sep="\t",
            compression="gzip",
            usecols=["gene_id"],
            dtype=str,
        )
        values.update(frame["gene_id"].astype(str))
    return len(values)


def read_failures(paths: Sequence[Path]) -> pd.DataFrame:
    frames = [
        pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
        for path in paths
    ]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def read_taxonomy_summary(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing taxonomy summary: {path}")
    summary = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
    for column in (
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
    ):
        if column in summary:
            summary[column] = pd.to_numeric(summary[column], errors="raise")
    return summary
