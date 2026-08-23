"""Resolve several compatible completed runs as one report input cohort."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import duckdb
import pandas as pd

from analytics.analyses.variant_summary_aggregation import (
    resolve_variant_aggregation_source,
)
from analytics.io.artifacts import (
    content_identity,
    path_metadata,
    write_text_atomic,
)
from analytics.io.run_inputs import (
    RunInputs,
    bulk_vep_manifest,
    read_json,
    resolve_run_inputs,
    validate_report_inputs,
)
from analytics.io.variant_source import (
    COHORT_VARIANT_SOURCE_KIND,
    COHORT_VARIANT_SOURCE_SCHEMA_VERSION,
    resolve_pre_vep_variant_source,
    sql_string,
    variant_source_sql,
)
from analytics.derivations.taxonomy import (
    TAXONOMY_SUMMARY_FIELDS,
    build_taxonomy_summary_rows,
    load_taxonomy_profiles,
)


COHORT_MANIFEST_SCHEMA_VERSION = 1
COHORT_CONTRACT_VERSION = 1
COHORT_RESOLVED_MANIFEST = "cohort.resolved.json"


@dataclass(frozen=True)
class CohortMember:
    label: str
    inputs: RunInputs
    requested_gene_ids: frozenset[str]
    target_gene_ids: frozenset[str]
    fingerprint: str
    scientific_files: tuple[dict[str, object], ...]
    root_manifest: dict[str, object]
    fetch_manifest: dict[str, object]
    alignment_manifest: dict[str, object]
    annotation_manifest: dict[str, object]
    vep_manifest: dict[str, object]


def resolve_cohort_inputs(
    manifest_path: Path,
    *,
    cohort_root: Path | None,
    clinvar_vcf: Path,
) -> RunInputs:
    """Validate, fingerprint, and expose a run set through the single-run contract."""

    manifest_path = manifest_path.expanduser().resolve()
    manifest = _read_cohort_manifest(manifest_path)
    declared = _declared_members(manifest, manifest_path.parent)
    members = [_load_member(label, run_dir) for label, run_dir in declared]
    _require_distinct_runs(members)
    _require_disjoint_requested_genes(members)
    compatibility = _validate_compatibility(members, clinvar_vcf.expanduser().resolve())

    members = sorted(members, key=lambda member: member.fingerprint)
    cohort_contract = {
        "contract_version": COHORT_CONTRACT_VERSION,
        "member_fingerprints": [member.fingerprint for member in members],
        "gene_overlap_policy": "reject",
    }
    cohort_id = hashlib.sha256(_canonical_json(cohort_contract)).hexdigest()[:24]
    root = (
        cohort_root.expanduser().resolve()
        if cohort_root is not None
        else manifest_path.parent / "cohorts"
    )
    cohort_dir = root / cohort_id
    inputs_dir = cohort_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)

    variant_descriptor = inputs_dir / "variant_annotations.cohort.json"
    _write_json_stable(
        variant_descriptor,
        {
            "schema_version": COHORT_VARIANT_SOURCE_SCHEMA_VERSION,
            "kind": COHORT_VARIANT_SOURCE_KIND,
            "members": [
                {
                    "fingerprint": member.fingerprint,
                    "path": str(member.inputs.variant_annotations_tsv),
                }
                for member in members
            ],
        },
    )
    _validate_shared_allele_evidence(variant_descriptor)

    genes_tsv = inputs_dir / "genes.tsv.gz"
    target_features_tsv = inputs_dir / "target_features.tsv.gz"
    coverage_tsv = inputs_dir / "feature_coverage.tsv.gz"
    support_tsv = inputs_dir / "variant_strategy_support.tsv.gz"
    ortholog_evidence_tsv = inputs_dir / "ortholog_evidence_summary.tsv.gz"
    failures_tsv = inputs_dir / "annotation_failures.tsv.gz"
    strategy_summary_tsv = inputs_dir / "strategy_summary.tsv.gz"
    taxonomy_summary_tsv = inputs_dir / "taxonomy_summary.tsv.gz"
    target_sequences_dir = inputs_dir / "target_sequences"
    _concatenate_tsv([member.inputs.genes_tsv for member in members], genes_tsv)
    _concatenate_tsv(
        [member.inputs.target_features_tsv for member in members],
        target_features_tsv,
    )
    _concatenate_tsv(
        [member.inputs.feature_coverage_tsv for member in members],
        coverage_tsv,
    )
    _concatenate_tsv(
        [member.inputs.variant_strategy_support_tsv for member in members],
        support_tsv,
    )
    _concatenate_tsv(
        [member.inputs.ortholog_evidence_summary_tsv for member in members],
        ortholog_evidence_tsv,
    )
    _concatenate_tsv(
        [member.inputs.annotation_failures_tsv for member in members],
        failures_tsv,
    )
    _aggregate_strategy_summaries(
        [member.inputs.strategy_summary_tsv for member in members],
        strategy_summary_tsv,
    )
    _aggregate_taxonomy_summary(members, taxonomy_summary_tsv)
    _link_target_sequences(members, target_sequences_dir)

    fetch_manifest_path = inputs_dir / "fetch_manifest.json"
    alignment_manifest_path = inputs_dir / "alignment_manifest.json"
    annotation_manifest_path = inputs_dir / "annotation_manifest.json"
    vep_manifest_path = cohort_dir / "analytics" / "vep_consequences" / "manifest.json"
    _write_json_stable(fetch_manifest_path, _cohort_fetch_manifest(members, compatibility))
    _write_json_stable(
        alignment_manifest_path,
        _cohort_alignment_manifest(members, compatibility),
    )
    _write_json_stable(
        annotation_manifest_path,
        _cohort_annotation_manifest(members, compatibility),
    )
    _write_json_stable(vep_manifest_path, _cohort_vep_manifest(members))

    resolved_manifest_path = cohort_dir / COHORT_RESOLVED_MANIFEST
    resolved_manifest = {
        "schema_version": COHORT_MANIFEST_SCHEMA_VERSION,
        "status": "ready",
        "cohort_id": cohort_id,
        "source_manifest": str(manifest_path),
        "contract": cohort_contract,
        "compatibility": compatibility,
        "members": [
            {
                "label": member.label,
                "run_dir": str(member.inputs.run_dir),
                "fingerprint": member.fingerprint,
                "requested_gene_count": len(member.requested_gene_ids),
                "target_gene_count": len(member.target_gene_ids),
                "scientific_files": list(member.scientific_files),
            }
            for member in members
        ],
        "limitations": {
            "gene_overlap": (
                "Rejected because current durable run aggregates cannot be subset by gene."
            ),
        },
    }
    _write_json_stable(resolved_manifest_path, resolved_manifest)

    inputs = RunInputs(
        run_dir=cohort_dir,
        fetch_manifest_json=fetch_manifest_path,
        genes_tsv=genes_tsv,
        target_features_tsv=target_features_tsv,
        target_sequences_dir=target_sequences_dir,
        variant_annotations_tsv=variant_descriptor,
        variant_strategy_support_tsv=support_tsv,
        ortholog_evidence_summary_tsv=ortholog_evidence_tsv,
        annotation_manifest_json=annotation_manifest_path,
        annotation_failures_tsv=failures_tsv,
        feature_coverage_tsv=coverage_tsv,
        alignment_manifest_json=alignment_manifest_path,
        strategy_summary_tsv=strategy_summary_tsv,
        taxonomy_summary_tsv=taxonomy_summary_tsv,
        cohort_manifest_json=resolved_manifest_path,
        source_run_dirs=tuple(member.inputs.run_dir for member in members),
        cohort_id=cohort_id,
    )
    validate_report_inputs(inputs)
    return inputs


def _read_cohort_manifest(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(f"Cohort manifest does not exist: {path}")
    try:
        payload = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid cohort manifest JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("Cohort manifest must be a JSON object")
    if payload.get("schema_version") != COHORT_MANIFEST_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported cohort manifest schema_version: "
            f"{payload.get('schema_version')!r}"
        )
    unknown = set(payload) - {"schema_version", "runs"}
    if unknown:
        raise ValueError(
            "Cohort manifest has unknown fields: " + ", ".join(sorted(unknown))
        )
    return payload


def _declared_members(
    manifest: dict[str, object],
    manifest_dir: Path,
) -> list[tuple[str, Path]]:
    raw_runs = manifest.get("runs")
    if not isinstance(raw_runs, list) or not raw_runs:
        raise ValueError("Cohort manifest must contain a non-empty runs array")
    members = []
    labels = set()
    for index, raw in enumerate(raw_runs, start=1):
        if not isinstance(raw, dict):
            raise ValueError(f"Cohort run #{index} must be a JSON object")
        unknown = set(raw) - {"label", "run_dir"}
        if unknown:
            raise ValueError(
                f"Cohort run #{index} has unknown fields: {', '.join(sorted(unknown))}"
            )
        raw_path = str(raw.get("run_dir") or "").strip()
        if not raw_path:
            raise ValueError(f"Cohort run #{index} is missing run_dir")
        path = Path(raw_path).expanduser()
        if not path.is_absolute():
            path = manifest_dir / path
        path = path.resolve()
        label = str(raw.get("label") or path.name).strip()
        if not label or label in labels:
            raise ValueError(f"Cohort run labels must be non-empty and unique: {label!r}")
        labels.add(label)
        members.append((label, path))
    return members


def _load_member(label: str, run_dir: Path) -> CohortMember:
    inputs = resolve_run_inputs(run_dir)
    validate_report_inputs(inputs)
    root_manifest = _required_json(run_dir / "run_manifest.json")
    fetch_manifest = _required_json(inputs.fetch_manifest_json)
    alignment_manifest = _required_json(inputs.alignment_manifest_json)
    annotation_manifest = _required_json(inputs.annotation_manifest_json)
    vep_manifest = bulk_vep_manifest(inputs)
    _validate_completed_run(root_manifest, run_dir)
    requested_gene_ids = frozenset(_read_requested_gene_ids(run_dir / "fetch" / "input.ids.tsv"))
    target_gene_ids = frozenset(_read_gene_ids(inputs.genes_tsv))
    if not target_gene_ids.issubset(requested_gene_ids):
        raise ValueError(f"Target genes fall outside accepted input IDs: {run_dir}")
    manifest_gene_ids = {str(value) for value in alignment_manifest.get("gene_ids", [])}
    if manifest_gene_ids != target_gene_ids:
        raise ValueError(
            f"Alignment manifest gene_ids do not match fetch/genes.tsv.gz: {run_dir}"
        )
    scientific_files = tuple(_scientific_file_records(inputs, run_dir))
    fingerprint_payload = {
        "contract_version": COHORT_CONTRACT_VERSION,
        "files": scientific_files,
    }
    fingerprint = hashlib.sha256(_canonical_json(fingerprint_payload)).hexdigest()
    return CohortMember(
        label=label,
        inputs=inputs,
        requested_gene_ids=requested_gene_ids,
        target_gene_ids=target_gene_ids,
        fingerprint=fingerprint,
        scientific_files=scientific_files,
        root_manifest=root_manifest,
        fetch_manifest=fetch_manifest,
        alignment_manifest=alignment_manifest,
        annotation_manifest=annotation_manifest,
        vep_manifest=vep_manifest,
    )


def _validate_completed_run(manifest: dict[str, object], run_dir: Path) -> None:
    if manifest.get("pipeline") != "gaph_v2":
        raise ValueError(f"Not a gaph_v2 run: {run_dir}")
    if manifest.get("schema_version") != 2:
        raise ValueError(f"Unsupported run manifest schema: {run_dir}")
    if manifest.get("status") != "complete" or manifest.get("success") is not True:
        raise ValueError(f"Run is not successfully complete: {run_dir}")
    if manifest.get("git_dirty") is not False:
        raise ValueError(f"Cohort reports require a clean pipeline revision: {run_dir}")


def _require_distinct_runs(members: list[CohortMember]) -> None:
    paths = [member.inputs.run_dir for member in members]
    if len(set(paths)) != len(paths):
        raise ValueError("Cohort manifest repeats the same run directory")
    fingerprints = [member.fingerprint for member in members]
    if len(set(fingerprints)) != len(fingerprints):
        raise ValueError("Cohort manifest repeats the same scientific run contents")


def _require_disjoint_requested_genes(members: list[CohortMember]) -> None:
    owners: dict[str, list[str]] = {}
    for member in members:
        for gene_id in member.requested_gene_ids:
            owners.setdefault(gene_id, []).append(member.label)
    conflicts = {gene: labels for gene, labels in owners.items() if len(labels) > 1}
    if conflicts:
        ordered = sorted(
            conflicts.items(), key=lambda item: _gene_sort_key(item[0])
        )
        examples = "; ".join(
            f"{gene}: {', '.join(labels)}" for gene, labels in ordered[:20]
        )
        suffix = "" if len(conflicts) <= 20 else f"; and {len(conflicts) - 20} more"
        raise ValueError(
            "Cohort runs contain overlapping accepted gene IDs. Current durable run "
            "aggregates cannot select one copy without bias; choose non-overlapping runs. "
            f"Conflicts: {examples}{suffix}"
        )


def _validate_compatibility(
    members: list[CohortMember],
    clinvar_vcf: Path,
) -> dict[str, object]:
    if not clinvar_vcf.is_file():
        raise FileNotFoundError(f"ClinVar VCF does not exist: {clinvar_vcf}")
    clinvar_tbi = Path(f"{clinvar_vcf}.tbi")
    if not clinvar_tbi.is_file():
        raise FileNotFoundError(f"ClinVar index does not exist: {clinvar_tbi}")
    report_clinvar = {
        "vcf": content_identity(clinvar_vcf),
        "tbi": content_identity(clinvar_tbi),
    }
    contracts = [_member_compatibility(member) for member in members]
    baseline = contracts[0]
    for member, observed in zip(members[1:], contracts[1:]):
        differing = [key for key in baseline if observed.get(key) != baseline.get(key)]
        if differing:
            raise ValueError(
                f"Incompatible cohort run {member.label!r}; contracts differ: "
                + ", ".join(differing)
            )
    missing = [
        key
        for key, value in baseline.items()
        if value is None or value == "" or value == [] or value == {}
    ]
    if missing:
        raise ValueError(
            "Cohort runs lack required modern compatibility provenance: "
            + ", ".join(missing)
        )
    expected_constants = {
        "target_assembly_accession": "GCF_000001405.40",
        "target_assembly_name": "GRCh38.p14",
        "ortholog_scope": "all",
        "alignment_event_mode": "compact_support",
        "event_ortholog_support_format": "event_group_id_v2",
        "annotation_schema": "normalized_annotation_evidence_v2",
        "gnomad_dataset": "gnomad_r4",
    }
    incompatible_constants = [
        key
        for key, expected in expected_constants.items()
        if str(baseline.get(key)) != expected
    ]
    if str(baseline.get("target_tax_id")) != "9606":
        incompatible_constants.append("target_tax_id")
    if incompatible_constants:
        raise ValueError(
            "Cohort runs do not use the current scientific constants: "
            + ", ".join(incompatible_constants)
        )
    for member in members:
        declared = _declared_clinvar_identity(member)
        if declared != report_clinvar:
            raise ValueError(
                f"Run {member.label!r} was annotated with different ClinVar VCF/index contents"
            )
    return {**baseline, "clinvar": report_clinvar}


def _member_compatibility(member: CohortMember) -> dict[str, object]:
    fetch = member.fetch_manifest
    alignment = member.alignment_manifest
    annotation = member.annotation_manifest
    vep = member.vep_manifest
    vep_config = vep.get("config")
    if not isinstance(vep_config, dict):
        vep_config = {}
    return {
        "pipeline_git_commit": member.root_manifest.get("git_commit"),
        "target_assembly_accession": fetch.get("target_assembly_accession"),
        "target_assembly_name": fetch.get("target_assembly_name"),
        "target_tax_id": fetch.get("target_tax_id"),
        "target_annotation_gff3_sha256": fetch.get("target_annotation_gff3_sha256"),
        "ortholog_scope": fetch.get("ortholog_scope"),
        "datasets_versions": fetch.get("datasets_versions"),
        "strategies": alignment.get("strategies"),
        "strategy_parameters": alignment.get("strategy_parameters"),
        "alignment_event_mode": alignment.get("alignment_event_mode"),
        "event_ortholog_support_format": alignment.get("event_ortholog_support_format"),
        "annotation_schema": annotation.get("schema"),
        "gnomad_api_url": annotation.get("gnomad_api_url"),
        "gnomad_dataset": annotation.get("gnomad_dataset"),
        "vep_backend": vep_config.get("backend"),
        "vep_release": vep_config.get("release"),
        "vep_columns": vep.get("columns"),
    }


def _declared_clinvar_identity(member: CohortMember) -> dict[str, object]:
    annotation = member.annotation_manifest
    identities = {}
    for key, manifest_key in (("vcf", "clinvar_vcf"), ("tbi", "clinvar_tbi")):
        metadata = annotation.get(manifest_key)
        if not isinstance(metadata, dict) or not str(metadata.get("path") or ""):
            raise ValueError(
                f"Run {member.label!r} annotation manifest lacks {manifest_key} provenance"
            )
        path = Path(str(metadata["path"])).expanduser()
        if not path.is_file():
            raise FileNotFoundError(
                f"Run {member.label!r} ClinVar provenance file is unavailable: {path}"
            )
        if path_metadata(path) != metadata:
            raise ValueError(
                f"Run {member.label!r} ClinVar provenance file changed after annotation: {path}"
            )
        identities[key] = content_identity(path)
    return identities


def _scientific_file_records(inputs: RunInputs, run_dir: Path) -> Iterable[dict[str, object]]:
    paths = [
        run_dir / "run_manifest.json",
        inputs.fetch_manifest_json,
        run_dir / "fetch" / "input.ids.tsv",
        inputs.genes_tsv,
        inputs.target_features_tsv,
        run_dir / "fetch" / "taxonomy.tsv.gz",
        run_dir / "fetch" / "orthologs.selected.tsv.gz",
        inputs.alignment_manifest_json,
        inputs.feature_coverage_tsv,
        inputs.strategy_summary_tsv,
        inputs.annotation_manifest_json,
        inputs.annotation_failures_tsv,
        inputs.variant_strategy_support_tsv,
        inputs.ortholog_evidence_summary_tsv,
        inputs.taxonomy_summary_tsv,
        inputs.variant_annotations_tsv.parent / "plan.json",
        inputs.variant_annotations_tsv.parent / "manifest.json",
    ]
    paths.extend(sorted(inputs.target_sequences_dir.glob("*.fa.gz")))
    source = resolve_variant_aggregation_source(inputs.variant_annotations_tsv)
    paths.extend(source.paths)
    unique = sorted(set(paths), key=lambda path: path.relative_to(run_dir).as_posix())
    for path in unique:
        if not path.is_file():
            raise FileNotFoundError(f"Missing scientific cohort input: {path}")
        identity = content_identity(path)
        yield {"path": path.relative_to(run_dir).as_posix(), **identity}


def _validate_shared_allele_evidence(variant_descriptor: Path) -> None:
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
    source = resolve_pre_vep_variant_source(
        variant_descriptor,
        required_columns=required,
    )
    candidate_checks = {
        "gnomad_af": "try_cast(nullif(gnomad_af, '') AS DOUBLE)",
        "gnomad_af_source": "nullif(gnomad_af_source, '')",
        "gnomad_csq": "nullif(gnomad_csq, '')",
        "clinvar_id": "nullif(clinvar_id, '')",
        "clinvar_allele_id": "nullif(clinvar_allele_id, '')",
        "clinvar_sig": "nullif(clinvar_sig, '')",
        "clinvar_revstat": "nullif(clinvar_revstat, '')",
        "clinvar_review_stars": "nullif(clinvar_review_stars, '')",
        "clinvar_review_stars_status": "nullif(clinvar_review_stars_status, '')",
        "clinvar_scv_count": "nullif(clinvar_scv_count, '')",
        "clinvar_hgvs": "nullif(clinvar_hgvs, '')",
        "clinvar_disease": "nullif(clinvar_disease, '')",
        "clinvar_variant_type": "nullif(clinvar_variant_type, '')",
    }
    checks = {
        field: expression
        for field, expression in candidate_checks.items()
        if field in source.columns
    }
    with tempfile.TemporaryDirectory(
        prefix=".cohort_preflight.", dir=variant_descriptor.parent
    ) as temporary:
        with duckdb.connect() as connection:
            connection.execute("SET preserve_insertion_order=false")
            connection.execute(f"SET temp_directory={sql_string(temporary)}")
            connection.execute(
                f"CREATE VIEW cohort_rows AS SELECT * FROM {variant_source_sql(source)}"
            )
            count_sql = ", ".join(
                f"count(DISTINCT {expression}) FILTER "
                f"(WHERE {expression} IS NOT NULL) AS {field}_count"
                for field, expression in checks.items()
            )
            conflict_sql = " OR ".join(
                f"{field}_count > 1" for field in checks
            )
            rows = connection.execute(
                "WITH evidence_counts AS (SELECT variant_key, "
                f"{count_sql} FROM cohort_rows WHERE variant_key <> '' "
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
        examples = ", ".join(f"{variant} ({field})" for field, variant in conflicts[:10])
        raise ValueError(
            "Cohort runs disagree on successful allele-level external evidence: " + examples
        )


def _concatenate_tsv(paths: list[Path], destination: Path) -> None:
    if not paths or any(not path.is_file() for path in paths):
        missing = next((path for path in paths if not path.is_file()), destination)
        raise FileNotFoundError(f"Missing cohort table: {missing}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(destination)
    try:
        expected_header: str | None = None
        with temporary.open("wb") as raw_output:
            with gzip.GzipFile(fileobj=raw_output, mode="wb", mtime=0) as output:
                for path in paths:
                    before = path.stat()
                    with gzip.open(path, "rt", newline="") as source:
                        header = source.readline()
                        if not header:
                            raise ValueError(f"Empty cohort table: {path}")
                        if expected_header is None:
                            expected_header = header
                            output.write(header.encode())
                        elif header != expected_header:
                            raise ValueError(f"Cohort table headers differ: {path}")
                        for line in source:
                            output.write(line.encode())
                    after = path.stat()
                    if (before.st_size, before.st_mtime_ns) != (
                        after.st_size,
                        after.st_mtime_ns,
                    ):
                        raise ValueError(f"Cohort source changed while reading: {path}")
        os.chmod(temporary, 0o644)
        _replace_if_changed(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _aggregate_strategy_summaries(paths: list[Path], destination: Path) -> None:
    snapshots = [(path, path.stat()) for path in paths]
    frames = [pd.read_csv(path, sep="\t", compression="gzip") for path in paths]
    for path, before in snapshots:
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (
            after.st_size,
            after.st_mtime_ns,
        ):
            raise ValueError(f"Cohort source changed while reading: {path}")
    columns = list(frames[0].columns)
    if any(list(frame.columns) != columns for frame in frames[1:]):
        raise ValueError("Cohort strategy summary columns differ")
    numeric = [column for column in columns if column != "strategy"]
    combined = pd.concat(frames, ignore_index=True)
    for column in numeric:
        combined[column] = pd.to_numeric(combined[column], errors="raise")
    grouped = combined.groupby("strategy", as_index=False, sort=True)[numeric].sum()
    _write_tsv_gzip_atomic(destination, grouped[columns])


def _aggregate_taxonomy_summary(
    members: list[CohortMember],
    destination: Path,
) -> None:
    profiles = {}
    selected_paths = []
    for member in members:
        taxonomy_path = member.inputs.run_dir / "fetch" / "taxonomy.tsv.gz"
        selected_path = member.inputs.run_dir / "fetch" / "orthologs.selected.tsv.gz"
        member_profiles = load_taxonomy_profiles(taxonomy_path)
        for tax_id, profile in member_profiles.items():
            existing = profiles.get(tax_id)
            if existing is not None and existing != profile:
                raise ValueError(
                    f"Cohort taxonomy disagrees for tax_id {tax_id}: {taxonomy_path}"
                )
            profiles[tax_id] = profile
        selected_paths.append(selected_path)

    def selected_rows() -> Iterable[dict[str, str]]:
        required = {"query_gene_id", "ortholog_gene_id", "tax_id"}
        for path in selected_paths:
            with gzip.open(path, "rt", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                missing = required - set(reader.fieldnames or [])
                if missing:
                    raise ValueError(
                        f"Selected ortholog table {path} is missing columns: "
                        + ", ".join(sorted(missing))
                    )
                yield from reader

    rows = build_taxonomy_summary_rows(selected_rows(), profiles)
    _write_tsv_gzip_atomic(
        destination,
        pd.DataFrame(rows, columns=TAXONOMY_SUMMARY_FIELDS),
    )


def _link_target_sequences(members: list[CohortMember], destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    expected = {}
    for member in members:
        for gene_id in member.target_gene_ids:
            source = member.inputs.target_sequences_dir / f"{gene_id}.fa.gz"
            if not source.is_file():
                raise FileNotFoundError(f"Missing target sequence for gene {gene_id}: {source}")
            expected[f"{gene_id}.fa.gz"] = source.resolve()
    for name, source in expected.items():
        link = destination / name
        if link.is_symlink() and link.resolve() == source:
            continue
        if link.is_symlink():
            link.unlink()
        elif link.exists():
            raise ValueError(f"Unexpected cohort target-sequence entry: {link}")
        link.symlink_to(source)
    unexpected = sorted(path.name for path in destination.iterdir() if path.name not in expected)
    if unexpected:
        raise ValueError(
            "Unexpected files in cohort target-sequence directory: " + ", ".join(unexpected)
        )


def _cohort_fetch_manifest(
    members: list[CohortMember], compatibility: dict[str, object]
) -> dict[str, object]:
    return {
        "stage": "fetch",
        "status": "complete",
        "cohort": True,
        "unique_gene_count": sum(len(member.requested_gene_ids) for member in members),
        "target_gene_count": sum(len(member.target_gene_ids) for member in members),
        **{
            key: compatibility[key]
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


def _cohort_alignment_manifest(
    members: list[CohortMember], compatibility: dict[str, object]
) -> dict[str, object]:
    return {
        "stage": "alignment",
        "cohort": True,
        "gene_ids": sorted(
            {gene for member in members for gene in member.target_gene_ids},
            key=_gene_sort_key,
        ),
        "gene_count": sum(len(member.target_gene_ids) for member in members),
        "strategies": compatibility["strategies"],
        "strategy_parameters": compatibility["strategy_parameters"],
        "alignment_event_mode": compatibility["alignment_event_mode"],
        "event_ortholog_support_format": compatibility["event_ortholog_support_format"],
    }


def _cohort_annotation_manifest(
    members: list[CohortMember], compatibility: dict[str, object]
) -> dict[str, object]:
    manifests = [member.annotation_manifest for member in members]
    return {
        "stage": "annotation",
        "schema": compatibility["annotation_schema"],
        "cohort": True,
        "gnomad_api_url": compatibility["gnomad_api_url"],
        "gnomad_dataset": compatibility["gnomad_dataset"],
        "failure_count": _sum_int(manifests, "failure_count"),
        "gnomad_region_failure_count": _sum_int(manifests, "gnomad_region_failure_count"),
        "event_key_status_counts": _sum_counters(manifests, "event_key_status_counts"),
    }


def _cohort_vep_manifest(members: list[CohortMember]) -> dict[str, object]:
    manifests = [member.vep_manifest for member in members]
    return {
        "schema_version": manifests[0].get("schema_version"),
        "status": "complete",
        "cohort": True,
        "config": manifests[0].get("config"),
        "row_count": _sum_int(manifests, "row_count"),
        "status_counts": _sum_counters(manifests, "status_counts"),
        "columns": manifests[0].get("columns"),
        "members": [member.fingerprint for member in members],
    }


def _sum_int(manifests: list[dict[str, object]], key: str) -> int:
    return sum(int(manifest.get(key, 0) or 0) for manifest in manifests)


def _sum_counters(manifests: list[dict[str, object]], key: str) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for manifest in manifests:
        values = manifest.get(key, {})
        if isinstance(values, dict):
            counts.update({str(name): int(value) for name, value in values.items()})
    return dict(sorted(counts.items()))


def _read_requested_gene_ids(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing normalized input IDs: {path}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"gene_id", "accepted"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Normalized input IDs have an incompatible header: {path}")
    return frame.loc[frame["accepted"].eq("true"), "gene_id"].astype(str).tolist()


def _read_gene_ids(path: Path) -> list[str]:
    frame = pd.read_csv(path, sep="\t", compression="gzip", usecols=["gene_id"], dtype=str)
    gene_ids = frame["gene_id"].astype(str).tolist()
    if len(gene_ids) != len(set(gene_ids)):
        raise ValueError(f"Run contains duplicate target gene IDs: {path}")
    return gene_ids


def _required_json(path: Path) -> dict[str, object]:
    value = read_json(path)
    if not value:
        raise FileNotFoundError(f"Missing or empty JSON manifest: {path}")
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _write_tsv_gzip_atomic(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = _temporary_path(path)
    try:
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            lineterminator="\n",
            compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
        )
        os.chmod(temporary, 0o644)
        _replace_if_changed(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _temporary_path(destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        delete=False,
    )
    handle.close()
    return Path(handle.name)


def _canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def _write_json_stable(path: Path, value: Any) -> None:
    text = json.dumps(value, indent=2, sort_keys=True) + "\n"
    if path.is_file() and path.read_text() == text:
        return
    write_text_atomic(path, text)


def _replace_if_changed(temporary: Path, destination: Path) -> None:
    if destination.is_file() and content_identity(temporary) == content_identity(destination):
        temporary.unlink()
        return
    temporary.replace(destination)


def _gene_sort_key(value: str) -> tuple[int, int | str]:
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)
