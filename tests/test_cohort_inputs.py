from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from analytics.analyses.observed_variant_store import (
    build_or_load_observed_variant_store,
)
from analytics.analyses.variant_summary import build_variant_summary
from analytics.analyses.variant_summary_aggregation import (
    resolve_variant_aggregation_source,
)
from analytics.io.artifacts import file_identity, path_metadata
from analytics.io.cohort_inputs import resolve_cohort_inputs


BASE_COLUMNS = [
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
]
VEP_COLUMNS = [
    "vep_status",
    "vep_primary_consequence",
    "vep_consequence_terms",
]


def _write_gzip_table(path: Path, rows: list[dict], columns: list[str], *, header=True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=columns).fillna("").to_csv(
        path,
        sep="\t",
        index=False,
        header=header,
        lineterminator="\n",
        compression={"method": "gzip", "mtime": 0},
    )


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _make_run(
    root: Path,
    *,
    gene_id: str,
    clinvar_vcf: Path,
    strategies: list[str] | None = None,
    variant_key: str = "1:100:A>G",
    gnomad_af: str = "0.01",
    gnomad_failure: bool = False,
) -> Path:
    strategies = strategies or ["s1"]
    run = root
    fetch = run / "fetch"
    alignment = run / "alignment"
    annotation = run / "annotation"
    vep = run / "analytics" / "vep_consequences"
    targets = fetch / "sequences" / "targets"
    targets.mkdir(parents=True)
    alignment.mkdir(parents=True)
    annotation.mkdir(parents=True)
    vep.mkdir(parents=True)

    _write_json(
        run / "run_manifest.json",
        {
            "schema_version": 1,
            "pipeline": "gaph_v2",
            "status": "complete",
            "success": True,
            "stage": "all",
            "git_commit": "a" * 40,
            "git_dirty": False,
        },
    )
    (fetch / "input.ids.tsv").write_text(
        "input_position\tline_number\traw_value\tgene_id\taccepted\taccepted_index\tduplicate_of_index\n"
        f"1\t1\t{gene_id}\t{gene_id}\ttrue\t1\t\n"
    )
    _write_gzip_table(
        fetch / "genes.tsv.gz",
        [
            {
                "gene_id": gene_id,
                "chromosome": "1",
                "begin": "1",
                "end": "10",
                "sequence_length": "10",
            }
        ],
        ["gene_id", "chromosome", "begin", "end", "sequence_length"],
    )
    _write_gzip_table(
        fetch / "target_features.tsv.gz",
        [
            {
                "gene_id": gene_id,
                "feature_type": "gene",
                "target_start0": "0",
                "target_end0": "10",
            }
        ],
        ["gene_id", "feature_type", "target_start0", "target_end0"],
    )
    with gzip.open(targets / f"{gene_id}.fa.gz", "wt") as handle:
        handle.write(f">{gene_id}\nAAAAAAAAAA\n")
    _write_json(
        fetch / "manifest.json",
        {
            "stage": "fetch",
            "status": "complete",
            "unique_gene_count": 1,
            "target_gene_count": 1,
            "target_assembly_accession": "GCF_000001405.40",
            "target_assembly_name": "GRCh38.p14",
            "target_tax_id": "9606",
            "target_annotation_gff3_sha256": "b" * 64,
            "ortholog_scope": "all",
            "datasets_versions": ["datasets 18"],
        },
    )

    strategy_text = ",".join(strategies)
    _write_gzip_table(
        alignment / "feature_coverage.tsv.gz",
        [
            {
                "gene_id": gene_id,
                "strategy": strategy,
                "feature_type": "gene",
                "length_bp": 10,
                "covered_bases": 10,
            }
            for strategy in strategies
        ],
        ["gene_id", "strategy", "feature_type", "length_bp", "covered_bases"],
    )
    _write_gzip_table(
        alignment / "strategy_summary.tsv.gz",
        [
            {
                "strategy": strategy,
                "gene_count": 1,
                "summary_row_count": 2,
                "aligned_summary_row_count": 2,
                "event_count": 1,
                "aligned_target_bp": 10,
            }
            for strategy in strategies
        ],
        [
            "strategy",
            "gene_count",
            "summary_row_count",
            "aligned_summary_row_count",
            "event_count",
            "aligned_target_bp",
        ],
    )
    _write_json(
        alignment / "manifest.json",
        {
            "stage": "alignment",
            "gene_ids": [gene_id],
            "strategies": strategies,
            "strategy_parameters": {strategy: {"preset": "fixed"} for strategy in strategies},
            "alignment_event_mode": "compact_support",
            "event_ortholog_support_format": "event_group_id_v1",
        },
    )

    base_row = {
        "variant_key": variant_key,
        "gene_id": gene_id,
        "event_type": "snv",
        "ref": "A",
        "alt": "G",
        "lookup_status": "ok",
        "strategies": strategy_text,
        "support_row_count": "2",
        "support_ortholog_count": "2",
        "clinvar_id": "VCV1",
        "clinvar_sig": "Pathogenic",
        "clinvar_review_stars": "2",
        "clinvar_scv_count": "3",
        "gnomad_af": gnomad_af,
    }
    source = annotation / "variant_annotations.tsv.gz"
    _write_gzip_table(source, [base_row], BASE_COLUMNS)
    _write_gzip_table(
        annotation / "variant_strategy_support.tsv.gz",
        [
            {
                "variant_key": variant_key,
                "gene_id": gene_id,
                "strategy": strategy,
                "alt_support_row_count": 2,
                "alt_support_ortholog_count": 2,
                "alt_support_genus_count": 1,
            }
            for strategy in strategies
        ],
        [
            "variant_key",
            "gene_id",
            "strategy",
            "alt_support_row_count",
            "alt_support_ortholog_count",
            "alt_support_genus_count",
        ],
    )
    _write_gzip_table(
        annotation / "ortholog_evidence_summary.tsv.gz",
        [
            {
                "strategy": strategy,
                "target_context": "cds",
                "taxonomic_scope": "all",
                "evidence_unit": "ortholog",
                "site_aligned_count": 2,
                "alt_support_count": 2,
                "gnomad_found_count": 1,
                "gnomad_not_found_count": 0,
                "gnomad_lookup_failed_count": 0,
            }
            for strategy in strategies
        ],
        [
            "strategy",
            "target_context",
            "taxonomic_scope",
            "evidence_unit",
            "site_aligned_count",
            "alt_support_count",
            "gnomad_found_count",
            "gnomad_not_found_count",
            "gnomad_lookup_failed_count",
        ],
    )
    _write_gzip_table(
        annotation / "failures.tsv.gz",
        [
            {
                "source": "gnomad",
                "scope": "region",
                "chrom": "1",
                "start": "100",
                "end": "100",
                "failure_type": "request_failed",
                "message": "test failure",
            }
        ]
        if gnomad_failure
        else [],
        ["source", "scope", "chrom", "start", "end", "failure_type", "message"],
    )
    clinvar_tbi = Path(f"{clinvar_vcf}.tbi")
    _write_json(
        annotation / "manifest.json",
        {
            "output_mode": "unique_variant_context",
            "gnomad_api_url": "https://gnomad.test/graphql",
            "gnomad_dataset": "gnomad_r4",
            "clinvar_vcf": path_metadata(clinvar_vcf),
            "clinvar_tbi": path_metadata(clinvar_tbi),
            "failure_count": int(gnomad_failure),
            "gnomad_region_failure_count": int(gnomad_failure),
            "event_key_status_counts": {"ok": 1},
        },
    )

    input_partition = vep / "input" / "p1.tsv.gz"
    output_partition = vep / "partitions" / "p1.tsv.gz"
    enriched = {
        **base_row,
        "vep_status": "ok",
        "vep_primary_consequence": "missense_variant",
        "vep_consequence_terms": "missense_variant",
    }
    _write_gzip_table(input_partition, [base_row], BASE_COLUMNS)
    _write_gzip_table(output_partition, [enriched], [*BASE_COLUMNS, *VEP_COLUMNS], header=False)
    entry = {
        "partition_id": "p1",
        "path": "input/p1.tsv.gz",
        "row_count": 1,
        "file": file_identity(input_partition),
    }
    _write_json(
        vep / "plan.json",
        {
            "schema_version": 1,
            "status": "complete",
            "source": path_metadata(source),
            "row_count": 1,
            "input_columns": BASE_COLUMNS,
            "output_columns": [*BASE_COLUMNS, *VEP_COLUMNS],
            "partitions": [entry],
        },
    )
    config = {"backend": "local", "release": "116"}
    _write_json(
        vep / "partitions" / "p1.json",
        {
            "status": "complete",
            "input": entry,
            "row_count": 1,
            "config": config,
            "output_columns": [*BASE_COLUMNS, *VEP_COLUMNS],
            "output": file_identity(output_partition),
        },
    )
    final_output = vep / "variant_annotations.vep.tsv.gz"
    _write_gzip_table(final_output, [enriched], [*BASE_COLUMNS, *VEP_COLUMNS])
    _write_json(
        vep / "manifest.json",
        {
            "schema_version": 1,
            "status": "complete",
            "source": path_metadata(source),
            "config": config,
            "row_count": 1,
            "status_counts": {"ok": 1},
            "columns": [*BASE_COLUMNS, *VEP_COLUMNS],
            "output": file_identity(final_output),
        },
    )
    return run


def _manifest(path: Path, runs: list[tuple[str, Path]]) -> Path:
    _write_json(
        path,
        {
            "schema_version": 1,
            "runs": [
                {"label": label, "run_dir": str(run_dir)} for label, run_dir in runs
            ],
        },
    )
    return path


def test_cohort_resolves_disjoint_runs_as_stable_virtual_union(tmp_path: Path) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run_a = _make_run(tmp_path / "run-a", gene_id="1", clinvar_vcf=clinvar)
    run_b = _make_run(tmp_path / "run-b", gene_id="2", clinvar_vcf=clinvar)
    cohort_root = tmp_path / "cohort-output"

    first = resolve_cohort_inputs(
        _manifest(tmp_path / "cohort-a.json", [("a", run_a), ("b", run_b)]),
        cohort_root=cohort_root,
        clinvar_vcf=clinvar,
    )
    descriptor_mtime = first.variant_annotations_tsv.stat().st_mtime_ns
    reordered = resolve_cohort_inputs(
        _manifest(tmp_path / "cohort-b.json", [("b", run_b), ("a", run_a)]),
        cohort_root=cohort_root,
        clinvar_vcf=clinvar,
    )

    assert first.cohort_id == reordered.cohort_id
    assert first.run_dir == reordered.run_dir
    assert descriptor_mtime == reordered.variant_annotations_tsv.stat().st_mtime_ns
    assert sorted(
        pd.read_csv(first.genes_tsv, sep="\t")["gene_id"].astype(str).tolist()
    ) == ["1", "2"]
    source = resolve_variant_aggregation_source(first.variant_annotations_tsv)
    assert source.row_count == 2
    assert len(source.paths) == 2
    strategy_summary = pd.read_csv(first.strategy_summary_tsv, sep="\t")
    assert strategy_summary.loc[0, "gene_count"] == 2
    assert strategy_summary.loc[0, "summary_row_count"] == 4

    observed = build_or_load_observed_variant_store(
        variant_annotations_tsv=first.variant_annotations_tsv,
        analytics_dir=first.run_dir / "analytics",
        strategies=["s1"],
    )
    assert observed.manifest["source_row_count"] == 2
    assert observed.manifest["allele_gene_count"] == 2
    assert observed.manifest["allele_count"] == 1


def test_cohort_rejects_overlapping_accepted_gene_ids(tmp_path: Path) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run_a = _make_run(tmp_path / "run-a", gene_id="1", clinvar_vcf=clinvar)
    run_b = _make_run(tmp_path / "run-b", gene_id="1", clinvar_vcf=clinvar)

    with pytest.raises(ValueError, match="overlapping accepted gene IDs"):
        resolve_cohort_inputs(
            _manifest(tmp_path / "cohort.json", [("a", run_a), ("b", run_b)]),
            cohort_root=tmp_path / "cohorts",
            clinvar_vcf=clinvar,
        )


def test_cohort_rejects_incompatible_strategy_contracts(tmp_path: Path) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run_a = _make_run(tmp_path / "run-a", gene_id="1", clinvar_vcf=clinvar)
    run_b = _make_run(
        tmp_path / "run-b",
        gene_id="2",
        clinvar_vcf=clinvar,
        strategies=["s1", "s2"],
    )

    with pytest.raises(ValueError, match="contracts differ: strategies, strategy_parameters"):
        resolve_cohort_inputs(
            _manifest(tmp_path / "cohort.json", [("a", run_a), ("b", run_b)]),
            cohort_root=tmp_path / "cohorts",
            clinvar_vcf=clinvar,
        )


def test_cohort_rejects_conflicting_successful_allele_evidence(tmp_path: Path) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run_a = _make_run(
        tmp_path / "run-a", gene_id="1", clinvar_vcf=clinvar, gnomad_af="0.01"
    )
    run_b = _make_run(
        tmp_path / "run-b", gene_id="2", clinvar_vcf=clinvar, gnomad_af="0.02"
    )

    with pytest.raises(ValueError, match="disagree on successful allele-level"):
        resolve_cohort_inputs(
            _manifest(tmp_path / "cohort.json", [("a", run_a), ("b", run_b)]),
            cohort_root=tmp_path / "cohorts",
            clinvar_vcf=clinvar,
        )


def test_cohort_successful_gnomad_evidence_supersedes_failed_duplicate(
    tmp_path: Path,
) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run_a = _make_run(
        tmp_path / "run-a",
        gene_id="1",
        clinvar_vcf=clinvar,
        gnomad_af="",
        gnomad_failure=True,
    )
    run_b = _make_run(
        tmp_path / "run-b",
        gene_id="2",
        clinvar_vcf=clinvar,
        gnomad_af="0.01",
    )
    inputs = resolve_cohort_inputs(
        _manifest(tmp_path / "cohort.json", [("a", run_a), ("b", run_b)]),
        cohort_root=tmp_path / "cohorts",
        clinvar_vcf=clinvar,
    )

    summary = build_variant_summary(
        inputs.variant_annotations_tsv,
        inputs.run_dir / "analytics",
        strategy_label=str,
        target_features_path=inputs.target_features_tsv,
        genes_path=inputs.genes_tsv,
        annotation_failures_path=inputs.annotation_failures_tsv,
        variant_strategy_support_path=inputs.variant_strategy_support_tsv,
        ortholog_evidence_summary_path=inputs.ortholog_evidence_summary_tsv,
    )

    assert summary.unique_variant_count == 1
    assert summary.gnomad_found == 1
    assert summary.gnomad_lookup_failed == 0
