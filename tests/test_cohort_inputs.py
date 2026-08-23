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
from analytics.io.artifacts import path_metadata
from analytics.io.cohort_inputs import resolve_cohort_inputs
from bin.alignment_table_schema import SEGMENT_FIELDS, SUMMARY_FIELDS
from bin.fetch_taxonomy import TAXONOMY_FIELDS
from bin.finalize_annotation_partitions import EVENT_VARIANT_MAP_FIELDS
from bin.merge_alignment_results import (
    COMPACT_EVENT_FIELDS,
    EVENT_ORTHOLOG_SUPPORT_FIELDS,
)


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
    targets = fetch / "sequences" / "targets"
    targets.mkdir(parents=True)
    alignment.mkdir(parents=True)
    annotation.mkdir(parents=True)

    _write_json(
        run / "run_manifest.json",
        {
            "schema_version": 2,
            "pipeline": "gaph_v2",
            "status": "complete",
            "success": True,
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
                "feature_id": f"gene_{gene_id}",
                "genomic_accession": "NC_000001.11",
                "genomic_start1": "1",
                "genomic_end1": "10",
                "target_start0": "0",
                "target_end0": "10",
                "length_bp": "10",
                "strand": "+",
            }
        ],
        [
            "gene_id",
            "feature_type",
            "feature_id",
            "genomic_accession",
            "genomic_start1",
            "genomic_end1",
            "target_start0",
            "target_end0",
            "length_bp",
            "strand",
        ],
    )
    ortholog_ids = [f"{gene_id}_o1", f"{gene_id}_o2"]
    _write_gzip_table(
        fetch / "orthologs.selected.tsv.gz",
        [
            {
                "query_gene_id": gene_id,
                "ortholog_gene_id": ortholog_gene_id,
                "tax_id": "9598",
            }
            for ortholog_gene_id in ortholog_ids
        ],
        ["query_gene_id", "ortholog_gene_id", "tax_id"],
    )
    _write_gzip_table(
        fetch / "taxonomy.tsv.gz",
        [
            {
                "tax_id": "9598",
                "taxonomy_status": "resolved",
                "lineage_tax_ids": "2759,33208,7742,32523,32524,40674,9443,9598",
                "species_id": "9598",
                "genus_id": "9596",
                "family_id": "9604",
                "order_id": "9443",
            }
        ],
        TAXONOMY_FIELDS,
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
    partition_id = "partition_000001"
    evidence = alignment / "evidence" / "partitions" / partition_id
    map_partition = annotation / "event_variant_map" / "partitions" / partition_id
    _write_gzip_table(
        evidence / "ortholog_alignment_summary.tsv.gz",
        [
            {
                "gene_id": gene_id,
                "ortholog_gene_id": ortholog_gene_id,
                "tax_id": "9598",
                "taxname": "Pan troglodytes",
                "strategy": strategy,
                "tool": "test",
                "preset": "fixed",
                "status": "aligned",
                "target_length": 10,
                "query_length": 10,
                "segment_count": 1,
                "primary_segment_count": 1,
                "secondary_segment_count": 0,
                "aligned_target_bp": 10,
                "aligned_query_bp": 10,
                "target_coverage": "1.000000",
                "query_coverage": "1.000000",
                "best_identity": "1.000000",
                "mean_identity": "1.000000",
                "event_count": 1,
                "qc_flags": "",
            }
            for strategy in strategies
            for ortholog_gene_id in ortholog_ids
        ],
        SUMMARY_FIELDS,
    )
    _write_gzip_table(
        evidence / "alignment_segments.tsv.gz",
        [
            {
                "gene_id": gene_id,
                "ortholog_gene_id": ortholog_gene_id,
                "tax_id": "9598",
                "taxname": "Pan troglodytes",
                "strategy": strategy,
                "tool": "test",
                "preset": "fixed",
                "sequence_id": ortholog_gene_id,
                "target_id": f"target_{gene_id}",
                "query_id": ortholog_gene_id,
                "target_start0": 0,
                "target_end0": 10,
                "query_start0": 0,
                "query_end0": 10,
                "strand": "+",
                "matches": 10,
                "block_length": 10,
                "identity": "1.000000",
                "mapq": 60,
                "is_primary": "true",
                "divergence": "0.000000",
                "gap_compressed_divergence": "0.000000",
                "native_record_id": f"record_{strategy}_{ortholog_gene_id}",
                "qc_flags": "",
            }
            for strategy in strategies
            for ortholog_gene_id in ortholog_ids
        ],
        SEGMENT_FIELDS,
    )
    event_rows = []
    support_rows = []
    map_rows = []
    for event_group_id, strategy in enumerate(strategies, start=1):
        event_rows.append(
            {
                "event_group_id": event_group_id,
                "gene_id": gene_id,
                "event_type": "snv",
                "target_start0": 0,
                "target_end0": 1,
                "genomic_accession": "NC_000001.11",
                "genomic_start1": 100,
                "genomic_end1": 100,
                "ref": "A",
                "alt": "G",
                "strategy": strategy,
                "support_row_count": 2,
                "support_ortholog_count": 2,
                "qc_flags": "",
            }
        )
        support_rows.extend(
            {
                "event_group_id": event_group_id,
                "ortholog_gene_id": ortholog_gene_id,
                "tax_id": "9598",
                "taxname": "Pan troglodytes",
                "mapq": 60,
                "native_alignment_type": "primary",
                "support_row_count": 1,
            }
            for ortholog_gene_id in ortholog_ids
        )
        map_rows.append(
            {
                "event_group_id": event_group_id,
                "variant_key": variant_key,
                "normalization_status": "ok",
            }
        )
    _write_gzip_table(
        evidence / "alignment_events.tsv.gz", event_rows, COMPACT_EVENT_FIELDS
    )
    _write_gzip_table(
        evidence / "event_ortholog_support.tsv.gz",
        support_rows,
        EVENT_ORTHOLOG_SUPPORT_FIELDS,
    )
    _write_gzip_table(
        map_partition / "event_variant_map.tsv.gz",
        map_rows,
        EVENT_VARIANT_MAP_FIELDS,
    )
    _write_json(
        alignment / "manifest.json",
        {
            "stage": "alignment",
            "schema": "normalized_alignment_evidence_v2",
            "gene_ids": [gene_id],
            "strategies": strategies,
            "strategy_parameters": {strategy: {"preset": "fixed"} for strategy in strategies},
            "alignment_event_mode": "compact_support",
            "event_ortholog_support_format": "event_group_id_v2",
            "normalized_evidence": {
                "layout": "partitioned",
                "format": "tsv_gzip_v1",
                "path": "evidence/partitions",
                "partition_count": 1,
                "partition_files": [
                    "manifest.json",
                    "ortholog_alignment_summary.tsv.gz",
                    "alignment_segments.tsv.gz",
                    "alignment_events.tsv.gz",
                    "event_ortholog_support.tsv.gz",
                ],
                "event_group_id_scope": "partition",
            },
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
            "stage": "annotation",
            "schema": "normalized_annotation_evidence_v3",
            "partition_ids": [partition_id],
            "event_variant_map": {
                "layout": "partitioned",
                "format": "tsv_gzip_v1",
                "path": "event_variant_map/partitions",
                "partition_count": 1,
                "row_count": len(map_rows),
                "fields": EVENT_VARIANT_MAP_FIELDS,
                "event_group_id_scope": "partition",
            },
            "gnomad_api_url": "https://gnomad.test/graphql",
            "gnomad_dataset": "gnomad_r4",
            "clinvar_vcf": path_metadata(clinvar_vcf),
            "clinvar_tbi": path_metadata(clinvar_tbi),
            "failure_count": int(gnomad_failure),
            "gnomad_region_failure_count": int(gnomad_failure),
            "event_key_status_counts": {"ok": 1},
        },
    )
    enriched = {
        **base_row,
        "vep_status": "ok",
        "vep_primary_consequence": "missense_variant",
        "vep_consequence_terms": "missense_variant",
    }
    config = {"backend": "local", "release": "116"}
    dataset_dir = annotation / "variant_annotations"
    shard_relative = Path("partitions") / partition_id / "shard_000001.tsv.gz"
    shard = dataset_dir / shard_relative
    fields = [*BASE_COLUMNS, *VEP_COLUMNS]
    _write_gzip_table(shard, [enriched], fields)
    variant_annotations = {
        "schema": "gaph_variant_annotation_dataset_v1",
        "status": "complete",
        "layout": "partitioned",
        "format": "tsv_gzip_v1",
        "path": "variant_annotations/manifest.json",
        "partition_count": 1,
        "shard_count": 1,
        "row_count": 1,
        "fields": fields,
        "vep_config": config,
        "vep_status_counts": {"ok": 1},
        "partitions": [
            {
                "partition_id": partition_id,
                "shard_count": 1,
                "row_count": 1,
                "shards": [
                    {
                        "shard_id": "shard_000001",
                        "path": str(shard_relative),
                        "row_count": 1,
                        "size_bytes": shard.stat().st_size,
                    }
                ],
            }
        ],
    }
    _write_json(dataset_dir / "manifest.json", variant_annotations)
    annotation_manifest = json.loads((annotation / "manifest.json").read_text())
    annotation_manifest["variant_annotations"] = variant_annotations
    _write_json(annotation / "manifest.json", annotation_manifest)
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
    descriptor_mtime = first.variant_annotations_source.stat().st_mtime_ns
    reordered = resolve_cohort_inputs(
        _manifest(tmp_path / "cohort-b.json", [("b", run_b), ("a", run_a)]),
        cohort_root=cohort_root,
        clinvar_vcf=clinvar,
    )

    assert first.cohort_id == reordered.cohort_id
    assert first.run_dir == reordered.run_dir
    assert descriptor_mtime == reordered.variant_annotations_source.stat().st_mtime_ns
    assert sorted(
        pd.read_csv(first.genes_tsv, sep="\t")["gene_id"].astype(str).tolist()
    ) == ["1", "2"]
    source = resolve_variant_aggregation_source(first.variant_annotations_source)
    assert source.row_count == 2
    assert len(source.paths) == 2
    strategy_summary = pd.read_csv(first.strategy_summary_tsv, sep="\t")
    assert strategy_summary.loc[0, "gene_count"] == 2
    assert strategy_summary.loc[0, "summary_row_count"] == 4
    taxonomy_summary = pd.read_csv(first.taxonomy_summary_tsv, sep="\t")
    all_orthologs = taxonomy_summary[
        taxonomy_summary["taxonomic_scope"].eq("all")
        & taxonomy_summary["evidence_unit"].eq("ortholog")
    ].iloc[0]
    assert int(all_orthologs["gene_count"]) == 2
    assert int(all_orthologs["ortholog_count"]) == 4
    assert int(all_orthologs["taxon_count"]) == 1
    assert int(all_orthologs["unit_count"]) == 4
    assert float(all_orthologs["orthologs_per_gene_median"]) == 2.0
    resolved = json.loads(first.cohort_manifest_json.read_text())
    assert "taxonomy_summary" not in resolved["limitations"]
    scientific_paths = {
        record["path"]
        for member in resolved["members"]
        for record in member["scientific_files"]
    }
    assert "analytics/alignment_aggregates/strategy_summary.tsv.gz" in scientific_paths
    assert "analytics/alignment_aggregates/feature_coverage.tsv.gz" in scientific_paths
    assert "analytics/annotation_support/variant_strategy_support.tsv.gz" in scientific_paths
    assert "analytics/annotation_support/ortholog_evidence_summary.tsv.gz" in scientific_paths
    assert "analytics/taxonomy_summary/taxonomy_summary.tsv.gz" in scientific_paths
    assert "fetch/taxonomy.tsv.gz" in scientific_paths
    assert "fetch/orthologs.selected.tsv.gz" in scientific_paths
    assert "annotation/variant_annotations/manifest.json" in scientific_paths
    assert (
        "annotation/variant_annotations/partitions/partition_000001/shard_000001.tsv.gz"
        in scientific_paths
    )
    assert not (first.run_dir / "analytics" / "vep_consequences").exists()
    descriptor = json.loads(first.variant_annotations_source.read_text())
    assert {Path(member["path"]) for member in descriptor["members"]} == {
        run_a / "annotation" / "variant_annotations" / "manifest.json",
        run_b / "annotation" / "variant_annotations" / "manifest.json",
    }
    annotation_manifest = json.loads(first.annotation_manifest_json.read_text())
    variant_contract = annotation_manifest["variant_annotations"]
    assert variant_contract["schema"] == "gaph_variant_annotation_dataset_v1"
    assert variant_contract["layout"] == "cohort_partitioned"
    assert variant_contract["row_count"] == 2
    assert variant_contract["shard_count"] == 2
    assert variant_contract["partition_count"] == 2
    assert variant_contract["vep_config"] == {"backend": "local", "release": "116"}
    assert variant_contract["vep_status_counts"] == {"ok": 2}
    assert "alignment/strategy_summary.tsv.gz" not in scientific_paths
    assert "annotation/variant_strategy_support.tsv.gz" not in scientific_paths

    observed = build_or_load_observed_variant_store(
        variant_annotations_source=first.variant_annotations_source,
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


def test_cohort_rejects_member_without_normalized_alignment_evidence(
    tmp_path: Path,
) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run = _make_run(tmp_path / "run", gene_id="1", clinvar_vcf=clinvar)
    evidence = run / "alignment" / "evidence"
    for path in sorted(evidence.rglob("*"), reverse=True):
        if path.is_file():
            path.unlink()
        else:
            path.rmdir()
    evidence.rmdir()

    with pytest.raises(FileNotFoundError, match="Missing normalized alignment evidence"):
        resolve_cohort_inputs(
            _manifest(tmp_path / "cohort.json", [("old", run)]),
            cohort_root=tmp_path / "cohorts",
            clinvar_vcf=clinvar,
        )


def test_cohort_rejects_member_without_pipeline_variant_dataset(tmp_path: Path) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run = _make_run(tmp_path / "run", gene_id="1", clinvar_vcf=clinvar)
    dataset_manifest = run / "annotation" / "variant_annotations" / "manifest.json"
    dataset_manifest.unlink()
    legacy = run / "analytics" / "vep_consequences" / "variant_annotations.vep.tsv.gz"
    _write_gzip_table(legacy, [], [*BASE_COLUMNS, *VEP_COLUMNS])

    with pytest.raises(
        FileNotFoundError,
        match="Missing pipeline variant-annotation dataset manifest",
    ):
        resolve_cohort_inputs(
            _manifest(tmp_path / "cohort.json", [("old", run)]),
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


def test_cohort_rejects_incompatible_variant_annotation_config(tmp_path: Path) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run_a = _make_run(tmp_path / "run-a", gene_id="1", clinvar_vcf=clinvar)
    run_b = _make_run(tmp_path / "run-b", gene_id="2", clinvar_vcf=clinvar)
    dataset_path = run_b / "annotation" / "variant_annotations" / "manifest.json"
    dataset = json.loads(dataset_path.read_text())
    dataset["vep_config"]["release"] = "115"
    _write_json(dataset_path, dataset)
    annotation_path = run_b / "annotation" / "manifest.json"
    annotation = json.loads(annotation_path.read_text())
    annotation["variant_annotations"] = dataset
    _write_json(annotation_path, annotation)

    with pytest.raises(ValueError, match="contracts differ: vep_release"):
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
        inputs.variant_annotations_source,
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
