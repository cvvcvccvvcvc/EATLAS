from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from analytics.io import run_inputs as run_inputs_module
from analytics.io.artifacts import content_identity
from analytics.io.performance import PerformanceProfile
from analytics.io.run_inputs import (
    build_analysis_inputs,
    resolve_analysis_workspace,
    resolve_report_html,
    resolve_source_runs,
)
from analytics.analyses.variant_summary_aggregation import (
    resolve_variant_aggregation_source,
)
from genomics.taxonomy import TAXONOMY_FIELDS


VARIANT_FIELDS = [
    "variant_key",
    "gene_id",
    "event_type",
    "ref",
    "alt",
    "lookup_status",
    "strategies",
    "gnomad_af",
    "clinvar_id",
    "clinvar_sig",
    "clinvar_review_stars",
    "clinvar_scv_count",
    "vep_status",
    "vep_primary_consequence",
]


def _json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _table(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=fields).fillna("").to_csv(
        path,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
        lineterminator="\n",
    )


def _make_source_run(
    path: Path,
    *,
    gene_id: str,
    clinvar_vcf: Path,
    variant_key: str | None = None,
    gnomad_af: str = "",
) -> Path:
    fetch = path / "fetch"
    alignment = path / "alignment"
    annotation = path / "annotation"
    targets = fetch / "sequences" / "targets"
    targets.mkdir(parents=True)
    alignment.mkdir(parents=True)
    annotation.mkdir(parents=True)
    (fetch / "input.ids.tsv").write_text(
        "input_position\tline_number\traw_value\tgene_id\taccepted\t"
        "accepted_index\tduplicate_of_index\n"
        f"1\t1\t{gene_id}\t{gene_id}\ttrue\t1\t\n"
    )
    _table(
        fetch / "genes.tsv.gz",
        [
            "gene_id",
            "genomic_accession",
            "chromosome",
            "begin",
            "end",
            "sequence_length",
        ],
        [
            {
                "gene_id": gene_id,
                "genomic_accession": "NC_000001.11",
                "chromosome": "1",
                "begin": 1,
                "end": 10,
                "sequence_length": 10,
            }
        ],
    )
    _table(
        fetch / "target_features.tsv.gz",
        ["gene_id", "feature_type", "target_start0", "target_end0"],
        [
            {
                "gene_id": gene_id,
                "feature_type": "gene",
                "target_start0": 0,
                "target_end0": 10,
            }
        ],
    )
    _table(
        fetch / "orthologs.selected.tsv.gz",
        ["query_gene_id", "ortholog_gene_id", "tax_id"],
        [
            {
                "query_gene_id": gene_id,
                "ortholog_gene_id": f"{gene_id}_o1",
                "tax_id": "9598",
            }
        ],
    )
    taxonomy = {field: "" for field in TAXONOMY_FIELDS}
    taxonomy.update(
        {
            "tax_id": "9598",
            "taxonomy_status": "resolved",
            "lineage_tax_ids": "2759,33208,7742,32523,32524,40674,9443,9598",
            "species_id": "9598",
            "genus_id": "9596",
            "family_id": "9604",
            "order_id": "9443",
        }
    )
    _table(fetch / "taxonomy.tsv.gz", TAXONOMY_FIELDS, [taxonomy])
    (targets / f"{gene_id}.fa.gz").write_bytes(b"target")
    _table(
        annotation / "failures.tsv.gz",
        ["source", "scope", "chrom", "start", "end"],
        [],
    )

    shard = (
        annotation
        / "variant_annotations"
        / "partitions"
        / "partition_000001"
        / "shard_000001.tsv.gz"
    )
    _table(
        shard,
        VARIANT_FIELDS,
        [
            {
                "variant_key": variant_key or f"1:{100 + int(gene_id)}:A>G",
                "gene_id": gene_id,
                "event_type": "snv",
                "ref": "A",
                "alt": "G",
                "lookup_status": "ok",
                "strategies": "s1",
                "gnomad_af": gnomad_af,
                "clinvar_id": "",
                "clinvar_sig": "",
                "clinvar_review_stars": "",
                "clinvar_scv_count": "",
                "vep_status": "ok",
                "vep_primary_consequence": "missense_variant",
            }
        ],
    )
    descriptor = {
        "schema": "gaph_variant_annotation_dataset_v1",
        "status": "complete",
        "layout": "partitioned",
        "format": "tsv_gzip_v1",
        "path": "variant_annotations/manifest.json",
        "partition_count": 1,
        "shard_count": 1,
        "row_count": 1,
        "fields": VARIANT_FIELDS,
        "vep_config": {"backend": "local", "release": "116"},
        "vep_status_counts": {"ok": 1},
        "partitions": [
            {
                "partition_id": "partition_000001",
                "shard_count": 1,
                "row_count": 1,
                "shards": [
                    {
                        "shard_id": "shard_000001",
                        "path": "partitions/partition_000001/shard_000001.tsv.gz",
                        "row_count": 1,
                        "size_bytes": shard.stat().st_size,
                    }
                ],
            }
        ],
    }
    _json(annotation / "variant_annotations" / "manifest.json", descriptor)
    _json(
        path / "run_manifest.json",
        {
            "schema_version": 2,
            "pipeline": "gaph_v2",
            "status": "complete",
            "success": True,
            "git_commit": "a" * 40,
            "git_dirty": False,
            "session_id": path.name,
        },
    )
    _json(
        fetch / "manifest.json",
        {
            "stage": "fetch",
            "status": "complete",
            "unique_gene_count": 1,
            "target_assembly_accession": "GCF_000001405.40",
            "target_assembly_name": "GRCh38.p14",
            "target_tax_id": "9606",
            "target_annotation_gff3_sha256": "gff3",
            "ortholog_scope": "all",
            "datasets_versions": {"datasets": "18"},
        },
    )
    _json(
        alignment / "manifest.json",
        {
            "stage": "alignment",
            "gene_ids": [gene_id],
            "strategies": ["s1"],
            "strategy_parameters": {"s1": {}},
            "alignment_event_mode": "compact_support",
            "event_ortholog_support_format": "event_group_id_v2",
        },
    )
    _json(
        annotation / "manifest.json",
        {
            "stage": "annotation",
            "schema": "normalized_annotation_evidence_v5",
            "gnomad_api_url": "https://example.invalid",
            "gnomad_dataset": "gnomad_r4",
            "gnomad_observation_window": {
                "started_at_utc": f"2026-03-0{gene_id}T00:00:00+00:00",
                "finished_at_utc": f"2026-03-0{gene_id}T00:00:01+00:00",
            },
            "clinvar_vcf": content_identity(clinvar_vcf),
            "clinvar_tbi": content_identity(Path(f"{clinvar_vcf}.tbi")),
            "variant_annotations": descriptor,
        },
    )
    return path


def _snapshot(path: Path) -> dict[str, tuple[int, int]]:
    return {
        item.relative_to(path).as_posix(): (item.stat().st_size, item.stat().st_mtime_ns)
        for item in path.rglob("*")
        if item.is_file()
    }


def _stub_derived_builders(monkeypatch: pytest.MonkeyPatch) -> None:
    def alignment(_run: Path, *, analytics_dir: Path):
        summary = analytics_dir / "alignment_aggregates" / "strategy_summary.tsv.gz"
        coverage = analytics_dir / "alignment_aggregates" / "feature_coverage.tsv.gz"
        _table(
            summary,
            [
                "strategy",
                "gene_count",
                "summary_row_count",
                "aligned_summary_row_count",
                "event_count",
            ],
            [{"strategy": "s1", "gene_count": 1}],
        )
        _table(coverage, ["gene_id", "strategy", "feature_type"], [])
        return SimpleNamespace(
            strategy_summary_tsv=summary,
            feature_coverage_tsv=coverage,
        )

    def annotation(
        _run: Path,
        *,
        analytics_dir: Path,
        workers: int = 1,
        progress_callback=None,
    ):
        assert workers >= 1
        support = analytics_dir / "annotation_support" / "variant_strategy_support.tsv.gz"
        evidence = analytics_dir / "annotation_support" / "ortholog_evidence_summary.tsv.gz"
        _table(
            support,
            [
                "variant_key",
                "gene_id",
                "strategy",
                "alt_support_row_count",
                "alt_support_ortholog_count",
                "alt_support_family_count",
            ],
            [],
        )
        _table(
            evidence,
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
            [],
        )
        if progress_callback is not None:
            progress_callback(
                {
                    "partition_total": 2,
                    "partition_completed": 2,
                    "partition_built": 2,
                }
            )
        return SimpleNamespace(
            variant_strategy_support_tsv=support,
            ortholog_evidence_summary_tsv=evidence,
        )

    def taxonomy(_run: Path, *, analytics_dir: Path):
        output = analytics_dir / "taxonomy_summary" / "taxonomy_summary.tsv.gz"
        _table(
            output,
            [
                "taxonomic_scope",
                "evidence_unit",
                "gene_count",
                "ortholog_count",
                "taxon_count",
                "unit_count",
                "orthologs_per_gene_median",
                "units_per_gene_median",
            ],
            [],
        )
        return output

    monkeypatch.setattr(run_inputs_module, "resolve_alignment_aggregate_paths", alignment)
    monkeypatch.setattr(run_inputs_module, "resolve_annotation_support_paths", annotation)
    monkeypatch.setattr(run_inputs_module, "resolve_taxonomy_summary_path", taxonomy)


def test_analysis_workspace_is_order_independent_and_source_runs_are_read_only(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run_a = _make_source_run(tmp_path / "run-a", gene_id="1", clinvar_vcf=clinvar)
    run_b = _make_source_run(tmp_path / "run-b", gene_id="2", clinvar_vcf=clinvar)
    before = {run: _snapshot(run) for run in (run_a, run_b)}
    _stub_derived_builders(monkeypatch)

    first_sources = resolve_source_runs([run_a, run_b], clinvar_vcf=clinvar)
    second_sources = resolve_source_runs([run_b, run_a], clinvar_vcf=clinvar)
    first = build_analysis_inputs(
        first_sources,
        analytics_root=tmp_path / "analytics",
        scientific_config={"test": True},
    )
    second = build_analysis_inputs(
        second_sources,
        analytics_root=tmp_path / "analytics",
        scientific_config={"test": True},
    )

    assert first.analysis_id == second.analysis_id
    assert first.analysis_dir == second.analysis_dir
    assert first.analysis_dir.is_relative_to(tmp_path / "analytics")
    assert len(resolve_variant_aggregation_source(first.variant_annotation_sources).paths) == 2
    assert all(not (run / "analytics").exists() for run in (run_a, run_b))
    assert {run: _snapshot(run) for run in (run_a, run_b)} == before
    manifest = json.loads(first.analysis_manifest_json.read_text())
    assert manifest["status"] == "ready"
    assert len(manifest["sources"]) == 2
    assert first.annotation_manifest["gnomad_observation_window"] == {
        "started_at_utc": "2026-03-01T00:00:00+00:00",
        "finished_at_utc": "2026-03-02T00:00:01+00:00",
    }


def test_analysis_input_preparation_is_profiled_before_cache_builds(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run = _make_source_run(tmp_path / "run", gene_id="1", clinvar_vcf=clinvar)
    _stub_derived_builders(monkeypatch)
    sources = resolve_source_runs([run], clinvar_vcf=clinvar)
    workspace = resolve_analysis_workspace(
        sources,
        analytics_root=tmp_path / "analytics",
        scientific_config={"test": True},
    )
    report_path = resolve_report_html(workspace, "report")
    profile_path = workspace.analysis_dir / "performance" / "report.json"
    profile = PerformanceProfile(
        profile_path,
        analysis_dir=workspace.analysis_dir,
        analysis_id=workspace.analysis_id,
        report_path=report_path,
        source_run_dirs=[run],
    )
    assert profile_path.is_file()

    with profile.stage("Prepare analysis inputs"):
        build_analysis_inputs(
            sources,
            analytics_root=tmp_path / "analytics",
            scientific_config={"test": True},
            annotation_support_workers=2,
            workspace=workspace,
            performance_profile=profile,
        )

    payload = json.loads(profile_path.read_text())
    annotation_stage = next(
        stage for stage in payload["stages"] if stage["name"].startswith("Annotation support")
    )
    assert annotation_stage["parent_id"] == payload["stages"][0]["id"]
    assert annotation_stage["metrics"]["partition_completed"] == 2


def test_analysis_rejects_overlapping_source_genes(tmp_path: Path) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run_a = _make_source_run(tmp_path / "run-a", gene_id="1", clinvar_vcf=clinvar)
    run_b = _make_source_run(tmp_path / "run-b", gene_id="1", clinvar_vcf=clinvar)

    with pytest.raises(ValueError, match="overlapping accepted Gene IDs"):
        resolve_source_runs([run_a, run_b], clinvar_vcf=clinvar)


def test_analysis_rejects_incompatible_source_contracts(tmp_path: Path) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run_a = _make_source_run(tmp_path / "run-a", gene_id="1", clinvar_vcf=clinvar)
    run_b = _make_source_run(tmp_path / "run-b", gene_id="2", clinvar_vcf=clinvar)
    manifest_path = run_b / "alignment" / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["strategies"] = ["s1", "s2"]
    manifest["strategy_parameters"] = {"s1": {}, "s2": {}}
    _json(manifest_path, manifest)

    with pytest.raises(ValueError, match="contracts differ: strategies"):
        resolve_source_runs([run_a, run_b], clinvar_vcf=clinvar)


def test_analysis_accepts_moved_clinvar_with_identical_contents(tmp_path: Path) -> None:
    original = tmp_path / "original" / "clinvar.vcf.gz"
    original.parent.mkdir()
    original.write_bytes(b"clinvar")
    Path(f"{original}.tbi").write_bytes(b"index")
    run = _make_source_run(tmp_path / "run", gene_id="1", clinvar_vcf=original)

    moved = tmp_path / "moved" / "clinvar.vcf.gz"
    moved.parent.mkdir()
    moved.write_bytes(original.read_bytes())
    Path(f"{moved}.tbi").write_bytes(Path(f"{original}.tbi").read_bytes())
    original.unlink()
    Path(f"{original}.tbi").unlink()

    assert resolve_source_runs([run], clinvar_vcf=moved)[0].run_dir == run


def test_analysis_rejects_different_clinvar_contents(tmp_path: Path) -> None:
    original = tmp_path / "original.vcf.gz"
    original.write_bytes(b"clinvar")
    Path(f"{original}.tbi").write_bytes(b"index")
    run = _make_source_run(tmp_path / "run", gene_id="1", clinvar_vcf=original)
    different = tmp_path / "different.vcf.gz"
    different.write_bytes(b"different")
    Path(f"{different}.tbi").write_bytes(b"index")

    with pytest.raises(ValueError, match="different ClinVar contents"):
        resolve_source_runs([run], clinvar_vcf=different)


def test_analysis_rejects_previous_annotation_schema(tmp_path: Path) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run = _make_source_run(tmp_path / "run", gene_id="1", clinvar_vcf=clinvar)
    manifest_path = run / "annotation" / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["schema"] = "normalized_annotation_evidence_v4"
    _json(manifest_path, manifest)

    with pytest.raises(ValueError, match="Unsupported pipeline annotation contract"):
        resolve_source_runs([run], clinvar_vcf=clinvar)


def test_analysis_rejects_conflicting_shared_allele_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run_a = _make_source_run(
        tmp_path / "run-a",
        gene_id="1",
        clinvar_vcf=clinvar,
        variant_key="1:100:A>G",
        gnomad_af="0.01",
    )
    run_b = _make_source_run(
        tmp_path / "run-b",
        gene_id="2",
        clinvar_vcf=clinvar,
        variant_key="1:100:A>G",
        gnomad_af="0.02",
    )
    _stub_derived_builders(monkeypatch)
    sources = resolve_source_runs([run_a, run_b], clinvar_vcf=clinvar)

    with pytest.raises(ValueError, match="disagree on successful allele-level"):
        build_analysis_inputs(
            sources,
            analytics_root=tmp_path / "analytics",
            scientific_config={"test": True},
        )


def test_analysis_root_must_be_outside_source_run(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run = _make_source_run(tmp_path / "run", gene_id="1", clinvar_vcf=clinvar)
    _stub_derived_builders(monkeypatch)
    sources = resolve_source_runs([run], clinvar_vcf=clinvar)

    with pytest.raises(ValueError, match="outside every immutable source run"):
        build_analysis_inputs(
            sources,
            analytics_root=run / "analytics",
            scientific_config={"test": True},
        )

    assert not (run / "analytics").exists()


def test_analysis_has_no_legacy_variant_fallback(tmp_path: Path) -> None:
    clinvar = tmp_path / "clinvar.vcf.gz"
    clinvar.write_bytes(b"clinvar")
    Path(f"{clinvar}.tbi").write_bytes(b"index")
    run = _make_source_run(tmp_path / "run", gene_id="1", clinvar_vcf=clinvar)
    (run / "annotation" / "variant_annotations" / "manifest.json").unlink()
    legacy = run / "analytics" / "variant_annotations.tsv.gz"
    _table(legacy, VARIANT_FIELDS, [])

    with pytest.raises(FileNotFoundError, match="dataset manifest"):
        resolve_source_runs([run], clinvar_vcf=clinvar)
