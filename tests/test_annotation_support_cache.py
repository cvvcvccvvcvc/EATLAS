from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from analytics.io import annotation_support as annotation_support_module
from analytics.io import run_inputs as run_inputs_module
from analytics.io.alignment_aggregates import AlignmentAggregatePaths
from analytics.io.annotation_support import (
    AnnotationSupportPaths,
    build_or_load_annotation_support,
    resolve_annotation_support_paths,
)
from analytics.derivations.ortholog_evidence import write_ortholog_evidence_summary
from analytics.derivations.support import (
    ORTHOLOG_EVIDENCE_FIELDS,
    VARIANT_STRATEGY_SUPPORT_FIELDS,
    merge_ortholog_evidence,
)
from analytics.derivations.taxonomy import (
    COUNT_KEYS,
    count_member_groups,
    load_taxonomy_profiles,
)
from bin.alignment_table_schema import SEGMENT_FIELDS
from bin.annotate_events import (
    EVENT_VARIANT_MAP_FIELDS,
    FAILURE_FIELDS,
)
from bin.merge_alignment_results import (
    COMPACT_EVENT_FIELDS,
    EVENT_ORTHOLOG_SUPPORT_FIELDS,
)
from bin.fetch_taxonomy import TAXONOMY_FIELDS


BEDTOOLS_AVAILABLE = shutil.which("bedtools") is not None
TARGET_FEATURE_FIELDS = [
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
]
SOURCE_ANNOTATION_FIELDS = [
    "variant_key",
    "gene_id",
    "lookup_status",
    "gnomad_af",
]


def _write_tsv_gz(
    path: Path,
    fields: list[str],
    rows: list[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _read_rows(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_text(path: Path) -> str:
    with gzip.open(path, "rt", newline="") as handle:
        return handle.read()


def _event(
    event_group_id: int,
    event_type: str,
    target_start0: int,
    ref: str,
    alt: str,
) -> dict[str, object]:
    return {
        "event_group_id": event_group_id,
        "gene_id": "1",
        "event_type": event_type,
        "target_start0": target_start0,
        "target_end0": target_start0 + max(1, len(ref)),
        "genomic_accession": "NC_000001.11",
        "genomic_start1": 100 + target_start0,
        "genomic_end1": 100 + target_start0 + max(0, len(ref) - 1),
        "ref": ref,
        "alt": alt,
        "strategy": "s1",
        "qc_flags": "",
    }


def _support(
    event_group_id: int,
    ortholog_gene_id: str,
    tax_id: str,
    support_row_count: int = 1,
) -> dict[str, object]:
    return {
        "event_group_id": event_group_id,
        "ortholog_gene_id": ortholog_gene_id,
        "tax_id": tax_id,
        "mapq": 60,
        "native_alignment_type": "primary",
        "support_row_count": support_row_count,
    }


def _segment(ortholog_gene_id: str, tax_id: str, start0: int = 0, end0: int = 30):
    return {
        "gene_id": "1",
        "ortholog_gene_id": ortholog_gene_id,
        "tax_id": tax_id,
        "taxname": f"taxon_{tax_id}",
        "strategy": "s1",
        "tool": "test",
        "preset": "test",
        "sequence_id": f"ortholog_{ortholog_gene_id}",
        "target_id": "target_1",
        "query_id": f"ortholog_{ortholog_gene_id}",
        "target_start0": start0,
        "target_end0": end0,
        "query_start0": start0,
        "query_end0": end0,
        "strand": "+",
        "matches": end0 - start0,
        "block_length": end0 - start0,
        "identity": "1.000000",
        "mapq": 60,
        "is_primary": "true",
        "divergence": "0.000000",
        "gap_compressed_divergence": "0.000000",
        "native_record_id": f"record_{ortholog_gene_id}_{start0}",
        "qc_flags": "",
    }


def _write_source_contract(run_dir: Path) -> dict[str, object]:
    partition_id = "partition_000001"
    partition = run_dir / "alignment" / "evidence" / "partitions" / partition_id
    map_root = run_dir / "annotation" / "event_variant_map" / "partitions"
    map_dir = map_root / partition_id

    events = [
        _event(1, "snv", 2, "A", "G"),
        _event(2, "snv", 3, "C", "T"),
        _event(3, "snv", 4, "G", "A"),
        _event(4, "del", 5, "A", ""),
        _event(5, "del", 6, "A", ""),
        _event(6, "snv", 7, "N", "A"),
    ]
    _write_tsv_gz(partition / "alignment_events.tsv.gz", COMPACT_EVENT_FIELDS, events)
    _write_tsv_gz(
        partition / "event_ortholog_support.tsv.gz",
        EVENT_ORTHOLOG_SUPPORT_FIELDS,
        [
            _support(1, "o1", "9598", 2),
            _support(1, "o2", "10090"),
            _support(2, "o1", "9598"),
            _support(3, "o2", "10090"),
            _support(3, "o3", "999999"),
            _support(4, "o1", "9598"),
            _support(4, "o2", "10090"),
            _support(5, "o1", "9598", 2),
            _support(6, "o1", "9598"),
        ],
    )
    _write_tsv_gz(
        partition / "alignment_segments.tsv.gz",
        SEGMENT_FIELDS,
        [
            _segment("o1", "9598", 0, 20),
            _segment("o1", "9598", 10, 30),
            _segment("o2", "10090"),
            _segment("o3", "999999"),
        ],
    )
    _write_tsv_gz(
        map_dir / "event_variant_map.tsv.gz",
        EVENT_VARIANT_MAP_FIELDS,
        [
            {"event_group_id": 1, "variant_key": "1:102:A>G", "normalization_status": "ok"},
            {"event_group_id": 2, "variant_key": "1:103:C>T", "normalization_status": "ok"},
            {"event_group_id": 3, "variant_key": "1:104:G>A", "normalization_status": "ok"},
            {"event_group_id": 4, "variant_key": "1:105:AA>A", "normalization_status": "ok"},
            {"event_group_id": 5, "variant_key": "1:105:AA>A", "normalization_status": "ok"},
            {
                "event_group_id": 6,
                "variant_key": "",
                "normalization_status": "non_concrete_allele",
            },
        ],
    )

    fetch = run_dir / "fetch"
    taxonomy = fetch / "taxonomy.tsv.gz"
    _write_tsv_gz(
        taxonomy,
        TAXONOMY_FIELDS,
        [
            {
                "tax_id": "9598",
                "taxonomy_status": "resolved",
                "lineage_tax_ids": "2759,33208,7742,32523,32524,40674,9443",
                "species_id": "9598",
                "genus_id": "9596",
                "family_id": "9604",
                "order_id": "9443",
            },
            {
                "tax_id": "10090",
                "taxonomy_status": "resolved",
                "lineage_tax_ids": "2759,33208,7742,32523,32524,40674",
                "species_id": "10090",
                "genus_id": "9596",
                "family_id": "10066",
                "order_id": "9989",
            },
            {
                "tax_id": "999999",
                "taxonomy_status": "not_returned",
                "lineage_tax_ids": "",
                "species_id": "",
                "genus_id": "",
                "family_id": "",
                "order_id": "",
            },
        ],
    )
    target_features = fetch / "target_features.tsv.gz"
    _write_tsv_gz(
        target_features,
        TARGET_FEATURE_FIELDS,
        [
            {
                "gene_id": "1",
                "feature_type": "gene",
                "feature_id": "gene_1",
                "genomic_accession": "NC_000001.11",
                "genomic_start1": 100,
                "genomic_end1": 129,
                "target_start0": 0,
                "target_end0": 30,
                "length_bp": 30,
                "strand": "+",
            },
            {
                "gene_id": "1",
                "feature_type": "cds",
                "feature_id": "cds_1",
                "genomic_accession": "NC_000001.11",
                "genomic_start1": 100,
                "genomic_end1": 103,
                "target_start0": 0,
                "target_end0": 4,
                "length_bp": 4,
                "strand": "+",
            },
            {
                "gene_id": "1",
                "feature_type": "utr",
                "feature_id": "utr_1",
                "genomic_accession": "NC_000001.11",
                "genomic_start1": 104,
                "genomic_end1": 106,
                "target_start0": 4,
                "target_end0": 7,
                "length_bp": 3,
                "strand": "+",
            },
        ],
    )

    annotation = run_dir / "annotation"
    source_annotations = annotation / "variant_annotations.tsv.gz"
    _write_tsv_gz(
        source_annotations,
        SOURCE_ANNOTATION_FIELDS,
        [
            {"variant_key": "1:102:A>G", "gene_id": "1", "lookup_status": "ok", "gnomad_af": "0"},
            {"variant_key": "1:103:C>T", "gene_id": "1", "lookup_status": "ok", "gnomad_af": ""},
            {"variant_key": "1:104:G>A", "gene_id": "1", "lookup_status": "ok", "gnomad_af": ""},
            {"variant_key": "1:105:AA>A", "gene_id": "1", "lookup_status": "ok", "gnomad_af": ""},
        ],
    )
    failures = annotation / "failures.tsv.gz"
    _write_tsv_gz(
        failures,
        FAILURE_FIELDS,
        [
            {
                "source": "gnomad",
                "scope": "region",
                "chrom": "1",
                "start": 104,
                "end": 104,
                "failure_type": "timeout",
                "message": "test failure",
            }
        ],
    )

    alignment_manifest = run_dir / "alignment" / "manifest.json"
    alignment_manifest.write_text(
        json.dumps(
            {
                "stage": "alignment",
                "schema": "normalized_alignment_evidence_v2",
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
                }
            }
        )
        + "\n"
    )
    annotation_manifest = annotation / "manifest.json"
    annotation_manifest.write_text(
        json.dumps(
            {
                "stage": "annotation",
                "schema": "normalized_annotation_evidence_v4",
                "partition_ids": [partition_id],
                "event_variant_map": {
                    "layout": "partitioned",
                    "format": "tsv_gzip_v1",
                    "path": "event_variant_map/partitions",
                    "partition_count": 1,
                    "row_count": len(events),
                    "fields": EVENT_VARIANT_MAP_FIELDS,
                    "event_group_id_scope": "partition",
                },
            }
        )
        + "\n"
    )
    return {
        "partition": partition,
        "map_root": map_root,
        "taxonomy": taxonomy,
        "target_features": target_features,
        "source_annotations": source_annotations,
        "failures": failures,
        "alignment_manifest": alignment_manifest,
        "annotation_manifest": annotation_manifest,
    }


@pytest.mark.skipif(not BEDTOOLS_AVAILABLE, reason="bedtools is not installed")
def test_annotation_support_cache_reproduces_current_report_contract(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    analytics_dir = tmp_path / "analytics"
    contract = _write_source_contract(run_dir)

    outputs = build_or_load_annotation_support(
        partition_dirs=[contract["partition"]],
        map_root=contract["map_root"],
        taxonomy=contract["taxonomy"],
        target_features=contract["target_features"],
        variant_annotations_source=contract["source_annotations"],
        failures=contract["failures"],
        alignment_manifest=contract["alignment_manifest"],
        annotation_manifest=contract["annotation_manifest"],
        analytics_dir=analytics_dir,
    )
    cache_manifest = json.loads(
        (analytics_dir / "annotation_support" / "manifest.json").read_text()
    )
    assert cache_manifest["schema_version"] == 6
    assert cache_manifest["exact_support"] == {
        "canonical_variant_count": 4,
        "canonical_support_edge_count": 8,
        "exact_support_edge_count": 7,
    }
    assert cache_manifest["timings_seconds"]["total"] >= 0
    assert cache_manifest["timings_seconds"]["aggregate_exact_support"] >= 0

    support_rows = _read_rows(outputs.variant_strategy_support_tsv)
    assert list(support_rows[0]) == VARIANT_STRATEGY_SUPPORT_FIELDS
    assert support_rows == [
        {
            "variant_key": "1:102:A>G",
            "gene_id": "1",
            "strategy": "s1",
            "alt_support_row_count": "3",
            "alt_support_ortholog_count": "2",
            "alt_support_genus_count": "1",
            "site_aligned_ortholog_count": "3",
        },
        {
            "variant_key": "1:103:C>T",
            "gene_id": "1",
            "strategy": "s1",
            "alt_support_row_count": "1",
            "alt_support_ortholog_count": "1",
            "alt_support_genus_count": "1",
            "site_aligned_ortholog_count": "3",
        },
        {
            "variant_key": "1:104:G>A",
            "gene_id": "1",
            "strategy": "s1",
            "alt_support_row_count": "2",
            "alt_support_ortholog_count": "2",
            "alt_support_genus_count": "2",
            "site_aligned_ortholog_count": "3",
        },
        {
            "variant_key": "1:105:AA>A",
            "gene_id": "1",
            "strategy": "s1",
            "alt_support_row_count": "4",
            "alt_support_ortholog_count": "2",
            "alt_support_genus_count": "",
            "site_aligned_ortholog_count": "",
        },
    ]

    evidence_rows = _read_rows(outputs.ortholog_evidence_summary_tsv)
    assert list(evidence_rows[0]) == ORTHOLOG_EVIDENCE_FIELDS
    by_key = {
        (
            row["target_context"],
            row["taxonomic_scope"],
            row["evidence_unit"],
            row["site_aligned_count"],
            row["alt_support_count"],
        ): row
        for row in evidence_rows
    }
    assert by_key[("cds", "all", "ortholog", "3", "2")][
        "gnomad_found_count"
    ] == "1"
    assert by_key[("cds", "all", "ortholog", "3", "1")][
        "gnomad_not_found_count"
    ] == "1"
    assert by_key[("utr", "all", "ortholog", "3", "2")][
        "gnomad_lookup_failed_count"
    ] == "1"
    primate = by_key[("cds", "primates", "genus", "1", "1")]
    assert primate["gnomad_found_count"] == "1"
    assert primate["gnomad_not_found_count"] == "1"
    assert (
        "utr",
        "primates",
        "ortholog",
        "1",
        "0",
    ) in by_key
    assert all(
        int(row["gnomad_found_count"])
        + int(row["gnomad_not_found_count"])
        + int(row["gnomad_lookup_failed_count"])
        > 0
        for row in evidence_rows
    )

    expected_depth = tmp_path / "expected_snv_taxonomic_depth.tsv.gz"
    expected_alt = tmp_path / "expected_snv_alt_taxonomic_support.tsv.gz"
    profiles = load_taxonomy_profiles(contract["taxonomy"])
    site_counts = count_member_groups(
        [("o1", "9598"), ("o2", "10090"), ("o3", "999999")],
        profiles,
    )
    _write_tsv_gz(
        expected_depth,
        ["gene_id", "strategy", "target_start0", *COUNT_KEYS],
        [
            {
                "gene_id": "1",
                "strategy": "s1",
                "target_start0": position,
                **site_counts,
            }
            for position in (2, 3, 4)
        ],
    )
    alt_members = {
        2: [("o1", "9598"), ("o2", "10090")],
        3: [("o1", "9598")],
        4: [("o2", "10090"), ("o3", "999999")],
    }
    _write_tsv_gz(
        expected_alt,
        ["gene_id", "strategy", "target_start0", "ref", "alt", *COUNT_KEYS],
        [
            {
                "gene_id": "1",
                "strategy": "s1",
                "target_start0": position,
                "ref": ref,
                "alt": alt,
                **count_member_groups(alt_members[position], profiles),
            }
            for position, ref, alt in ((2, "A", "G"), (3, "C", "T"), (4, "G", "A"))
        ],
    )
    expected_partition = tmp_path / "expected_partition"
    expected_partition.mkdir()
    expected_partition_evidence = expected_partition / "ortholog_evidence_summary.tsv.gz"
    expected_partition_count = write_ortholog_evidence_summary(
        expected_depth,
        expected_alt,
        [contract["target_features"]],
        {
            ("1", 2, "A", "G"): "found",
            ("1", 3, "C", "T"): "not_found",
            ("1", 4, "G", "A"): "lookup_failed",
        },
        expected_partition_evidence,
    )
    expected_evidence = tmp_path / "expected_ortholog_evidence_summary.tsv.gz"
    merge_ortholog_evidence(
        [
            (
                expected_partition,
                {"ortholog_evidence_summary_count": expected_partition_count},
            )
        ],
        expected_evidence,
    )
    assert _read_text(outputs.ortholog_evidence_summary_tsv) == _read_text(
        expected_evidence
    )


@pytest.mark.skipif(not BEDTOOLS_AVAILABLE, reason="bedtools is not installed")
def test_annotation_support_cache_hit_and_failure_invalidation(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    analytics_dir = tmp_path / "analytics"
    contract = _write_source_contract(run_dir)
    kwargs = {
        "partition_dirs": [contract["partition"]],
        "map_root": contract["map_root"],
        "taxonomy": contract["taxonomy"],
        "target_features": contract["target_features"],
        "variant_annotations_source": contract["source_annotations"],
        "failures": contract["failures"],
        "alignment_manifest": contract["alignment_manifest"],
        "annotation_manifest": contract["annotation_manifest"],
        "analytics_dir": analytics_dir,
    }
    first = build_or_load_annotation_support(**kwargs)
    cache_manifest = analytics_dir / "annotation_support" / "manifest.json"
    first_manifest = json.loads(cache_manifest.read_text())
    mtimes = (
        first.variant_strategy_support_tsv.stat().st_mtime_ns,
        first.ortholog_evidence_summary_tsv.stat().st_mtime_ns,
    )

    second = build_or_load_annotation_support(**kwargs)
    assert second == first
    assert (
        second.variant_strategy_support_tsv.stat().st_mtime_ns,
        second.ortholog_evidence_summary_tsv.stat().st_mtime_ns,
    ) == mtimes

    _write_tsv_gz(contract["failures"], FAILURE_FIELDS, [])
    build_or_load_annotation_support(**kwargs)
    second_manifest = json.loads(cache_manifest.read_text())
    assert second_manifest["fingerprint"] != first_manifest["fingerprint"]
    assert second_manifest["inputs"] != first_manifest["inputs"]


@pytest.mark.skipif(not BEDTOOLS_AVAILABLE, reason="bedtools is not installed")
def test_annotation_support_resumes_completed_partitions(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    run_dir = tmp_path / "run"
    analytics_dir = tmp_path / "analytics"
    contract = _write_source_contract(run_dir)
    first_partition = contract["partition"]
    second_partition = first_partition.parent / "partition_000002"
    shutil.copytree(first_partition, second_partition)
    first_map = contract["map_root"] / first_partition.name
    second_map = contract["map_root"] / second_partition.name
    shutil.copytree(first_map, second_map)

    alignment_manifest = json.loads(contract["alignment_manifest"].read_text())
    alignment_manifest["normalized_evidence"]["partition_count"] = 2
    contract["alignment_manifest"].write_text(json.dumps(alignment_manifest) + "\n")
    annotation_manifest = json.loads(contract["annotation_manifest"].read_text())
    annotation_manifest["partition_ids"] = [first_partition.name, second_partition.name]
    annotation_manifest["event_variant_map"]["partition_count"] = 2
    annotation_manifest["event_variant_map"]["row_count"] *= 2
    contract["annotation_manifest"].write_text(json.dumps(annotation_manifest) + "\n")

    dataset_root = run_dir / "annotation" / "variant_annotations"
    dataset_partitions = []
    for partition in (first_partition, second_partition):
        shard = (
            dataset_root
            / "partitions"
            / partition.name
            / "shard_000001.tsv.gz"
        )
        shard.parent.mkdir(parents=True)
        shutil.copyfile(contract["source_annotations"], shard)
        dataset_partitions.append(
            {
                "partition_id": partition.name,
                "shard_count": 1,
                "row_count": 4,
                "shards": [
                    {
                        "shard_id": "shard_000001",
                        "path": str(shard.relative_to(dataset_root)),
                        "row_count": 4,
                        "size_bytes": shard.stat().st_size,
                    }
                ],
            }
        )
    variant_manifest = dataset_root / "manifest.json"
    variant_manifest.write_text(
        json.dumps(
            {
                "schema": "gaph_variant_annotation_dataset_v1",
                "status": "complete",
                "layout": "partitioned",
                "format": "tsv_gzip_v1",
                "partition_count": 2,
                "shard_count": 2,
                "row_count": 8,
                "fields": SOURCE_ANNOTATION_FIELDS,
                "partitions": dataset_partitions,
            }
        )
        + "\n"
    )

    progress_updates: list[dict[str, object]] = []
    kwargs = {
        "partition_dirs": [first_partition, second_partition],
        "map_root": contract["map_root"],
        "taxonomy": contract["taxonomy"],
        "target_features": contract["target_features"],
        "variant_annotations_source": variant_manifest,
        "failures": contract["failures"],
        "alignment_manifest": contract["alignment_manifest"],
        "annotation_manifest": contract["annotation_manifest"],
        "analytics_dir": analytics_dir,
        "progress_callback": progress_updates.append,
    }
    original_collapse = annotation_support_module._collapse_partition

    def fail_second_partition(partition_id, *args, **inner_kwargs):
        if partition_id == second_partition.name:
            raise RuntimeError("injected partition failure")
        return original_collapse(partition_id, *args, **inner_kwargs)

    monkeypatch.setattr(
        annotation_support_module,
        "_collapse_partition",
        fail_second_partition,
    )
    with pytest.raises(RuntimeError, match="injected partition failure"):
        build_or_load_annotation_support(**kwargs)

    cache_dir = analytics_dir / "annotation_support"
    assert not (cache_dir / "manifest.json").exists()
    assert progress_updates[-1]["partition_completed"] == 1
    assert progress_updates[-1]["last_partition"] == first_partition.name
    completed = list((cache_dir / ".build").glob("*/partition_000001/*/manifest.json"))
    assert len(completed) == 1
    second_shard = (
        dataset_root
        / "partitions"
        / second_partition.name
        / "shard_000001.tsv.gz"
    )
    stat = second_shard.stat()
    os.utime(
        second_shard,
        ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000),
    )

    rebuilt_partitions: list[str] = []

    def record_collapse(partition_id, *args, **inner_kwargs):
        rebuilt_partitions.append(partition_id)
        return original_collapse(partition_id, *args, **inner_kwargs)

    monkeypatch.setattr(
        annotation_support_module,
        "_collapse_partition",
        record_collapse,
    )
    progress_updates.clear()
    outputs = build_or_load_annotation_support(**kwargs)

    assert rebuilt_partitions == [second_partition.name]
    assert progress_updates[0]["partition_reused"] == 1
    assert progress_updates[-1]["partition_completed"] == 2
    assert len(_read_rows(outputs.variant_strategy_support_tsv)) == 8
    final_manifest = json.loads((cache_dir / "manifest.json").read_text())
    assert final_manifest["partitions"] == {
        "total": 2,
        "built": 1,
        "reused": 1,
        "workers": 1,
    }
    assert not (cache_dir / ".build").exists()

    parallel_kwargs = {
        **kwargs,
        "analytics_dir": tmp_path / "parallel-analytics",
        "workers": 2,
    }
    parallel_outputs = build_or_load_annotation_support(**parallel_kwargs)
    assert _read_text(parallel_outputs.variant_strategy_support_tsv) == _read_text(
        outputs.variant_strategy_support_tsv
    )
    assert _read_text(parallel_outputs.ortholog_evidence_summary_tsv) == _read_text(
        outputs.ortholog_evidence_summary_tsv
    )


def test_annotation_support_resolution_rejects_missing_or_incomplete_contract(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "run"
    annotation = run_dir / "annotation"
    annotation.mkdir(parents=True)
    with pytest.raises(FileNotFoundError, match="Missing annotation manifest"):
        resolve_annotation_support_paths(run_dir, analytics_dir=tmp_path / "analytics")

    (annotation / "manifest.json").write_text(
        json.dumps(
            {
                "stage": "annotation",
                "schema": "normalized_annotation_evidence_v4",
            }
        )
        + "\n"
    )
    with pytest.raises(ValueError, match="does not declare the required partitioned"):
        resolve_annotation_support_paths(run_dir, analytics_dir=tmp_path / "analytics")

    (annotation / "manifest.json").write_text(
        json.dumps(
            {
                "stage": "annotation",
                "schema": "normalized_annotation_evidence_v4",
                "partition_ids": ["partition_000001"],
                "event_variant_map": {
                    "layout": "partitioned",
                    "format": "tsv_gzip_v1",
                    "path": "event_variant_map/partitions",
                    "partition_count": 1,
                    "fields": EVENT_VARIANT_MAP_FIELDS,
                    "event_group_id_scope": "partition",
                },
            }
        )
        + "\n"
    )
    with pytest.raises(FileNotFoundError, match="Incomplete analytics annotation-support"):
        resolve_annotation_support_paths(run_dir, analytics_dir=tmp_path / "analytics")


def test_annotation_support_import_does_not_require_pysam(tmp_path: Path) -> None:
    blocker = tmp_path / "blocked_import"
    blocker.mkdir()
    (blocker / "pysam.py").write_text("raise RuntimeError('pysam imported eagerly')\n")
    repository_root = Path(__file__).resolve().parents[1]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join((str(blocker), str(repository_root)))

    result = subprocess.run(
        [sys.executable, "-c", "import analytics.io.annotation_support"],
        cwd=repository_root,
        env=environment,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr
