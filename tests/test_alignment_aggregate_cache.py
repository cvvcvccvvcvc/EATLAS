from __future__ import annotations

import csv
import gzip
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from analytics.io import run_inputs as run_inputs_module
from analytics.derivations.alignment_summary import (
    concatenate_tsv_gz,
    merge_strategy_summaries,
    write_strategy_summary,
)
from analytics.derivations.feature_coverage import summarize_feature_coverage
from analytics.io.alignment_aggregates import (
    AlignmentAggregatePaths,
    build_or_load_alignment_aggregates,
    resolve_alignment_aggregate_paths,
)
from bin.alignment_table_schema import SEGMENT_FIELDS, SUMMARY_FIELDS


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


def _summary_row(
    gene_id: str,
    ortholog_gene_id: str,
    strategy: str,
    *,
    status: str = "aligned",
    event_count: int = 0,
    aligned_target_bp: int = 0,
) -> dict[str, object]:
    return {
        "gene_id": gene_id,
        "ortholog_gene_id": ortholog_gene_id,
        "tax_id": "9598",
        "taxname": "Pan troglodytes",
        "strategy": strategy,
        "tool": "test",
        "preset": "test",
        "status": status,
        "target_length": 10,
        "query_length": 10,
        "segment_count": int(status == "aligned"),
        "primary_segment_count": int(status == "aligned"),
        "secondary_segment_count": 0,
        "aligned_target_bp": aligned_target_bp,
        "aligned_query_bp": aligned_target_bp,
        "target_coverage": "0.800000" if aligned_target_bp else "0.000000",
        "query_coverage": "0.800000" if aligned_target_bp else "0.000000",
        "best_identity": "0.900000" if aligned_target_bp else "",
        "mean_identity": "0.900000" if aligned_target_bp else "",
        "event_count": event_count,
        "qc_flags": "",
    }


def _segment_row(
    gene_id: str,
    ortholog_gene_id: str,
    strategy: str,
    start0: int,
    end0: int,
) -> dict[str, object]:
    return {
        "gene_id": gene_id,
        "ortholog_gene_id": ortholog_gene_id,
        "tax_id": "9598",
        "taxname": "Pan troglodytes",
        "strategy": strategy,
        "tool": "test",
        "preset": "test",
        "sequence_id": f"ortholog_{ortholog_gene_id}",
        "target_id": f"target_{gene_id}",
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
        "native_record_id": "record_1",
        "qc_flags": "",
    }


def _write_evidence_run(run_dir: Path) -> tuple[list[Path], Path, Path]:
    alignment_dir = run_dir / "alignment"
    partitions_root = alignment_dir / "evidence" / "partitions"
    first = partitions_root / "partition_000001"
    second = partitions_root / "partition_000002"
    _write_tsv_gz(
        first / "ortholog_alignment_summary.tsv.gz",
        SUMMARY_FIELDS,
        [_summary_row("2", "201", "s1", event_count=2, aligned_target_bp=8)],
    )
    _write_tsv_gz(
        first / "alignment_segments.tsv.gz",
        SEGMENT_FIELDS,
        [_segment_row("2", "201", "s1", 0, 8)],
    )
    _write_tsv_gz(
        second / "ortholog_alignment_summary.tsv.gz",
        SUMMARY_FIELDS,
        [
            _summary_row("1", "101", "s1", event_count=1, aligned_target_bp=6),
            _summary_row("1", "102", "s2", status="no_alignment"),
        ],
    )
    _write_tsv_gz(
        second / "alignment_segments.tsv.gz",
        SEGMENT_FIELDS,
        [_segment_row("1", "101", "s1", 2, 8)],
    )

    target_features = run_dir / "fetch" / "target_features.tsv.gz"
    _write_tsv_gz(
        target_features,
        TARGET_FEATURE_FIELDS,
        [
            {
                "gene_id": gene_id,
                "feature_type": "gene",
                "feature_id": f"gene_{gene_id}",
                "genomic_accession": f"NC_{gene_id}",
                "genomic_start1": "101",
                "genomic_end1": "110",
                "target_start0": "0",
                "target_end0": "10",
                "length_bp": "10",
                "strand": "+",
            }
            for gene_id in ("1", "2")
        ],
    )
    alignment_manifest = alignment_dir / "manifest.json"
    alignment_manifest.write_text(
        json.dumps(
            {
                "stage": "alignment",
                "schema": "normalized_alignment_evidence_v2",
                "strategies": ["s1", "s2", "s_zero"],
                "normalized_evidence": {
                    "layout": "partitioned",
                    "format": "tsv_gzip_v1",
                    "path": "evidence/partitions",
                    "partition_count": 2,
                    "partition_files": [
                        "manifest.json",
                        "ortholog_alignment_summary.tsv.gz",
                        "alignment_segments.tsv.gz",
                        "alignment_events.tsv.gz",
                        "event_ortholog_support.tsv.gz",
                    ],
                    "event_group_id_scope": "partition",
                },
            }
        )
        + "\n"
    )
    return [first, second], target_features, alignment_manifest


def test_alignment_aggregate_rejects_noncanonical_manifest(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    _partitions, _target_features, alignment_manifest = _write_evidence_run(run_dir)
    manifest = json.loads(alignment_manifest.read_text())
    del manifest["schema"]
    alignment_manifest.write_text(json.dumps(manifest) + "\n")

    with pytest.raises(ValueError, match="unsupported schema"):
        resolve_alignment_aggregate_paths(run_dir, analytics_dir=tmp_path / "analytics")


def _read_text(path: Path) -> str:
    with gzip.open(path, "rt", newline="") as handle:
        return handle.read()


@pytest.mark.skipif(not BEDTOOLS_AVAILABLE, reason="bedtools is not installed")
def test_alignment_aggregate_cache_exactly_matches_current_builders(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "run"
    partitions, target_features, alignment_manifest = _write_evidence_run(run_dir)

    actual = resolve_alignment_aggregate_paths(
        run_dir,
        analytics_dir=tmp_path / "analytics",
    )

    expected_dir = tmp_path / "expected"
    expected_strategy = expected_dir / "strategy_summary.tsv.gz"
    expected_coverage = expected_dir / "feature_coverage.tsv.gz"
    partition_strategy_summaries = []
    for index, partition in enumerate(partitions, start=1):
        summary = expected_dir / f"strategy_{index:06d}.tsv.gz"
        write_strategy_summary(
            [partition / "ortholog_alignment_summary.tsv.gz"],
            summary,
            ["s1", "s2", "s_zero"],
        )
        partition_strategy_summaries.append(summary)
    merge_strategy_summaries(
        partition_strategy_summaries,
        expected_strategy,
        ["s1", "s2", "s_zero"],
    )
    partition_coverages = []
    for index, partition in enumerate(partitions, start=1):
        coverage = expected_dir / f"partition_{index:06d}.tsv.gz"
        summarize_feature_coverage(
            target_features,
            partition / "ortholog_alignment_summary.tsv.gz",
            partition / "alignment_segments.tsv.gz",
            coverage,
        )
        partition_coverages.append(coverage)
    concatenate_tsv_gz(partition_coverages, expected_coverage)

    assert _read_text(actual.strategy_summary_tsv) == _read_text(expected_strategy)
    assert _read_text(actual.feature_coverage_tsv) == _read_text(expected_coverage)
    with gzip.open(actual.strategy_summary_tsv, "rt", newline="") as handle:
        by_strategy = {
            row["strategy"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    assert by_strategy["s_zero"] == {
        "strategy": "s_zero",
        "summary_row_count": "0",
        "gene_count": "0",
        "aligned_summary_row_count": "0",
        "event_count": "0",
        "aligned_target_bp": "0",
    }
    assert json.loads(alignment_manifest.read_text())["strategies"] == [
        "s1",
        "s2",
        "s_zero",
    ]


@pytest.mark.skipif(not BEDTOOLS_AVAILABLE, reason="bedtools is not installed")
def test_alignment_aggregate_cache_hit_and_evidence_invalidation(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    analytics_dir = tmp_path / "analytics"
    partitions, target_features, alignment_manifest = _write_evidence_run(run_dir)
    first = build_or_load_alignment_aggregates(
        partition_dirs=partitions,
        target_features=target_features,
        alignment_manifest=alignment_manifest,
        analytics_dir=analytics_dir,
    )
    manifest_path = analytics_dir / "alignment_aggregates" / "manifest.json"
    first_manifest = json.loads(manifest_path.read_text())
    first_mtimes = (
        first.strategy_summary_tsv.stat().st_mtime_ns,
        first.feature_coverage_tsv.stat().st_mtime_ns,
    )

    second = build_or_load_alignment_aggregates(
        partition_dirs=partitions,
        target_features=target_features,
        alignment_manifest=alignment_manifest,
        analytics_dir=analytics_dir,
    )
    assert second == first
    assert (
        second.strategy_summary_tsv.stat().st_mtime_ns,
        second.feature_coverage_tsv.stat().st_mtime_ns,
    ) == first_mtimes

    _write_tsv_gz(
        partitions[0] / "alignment_segments.tsv.gz",
        SEGMENT_FIELDS,
        [_segment_row("2", "201", "s1", 0, 5)],
    )
    build_or_load_alignment_aggregates(
        partition_dirs=partitions,
        target_features=target_features,
        alignment_manifest=alignment_manifest,
        analytics_dir=analytics_dir,
    )
    second_manifest = json.loads(manifest_path.read_text())

    assert second_manifest["fingerprint"] != first_manifest["fingerprint"]
    assert second_manifest["inputs"] != first_manifest["inputs"]


def test_alignment_aggregate_resolution_requires_partitioned_evidence(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    with pytest.raises(FileNotFoundError, match="Missing normalized alignment evidence"):
        resolve_alignment_aggregate_paths(run_dir, analytics_dir=tmp_path / "analytics")

    (run_dir / "alignment" / "evidence" / "partitions").mkdir(parents=True)
    with pytest.raises(ValueError, match="contains no partitions"):
        resolve_alignment_aggregate_paths(run_dir, analytics_dir=tmp_path / "analytics")
