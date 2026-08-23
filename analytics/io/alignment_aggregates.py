"""Resolve analytics-owned aggregates from durable partitioned alignment evidence."""

from __future__ import annotations

import hashlib
import json
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

from analytics.io.artifacts import content_identity, file_identity, write_json_atomic


# The pipeline scripts are importable helpers, but still use sibling imports when
# executed as standalone staged files. Keep this compatibility bridge local to the
# shadow migration instead of duplicating their scientific algorithms.
_BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
_ADDED_BIN_PATH = str(_BIN_DIR) not in sys.path
if _ADDED_BIN_PATH:
    sys.path.insert(0, str(_BIN_DIR))
try:
    from feature_coverage import FEATURE_COVERAGE_FIELDS, summarize_feature_coverage
    from merge_alignment_results import (
        STRATEGY_SUMMARY_FIELDS,
        merge_strategy_summaries,
        merge_tsv_gz,
        write_strategy_summary,
    )
finally:
    if _ADDED_BIN_PATH:
        sys.path.remove(str(_BIN_DIR))


CACHE_SCHEMA_VERSION = 1
CACHE_DIRNAME = "alignment_aggregates"
STRATEGY_SUMMARY_FILENAME = "strategy_summary.tsv.gz"
FEATURE_COVERAGE_FILENAME = "feature_coverage.tsv.gz"
SUMMARY_FILENAME = "ortholog_alignment_summary.tsv.gz"
SEGMENTS_FILENAME = "alignment_segments.tsv.gz"


@dataclass(frozen=True)
class AlignmentAggregatePaths:
    strategy_summary_tsv: Path
    feature_coverage_tsv: Path


def resolve_alignment_aggregate_paths(run_dir: Path) -> AlignmentAggregatePaths:
    """Require durable partition evidence and expose analytics-owned aggregates."""

    alignment_dir = run_dir / "alignment"
    partitions_root = alignment_dir / "evidence" / "partitions"
    if not partitions_root.exists():
        raise FileNotFoundError(
            "Missing normalized alignment evidence required for analytics: "
            f"{partitions_root}"
        )
    if not partitions_root.is_dir():
        raise NotADirectoryError(
            f"Alignment evidence partitions path is not a directory: {partitions_root}"
        )

    partition_dirs = sorted(
        (path for path in partitions_root.iterdir() if path.is_dir()),
        key=lambda path: path.name,
    )
    if not partition_dirs:
        raise ValueError(f"Alignment evidence contains no partitions: {partitions_root}")
    for partition_dir in partition_dirs:
        for filename in (SUMMARY_FILENAME, SEGMENTS_FILENAME):
            path = partition_dir / filename
            if not path.is_file():
                raise FileNotFoundError(
                    f"Alignment evidence partition {partition_dir.name} is missing {filename}"
                )

    alignment_manifest = alignment_dir / "manifest.json"
    target_features = run_dir / "fetch" / "target_features.tsv.gz"
    for path in (alignment_manifest, target_features):
        if not path.is_file():
            raise FileNotFoundError(f"Missing alignment aggregate input: {path}")
    return build_or_load_alignment_aggregates(
        partition_dirs=partition_dirs,
        target_features=target_features,
        alignment_manifest=alignment_manifest,
        analytics_dir=run_dir / "analytics",
    )


def build_or_load_alignment_aggregates(
    *,
    partition_dirs: list[Path],
    target_features: Path,
    alignment_manifest: Path,
    analytics_dir: Path,
) -> AlignmentAggregatePaths:
    """Build report aggregates from normalized evidence partitions."""

    if not partition_dirs:
        raise ValueError("Alignment aggregates require at least one evidence partition")
    partition_dirs = sorted(partition_dirs, key=lambda path: path.name)
    partition_ids = [path.name for path in partition_dirs]
    if len(partition_ids) != len(set(partition_ids)):
        raise ValueError(f"Duplicate alignment evidence partition IDs: {partition_ids}")
    cache_dir = analytics_dir / CACHE_DIRNAME
    outputs = AlignmentAggregatePaths(
        strategy_summary_tsv=cache_dir / STRATEGY_SUMMARY_FILENAME,
        feature_coverage_tsv=cache_dir / FEATURE_COVERAGE_FILENAME,
    )
    manifest_path = cache_dir / "manifest.json"
    expected_strategies = _read_expected_strategies(
        alignment_manifest,
        partition_dirs,
    )
    inputs = _input_identities(
        partition_dirs,
        target_features,
        alignment_manifest,
        expected_strategies,
    )
    fingerprint = _fingerprint(inputs)
    if _cache_is_valid(manifest_path, outputs, inputs, fingerprint):
        return outputs

    cache_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".alignment_aggregates_",
        dir=cache_dir,
    ) as temporary_name:
        temporary_dir = Path(temporary_name)
        strategy_output = temporary_dir / STRATEGY_SUMMARY_FILENAME
        feature_output = temporary_dir / FEATURE_COVERAGE_FILENAME
        summary_row_count = 0
        partition_strategy_summaries: list[Path] = []
        partition_coverages: list[Path] = []
        for index, partition_dir in enumerate(partition_dirs, start=1):
            partition_strategy = temporary_dir / f"strategy_{index:06d}.tsv.gz"
            partition_summary_count, _ = write_strategy_summary(
                [partition_dir / SUMMARY_FILENAME],
                partition_strategy,
                expected_strategies,
            )
            summary_row_count += partition_summary_count
            partition_strategy_summaries.append(partition_strategy)
            partition_coverage = temporary_dir / f"coverage_{index:06d}.tsv.gz"
            summarize_feature_coverage(
                target_features,
                partition_dir / SUMMARY_FILENAME,
                partition_dir / SEGMENTS_FILENAME,
                partition_coverage,
            )
            partition_coverages.append(partition_coverage)
        strategy_count = merge_strategy_summaries(
            partition_strategy_summaries,
            strategy_output,
            expected_strategies,
        )
        feature_coverage_count = merge_tsv_gz(partition_coverages, feature_output)

        strategy_output.chmod(0o644)
        feature_output.chmod(0o644)
        strategy_output.replace(outputs.strategy_summary_tsv)
        feature_output.replace(outputs.feature_coverage_tsv)

    write_json_atomic(
        manifest_path,
        {
            "schema_version": CACHE_SCHEMA_VERSION,
            "status": "complete",
            "fingerprint": fingerprint,
            "inputs": inputs,
            "strategy_summary": {
                "columns": STRATEGY_SUMMARY_FIELDS,
                "source_row_count": summary_row_count,
                "row_count": strategy_count,
                "output": file_identity(outputs.strategy_summary_tsv),
            },
            "feature_coverage": {
                "columns": FEATURE_COVERAGE_FIELDS,
                "row_count": feature_coverage_count,
                "output": file_identity(outputs.feature_coverage_tsv),
            },
        },
    )
    return outputs


def _read_expected_strategies(path: Path, partition_dirs: list[Path]) -> list[str]:
    try:
        manifest = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid alignment manifest JSON: {path}") from exc
    if manifest.get("stage") != "alignment":
        raise ValueError(f"Alignment manifest has invalid stage: {path}")
    if manifest.get("schema") != "normalized_alignment_evidence_v1":
        raise ValueError(f"Alignment manifest has unsupported schema: {path}")
    evidence = manifest.get("normalized_evidence")
    expected_contract = {
        "layout": "partitioned",
        "format": "tsv_gzip_v1",
        "path": "evidence/partitions",
        "event_group_id_scope": "partition",
        "partition_files": [
            "manifest.json",
            "ortholog_alignment_summary.tsv.gz",
            "alignment_segments.tsv.gz",
            "alignment_events.tsv.gz",
            "event_ortholog_support.tsv.gz",
        ],
    }
    if not isinstance(evidence, dict) or any(
        evidence.get(field) != value for field, value in expected_contract.items()
    ):
        raise ValueError(f"Alignment manifest has invalid normalized_evidence: {path}")
    if evidence.get("partition_count") != len(partition_dirs):
        raise ValueError(
            "Alignment manifest partition_count does not match evidence directories: "
            f"{path}"
        )
    strategies = manifest.get("strategies")
    if (
        not isinstance(strategies, list)
        or not strategies
        or not all(isinstance(value, str) and value for value in strategies)
        or len(strategies) != len(set(strategies))
    ):
        raise ValueError(f"Alignment manifest has invalid strategies: {strategies!r}")
    return list(strategies)


def _input_identities(
    partition_dirs: list[Path],
    target_features: Path,
    alignment_manifest: Path,
    expected_strategies: list[str],
) -> dict[str, object]:
    return {
        "alignment_manifest": content_identity(alignment_manifest),
        "expected_strategies": expected_strategies,
        "target_features": content_identity(target_features),
        "partitions": [
            {
                "partition_id": partition_dir.name,
                "ortholog_alignment_summary": file_identity(
                    partition_dir / SUMMARY_FILENAME
                ),
                "alignment_segments": file_identity(
                    partition_dir / SEGMENTS_FILENAME
                ),
            }
            for partition_dir in partition_dirs
        ],
    }


def _fingerprint(inputs: dict[str, object]) -> str:
    payload = {
        "schema_version": CACHE_SCHEMA_VERSION,
        "inputs": inputs,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _cache_is_valid(
    manifest_path: Path,
    outputs: AlignmentAggregatePaths,
    inputs: dict[str, object],
    fingerprint: str,
) -> bool:
    if not manifest_path.is_file() or not all(
        path.is_file()
        for path in (outputs.strategy_summary_tsv, outputs.feature_coverage_tsv)
    ):
        return False
    try:
        manifest = json.loads(manifest_path.read_text())
        return (
            manifest.get("schema_version") == CACHE_SCHEMA_VERSION
            and manifest.get("status") == "complete"
            and manifest.get("fingerprint") == fingerprint
            and manifest.get("inputs") == inputs
            and manifest.get("strategy_summary", {}).get("columns")
            == STRATEGY_SUMMARY_FIELDS
            and manifest.get("strategy_summary", {}).get("output")
            == file_identity(outputs.strategy_summary_tsv)
            and manifest.get("feature_coverage", {}).get("columns")
            == FEATURE_COVERAGE_FIELDS
            and manifest.get("feature_coverage", {}).get("output")
            == file_identity(outputs.feature_coverage_tsv)
        )
    except (OSError, ValueError, json.JSONDecodeError):
        return False
