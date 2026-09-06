"""Stable identity and provenance for one analytics calculation runtime."""

from __future__ import annotations

import platform
import subprocess
from importlib import metadata
from pathlib import Path

from analytics.analyses import (
    basic_filtering,
    candidate_conservation,
    clinvar_conditions,
    clinvar_validation,
    conservation,
    conservation_validation,
    external_evidence,
    matched_control,
    observed_variant_store,
    variant_summary,
)
from analytics.io import alignment_aggregates, annotation_support, taxonomy_summary


ANALYSIS_SEMANTICS_VERSION = 2
COMPUTATION_DISTRIBUTIONS = (
    "duckdb",
    "numpy",
    "pandas",
    "pyBigWig",
    "scipy",
    "statsmodels",
)


def build_calculation_identity(
    *,
    firth_runtime: dict[str, str],
) -> dict[str, object]:
    """Describe everything intentionally allowed to reuse scientific caches."""

    return {
        "semantics_version": ANALYSIS_SEMANTICS_VERSION,
        "cache_versions": calculation_cache_versions(),
        "runtime": {
            "python": {
                "implementation": platform.python_implementation(),
                "version": platform.python_version(),
            },
            "packages": _distribution_versions(),
            "R": dict(firth_runtime),
        },
    }


def calculation_cache_versions() -> dict[str, object]:
    """Collect the owning version constants without duplicating their values."""

    return {
        "alignment_aggregates": alignment_aggregates.CACHE_SCHEMA_VERSION,
        "annotation_support": {
            "result": annotation_support.CACHE_SCHEMA_VERSION,
            "partition": annotation_support.PARTITION_CACHE_SCHEMA_VERSION,
        },
        "basic_filtering": basic_filtering.FILTER_SCORE_SCHEMA_VERSION,
        "candidate_conservation": candidate_conservation.CACHE_VERSION,
        "clinvar_conditions": clinvar_conditions.CACHE_SCHEMA_VERSION,
        "clinvar_validation": {
            "universe": clinvar_validation.CACHE_VERSION,
            "observed_membership": (
                clinvar_validation.OBSERVED_MEMBERSHIP_CACHE_VERSION
            ),
            "vep": clinvar_validation.VEP_CACHE_VERSION,
        },
        "conservation": conservation.CACHE_VERSION,
        "continuous_firth": conservation_validation.CONTINUOUS_CACHE_VERSION,
        "external_evidence": external_evidence.CACHE_VERSION,
        "matched_control": {
            "result": matched_control.CONTROL_VERSION,
            "focal": matched_control.FOCAL_CACHE_VERSION,
        },
        "observed_variant_store": observed_variant_store.STORE_SCHEMA_VERSION,
        "taxonomy_summary": taxonomy_summary.CACHE_SCHEMA_VERSION,
        "variant_summary": variant_summary.SUMMARY_CACHE_VERSION,
    }


def repository_provenance(project_dir: Path) -> dict[str, object]:
    """Record the checked-out revision without making it a cache invalidator."""

    root = project_dir.resolve()
    commit = _git(root, "rev-parse", "HEAD")
    status = _git(root, "status", "--porcelain", "--untracked-files=normal")
    return {
        "git_commit": commit,
        "git_dirty": bool(status),
    }


def _distribution_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    missing: list[str] = []
    for distribution in COMPUTATION_DISTRIBUTIONS:
        try:
            versions[distribution] = metadata.version(distribution)
        except metadata.PackageNotFoundError:
            missing.append(distribution)
    if missing:
        raise RuntimeError(
            "Analytics environment is missing required package metadata: "
            + ", ".join(missing)
        )
    return versions


def _git(root: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["git", *arguments],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or "unknown error"
        raise RuntimeError(f"Cannot record analytics code revision: {message}")
    return completed.stdout.strip()
