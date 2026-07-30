"""Orchestration for candidate and ClinVar conservation analyses."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from analytics.io.performance import PerformanceProfile, profile_stage
from analytics.io.run_inputs import RunInputs
from .candidate_conservation import CandidateConservation, build_candidate_conservation
from .conservation import DEFAULT_TRACK_NAMES, build_conservation_annotations, universe_rows
from .conservation_validation import (
    ConservationValidation,
    build_conservation_cohort,
    compute_conservation_validation,
)


@dataclass(frozen=True)
class ConservationAnalysis:
    annotations_path: Path
    manifest_path: Path
    manifest: dict
    validation: ConservationValidation
    candidate: CandidateConservation


def alignment_gene_ids_by_strategy(coverage: pd.DataFrame) -> dict[str, set[str]]:
    if coverage.empty or not {"strategy", "gene_id"}.issubset(coverage.columns):
        return {}
    return {
        str(strategy): set(group["gene_id"].astype(str))
        for strategy, group in coverage.groupby("strategy", sort=False)
    }


def build_conservation_analysis(
    *,
    inputs: RunInputs,
    validation,
    strategies: list[str],
    eligible_gene_ids_by_strategy: dict[str, set[str]],
    firth_workers: int = 1,
    performance_profile: PerformanceProfile | None = None,
) -> ConservationAnalysis:
    with profile_stage(performance_profile, "Candidate phyloP distributions"):
        candidate = build_candidate_conservation(
            variant_annotations_tsv=inputs.variant_annotations_tsv,
            analytics_dir=inputs.run_dir / "analytics",
            annotation_failures_tsv=inputs.annotation_failures_tsv,
            additional_rows=universe_rows(validation.universe),
            track_names=DEFAULT_TRACK_NAMES,
        )
    with profile_stage(performance_profile, "ClinVar phyloP annotations"):
        conservation = build_conservation_annotations(
            universe=validation.universe,
            universe_path=validation.universe_path,
            analytics_dir=inputs.run_dir / "analytics",
            track_names=DEFAULT_TRACK_NAMES,
            position_scores=candidate.position_scores,
        )
    with profile_stage(performance_profile, "ClinVar model cohort"):
        cohort = build_conservation_cohort(
            universe=validation.universe,
            conservation=conservation.annotations,
            genes_tsv=inputs.genes_tsv,
            target_features_tsv=inputs.target_features_tsv,
            consequence_column=validation.consequence_column,
        )
    results = compute_conservation_validation(
        cohort=cohort,
        observed_by_strategy_type=validation.observed_by_strategy_type,
        strategies=strategies,
        analytics_dir=inputs.run_dir / "analytics",
        eligible_gene_ids_by_strategy=eligible_gene_ids_by_strategy,
        firth_workers=firth_workers,
        performance_profile=performance_profile,
    )
    return ConservationAnalysis(
        annotations_path=conservation.annotations_path,
        manifest_path=conservation.manifest_path,
        manifest=conservation.manifest,
        validation=results,
        candidate=candidate.without_position_scores(),
    )
