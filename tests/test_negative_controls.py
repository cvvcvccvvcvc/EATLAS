from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from analytics.core.negative_controls import (
    _collapse_sorted_segments,
    _finalize_same_site_options,
    _matched_ecdf,
    _matched_summary,
    _read_disjoint_contexts,
)


def test_callable_blocks_count_each_species_once(tmp_path: Path) -> None:
    sorted_segments = tmp_path / "segments.sorted.tsv"
    sorted_segments.write_text(
        "s1\t1\tspecies_a\t0\t10\n"
        "s1\t1\tspecies_a\t5\t15\n"
        "s1\t1\tspecies_b\t5\t12\n"
    )
    output = tmp_path / "callable_blocks.tsv.gz"

    count = _collapse_sorted_segments(sorted_segments, output)

    assert count == 3
    frame = pd.read_csv(output, sep="\t", compression="gzip")
    assert frame.to_dict("records") == [
        {
            "strategy": "s1",
            "gene_id": 1,
            "target_start0": 0,
            "target_end0": 5,
            "callable_species": 1,
        },
        {
            "strategy": "s1",
            "gene_id": 1,
            "target_start0": 5,
            "target_end0": 12,
            "callable_species": 2,
        },
        {
            "strategy": "s1",
            "gene_id": 1,
            "target_start0": 12,
            "target_end0": 15,
            "callable_species": 1,
        },
    ]


def test_same_site_control_excludes_alts_observed_by_strategy() -> None:
    frame = pd.DataFrame(
        [
            {
                "focal_id": "f1",
                "strategy": "s1",
                "role": "observed",
                "option": 0,
                "variant_key": "1:10:A>G",
            },
            {
                "focal_id": "f1",
                "strategy": "s1",
                "role": "control",
                "option": 1,
                "variant_key": "1:10:A>C",
            },
            {
                "focal_id": "f1",
                "strategy": "s1",
                "role": "control",
                "option": 2,
                "variant_key": "1:10:A>T",
            },
        ]
    )

    result = _finalize_same_site_options(frame, {("1:10:A>C", "s1")})

    assert result["variant_key"].tolist() == ["1:10:A>G", "1:10:A>T"]
    assert result["option"].tolist() == [0, 1]


def test_contexts_keep_noncoding_exon_separate_from_other_sequence(tmp_path: Path) -> None:
    features = pd.DataFrame(
        [
            {"gene_id": "1", "feature_type": "gene", "target_start0": 0, "target_end0": 20},
            {"gene_id": "1", "feature_type": "exon", "target_start0": 2, "target_end0": 12},
            {"gene_id": "1", "feature_type": "cds", "target_start0": 4, "target_end0": 8},
            {"gene_id": "1", "feature_type": "intron", "target_start0": 12, "target_end0": 18},
        ]
    )
    path = tmp_path / "features.tsv.gz"
    features.to_csv(path, sep="\t", index=False, compression="gzip")

    contexts = _read_disjoint_contexts(path, {"1": {"length": 20}})["1"]

    assert contexts == [
        (0, 2, "other"),
        (2, 4, "other_exon"),
        (4, 8, "cds"),
        (8, 12, "other_exon"),
        (12, 18, "intron"),
        (18, 20, "other"),
    ]


def test_matched_summary_uses_paired_control_options_deterministically() -> None:
    rows = []
    for index, (observed, controls) in enumerate([(2.0, [0.0, 1.0]), (4.0, [1.0, 3.0])]):
        focal_id = f"f{index}"
        rows.append(
            {
                "focal_id": focal_id,
                "strategy": "s1",
                "role": "observed",
                "phyloP100way": observed,
            }
        )
        rows.extend(
            {
                "focal_id": focal_id,
                "strategy": "s1",
                "role": "control",
                "phyloP100way": value,
            }
            for value in controls
        )
    frame = pd.DataFrame(rows)

    first = _matched_summary(frame, ["strategy"], permutations=200, seed=7)
    second = _matched_summary(frame, ["strategy"], permutations=200, seed=7)

    pd.testing.assert_frame_equal(first, second)
    assert first.loc[0, "matched_focals"] == 2
    assert first.loc[0, "observed_median"] == 3.0
    assert np.isfinite(first.loc[0, "empirical_p"])


def test_matched_ecdf_weights_each_focal_equally() -> None:
    frame = pd.DataFrame(
        [
            {"focal_id": "f1", "strategy": "s1", "role": "observed", "phyloP100way": 0.0},
            {"focal_id": "f1", "strategy": "s1", "role": "control", "phyloP100way": 0.0},
            {"focal_id": "f2", "strategy": "s1", "role": "observed", "phyloP100way": 2.0},
            {"focal_id": "f2", "strategy": "s1", "role": "control", "phyloP100way": 2.0},
            {"focal_id": "f2", "strategy": "s1", "role": "control", "phyloP100way": 2.0},
            {"focal_id": "f2", "strategy": "s1", "role": "control", "phyloP100way": 2.0},
        ]
    )

    ecdf = _matched_ecdf(frame)
    control = ecdf[ecdf["set"] == "Matched callable"]
    midpoint = control.loc[(control["phyloP100way"] - 1.0).abs().idxmin()]

    assert midpoint["fraction_leq"] == 0.5
