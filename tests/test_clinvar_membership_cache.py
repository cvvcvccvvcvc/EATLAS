from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pandas as pd

from analytics.analyses import clinvar_validation


def test_observed_clinvar_memberships_are_cached_and_reused(
    tmp_path: Path,
    monkeypatch,
) -> None:
    universe_path = tmp_path / "clinvar_universe.tsv.gz"
    annotations_path = tmp_path / "variant_annotations.tsv.gz"
    universe = pd.DataFrame(
        [
            {"variant_key": "1:10:A>G", "variant_type": "snv"},
            {"variant_key": "1:20:A>AT", "variant_type": "indel"},
        ]
    )
    universe.to_csv(universe_path, sep="\t", index=False, compression="gzip")
    with gzip.open(annotations_path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["variant_key", "strategies"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            [
                {"variant_key": "1:10:A>G", "strategies": "s1,s2"},
                {"variant_key": "1:20:A>AT", "strategies": "s2"},
                {"variant_key": "1:30:C>T", "strategies": "s1"},
            ]
        )

    observed, manifest, output_path, manifest_path = (
        clinvar_validation.build_or_load_observed_keys_by_strategy_type(
            universe=universe,
            universe_path=universe_path,
            variant_annotations_tsv=annotations_path,
            strategies=["s1", "s2"],
            analytics_dir=tmp_path / "analytics",
        )
    )

    assert not manifest["cache_hit"]
    assert manifest["membership_count"] == 3
    assert observed[("s1", "snv")] == {"1:10:A>G"}
    assert observed[("s1", "indel")] == set()
    assert observed[("s2", "indel")] == {"1:20:A>AT"}
    assert output_path.exists()
    assert manifest_path.exists()

    monkeypatch.setattr(
        clinvar_validation,
        "collect_observed_keys_by_strategy_type",
        lambda **_kwargs: (_ for _ in ()).throw(
            AssertionError("annotation scan was not cached")
        ),
    )
    cached, cached_manifest, _output_path, _manifest_path = (
        clinvar_validation.build_or_load_observed_keys_by_strategy_type(
            universe=universe,
            universe_path=universe_path,
            variant_annotations_tsv=annotations_path,
            strategies=["s1", "s2"],
            analytics_dir=tmp_path / "analytics",
        )
    )

    assert cached_manifest["cache_hit"]
    assert cached == observed
