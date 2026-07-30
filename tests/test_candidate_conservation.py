from __future__ import annotations

import csv
import gzip
from pathlib import Path

from analytics.analyses import candidate_conservation as candidate
from analytics.analyses.conservation import PositionScores, parse_tracks


def test_candidate_conservation_deduplicates_memberships_and_reuses_cache(
    tmp_path: Path,
    monkeypatch,
) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    failures = tmp_path / "failures.tsv.gz"
    fields = sorted(candidate.REQUIRED_COLUMNS)
    rows = [
        {
            "variant_key": "1:1:A>G",
            "lookup_status": "ok",
            "strategies": "s1,s2",
            "gnomad_af": "0.1",
        },
        {
            "variant_key": "1:1:A>G",
            "lookup_status": "ok",
            "strategies": "s1,s2",
            "gnomad_af": "0.1",
        },
        {
            "variant_key": "1:2:C>T",
            "lookup_status": "ok",
            "strategies": "s1",
            "gnomad_af": "",
        },
        {
            "variant_key": "1:3:G>A",
            "lookup_status": "ok",
            "strategies": "s1",
            "gnomad_af": "",
        },
        {
            "variant_key": "raw",
            "lookup_status": "raw_no_context",
            "strategies": "s1",
            "gnomad_af": "",
        },
    ]
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with gzip.open(failures, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source", "scope", "chrom", "start", "end"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"source": "gnomad", "scope": "region", "chrom": "1", "start": 2, "end": 2})

    captured = {}
    track = parse_tracks("phyloP100way")[0]

    def fake_read_position_scores(**kwargs):
        captured["positions"] = kwargs["positions_by_chrom"]
        return PositionScores(
            track,
            {("chr1", 0): -1.0, ("chr1", 1): 2.0, ("chr1", 2): 3.0, ("chr1", 3): 0.5},
            {
                "track": track.name,
                "status": "complete",
                "unique_positions": 4,
                "annotated_positions": 4,
                "missing_positions": 0,
                "failed_block_count": 0,
            },
        )

    monkeypatch.setattr(candidate, "read_position_scores", fake_read_position_scores)
    result = candidate.build_candidate_conservation(
        variant_annotations_tsv=annotations,
        analytics_dir=tmp_path / "analytics",
        annotation_failures_tsv=failures,
        additional_rows=[{"chrom": "1", "pos": "4", "ref": "T", "alt": "C"}],
    )

    assert captured["positions"] == {"chr1": {0, 2, 3}}
    assert result.manifest["memberships"]["unique_usable_allele_count"] == 2
    assert result.manifest["memberships"]["strategy_variant_membership_count"] == 3
    medians = result.distributions[result.distributions["quantile"] == 0.5].set_index(
        ["strategy", "gnomad_status"]
    )["phyloP100way"]
    assert medians.loc[("s1", "found")] == -1.0
    assert medians.loc[("s1", "not_found")] == 3.0
    assert result.manifest["memberships"]["lookup_failed_allele_context_count"] == 2
    found_histogram = result.histograms[
        result.histograms["strategy"].eq("s1") & result.histograms["gnomad_status"].eq("found")
    ]
    not_found_histogram = result.histograms[
        result.histograms["strategy"].eq("s1") & result.histograms["gnomad_status"].eq("not_found")
    ]
    assert found_histogram[["bin_left", "bin_right"]].values.tolist() == (
        not_found_histogram[["bin_left", "bin_right"]].values.tolist()
    )
    assert found_histogram["fraction"].sum() == 1.0
    groups = {(row["strategy"], row["gnomad_status"]): row for row in result.manifest["groups"]}
    assert groups[("s1", "found")]["median"] == -1.0
    assert groups[("s1", "not_found")]["median"] == 3.0

    monkeypatch.setattr(candidate, "read_position_scores", lambda **_kwargs: (_ for _ in ()).throw(AssertionError()))
    cached = candidate.build_candidate_conservation(
        variant_annotations_tsv=annotations,
        analytics_dir=tmp_path / "analytics",
        annotation_failures_tsv=failures,
    )
    assert cached.position_scores is None
    assert len(cached.distributions) == len(result.distributions)
    assert len(cached.histograms) == len(result.histograms)
