from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

from analytics.analyses import candidate_conservation as candidate
from analytics.analyses.candidate_conservation_aggregation import (
    build_candidate_allele_store,
    resolve_candidate_aggregation_source,
)
from analytics.analyses.conservation import PositionScores, parse_tracks, track_identity


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
        {
            "variant_key": "1:4:A>A.",
            "lookup_status": "raw_no_context",
            "strategies": "s1",
            "gnomad_af": "0.2",
        },
        {
            "variant_key": "1:5:A>G:extra",
            "lookup_status": "ok",
            "strategies": "s1",
            "gnomad_af": "0.3",
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
    local_bigwig = tmp_path / "hg38.phyloP100way.bw"
    local_bigwig.write_bytes(b"test bigwig identity")
    track = parse_tracks("phyloP100way", phylop_bigwig=local_bigwig)[0]

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
        variant_annotations_source=annotations,
        analytics_dir=tmp_path / "analytics",
        annotation_failures_tsv=failures,
        additional_rows=[{"chrom": "1", "pos": "4", "ref": "T", "alt": "C"}],
        chunk_size=1,
        phylop_bigwig=local_bigwig,
    )

    assert captured["positions"] == {"chr1": {0, 2, 3}}
    assert result.manifest["candidate_scan"]["unsupported_allele_context_count"] == 2
    assert result.manifest["memberships"]["unique_usable_allele_count"] == 2
    assert result.manifest["memberships"]["strategy_variant_membership_count"] == 3
    assert result.manifest["score_materialization"] == {
        "attempted_allele_count": 2,
        "scored_allele_count": 2,
        "missing_score_allele_count": 0,
    }
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
        variant_annotations_source=annotations,
        analytics_dir=tmp_path / "analytics",
        annotation_failures_tsv=failures,
        phylop_bigwig=local_bigwig,
    )
    assert cached.position_scores is None
    assert len(cached.distributions) == len(result.distributions)
    assert len(cached.histograms) == len(result.histograms)
    assert result.manifest["aggregation"]["engine"] == "duckdb"
    assert result.manifest["inputs"]["track"] == track_identity(track)
    assert not list((tmp_path / "analytics").glob("*.sqlite3"))


def test_candidate_conservation_scores_indels_and_excludes_complex_alleles(
    tmp_path: Path,
    monkeypatch,
) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    fields = sorted(candidate.REQUIRED_COLUMNS)
    rows = [
        {
            "variant_key": "1:10:ATC>A",
            "lookup_status": "ok",
            "strategies": "s1",
            "gnomad_af": "",
        },
        {
            "variant_key": "1:20:A>AGG",
            "lookup_status": "ok",
            "strategies": "s1",
            "gnomad_af": "",
        },
        {
            "variant_key": "1:30:AC>GT",
            "lookup_status": "ok",
            "strategies": "s1",
            "gnomad_af": "",
        },
    ]
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    captured = {}
    track = parse_tracks("phyloP100way")[0]

    def fake_read_position_scores(**kwargs):
        captured["positions"] = kwargs["positions_by_chrom"]
        return PositionScores(
            track,
            {
                ("chr1", 10): 2.0,
                ("chr1", 11): 4.0,
                ("chr1", 19): 1.0,
                ("chr1", 20): 5.0,
            },
            {"track": track.name, "status": "complete", "failed_block_count": 0},
        )

    monkeypatch.setattr(candidate, "read_position_scores", fake_read_position_scores)
    result = candidate.build_candidate_conservation(
        variant_annotations_source=annotations,
        analytics_dir=tmp_path / "analytics",
    )

    assert captured["positions"] == {"chr1": {10, 11, 19, 20}}
    assert result.manifest["candidate_scan"]["unsupported_allele_context_count"] == 1
    assert result.manifest["memberships"]["unique_usable_allele_count"] == 2
    group = result.manifest["groups"][0]
    assert group["variant_count"] == 2
    assert group["scored_count"] == 2
    assert group["median"] == 3.0


def test_candidate_score_materialization_reuses_scores_and_tracks_missing_values(
    tmp_path: Path,
    monkeypatch,
) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    fields = sorted(candidate.REQUIRED_COLUMNS)
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            [
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
            ]
        )

    track = parse_tracks("phyloP100way")[0]
    monkeypatch.setattr(
        candidate,
        "read_position_scores",
        lambda **_kwargs: PositionScores(
            track,
            {("chr1", 0): 1.25},
            {"track": track.name, "status": "complete", "failed_block_count": 0},
        ),
    )

    result = candidate.build_candidate_conservation(
        variant_annotations_source=annotations,
        analytics_dir=tmp_path / "analytics",
        chunk_size=1,
    )

    assert result.manifest["score_materialization"] == {
        "attempted_allele_count": 2,
        "scored_allele_count": 1,
        "missing_score_allele_count": 1,
    }
    groups = {(row["strategy"], row["gnomad_status"]): row for row in result.manifest["groups"]}
    assert groups[("s1", "found")]["scored_count"] == 1
    assert groups[("s2", "found")]["scored_count"] == 1
    assert groups[("s1", "not_found")]["scored_count"] == 0


def test_candidate_source_reads_validated_pipeline_partitions(tmp_path: Path) -> None:
    artifact = tmp_path / "variant_annotations"
    inputs = artifact / "partitions"
    columns = sorted(candidate.REQUIRED_COLUMNS)
    entries = []
    for index, variant_key in enumerate(["1:1:A>G", "1:2:C>T"], start=1):
        partition_id = f"partition_{index:06d}"
        path = inputs / partition_id / "shard_000001.tsv.gz"
        path.parent.mkdir(parents=True)
        with gzip.open(path, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(
                {
                    "variant_key": variant_key,
                    "lookup_status": "ok",
                    "strategies": "s1",
                    "gnomad_af": "",
                }
            )
        entries.append({
            "partition_id": partition_id,
            "shard_count": 1,
            "row_count": 1,
            "shards": [{
                "shard_id": "shard_000001",
                "path": f"partitions/{partition_id}/shard_000001.tsv.gz",
                "row_count": 1,
                "size_bytes": path.stat().st_size,
            }],
        })
    manifest = artifact / "manifest.json"
    manifest.write_text(json.dumps({
        "schema": "gaph_variant_annotation_dataset_v1",
        "status": "complete",
        "layout": "partitioned",
        "format": "tsv_gzip_v1",
        "partition_count": 2,
        "shard_count": 2,
        "row_count": 2,
        "fields": columns,
        "partitions": entries,
    }))

    source = resolve_candidate_aggregation_source(manifest)

    assert source.mode == "partitioned"
    assert source.header is True
    assert source.row_count == 2
    assert source.paths == tuple(
        (inputs / f"partition_{index:06d}" / "shard_000001.tsv.gz").resolve()
        for index in (1, 2)
    )

    store = build_candidate_allele_store(
        variant_annotations_source=manifest,
        strategies=["s1"],
        annotation_failures_path=None,
        temp_dir=tmp_path / "duckdb_tmp",
    )
    try:
        assert store.summary()["variant_context_row_count"] == 2
        assert store.summary()["unique_usable_allele_count"] == 2
    finally:
        store.close()


def test_candidate_store_keeps_gnomad_status_per_strategy(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    columns = sorted(candidate.REQUIRED_COLUMNS)
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            [
                {
                    "variant_key": "chr1:8:A>G",
                    "lookup_status": "ok",
                    "strategies": "found_strategy",
                    "gnomad_af": "0.1",
                },
                {
                    "variant_key": "chr1:8:A>G",
                    "lookup_status": "ok",
                    "strategies": "not_found_strategy",
                    "gnomad_af": "",
                },
            ]
        )

    store = build_candidate_allele_store(
        variant_annotations_source=annotations,
        strategies=["found_strategy", "not_found_strategy"],
        annotation_failures_path=None,
        temp_dir=tmp_path / "duckdb_tmp",
    )
    try:
        groups = {
            (row.strategy, row.gnomad_status): int(row.variant_count)
            for row in store.group_counts().itertuples(index=False)
        }
        assert groups == {
            ("found_strategy", "found"): 1,
            ("not_found_strategy", "not_found"): 1,
        }
        assert store.summary()["gnomad_status_conflict_membership_count"] == 0
    finally:
        store.close()
