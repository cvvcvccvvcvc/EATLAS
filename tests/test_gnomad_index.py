from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from genomics.gnomad_cache import GnomadRegionCache
from genomics.gnomad_index import FRAGMENT_TILE_COUNT, GnomadAlleleIndex


def _variant(pos: int) -> dict:
    return {
        "variant_id": f"1-{pos}-A-G",
        "chrom": "1",
        "pos": pos,
        "ref": "A",
        "alt": "G",
        "consequence": "intron_variant",
        "exome": {"af": 0.03},
        "genome": {"af": 0.02},
        "joint": {"an": 100, "ac": [4]},
    }


def _requests() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "variant_key": "1:100:A>G",
                "chrom": "1",
                "pos": 100,
                "ref": "A",
                "alt": "G",
            },
            {
                "variant_key": "1:101:A>T",
                "chrom": "1",
                "pos": 101,
                "ref": "A",
                "alt": "T",
            },
        ]
    )


def test_allele_index_builds_from_json_and_reuses_parquet(tmp_path: Path) -> None:
    region_cache = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: [_variant(100)],
    )
    region_cache.fetch_region("1", 100, 101)
    index = GnomadAlleleIndex(tmp_path, region_cache=region_cache)

    first, unresolved, first_summary = index.lookup(_requests())
    first = first.set_index("variant_key")

    assert unresolved.empty
    assert first_summary["tile_build_count"] == 1
    assert first_summary["fragment_build_count"] == 1
    assert first_summary["indexed_variant_count"] == 1
    assert first_summary["raw_tile_missing_count"] == 0
    assert first_summary["observation_window"] is not None
    assert first.loc["1:100:A>G", "gnomad_status"] == "ok"
    assert bool(first.loc["1:100:A>G", "gnomad_found"])
    assert first.loc["1:100:A>G", "gnomad_af"] == 0.04
    assert not bool(first.loc["1:101:A>T", "gnomad_found"])

    warm_region_cache = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("warm allele index must not use the network")
        ),
    )
    warm_index = GnomadAlleleIndex(tmp_path, region_cache=warm_region_cache)
    second, second_unresolved, second_summary = warm_index.lookup(_requests())

    pd.testing.assert_frame_equal(
        first.reset_index().sort_values("variant_key").reset_index(drop=True),
        second.sort_values("variant_key").reset_index(drop=True),
        check_dtype=False,
    )
    assert second_unresolved.empty
    assert second_summary["tile_hit_count"] == 1
    assert second_summary["tile_build_count"] == 0
    assert second_summary["observation_window"] == first_summary["observation_window"]
    assert warm_region_cache.snapshot()["fetch_batch_count"] == 0
    assert warm_region_cache.snapshot()["tile_hit_count"] == 0


def test_empty_complete_tile_proves_allele_absence(tmp_path: Path) -> None:
    region_cache = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: [],
    )
    region_cache.fetch_region("1", 100, 101)
    index = GnomadAlleleIndex(tmp_path, region_cache=region_cache)

    evidence, unresolved, _summary = index.lookup(_requests())

    assert unresolved.empty
    assert evidence["gnomad_status"].eq("ok").all()
    assert not evidence["gnomad_found"].any()


def test_incremental_chromosome_fragments_do_not_duplicate_coverage(
    tmp_path: Path,
) -> None:
    def fetch(_chrom, start, _end, **_kwargs):
        position = 100 if start == 1 else 30_000
        return [_variant(position)]

    region_cache = GnomadRegionCache(tmp_path, fetcher=fetch)
    region_cache.fetch_region("1", 100, 100)
    index = GnomadAlleleIndex(tmp_path, region_cache=region_cache)
    index.lookup(_requests().iloc[[0]])

    region_cache.fetch_region("1", 30_000, 30_000)
    expanded = pd.concat(
        [
            _requests().iloc[[0]],
            pd.DataFrame(
                [
                    {
                        "variant_key": "1:30000:A>G",
                        "chrom": "1",
                        "pos": 30_000,
                        "ref": "A",
                        "alt": "G",
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    evidence, unresolved, summary = index.lookup(expanded)

    assert unresolved.empty
    assert evidence["gnomad_found"].all()
    assert summary["requested_tile_count"] == 2
    assert summary["tile_hit_count"] == 1
    assert summary["tile_build_count"] == 1
    assert len(index.fragment_paths()) == 2


def test_large_tile_set_is_written_in_bounded_fragments(tmp_path: Path) -> None:
    positions = [
        index * 25_000 + 1
        for index in range(FRAGMENT_TILE_COUNT + 1)
    ]
    region_cache = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: [_variant(pos) for pos in positions],
    )
    region_cache.fetch_region("1", positions[0], positions[-1])
    requests = pd.DataFrame(
        [
            {
                "variant_key": f"1:{pos}:A>G",
                "chrom": "1",
                "pos": pos,
                "ref": "A",
                "alt": "G",
            }
            for pos in positions
        ]
    )
    index = GnomadAlleleIndex(tmp_path, region_cache=region_cache)

    evidence, unresolved, summary = index.lookup(requests)

    assert unresolved.empty
    assert evidence["gnomad_found"].all()
    assert summary["tile_build_count"] == len(positions)
    assert summary["fragment_count"] == 2
    assert len(index.fragment_paths()) == 2


def test_corrupt_index_fails_instead_of_claiming_absence(tmp_path: Path) -> None:
    region_cache = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: [_variant(100)],
    )
    region_cache.fetch_region("1", 100, 101)
    index = GnomadAlleleIndex(tmp_path, region_cache=region_cache)
    index.lookup(_requests())
    fragment_path = index.fragment_paths()[0]
    fragment_path.write_bytes(b"not parquet")

    with pytest.raises(ValueError, match="gnomAD allele-index"):
        index.lookup(_requests())


def test_modified_coverage_manifest_fails_instead_of_claiming_absence(
    tmp_path: Path,
) -> None:
    region_cache = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: [_variant(100)],
    )
    region_cache.fetch_region("1", 100, 101)
    index = GnomadAlleleIndex(tmp_path, region_cache=region_cache)
    index.lookup(_requests())
    manifest_path = next(index.namespace_dir.glob("chrom=*/fragment-*.json"))
    manifest = json.loads(manifest_path.read_text())
    manifest["tiles"][0]["observation_window"]["finished_at_utc"] = (
        "2099-01-01T00:00:00+00:00"
    )
    manifest_path.write_text(json.dumps(manifest))

    with pytest.raises(ValueError, match="manifest identity changed"):
        index.lookup(_requests())
