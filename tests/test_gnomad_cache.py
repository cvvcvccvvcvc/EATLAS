from __future__ import annotations

from pathlib import Path
from urllib.error import HTTPError

from bin.gnomad_cache import GnomadRegionCache, tiles_for_region


def variant(pos: int) -> dict:
    return {
        "variant_id": f"1-{pos}-A-G",
        "chrom": "1",
        "pos": pos,
        "ref": "A",
        "alt": "G",
    }


def test_tile_boundaries_are_one_based_and_inclusive() -> None:
    assert tiles_for_region("chr1", 25_000, 25_001) == [
        tiles_for_region("1", 1, 1)[0],
        tiles_for_region("1", 25_001, 25_001)[0],
    ]


def test_cold_grouped_fetch_populates_tiles_and_warm_read_uses_cache(tmp_path: Path) -> None:
    calls = []

    def fetch(chrom, start, end, *, max_attempts):
        calls.append((chrom, start, end, max_attempts))
        return [variant(10_000), variant(30_000), variant(60_000)]

    cache = GnomadRegionCache(tmp_path, fetcher=fetch)
    records = cache.fetch_region("1", 10_000, 60_000)

    assert [record["pos"] for record in records] == [10_000, 30_000, 60_000]
    assert calls == [("1", 1, 75_000, 2)]
    assert cache.snapshot()["tile_write_count"] == 3
    assert len(list(tmp_path.rglob("*.json.gz"))) == 3

    def unexpected_fetch(*_args, **_kwargs):
        raise AssertionError("warm cache must not use the network")

    warm_cache = GnomadRegionCache(tmp_path, fetcher=unexpected_fetch)
    warm_records = warm_cache.fetch_region("1", 10_000, 60_000)

    assert warm_records == records
    assert warm_cache.snapshot()["tile_hit_count"] == 3
    assert warm_cache.snapshot()["fetch_batch_count"] == 0


def test_timeout_splits_group_until_tiles_succeed(tmp_path: Path) -> None:
    calls = []

    def fetch(chrom, start, end, *, max_attempts):
        calls.append((chrom, start, end, max_attempts))
        if end - start + 1 > 25_000:
            raise TimeoutError("timed out")
        return [variant(start)]

    cache = GnomadRegionCache(tmp_path, fetcher=fetch)
    records = cache.fetch_region("1", 1, 50_000)

    assert [record["pos"] for record in records] == [1, 25_001]
    assert calls == [
        ("1", 1, 50_000, 2),
        ("1", 1, 25_000, 10),
        ("1", 25_001, 50_000, 10),
    ]
    assert cache.snapshot()["split_count"] == 1
    assert cache.snapshot()["tile_write_count"] == 2


def test_transient_http_error_splits_group_instead_of_retrying_large_region(
    tmp_path: Path,
) -> None:
    calls = []

    def fetch(chrom, start, end, *, max_attempts):
        calls.append((chrom, start, end, max_attempts))
        if end - start + 1 > 25_000:
            raise HTTPError("", 500, "Internal Server Error", None, None)
        return [variant(start)]

    cache = GnomadRegionCache(tmp_path, fetcher=fetch)
    records = cache.fetch_region("1", 1, 50_000)

    assert [record["pos"] for record in records] == [1, 25_001]
    assert calls == [
        ("1", 1, 50_000, 2),
        ("1", 1, 25_000, 10),
        ("1", 25_001, 50_000, 10),
    ]
    assert cache.snapshot()["split_count"] == 1
    assert cache.snapshot()["tile_write_count"] == 2


def test_empty_and_corrupt_tiles_are_handled_safely(tmp_path: Path) -> None:
    empty_cache = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: [],
    )
    assert empty_cache.fetch_region("1", 1, 25_000) == []

    warm_cache = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("empty tiles must be cached")
        ),
    )
    assert warm_cache.fetch_region("1", 1, 25_000) == []

    tile_path = next(tmp_path.rglob("*.json.gz"))
    tile_path.write_bytes(b"not gzip")
    repaired_cache = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: [variant(100)],
    )
    assert [record["pos"] for record in repaired_cache.fetch_region("1", 1, 25_000)] == [100]
    assert repaired_cache.snapshot()["corrupt_tile_count"] == 1
    assert repaired_cache.snapshot()["tile_write_count"] == 1


def test_disabled_cache_preserves_direct_fetch_behavior() -> None:
    calls = []

    def fetch(chrom, start, end, *, max_attempts):
        calls.append((chrom, start, end, max_attempts))
        return [variant(100)]

    cache = GnomadRegionCache(None, fetcher=fetch)

    assert cache.fetch_region("1", 90, 110) == [variant(100)]
    assert calls == [("1", 90, 110, 10)]
    assert cache.snapshot()["enabled"] is False
