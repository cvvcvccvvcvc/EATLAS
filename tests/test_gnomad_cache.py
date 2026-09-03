from __future__ import annotations

import gzip
import json
import multiprocessing
from pathlib import Path
from urllib.error import HTTPError

import pytest

from genomics import gnomad_cache as gnomad_cache_module
from genomics.gnomad_cache import GnomadRegionCache, tiles_for_region


def variant(pos: int) -> dict:
    return {
        "variant_id": f"1-{pos}-A-G",
        "chrom": "1",
        "pos": pos,
        "ref": "A",
        "alt": "G",
    }


class _SynchronizedRegionCache(GnomadRegionCache):
    def __init__(self, *args, initial_tile_count: int, barrier, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._initial_tile_count = initial_tile_count
        self._initial_reads = 0
        self._barrier = barrier

    def read_cached_tile(self, tile):
        result = super().read_cached_tile(tile)
        self._initial_reads += 1
        if self._initial_reads == self._initial_tile_count:
            self._barrier.wait(timeout=20)
        return result


def _concurrent_fetch_worker(
    cache_dir: str,
    start: int,
    end: int,
    initial_tile_count: int,
    barrier,
    result_queue,
) -> None:
    calls: list[tuple[int, int]] = []

    def fetch(_chrom, fetch_start, fetch_end, *, max_attempts):
        calls.append((fetch_start, fetch_end))
        return [
            variant(position)
            for position in (100, 25_100, 50_100)
            if fetch_start <= position <= fetch_end
        ]

    try:
        cache = _SynchronizedRegionCache(
            cache_dir,
            fetcher=fetch,
            initial_tile_count=initial_tile_count,
            barrier=barrier,
        )
        records = cache.fetch_region("1", start, end)
        result_queue.put(("ok", calls, [record["pos"] for record in records]))
    except Exception as exc:
        result_queue.put(("error", type(exc).__name__, str(exc)))


def _concurrent_fetches(
    cache_dir: Path,
    regions: list[tuple[int, int]],
) -> list[tuple]:
    context = multiprocessing.get_context("spawn")
    barrier = context.Barrier(len(regions))
    result_queue = context.Queue()
    processes = [
        context.Process(
            target=_concurrent_fetch_worker,
            args=(
                str(cache_dir),
                start,
                end,
                len(tiles_for_region("1", start, end)),
                barrier,
                result_queue,
            ),
        )
        for start, end in regions
    ]
    for process in processes:
        process.start()
    for process in processes:
        process.join(timeout=30)
        if process.is_alive():
            process.terminate()
            process.join()
            pytest.fail("Concurrent gnomAD cache fetch did not finish")
        assert process.exitcode == 0
    return [result_queue.get(timeout=5) for _ in processes]


def test_tile_boundaries_are_one_based_and_inclusive() -> None:
    assert tiles_for_region("chr1", 25_000, 25_001) == [
        tiles_for_region("1", 1, 1)[0],
        tiles_for_region("1", 25_001, 25_001)[0],
    ]


def test_cold_grouped_fetch_populates_tiles_and_warm_read_uses_cache(
    tmp_path: Path,
    monkeypatch,
) -> None:
    calls = []
    timestamps = iter(["2026-03-01T10:00:00+00:00", "2026-03-01T10:00:02+00:00"])
    monkeypatch.setattr(gnomad_cache_module, "utc_now", lambda: next(timestamps))

    def fetch(chrom, start, end, *, max_attempts):
        calls.append((chrom, start, end, max_attempts))
        return [variant(10_000), variant(30_000), variant(60_000)]

    cache = GnomadRegionCache(tmp_path, fetcher=fetch)
    records = cache.fetch_region("1", 10_000, 60_000)

    assert [record["pos"] for record in records] == [10_000, 30_000, 60_000]
    assert calls == [("1", 1, 75_000, 2)]
    assert cache.snapshot()["tile_write_count"] == 3
    expected_window = {
        "started_at_utc": "2026-03-01T10:00:00+00:00",
        "finished_at_utc": "2026-03-01T10:00:02+00:00",
    }
    assert cache.snapshot()["observation_window"] == expected_window
    assert len(list(tmp_path.rglob("*.json.gz"))) == 3

    def unexpected_fetch(*_args, **_kwargs):
        raise AssertionError("warm cache must not use the network")

    warm_cache = GnomadRegionCache(tmp_path, fetcher=unexpected_fetch)
    warm_records = warm_cache.fetch_region("1", 10_000, 60_000)

    assert warm_records == records
    assert warm_cache.snapshot()["tile_hit_count"] == 3
    assert warm_cache.snapshot()["fetch_batch_count"] == 0
    assert warm_cache.snapshot()["observation_window"] == expected_window


def test_concurrent_cold_tile_is_fetched_once(tmp_path: Path) -> None:
    results = _concurrent_fetches(tmp_path / "cache", [(1, 25_000), (1, 25_000)])

    assert all(result[0] == "ok" for result in results)
    assert sum(len(result[1]) for result in results) == 1
    assert all(result[2] == [100] for result in results)
    assert len(list((tmp_path / "cache").rglob("*.json.gz"))) == 1


def test_concurrent_overlapping_groups_fetch_disjoint_regions(tmp_path: Path) -> None:
    results = _concurrent_fetches(
        tmp_path / "cache",
        [(1, 50_000), (25_001, 75_000)],
    )

    assert all(result[0] == "ok" for result in results)
    calls = sorted(call for result in results for call in result[1])
    assert len(calls) == 2
    assert calls[0][0] == 1
    assert calls[0][1] + 1 == calls[1][0]
    assert calls[1][1] == 75_000
    assert sorted(result[2] for result in results) == [
        [100, 25_100],
        [25_100, 50_100],
    ]
    assert len(list((tmp_path / "cache").rglob("*.json.gz"))) == 3


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


def test_disabled_cache_preserves_direct_fetch_behavior(monkeypatch) -> None:
    calls = []
    timestamps = iter(["2026-03-02T10:00:00+00:00", "2026-03-02T10:00:01+00:00"])
    monkeypatch.setattr(gnomad_cache_module, "utc_now", lambda: next(timestamps))

    def fetch(chrom, start, end, *, max_attempts):
        calls.append((chrom, start, end, max_attempts))
        return [variant(100)]

    cache = GnomadRegionCache(None, fetcher=fetch)

    assert cache.fetch_region("1", 90, 110) == [variant(100)]
    assert calls == [("1", 90, 110, 10)]
    assert cache.snapshot()["enabled"] is False
    assert cache.snapshot()["observation_window"] == {
        "started_at_utc": "2026-03-02T10:00:00+00:00",
        "finished_at_utc": "2026-03-02T10:00:01+00:00",
    }


def test_mixed_warm_and_fresh_tiles_merge_observation_windows(
    tmp_path: Path,
    monkeypatch,
) -> None:
    timestamps = iter(
        [
            "2026-03-01T10:00:00+00:00",
            "2026-03-01T10:00:01+00:00",
            "2026-03-02T11:00:00+00:00",
            "2026-03-02T11:00:01+00:00",
        ]
    )
    monkeypatch.setattr(gnomad_cache_module, "utc_now", lambda: next(timestamps))
    GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: [variant(100)],
    ).fetch_region("1", 1, 25_000)

    calls = []
    cache = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: calls.append(True) or [variant(30_000)],
    )
    records = cache.fetch_region("1", 1, 50_000)

    assert [record["pos"] for record in records] == [100, 30_000]
    assert calls == [True]
    assert cache.snapshot()["observation_window"] == {
        "started_at_utc": "2026-03-01T10:00:00+00:00",
        "finished_at_utc": "2026-03-02T11:00:01+00:00",
    }


def test_tile_without_observation_window_is_refetched(tmp_path: Path) -> None:
    initial = GnomadRegionCache(tmp_path, fetcher=lambda *_args, **_kwargs: [])
    initial.fetch_region("1", 1, 25_000)
    tile_path = next(tmp_path.rglob("*.json.gz"))

    with gzip.open(tile_path, "rt", encoding="utf-8") as handle:
        payload = json.load(handle)
    del payload["observation_window"]
    with gzip.open(tile_path, "wt", encoding="utf-8") as handle:
        json.dump(payload, handle)

    repaired = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: [variant(100)],
    )
    assert [item["pos"] for item in repaired.fetch_region("1", 1, 25_000)] == [100]
    assert repaired.snapshot()["corrupt_tile_count"] == 1


def test_schema_v1_cache_is_not_reused(tmp_path: Path) -> None:
    legacy = (
        tmp_path
        / "gnomad_r4"
        / "GRCh38"
        / "schema_v1"
        / "tiles_25000bp"
        / "1"
        / "000000000001-000000025000.json.gz"
    )
    legacy.parent.mkdir(parents=True)
    legacy.write_bytes(b"legacy")
    calls = []
    cache = GnomadRegionCache(
        tmp_path,
        fetcher=lambda *_args, **_kwargs: calls.append(True) or [],
    )

    assert cache.fetch_region("1", 1, 25_000) == []
    assert calls == [True]
    assert "schema_v2" in str(cache.namespace_dir)
