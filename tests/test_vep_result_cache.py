from __future__ import annotations

import multiprocessing
from pathlib import Path

import pandas as pd
import pytest

from analytics.vep.result_cache import VepResultCache


class _SynchronizedLookupCache(VepResultCache):
    def __init__(self, *args, barrier, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._barrier = barrier
        self._initial_lookup_complete = False

    def lookup(
        self,
        requests: pd.DataFrame,
    ) -> tuple[pd.DataFrame, dict[str, object]]:
        result = super().lookup(requests)
        if not self._initial_lookup_complete:
            self._initial_lookup_complete = True
            self._barrier.wait(timeout=20)
        return result


def _publish_worker(
    cache_dir: str,
    records: list[dict[str, object]],
    barrier,
    result_queue,
) -> None:
    try:
        cache = _SynchronizedLookupCache(
            cache_dir,
            config={"backend": "local", "release": "116"},
            barrier=barrier,
        )
        result_queue.put(("ok", cache.publish(pd.DataFrame.from_records(records))))
    except Exception as exc:
        result_queue.put(("error", type(exc).__name__, str(exc)))


def _concurrent_publish(
    cache_dir: Path,
    frames: list[pd.DataFrame],
) -> list[tuple]:
    context = multiprocessing.get_context("spawn")
    barrier = context.Barrier(len(frames))
    result_queue = context.Queue()
    processes = [
        context.Process(
            target=_publish_worker,
            args=(
                str(cache_dir),
                frame.to_dict(orient="records"),
                barrier,
                result_queue,
            ),
        )
        for frame in frames
    ]
    for process in processes:
        process.start()
    for process in processes:
        process.join(timeout=30)
        if process.is_alive():
            process.terminate()
            process.join()
            pytest.fail("Concurrent VEP cache publish did not finish")
        assert process.exitcode == 0
    return [result_queue.get(timeout=5) for _ in processes]


def _annotations() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "variant_key": "1:10:A>G",
                "gene_id": "1",
                "status": "ok",
                "primary_consequence": "missense_variant",
                "consequence_terms": "missense_variant",
                "transcript_id": "NM_1",
                "mane_select": "NM_1",
                "canonical": True,
                "impact": "MODERATE",
                "variant_class": "SNV",
            },
            {
                "variant_key": "1:1000010:C>T",
                "gene_id": "2",
                "status": "no_target_gene",
                "primary_consequence": "",
                "consequence_terms": "",
                "transcript_id": "",
                "mane_select": "",
                "canonical": False,
                "impact": "",
                "variant_class": "",
            },
            {
                "variant_key": "2:20:G>A",
                "gene_id": "3",
                "status": "no_response",
                "primary_consequence": "",
                "consequence_terms": "",
                "transcript_id": "",
                "mane_select": "",
                "canonical": False,
                "impact": "",
                "variant_class": "",
            },
        ]
    )


def _requests() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "variant_key": "1:10:A>G",
                "gene_id": "1",
                "chrom": "1",
                "pos": 10,
                "ref": "A",
                "alt": "G",
            },
            {
                "variant_key": "1:1000010:C>T",
                "gene_id": "2",
                "chrom": "1",
                "pos": 1000010,
                "ref": "C",
                "alt": "T",
            },
            {
                "variant_key": "2:20:G>A",
                "gene_id": "3",
                "chrom": "2",
                "pos": 20,
                "ref": "G",
                "alt": "A",
            },
        ]
    )


def test_cache_publishes_terminal_results_and_reuses_them(tmp_path: Path) -> None:
    cache = VepResultCache(
        tmp_path / "cache",
        config={"backend": "local", "release": "116", "refseq": True},
    )

    first = cache.publish(_annotations())
    hits, lookup = cache.lookup(_requests())
    second = cache.publish(_annotations())

    assert first["accepted_count"] == 2
    assert first["skipped_count"] == 1
    assert first["published_count"] == 2
    assert first["fragment_count"] == 2
    assert lookup["hit_count"] == 2
    assert lookup["miss_count"] == 1
    assert hits["variant_key"].tolist() == ["1:1000010:C>T", "1:10:A>G"]
    assert second["existing_count"] == 2
    assert second["published_count"] == 0
    assert len(list(cache.namespace_dir.rglob("part-*.parquet"))) == 2


def test_cache_rejects_conflicting_result(tmp_path: Path) -> None:
    cache = VepResultCache(
        tmp_path / "cache",
        config={"backend": "local", "release": "116"},
    )
    rows = _annotations().iloc[[0]].copy()
    cache.publish(rows)
    rows.loc[:, "primary_consequence"] = "stop_gained"
    rows.loc[:, "consequence_terms"] = "stop_gained"

    with pytest.raises(ValueError, match="Conflicting VEP result"):
        cache.publish(rows)


def test_lookup_ignores_unpublished_temporary_files(tmp_path: Path) -> None:
    cache = VepResultCache(
        tmp_path / "cache",
        config={"backend": "local", "release": "116"},
    )
    cache.publish(_annotations().iloc[[0]])
    tile_dir = next(path.parent for path in cache.namespace_dir.rglob("part-*.parquet"))
    (tile_dir / ".interrupted.parquet.tmp").write_bytes(b"incomplete")

    hits, summary = cache.lookup(_requests().iloc[[0]])

    assert summary["hit_count"] == 1
    assert hits.loc[0, "primary_consequence"] == "missense_variant"


def test_cache_validates_request_coordinates(tmp_path: Path) -> None:
    cache = VepResultCache(
        tmp_path / "cache",
        config={"backend": "local", "release": "116"},
    )
    requests = _requests().iloc[[0]].copy()
    requests.loc[:, "pos"] = 11

    with pytest.raises(ValueError, match="coordinates do not match"):
        cache.lookup(requests)


def test_tile_size_creates_an_independent_namespace(tmp_path: Path) -> None:
    one_mb = VepResultCache(
        tmp_path / "cache",
        config={"backend": "local", "release": "116"},
        tile_size_bp=1_000_000,
    )
    five_mb = VepResultCache(
        tmp_path / "cache",
        config={"backend": "local", "release": "116"},
        tile_size_bp=5_000_000,
    )

    one_mb.publish(_annotations().iloc[[0]])
    hits, summary = five_mb.lookup(_requests().iloc[[0]])

    assert hits.empty
    assert summary["miss_count"] == 1
    assert one_mb.namespace_dir != five_mb.namespace_dir


def test_concurrent_identical_publications_create_one_fragment(tmp_path: Path) -> None:
    rows = _annotations().iloc[[0]]

    results = _concurrent_publish(tmp_path / "cache", [rows, rows])

    assert [result[0] for result in results].count("ok") == 2
    summaries = [result[1] for result in results]
    assert sorted(summary["published_count"] for summary in summaries) == [0, 1]
    assert sorted(summary["existing_count"] for summary in summaries) == [0, 1]
    cache = VepResultCache(
        tmp_path / "cache",
        config={"backend": "local", "release": "116"},
    )
    hits, summary = cache.lookup(_requests().iloc[[0]])
    assert summary["hit_count"] == 1
    assert hits.loc[0, "primary_consequence"] == "missense_variant"
    assert len(list(cache.namespace_dir.rglob("part-*.parquet"))) == 1
    assert not list(cache.namespace_dir.rglob("*.tmp"))


def test_concurrent_conflicting_publication_keeps_cache_consistent(
    tmp_path: Path,
) -> None:
    first = _annotations().iloc[[0]].copy()
    conflicting = first.copy()
    conflicting.loc[:, "primary_consequence"] = "stop_gained"
    conflicting.loc[:, "consequence_terms"] = "stop_gained"

    results = _concurrent_publish(tmp_path / "cache", [first, conflicting])

    assert [result[0] for result in results].count("ok") == 1
    errors = [result for result in results if result[0] == "error"]
    assert len(errors) == 1
    assert errors[0][1:] == (
        "ValueError",
        "Conflicting VEP result for cached key 1:10:A>G / gene 1",
    )
    cache = VepResultCache(
        tmp_path / "cache",
        config={"backend": "local", "release": "116"},
    )
    hits, summary = cache.lookup(_requests().iloc[[0]])
    assert summary["hit_count"] == 1
    assert hits.loc[0, "primary_consequence"] in {
        "missense_variant",
        "stop_gained",
    }
    assert len(list(cache.namespace_dir.rglob("part-*.parquet"))) == 1
    assert not list(cache.namespace_dir.rglob("*.tmp"))
