from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path

import pandas as pd
import pytest

from analytics.analyses import matched_control as controls
from analytics.analyses.observed_variant_store import (
    ALLELE_COLUMNS,
    ALLELE_GENE_COLUMNS,
    build_or_load_observed_variant_store,
)
from analytics.analyses.target_context import context_at
from genomics.variants import changed_target_position, parse_variant_key


FIELDS = [
    "variant_key",
    "gene_id",
    "event_type",
    "ref",
    "alt",
    "lookup_status",
    "strategies",
]


def _write_annotations(path: Path, rows: list[dict[str, str]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _expected_md5_sample(
    path: Path,
    contexts: dict[str, list[tuple[int, int, str]]],
    genes: dict[str, dict[str, object]],
    strategies: list[str],
    limit: int,
    seed: int,
) -> pd.DataFrame:
    strategy_set = set(strategies)
    ranked: dict[str, list[tuple[str, str, dict[str, object]]]] = {
        strategy: [] for strategy in strategies
    }
    frame = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
    frame = frame[
        frame["event_type"].eq("snv")
        & frame["ref"].str.len().eq(1)
        & frame["alt"].str.len().eq(1)
        & frame["ref"].str.upper().isin({"A", "C", "G", "T"})
        & frame["alt"].str.upper().isin({"A", "C", "G", "T"})
        & frame["lookup_status"].eq("ok")
    ]
    for row in frame.itertuples(index=False):
        gene_id = str(row.gene_id)
        parsed = parse_variant_key(row.variant_key)
        gene = genes.get(gene_id)
        if parsed is None or gene is None:
            continue
        chrom, pos, ref, alt = parsed
        target_pos = changed_target_position(parsed, int(gene["begin"]))
        record_base = {
            "gene_id": gene_id,
            "variant_key": str(row.variant_key),
            "target_pos": target_pos,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "context": context_at(contexts.get(gene_id, []), target_pos),
        }
        for strategy in (
            item.strip() for item in str(row.strategies).split(",") if item.strip()
        ):
            if strategy_set and strategy not in strategy_set:
                continue
            record = {**record_base, "strategy": strategy}
            token = f"{gene_id}:{record['variant_key']}"
            payload = f"{seed}|{strategy}|{token}".encode()
            rank = hashlib.md5(payload, usedforsecurity=False).hexdigest()
            ranked[strategy].append((rank, token, record))

    rows = [
        record
        for strategy in strategies
        for _rank, _token, record in sorted(
            ranked[strategy],
            key=lambda item: (item[0], item[1]),
        )[:limit]
    ]
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result = result.sort_values(
        ["strategy", "gene_id", "variant_key"],
        kind="mergesort",
    ).reset_index(drop=True)
    result.insert(0, "focal_id", [f"focal_{index:09d}" for index in range(len(result))])
    return result


def test_observed_store_reuses_cache_and_queries_strategy_memberships(
    tmp_path: Path,
) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    rows = [
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "event_type": "snv",
            "ref": "A",
            "alt": "G",
            "lookup_status": "ok",
            "strategies": "s1,s2",
        },
        {
            "variant_key": "1:101:C>T",
            "gene_id": "1",
            "event_type": "snv",
            "ref": "C",
            "alt": "T",
            "lookup_status": "ok",
            "strategies": "s2",
        },
        {
            "variant_key": "1:100:A>G",
            "gene_id": "2",
            "event_type": "snv",
            "ref": "A",
            "alt": "G",
            "lookup_status": "ok",
            "strategies": "s3",
        },
        {
            "variant_key": "",
            "gene_id": "2",
            "event_type": "snv",
            "ref": "N",
            "alt": "G",
            "lookup_status": "non_concrete_allele",
            "strategies": "s1",
        },
    ]
    _write_annotations(annotations, rows)

    store = build_or_load_observed_variant_store(
        variant_annotations_tsv=annotations,
        analytics_dir=tmp_path / "analytics",
        strategies=["s1", "s2", "s3"],
    )

    assert store.cache_hit is False
    assert store.strategies == ("s1", "s2", "s3")
    assert store.manifest["source_row_count"] == 4
    assert store.manifest["allele_gene_count"] == 4
    assert store.manifest["allele_count"] == 2
    assert store.observed_strategy_keys(
        pd.Series(["1:100:A>G", "1:999:G>A"]),
        ["s1", "s2", "s3"],
    ) == {
        ("1:100:A>G", "s1"),
        ("1:100:A>G", "s2"),
        ("1:100:A>G", "s3"),
    }

    import duckdb

    with duckdb.connect() as connection:
        assert connection.read_parquet(str(store.allele_gene_path)).columns == ALLELE_GENE_COLUMNS
        assert connection.read_parquet(str(store.allele_path)).columns == ALLELE_COLUMNS

    cached = build_or_load_observed_variant_store(
        variant_annotations_tsv=annotations,
        analytics_dir=tmp_path / "analytics",
        strategies=["s1", "s2", "s3"],
    )
    assert cached.cache_hit is True


def test_observed_store_focal_sampling_matches_stable_md5_topk(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    rows = []
    for offset, strategies in enumerate(
        ["s1", "s2", "s1,s2", "s1", "s2", "s1,s2", "s1", "s2"],
    ):
        ref, alt = [("A", "G"), ("C", "T"), ("G", "A"), ("T", "C")][offset % 4]
        rows.append(
            {
                "variant_key": f"1:{100 + offset}:{ref}>{alt}",
                "gene_id": "1",
                "event_type": "snv",
                "ref": ref,
                "alt": alt,
                "lookup_status": "ok",
                "strategies": strategies,
            }
        )
    rows.extend(
        [
            {
                "variant_key": "1:120:A>AT",
                "gene_id": "1",
                "event_type": "ins",
                "ref": "A",
                "alt": "AT",
                "lookup_status": "ok",
                "strategies": "s1",
            },
            {
                "variant_key": "1:121:A>G",
                "gene_id": "1",
                "event_type": "snv",
                "ref": "A",
                "alt": "G",
                "lookup_status": "ref_mismatch",
                "strategies": "s1",
            },
            {
                "variant_key": "1:122:A>G:extra",
                "gene_id": "1",
                "event_type": "snv",
                "ref": "A",
                "alt": "G",
                "lookup_status": "ok",
                "strategies": "s1",
            },
            {
                "variant_key": "1:123:A>G",
                "gene_id": "missing_gene",
                "event_type": "snv",
                "ref": "A",
                "alt": "G",
                "lookup_status": "ok",
                "strategies": "s2",
            },
            {
                "variant_key": "chr:124:A>G",
                "gene_id": "1",
                "event_type": "snv",
                "ref": "A",
                "alt": "G",
                "lookup_status": "ok",
                "strategies": "s1",
            },
        ]
    )
    _write_annotations(annotations, rows)
    contexts = {"1": [(0, 4, "cds"), (4, 30, "intron")]}
    genes = {"1": {"begin": 100, "chrom": "1", "length": 30}}
    strategies = ["s1", "s2"]

    store = build_or_load_observed_variant_store(
        variant_annotations_tsv=annotations,
        analytics_dir=tmp_path / "analytics",
        strategies=strategies,
    )
    expected = _expected_md5_sample(
        annotations,
        contexts,
        genes,
        strategies,
        limit=3,
        seed=17,
    )
    observed = controls._sample_focal_snvs(
        store,
        contexts,
        genes,
        strategies,
        limit=3,
        seed=17,
    )

    pd.testing.assert_frame_equal(observed, expected)
    all_sampled_keys = {
        str(row[0])
        for chunk in store.iter_sampled_focal_rows(
            strategies,
            set(genes),
            limit=100,
            seed=17,
        )
        for row in chunk
    }
    assert "1:120:A>AT" not in all_sampled_keys
    assert "1:121:A>G" not in all_sampled_keys
    assert "1:122:A>G:extra" not in all_sampled_keys
    assert "1:123:A>G" not in all_sampled_keys
    assert "chr:124:A>G" not in all_sampled_keys
    pd.testing.assert_frame_equal(
        controls._sample_focal_snvs(
            store,
            contexts,
            genes,
            strategies,
            limit=3,
            seed=17,
        ),
        observed,
    )
    assert observed.groupby("strategy").size().to_dict() == {"s1": 3, "s2": 3}
    assert not observed.duplicated(["strategy", "gene_id", "variant_key"]).any()

    reordered_annotations = tmp_path / "variant_annotations.reordered.tsv.gz"
    _write_annotations(reordered_annotations, list(reversed(rows)))
    reordered_store = build_or_load_observed_variant_store(
        variant_annotations_tsv=reordered_annotations,
        analytics_dir=tmp_path / "reordered_analytics",
        strategies=strategies,
    )
    reordered = controls._sample_focal_snvs(
        reordered_store,
        contexts,
        genes,
        strategies,
        limit=3,
        seed=17,
    )
    pd.testing.assert_frame_equal(reordered, observed)

    changed_seed = controls._sample_focal_snvs(
        store,
        contexts,
        genes,
        strategies,
        limit=3,
        seed=18,
    )
    assert set(zip(changed_seed["strategy"], changed_seed["variant_key"])) != set(
        zip(observed["strategy"], observed["variant_key"])
    )


def test_observed_store_invalidates_when_source_changes(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    first_rows = [
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "event_type": "snv",
            "ref": "A",
            "alt": "G",
            "lookup_status": "ok",
            "strategies": "s1",
        }
    ]
    _write_annotations(annotations, first_rows)
    first = build_or_load_observed_variant_store(
        variant_annotations_tsv=annotations,
        analytics_dir=tmp_path / "analytics",
        strategies=["s1"],
    )
    assert first.manifest["allele_count"] == 1

    _write_annotations(
        annotations,
        [
            *first_rows,
            {
                "variant_key": "1:101:C>T",
                "gene_id": "1",
                "event_type": "snv",
                "ref": "C",
                "alt": "T",
                "lookup_status": "ok",
                "strategies": "s1",
            },
        ],
    )
    rebuilt = build_or_load_observed_variant_store(
        variant_annotations_tsv=annotations,
        analytics_dir=tmp_path / "analytics",
        strategies=["s1"],
    )

    assert rebuilt.cache_hit is False
    assert rebuilt.manifest["allele_count"] == 2


def test_observed_store_rejects_strategy_contract_mismatch(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    _write_annotations(
        annotations,
        [
            {
                "variant_key": "1:100:A>G",
                "gene_id": "1",
                "event_type": "snv",
                "ref": "A",
                "alt": "G",
                "lookup_status": "ok",
                "strategies": "s1,s2",
            }
        ],
    )

    with pytest.raises(ValueError, match="source strategies differ"):
        build_or_load_observed_variant_store(
            variant_annotations_tsv=annotations,
            analytics_dir=tmp_path / "analytics",
            strategies=["s1"],
        )
