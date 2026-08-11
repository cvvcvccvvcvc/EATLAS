from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest


duckdb = pytest.importorskip("duckdb")
PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "bin"))
sys.modules.setdefault("pysam", types.SimpleNamespace())

from annotate_events import (  # noqa: E402
    ExactSupportSpool,
    add_strategy_support,
    aggregate_exact_support,
)
from finalize_annotation_partitions import (  # noqa: E402
    merge_ortholog_support_dataset,
    sql_string,
)


def test_exact_support_collapses_normalized_event_collisions(tmp_path: Path) -> None:
    aggregate = {
        "variant_key": "1:100:AA>A",
        "gene_id": "1",
        "_variant_context_id": 1,
        "_exact_ortholog_count": 0,
        "_support_by_strategy": {},
    }
    for _raw_representation in range(2):
        add_strategy_support(
            aggregate,
            {
                "strategy": "s1",
                "support_row_count": "1",
                "support_ortholog_count": "1",
            },
        )
    edge = {
        "ortholog_gene_id": "101",
        "tax_id": "10090",
        "taxname": "Mus musculus",
        "support_row_count": "1",
    }
    spool = ExactSupportSpool(tmp_path)
    spool.add_group(aggregate, {"strategy": "s1"}, [edge])
    spool.add_group(aggregate, {"strategy": "s1"}, [edge])

    row_count = aggregate_exact_support(
        spool,
        [None, aggregate],
        tmp_path / "variant_ortholog_support",
    )

    rows = duckdb.connect().execute(
        "SELECT * FROM read_parquet(?)",
        [str(tmp_path / "variant_ortholog_support" / "*.parquet")],
    ).fetchall()
    assert row_count == 1
    assert rows == [("1:100:AA>A", "1", "s1", "101", "10090", "Mus musculus", 2)]
    support = aggregate["_support_by_strategy"]["s1"]
    assert (support.ortholog_count_hint, support.row_count) == (1, 2)
    assert aggregate["_exact_ortholog_count"] == 1


def test_raw_event_support_uses_the_same_exact_spool(tmp_path: Path) -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "_variant_context_id": 1,
        "_exact_ortholog_count": 0,
        "_support_by_strategy": {},
    }
    event = {
        "strategy": "s1",
        "ortholog_gene_id": "101",
        "tax_id": "10090",
        "taxname": "Mus musculus",
        "support_row_count": "1",
    }
    add_strategy_support(aggregate, event)
    spool = ExactSupportSpool(tmp_path)
    spool.add_group(aggregate, event, [event])

    row_count = aggregate_exact_support(
        spool,
        [None, aggregate],
        tmp_path / "variant_ortholog_support",
    )

    assert row_count == 1
    assert aggregate["_exact_ortholog_count"] == 1


def test_finalizer_copies_partition_parquet_without_row_rewrite(tmp_path: Path) -> None:
    partitions = []
    connection = duckdb.connect()
    connection.execute(
        """
        CREATE TABLE support (
            variant_key VARCHAR,
            gene_id VARCHAR,
            strategy VARCHAR,
            ortholog_gene_id VARCHAR,
            tax_id VARCHAR,
            taxname VARCHAR,
            support_row_count UBIGINT
        )
        """
    )
    for index, (variant_key, gene_id, ortholog_gene_id) in enumerate(
        [("1:1:A>G", "1", "101"), ("2:2:C>T", "2", "201")],
        start=1,
    ):
        partition = tmp_path / f"partition_{index:06d}"
        support_dir = partition / "variant_ortholog_support"
        support_dir.mkdir(parents=True)
        source = support_dir / "part-00000.parquet"
        connection.execute("DELETE FROM support")
        connection.execute(
            "INSERT INTO support VALUES (?, ?, 's1', ?, '10090', 'Mus musculus', 1)",
            [variant_key, gene_id, ortholog_gene_id],
        )
        connection.execute(
            f"COPY support TO {sql_string(source)} (FORMAT PARQUET, COMPRESSION ZSTD)"
        )
        partitions.append(
            (
                partition,
                {
                    "partition_id": partition.name,
                    "variant_ortholog_support_format": "parquet_dataset",
                    "variant_ortholog_support_file_count": 1,
                },
            )
        )

    output = tmp_path / "merged"
    row_count, file_count = merge_ortholog_support_dataset(partitions, output)

    assert (row_count, file_count) == (2, 2)
    assert len(list(output.glob("*.parquet"))) == 2
    assert connection.execute(
        "SELECT COUNT(*) FROM read_parquet(?)",
        [str(output / "*.parquet")],
    ).fetchone()[0] == 2
