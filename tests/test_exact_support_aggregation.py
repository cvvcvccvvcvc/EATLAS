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
        "mapq": "60",
        "native_alignment_type": "P",
        "support_row_count": "1",
    }
    spool = ExactSupportSpool(tmp_path)
    spool.add_group(aggregate, {"strategy": "s1"}, [edge])
    spool.add_group(
        aggregate,
        {"strategy": "s1"},
        [
            {
                **edge,
                "mapq": "20",
                "native_alignment_type": "S",
            }
        ],
    )

    row_count = aggregate_exact_support(
        spool,
        [None, aggregate],
    )

    assert row_count == 1
    assert not (tmp_path / "variant_ortholog_support").exists()
    support = aggregate["_support_by_strategy"]["s1"]
    assert (support.ortholog_count_hint, support.row_count) == (1, 2)
    assert aggregate["_exact_ortholog_count"] == 1


def test_compact_event_support_uses_the_same_exact_spool(tmp_path: Path) -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "_variant_context_id": 1,
        "_exact_ortholog_count": 0,
        "_support_by_strategy": {},
    }
    event = {
        "strategy": "s1",
        "support_row_count": "1",
        "support_ortholog_count": "1",
    }
    support_edge = {
        "ortholog_gene_id": "101",
        "tax_id": "10090",
        "taxname": "Mus musculus",
        "mapq": "",
        "native_alignment_type": "",
        "support_row_count": "1",
    }
    add_strategy_support(aggregate, event)
    spool = ExactSupportSpool(tmp_path)
    spool.add_group(aggregate, event, [support_edge])

    row_count = aggregate_exact_support(
        spool,
        [None, aggregate],
    )

    assert row_count == 1
    assert not (tmp_path / "variant_ortholog_support").exists()
    assert aggregate["_exact_ortholog_count"] == 1
