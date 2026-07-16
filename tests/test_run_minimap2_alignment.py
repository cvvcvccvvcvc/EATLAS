from __future__ import annotations

import sys
from pathlib import Path


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from run_minimap2_alignment import (  # noqa: E402
    EVENT_FIELDS,
    SEGMENT_FIELDS,
    empty_summary,
    parse_paf,
)


PAF_LINES = [
    "q1\t10\t0\t10\t+\ttarget_1\t10\t0\t10\t9\t10\t60\ttp:A:P\tcs:Z::4*ag:5",
    "q2\t10\t0\t10\t+\ttarget_1\t10\t0\t10\t9\t10\t60\ttp:A:P\tcs:Z::2*ct:7",
]


def parse_rows(path: Path, lines: list[str]):
    path.write_text("\n".join(lines) + "\n")
    metadata = {
        "q1": {"ortholog_gene_id": "101", "tax_id": "1", "taxname": "species 1"},
        "q2": {"ortholog_gene_id": "102", "tax_id": "2", "taxname": "species 2"},
    }
    summaries = {
        query_id: empty_summary("1", "minimap2_asm20", "asm20", meta, 10)
        for query_id, meta in metadata.items()
    }
    segments, events, _ = parse_paf(
        path,
        "1",
        "minimap2_asm20",
        "asm20",
        {"genomic_accession": "NC_1", "genomic_begin": "100"},
        metadata,
        summaries,
        1,
    )
    return (
        sorted(tuple(row[field] for field in SEGMENT_FIELDS) for row in segments),
        sorted(tuple(row[field] for field in EVENT_FIELDS) for row in events),
    )


def test_paf_identifiers_do_not_depend_on_record_order(tmp_path: Path) -> None:
    forward_segments, forward_events = parse_rows(tmp_path / "forward.paf", PAF_LINES)
    reverse_segments, reverse_events = parse_rows(tmp_path / "reverse.paf", PAF_LINES[::-1])

    assert forward_segments == reverse_segments
    assert forward_events == reverse_events
    assert len({row[7] for row in forward_events}) == len(forward_events)
    assert all(str(row[22]).startswith("paf:") for row in forward_segments)
