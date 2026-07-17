from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from run_ensembl_compara_maf_alignment import (  # noqa: E402
    AlignmentRow,
    EVENT_FIELDS,
    SEGMENT_FIELDS,
    MafSequence,
    TsvGzWriter,
    convert_pair,
    empty_summary,
    resolve_maf_dots,
    to_alignment_row,
)


def args() -> argparse.Namespace:
    return argparse.Namespace(strategy="maf", method="EPO_EXTENDED", species_set="test")


def read_rows(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_pair(
    tmp_path: Path,
    human_row: AlignmentRow,
    query_row: AlignmentRow,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, object]]:
    segment_path = tmp_path / "segments.tsv.gz"
    event_path = tmp_path / "events.tsv.gz"
    segment_writer = TsvGzWriter(segment_path, SEGMENT_FIELDS)
    event_writer = TsvGzWriter(event_path, EVENT_FIELDS)
    summary = empty_summary(args(), query_row, human_row.end1 - human_row.start1 + 1)
    summary["gene_id"] = "1"
    try:
        convert_pair(
            args(),
            "1",
            "NC_1",
            human_row.start1,
            human_row.end1,
            human_row,
            query_row,
            "record-1",
            summary,
            1,
            segment_writer,
            event_writer,
        )
    finally:
        segment_writer.close()
        event_writer.close()
    return read_rows(segment_path), read_rows(event_path), summary


def test_resolve_maf_dots_uses_reference_alignment_character() -> None:
    reference = AlignmentRow("human", "1", 101, 102, 1, "A-C-", "human.1")
    query = AlignmentRow("mouse", "1", 201, 203, 1, ".T..", "mouse.1")

    resolved = resolve_maf_dots(reference, query)

    assert resolved.seq == "ATC-"
    assert query.seq == ".T.."


def test_resolve_maf_dots_after_reverse_orientation() -> None:
    human = to_alignment_row(MafSequence("human.1", 10, 1, "-", 100, "A-"), True)
    query = to_alignment_row(MafSequence("mouse.1", 20, 1, "-", 100, ".C"), True)

    assert human.seq == "-T"
    assert query.seq == "G."
    assert resolve_maf_dots(human, query).seq == "GT"


def test_mixed_real_and_dot_insertion_emits_only_real_bases(tmp_path: Path) -> None:
    human = AlignmentRow("human", "1", 101, 101, 1, "---A", "human.1")
    query = AlignmentRow("mouse", "1", 201, 203, 1, "TG.A", "mouse.1")

    segments, events, summary = run_pair(tmp_path, human, resolve_maf_dots(human, query))

    assert [(row["event_type"], row["ref"], row["alt"]) for row in events] == [("ins", "", "TG")]
    assert segments[0]["query_start0"] == "202"
    assert segments[0]["query_end0"] == "203"
    assert summary["event_count"] == 1


def test_ambiguous_maf_indel_is_not_emitted(tmp_path: Path) -> None:
    human = AlignmentRow("human", "1", 101, 101, 1, "--A", "human.1")
    query = AlignmentRow("mouse", "1", 201, 203, 1, "RYA", "mouse.1")

    _segments, events, summary = run_pair(tmp_path, human, query)

    assert events == []
    assert summary["event_count"] == 0
    assert summary["qc_flags"] == {"ambiguous_base"}
