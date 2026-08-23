"""Tests for shared and chunk-based Ensembl Compara MAF normalization."""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import sys
from pathlib import Path
from urllib.error import HTTPError

import pytest


from bin import ensembl_compara_maf as maf
from bin import run_ensembl_compara_maf_chunk_alignment as maf_chunk
from bin.ensembl_compara_maf import (
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


def retry_args(retries: int = 8) -> argparse.Namespace:
    return argparse.Namespace(
        strategy="maf",
        method="EPO_EXTENDED",
        species_set="test",
        retries=retries,
        retry_base_seconds=5.0,
        retry_max_seconds=300.0,
        timeout=120.0,
    )


def test_chunk_cli_uses_fixed_epo_defaults(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_ensembl_compara_maf_chunk_alignment.py",
            "--chunk-task-dir",
            "task",
            "--outdir",
            "out",
        ],
    )

    parsed = maf_chunk.parse_args()

    assert parsed.strategy == maf.STRATEGY_NAME
    assert (parsed.release, parsed.species_set, parsed.method) == (
        maf.RELEASE,
        maf.SPECIES_SET,
        maf.METHOD,
    )
    assert (parsed.timeout, parsed.retries) == (
        maf.REQUEST_TIMEOUT_SECONDS,
        maf.DOWNLOAD_ATTEMPTS,
    )


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


def test_maf_retry_classification_and_backoff(monkeypatch: pytest.MonkeyPatch) -> None:
    transient = HTTPError("https://example.test/source.maf.gz", 503, "unavailable", None, None)
    missing = HTTPError("https://example.test/source.maf.gz", 404, "not found", None, None)
    monkeypatch.setattr(maf.random, "uniform", lambda _start, _end: 0.0)

    assert maf.retryable_maf_error(EOFError("truncated"))
    assert maf.retryable_maf_error(transient)
    assert not maf.retryable_maf_error(missing)
    assert maf.missing_maf_source_error(missing)
    assert not maf.retryable_maf_error(FileNotFoundError("missing.maf.gz"))
    assert maf.retry_sleep_seconds(retry_args(), 1) == 5.0
    assert maf.retry_sleep_seconds(retry_args(), 7) == 300.0


def test_chunk_scan_resumes_after_truncated_stream_without_duplicates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    gene = {
        "gene_id": "1",
        "human_src": "homo_sapiens.1",
        "genomic_accession": "NC_000001.11",
        "target_origin1": "1",
        "target_end1": "2",
        "target_length": "2",
    }
    blocks = [
        [
            MafSequence("homo_sapiens.1", 0, 1, "+", 100, "A"),
            MafSequence("mus_musculus.1", 0, 1, "+", 100, "G"),
        ],
        [
            MafSequence("homo_sapiens.1", 1, 1, "+", 100, "C"),
            MafSequence("mus_musculus.1", 1, 1, "+", 100, "T"),
        ],
    ]
    attempts = iter([(blocks[:1], EOFError("truncated")), (blocks, None)])
    converted: list[str] = []

    def iter_blocks(_handle):
        attempt_blocks, error = next(attempts)
        yield from attempt_blocks
        if error:
            raise error

    def convert(*call_args, **_kwargs):
        converted.append(call_args[7])
        return call_args[9] + 1

    monkeypatch.setattr(maf_chunk, "open_maf_text", lambda _source, _timeout: contextlib.nullcontext())
    monkeypatch.setattr(maf_chunk, "iter_maf_blocks", iter_blocks)
    monkeypatch.setattr(maf_chunk, "convert_pair", convert)
    monkeypatch.setattr(maf_chunk.time, "sleep", lambda _seconds: None)

    _event_id, used_blocks, row_count, failures = maf_chunk.scan_chunk_source(
        retry_args(retries=2),
        "source.maf.gz",
        maf_chunk.GeneIntervalIndex.build([gene]),
        ["1"],
        {},
        None,
        None,
    )

    assert converted == ["source.maf.gz:block1:row2", "source.maf.gz:block2:row2"]
    assert used_blocks == 2
    assert row_count == 2
    assert failures == []


def test_missing_chunk_source_is_not_retried(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = 0

    def missing(_source: str, _timeout: float):
        nonlocal calls
        calls += 1
        raise HTTPError("https://example.test/source.maf.gz", 404, "not found", None, None)

    monkeypatch.setattr(maf_chunk, "open_maf_text", missing)

    _event_id, used_blocks, row_count, failures = maf_chunk.scan_chunk_source(
        retry_args(),
        "https://example.test/source.maf.gz",
        maf_chunk.GeneIntervalIndex.build([]),
        ["1"],
        {},
        None,
        None,
    )

    assert calls == 1
    assert used_blocks == 0
    assert row_count == 0
    assert len(failures) == 1
    assert "failed after 1 attempts" in failures[0]["message"]
