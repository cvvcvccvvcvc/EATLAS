from __future__ import annotations

import sys

import pytest

from bin.run_minimap2_alignment import (
    EVENT_FIELDS,
    QuerySlice,
    SEGMENT_FIELDS,
    empty_summary,
    generate_long_pseudoreads,
    is_primary,
    iter_paf_records,
    parse_args,
    parse_paf,
    pseudoread_starts,
    select_pseudoread_backbone,
    validate_query_mode,
)


PAF_LINES = [
    "q1\t10\t0\t10\t+\ttarget_1\t10\t0\t10\t9\t10\t60\ttp:A:P\tcs:Z::4*ag:5",
    "q2\t10\t0\t10\t+\ttarget_1\t10\t0\t10\t9\t10\t60\ttp:A:P\tcs:Z::2*ct:7",
]


def minimap2_cli_args(preset: str) -> list[str]:
    return [
        "run_minimap2_alignment.py",
        "--task-dir",
        "task",
        "--source-target-fasta",
        "target.fa.gz",
        "--source-ortholog-fasta",
        "ortholog.fa.gz",
        "--outdir",
        "out",
        "--strategy",
        f"minimap2_{preset}",
        "--preset",
        preset,
        "--minimap2-bin",
        "minimap2",
    ]


@pytest.mark.parametrize("preset", ["asm10", "asm20"])
def test_cli_accepts_supported_fixed_presets(
    monkeypatch: pytest.MonkeyPatch,
    preset: str,
) -> None:
    monkeypatch.setattr(sys, "argv", minimap2_cli_args(preset))

    assert parse_args().preset == preset


def test_cli_rejects_other_presets(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", minimap2_cli_args("asm5"))

    with pytest.raises(SystemExit):
        parse_args()


def test_map_ont_requires_complete_pseudoread_parameters(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    args = minimap2_cli_args("map-ont") + [
        "--pseudoread-len",
        "30000",
        "--pseudoread-step",
        "15000",
    ]
    monkeypatch.setattr(sys, "argv", args)

    parsed = parse_args()
    validate_query_mode(parsed.preset, parsed.pseudoread_len, parsed.pseudoread_step)

    with pytest.raises(ValueError, match="requires pseudoread"):
        validate_query_mode("map-ont", 0, 0)
    with pytest.raises(ValueError, match="must not exceed"):
        validate_query_mode("map-ont", 1000, 1001)


def test_long_pseudoread_generation_is_gap_free_and_keeps_short_sequences(
    tmp_path: Path,
) -> None:
    source = tmp_path / "orthologs.fa"
    output = tmp_path / "long.fa"
    source.write_text(">short\n" + "A" * 5 + "\n>long\n" + "C" * 50 + "\n")
    metadata = {
        "short": {"sequence_length": "5"},
        "long": {"sequence_length": "50"},
    }

    generation = generate_long_pseudoreads(source, output, metadata, read_len=30, step=15)

    assert pseudoread_starts(5, 30, 15) == [0]
    assert pseudoread_starts(50, 30, 15) == [0, 15, 20]
    assert generation.total_reads == 4
    assert [
        (row.source_sequence_id, row.source_start0, row.source_end0)
        for row in generation.query_slices.values()
    ] == [
        ("short", 0, 5),
        ("long", 0, 30),
        ("long", 15, 45),
        ("long", 20, 50),
    ]


def parse_rows(path: Path, lines: list[str]):
    path.write_text("\n".join(lines) + "\n")
    metadata = {
        "q1": {"ortholog_gene_id": "101", "tax_id": "1", "taxname": "species 1", "sequence_length": "10"},
        "q2": {"ortholog_gene_id": "102", "tax_id": "2", "taxname": "species 2", "sequence_length": "10"},
    }
    summaries = {
        query_id: empty_summary("1", "minimap2_asm20", "asm20", meta, 10)
        for query_id, meta in metadata.items()
    }
    segments, events = parse_paf(
        path,
        "1",
        "minimap2_asm20",
        "asm20",
        {"genomic_accession": "NC_1", "genomic_begin": "100"},
        metadata,
        summaries,
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


def test_paf_rejects_malformed_record(tmp_path: Path) -> None:
    path = tmp_path / "malformed.paf"
    path.write_text("\nq1\t10\tbroken\n")

    with pytest.raises(
        ValueError,
        match=r"malformed\.paf at line 2.*observed 3",
    ):
        list(iter_paf_records(path))


def test_paf_rejects_unknown_query_id(tmp_path: Path) -> None:
    path = tmp_path / "unknown.paf"

    with pytest.raises(
        ValueError,
        match=r"unknown\.paf at line 1.*unknown query ID 'q3'",
    ):
        parse_rows(path, [PAF_LINES[0].replace("q1", "q3", 1)])


def test_paf_rejects_query_source_without_metadata(tmp_path: Path) -> None:
    path = tmp_path / "missing_source.paf"
    path.write_text(PAF_LINES[0] + "\n")

    with pytest.raises(
        ValueError,
        match=r"missing_source\.paf at line 1.*'missing'.*without metadata",
    ):
        parse_paf(
            path,
            "1",
            "minimap2_asm20",
            "asm20",
            {"genomic_accession": "NC_1", "genomic_begin": "100"},
            {},
            {},
            query_slices={
                "q1": QuerySlice("missing", 0, 10, 10, 1, False),
            },
        )


def test_paf_preserves_native_type_and_prefers_primary_duplicates(tmp_path: Path) -> None:
    path = tmp_path / "records.paf"
    path.write_text(
        "\n".join(
            [
                "q1\t10\t0\t10\t+\ttarget_1\t10\t0\t10\t9\t10\t20\ttp:A:S\tcs:Z::4*ag:5",
                "q1\t10\t0\t10\t+\ttarget_1\t10\t0\t10\t9\t10\t60\ttp:A:P\tcs:Z::4*ag:5",
                "q2\t10\t0\t10\t+\ttarget_1\t10\t0\t10\t9\t10\t20\ttp:A:S\tcs:Z::2*ct:7",
            ]
        )
        + "\n"
    )
    metadata = {
        "q1": {"ortholog_gene_id": "101", "tax_id": "1", "taxname": "species 1", "sequence_length": "10"},
        "q2": {"ortholog_gene_id": "102", "tax_id": "2", "taxname": "species 2", "sequence_length": "10"},
    }
    summaries = {
        query_id: empty_summary("1", "minimap2_asm20", "asm20", meta, 10)
        for query_id, meta in metadata.items()
    }

    _segments, events = parse_paf(
        path,
        "1",
        "minimap2_asm20",
        "asm20",
        {"genomic_accession": "NC_1", "genomic_begin": "100"},
        metadata,
        summaries,
    )

    assert len(events) == 2
    by_query = {row["query_id"]: row for row in events}
    assert by_query["q1"]["qc_flags"] == ""
    assert by_query["q1"]["mapq"] == 60
    assert by_query["q1"]["native_alignment_type"] == "P"
    assert by_query["q2"]["qc_flags"] == ""
    assert by_query["q2"]["mapq"] == 20
    assert by_query["q2"]["native_alignment_type"] == "S"
    assert summaries["q1"]["event_count"] == 1
    assert summaries["q2"]["event_count"] == 1


def test_paf_reports_actual_mapq_without_threshold_flag(tmp_path: Path) -> None:
    path = tmp_path / "low_mapq.paf"
    path.write_text(
        "q1\t10\t0\t10\t+\ttarget_1\t10\t0\t10\t9\t10\t0\ttp:A:S\tcs:Z::4*ag:5\n"
    )
    metadata = {
        "q1": {
            "ortholog_gene_id": "101",
            "tax_id": "1",
            "taxname": "species 1",
            "sequence_length": "10",
        }
    }
    summaries = {"q1": empty_summary("1", "minimap2_asm20", "asm20", metadata["q1"], 10)}

    segments, events = parse_paf(
        path,
        "1",
        "minimap2_asm20",
        "asm20",
        {"genomic_accession": "NC_1", "genomic_begin": "100"},
        metadata,
        summaries,
    )

    assert segments[0]["mapq"] == 0
    assert segments[0]["qc_flags"] == ""
    assert events[0]["mapq"] == 0
    assert events[0]["native_alignment_type"] == "S"
    assert events[0]["qc_flags"] == ""


def test_minimap_inversion_types_preserve_primary_semantics() -> None:
    assert is_primary({"tp": "I"}) is True
    assert is_primary({"tp": "i"}) is False


def test_long_pseudoread_coordinates_are_lifted_and_events_deduplicated(
    tmp_path: Path,
) -> None:
    path = tmp_path / "long.paf"
    path.write_text(
        "\n".join(
            [
                "read1\t30000\t100\t200\t+\ttarget_1\t50000\t10\t110\t99\t100\t60\ttp:A:P\tcs:Z::4*ag:95",
                "read2\t30000\t200\t300\t+\ttarget_1\t50000\t10\t110\t99\t100\t40\ttp:A:S\tcs:Z::4*ag:95",
            ]
        )
        + "\n"
    )
    metadata = {
        "ortholog_101": {
            "ortholog_gene_id": "101",
            "tax_id": "1",
            "taxname": "species 1",
            "sequence_length": "45000",
        }
    }
    query_slices = {
        "read1": QuerySlice("ortholog_101", 0, 30000, 45000, 1, True),
        "read2": QuerySlice("ortholog_101", 15000, 45000, 45000, 2, True),
    }
    summaries = {
        "ortholog_101": empty_summary(
            "1",
            "minimap2_map_ont_pseudoreads_30000_15000",
            "map-ont",
            metadata["ortholog_101"],
            50000,
        )
    }

    segments, events = parse_paf(
        path,
        "1",
        "minimap2_map_ont_pseudoreads_30000_15000",
        "map-ont",
        {"genomic_accession": "NC_1", "genomic_begin": "100"},
        metadata,
        summaries,
        query_slices,
    )

    assert [(row["query_start0"], row["query_end0"]) for row in segments] == [
        (100, 200),
        (15200, 15300),
    ]
    assert all(row["query_id"] == "ortholog_101" for row in segments)
    assert [row["mapq"] for row in segments] == [60, 40]
    assert len(events) == 1
    assert events[0]["qc_flags"] == "filtered_pseudoread"


def test_long_secondary_backbone_record_keeps_flagged_alt_support(
    tmp_path: Path,
) -> None:
    path = tmp_path / "long_secondary.paf"
    path.write_text(
        "read1\t30000\t0\t100\t+\ttarget_1\t50000\t10\t110\t99\t100\t1"
        "\ttp:A:S\tcs:Z::4*ag:95\n"
    )
    metadata = {
        "ortholog_101": {
            "ortholog_gene_id": "101",
            "tax_id": "1",
            "taxname": "species 1",
            "sequence_length": "30000",
        }
    }
    query_slices = {
        "read1": QuerySlice("ortholog_101", 0, 30000, 30000, 1, True),
    }
    summaries = {
        "ortholog_101": empty_summary(
            "1",
            "minimap2_map_ont_pseudoreads_30000_15000",
            "map-ont",
            metadata["ortholog_101"],
            50000,
        )
    }
    backbone = select_pseudoread_backbone(path, query_slices)

    segments, events = parse_paf(
        path,
        "1",
        "minimap2_map_ont_pseudoreads_30000_15000",
        "map-ont",
        {"genomic_accession": "NC_1", "genomic_begin": "100"},
        metadata,
        summaries,
        query_slices,
        backbone.accepted_record_ids,
    )

    assert segments[0]["mapq"] == 1
    assert segments[0]["is_primary"] == "false"
    assert segments[0]["qc_flags"] == "filtered_pseudoread"
    assert events[0]["mapq"] == 1
    assert events[0]["native_alignment_type"] == "S"
    assert events[0]["qc_flags"] == "filtered_pseudoread"


def test_pseudoread_backbone_keeps_dominant_monotonic_order(tmp_path: Path) -> None:
    path = tmp_path / "backbone.paf"
    path.write_text(
        "\n".join(
            [
                "read1\t30\t0\t30\t+\ttarget\t100\t10\t40\t30\t30\t60\ttp:A:P",
                "read3\t30\t0\t30\t+\ttarget\t100\t20\t50\t30\t30\t60\ttp:A:P",
                "read2\t30\t0\t30\t+\ttarget\t100\t30\t60\t30\t30\t60\ttp:A:P",
                "reverse\t30\t0\t30\t-\ttarget\t100\t40\t70\t30\t30\t60\ttp:A:S",
            ]
        )
        + "\n"
    )
    slices = {
        "read1": QuerySlice("ortholog_101", 0, 30, 60, 1, True),
        "read2": QuerySlice("ortholog_101", 15, 45, 60, 2, True),
        "read3": QuerySlice("ortholog_101", 30, 60, 60, 3, True),
        "reverse": QuerySlice("ortholog_101", 30, 60, 60, 3, True),
    }

    selected = select_pseudoread_backbone(path, slices)
    retained_query_names = {
        fields[0]
        for fields, native_record_id, _line_number in iter_paf_records(path)
        if native_record_id in selected.accepted_record_ids
    }

    assert retained_query_names == {"read1", "read2"}
    assert selected.input_alignment_count == 4
    assert selected.after_strand_count == 3
    assert selected.retained_alignment_count == 2
