from __future__ import annotations

from pathlib import Path

import pysam


from bin.run_nucmer_alignment import (
    empty_summary,
    parse_sam,
    query_interval,
    sam_alignment_type,
)


TARGET_SEQ = "AACCGGTTAACCGGTT"


def test_parse_sam_emits_contiguous_indels_and_deduplicates_alignments(tmp_path: Path) -> None:
    sam_path = tmp_path / "nucmer.sam"
    query_seq = "AACCTAGGTTCAGTT"
    ambiguous_seq = "AACCN" + TARGET_SEQ[5:]
    sam_path.write_text(
        "\n".join(
            [
                "@HD\tVN:1.4\tSO:unsorted",
                f"@SQ\tSN:target_1\tLN:{len(TARGET_SEQ)}",
                (
                    f"query_1\t0\ttarget_1\t1\t10\t4M2I4M3D5M\t*\t0\t0\t{query_seq}\t*"
                    "\tNM:i:6"
                ),
                (
                    f"query_1\t2048\ttarget_1\t1\t10\t4M2I4M3D5M\t*\t0\t0\t{query_seq}\t*"
                    "\tNM:i:6"
                ),
                (
                    f"query_2\t0\ttarget_1\t1\t10\t16M\t*\t0\t0\t{ambiguous_seq}\t*"
                    "\tNM:i:1"
                ),
            ]
        )
        + "\n"
    )
    metadata = {
        "query_1": {
            "ortholog_gene_id": "101",
            "tax_id": "10090",
            "taxname": "Mus musculus",
            "sequence_length": str(len(query_seq)),
        },
        "query_2": {
            "ortholog_gene_id": "102",
            "tax_id": "10116",
            "taxname": "Rattus norvegicus",
            "sequence_length": str(len(ambiguous_seq)),
        },
    }
    summaries = {
        query_id: empty_summary("1", meta, len(TARGET_SEQ))
        for query_id, meta in metadata.items()
    }

    segments, events, ambiguous_count = parse_sam(
        sam_path,
        "1",
        {"genomic_accession": "NC_000001.11", "genomic_begin": "100"},
        TARGET_SEQ,
        metadata,
        summaries,
    )

    assert len(segments) == 3
    assert [
        (
            row["event_type"],
            row["target_start0"],
            row["target_end0"],
            row["ref"],
            row["alt"],
        )
        for row in events
    ] == [
        ("ins", 4, 4, "", "TA"),
        ("del", 8, 11, "AAC", ""),
        ("snv", 12, 13, "G", "A"),
    ]
    assert all(row["native_record_id"] == 1 for row in events)
    assert all(row["mapq"] == 10 for row in events)
    assert all(row["native_alignment_type"] == "primary" for row in events)
    assert all(row["qc_flags"] == "unfiltered_nucmer" for row in events)
    assert all(row["mapq"] == 10 for row in segments)
    assert summaries["query_1"]["event_count"] == 3
    assert summaries["query_1"]["primary_segment_count"] == 1
    assert summaries["query_1"]["secondary_segment_count"] == 1
    assert summaries["query_2"]["event_count"] == 0
    assert summaries["query_2"]["qc_flags"] == {
        "ambiguous_event_allele",
        "unfiltered_nucmer",
    }
    assert ambiguous_count == 1


def test_sam_alignment_type_preserves_combined_flags() -> None:
    read = pysam.AlignedSegment()
    read.flag = 0
    assert sam_alignment_type(read) == "primary"
    read.flag = 256
    assert sam_alignment_type(read) == "secondary"
    read.flag = 2048
    assert sam_alignment_type(read) == "supplementary"
    read.flag = 256 | 2048
    assert sam_alignment_type(read) == "secondary_supplementary"


def test_query_interval_converts_reverse_hard_clipping_to_forward_coordinates() -> None:
    read = pysam.AlignedSegment()
    read.query_name = "query_1"
    read.query_sequence = "A" * 10
    read.flag = 16
    read.cigar = ((5, 5), (0, 10), (5, 3))

    assert query_interval(read) == (3, 13, 18)
