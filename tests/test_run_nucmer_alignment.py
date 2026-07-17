from __future__ import annotations

import sys
from pathlib import Path


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from run_nucmer_alignment import empty_summary, parse_snps  # noqa: E402


def test_parse_snps_excludes_ambiguous_event_alleles(tmp_path: Path) -> None:
    snps_path = tmp_path / "nucmer.snps"
    snps_path.write_text(
        "\n".join(
            [
                "1\tG\tY\t0\t0\tquery_1",
                "2\tA\tT\t0\t0\tquery_1",
                "3\t.\tC\t0\t0\tquery_1",
                "4\tT\t.\t0\t0\tquery_1",
                "5\tT\tR\t0\t0\tquery_1",
                "6\tA\tN\t0\t0\tquery_1",
            ]
        )
        + "\n"
    )
    meta = {
        "ortholog_gene_id": "101",
        "tax_id": "10090",
        "taxname": "Mus musculus",
        "sequence_length": "10",
    }
    summaries = {"query_1": empty_summary("1", meta, 10)}

    events, ambiguous_count = parse_snps(
        snps_path,
        "1",
        {"genomic_accession": "NC_000001.11", "genomic_begin": "100"},
        {"query_1": meta},
        {"query_1": []},
        summaries,
    )

    assert [row["event_id"] for row in events] == [2, 3, 4]
    assert [
        (row["event_type"], row["ref"], row["alt"])
        for row in events
    ] == [
        ("snv", "A", "T"),
        ("ins", "", "C"),
        ("del", "T", ""),
    ]
    assert ambiguous_count == 3
    assert summaries["query_1"]["event_count"] == 3
    assert summaries["query_1"]["qc_flags"] == {
        "ambiguous_event_allele",
        "unfiltered_nucmer",
    }
