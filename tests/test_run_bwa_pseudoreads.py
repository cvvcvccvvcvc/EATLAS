from __future__ import annotations

import sys
from pathlib import Path


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from bam_filtering_v1 import expected_pseudoreads, pseudoread_starts  # noqa: E402
from run_bwa_pseudoreads import generate_pseudoreads  # noqa: E402


def test_pseudoread_starts_include_final_window() -> None:
    assert pseudoread_starts(109, read_len=75, step=35) == [0, 34]
    assert pseudoread_starts(110, read_len=75, step=35) == [0, 35]
    assert pseudoread_starts(111, read_len=75, step=35) == [0, 35, 36]


def test_short_sequence_count_matches_generation() -> None:
    assert pseudoread_starts(19, read_len=75, step=35) == []
    assert pseudoread_starts(20, read_len=75, step=35) == [0]
    assert expected_pseudoreads(20, read_len=75, step=35) == 1


def test_generate_pseudoreads_uses_endpoint_inclusive_starts(tmp_path: Path) -> None:
    fasta = tmp_path / "orthologs.fa"
    fastq = tmp_path / "pseudoreads.fastq"
    fasta.write_text(">ortholog_1\n" + "A" * 109 + "\n")

    count = generate_pseudoreads(
        fasta,
        fastq,
        read_len=75,
        step=35,
        phred=30,
    )

    headers = fastq.read_text().splitlines()[::4]
    assert count == expected_pseudoreads(109, read_len=75, step=35)
    assert headers == [
        "@ortholog_1_pseudo_1_1-75",
        "@ortholog_1_pseudo_2_35-109",
    ]
