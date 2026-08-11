from __future__ import annotations

import sys
import subprocess
from pathlib import Path

import pytest


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from bam_filtering_v1 import expected_pseudoreads, pseudoread_starts  # noqa: E402
import run_bwa_pseudoreads as bwa_runner  # noqa: E402
from run_bwa_pseudoreads import generate_pseudoreads  # noqa: E402


def test_bwa_cli_uses_fixed_pseudoread_defaults(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_bwa_pseudoreads.py",
            "--task-dir",
            "task",
            "--source-target-fasta",
            "target.fa.gz",
            "--source-ortholog-fasta",
            "ortholog.fa.gz",
            "--outdir",
            "out",
        ],
    )

    args = bwa_runner.parse_args()

    assert (
        args.pseudoread_len,
        args.pseudoread_step,
        args.pseudoread_phred,
    ) == (150, 75, 30)


def test_pseudoread_starts_include_final_window() -> None:
    assert pseudoread_starts(224, read_len=150, step=75) == [0, 74]
    assert pseudoread_starts(225, read_len=150, step=75) == [0, 75]
    assert pseudoread_starts(226, read_len=150, step=75) == [0, 75, 76]


def test_short_sequence_count_matches_generation() -> None:
    assert pseudoread_starts(19, read_len=150, step=75) == []
    assert pseudoread_starts(20, read_len=150, step=75) == [0]
    assert expected_pseudoreads(20, read_len=150, step=75) == 1


def test_generate_pseudoreads_uses_endpoint_inclusive_starts(tmp_path: Path) -> None:
    fasta = tmp_path / "orthologs.fa"
    fastq = tmp_path / "pseudoreads.fastq"
    fasta.write_text(">ortholog_1\n" + "A" * 224 + "\n")

    generation = generate_pseudoreads(
        fasta,
        fastq,
        read_len=150,
        step=75,
        phred=30,
    )

    headers = fastq.read_text().splitlines()[::4]
    assert generation.total_reads == expected_pseudoreads(224, read_len=150, step=75)
    assert generation.query_lengths == {"1": 224}
    assert generation.generated_counts == {"ortholog_1": 2}
    assert headers == [
        "@ortholog_1_pseudo_1_1-150",
        "@ortholog_1_pseudo_2_75-224",
    ]


def test_bwa_pipeline_uses_declared_cpu_budget_without_samtools_view(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    popen_calls: list[list[str]] = []
    checked_calls: list[list[str]] = []

    class FakeStdout:
        def close(self) -> None:
            pass

    class FakeProcess:
        def __init__(self, cmd, **kwargs):
            popen_calls.append(cmd)
            self.stdout = FakeStdout() if kwargs.get("stdout") == subprocess.PIPE else None

        def wait(self) -> int:
            return 0

    monkeypatch.setattr(bwa_runner.subprocess, "Popen", FakeProcess)
    monkeypatch.setattr(
        bwa_runner,
        "run_checked",
        lambda cmd, **_kwargs: checked_calls.append(cmd),
    )

    bwa_threads = bwa_runner.run_bwa_mem_pipeline(
        "bwa",
        "samtools",
        tmp_path / "target.fa",
        tmp_path / "reads.fastq",
        tmp_path / "sorted.bam",
        threads=3,
    )

    assert bwa_threads == 2
    assert popen_calls == [
        [
            "bwa",
            "mem",
            "-t",
            "2",
            str(tmp_path / "target.fa"),
            str(tmp_path / "reads.fastq"),
        ],
        ["samtools", "sort", "-o", str(tmp_path / "sorted.bam"), "-"],
    ]
    assert checked_calls == [
        ["bwa", "index", str(tmp_path / "target.fa")],
        ["samtools", "index", str(tmp_path / "sorted.bam")],
    ]


def test_bwa_pipeline_rejects_single_cpu_budget(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="at least 2"):
        bwa_runner.run_bwa_mem_pipeline(
            "bwa",
            "samtools",
            tmp_path / "target.fa",
            tmp_path / "reads.fastq",
            tmp_path / "sorted.bam",
            threads=1,
        )


def test_scan_bam_deduplicates_event_support_by_ortholog(tmp_path: Path) -> None:
    bam_path = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "target", "LN": 100}]}

    def make_read(name: str, *, secondary: bool = False) -> bwa_runner.pysam.AlignedSegment:
        read = bwa_runner.pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = "A" * 20 + "G" + "A" * 9
        read.flag = 256 if secondary else 0
        read.reference_id = 0
        read.reference_start = 0
        read.mapping_quality = 60
        read.cigar = ((0, 30),)
        read.query_qualities = bwa_runner.pysam.qualitystring_to_array("I" * 30)
        read.set_tag("NM", 1)
        return read

    with bwa_runner.pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        bam.write(make_read("ortholog_101_pseudo_0_1-30", secondary=True))
        bam.write(make_read("ortholog_101_pseudo_2_1-30"))
        bam.write(make_read("ortholog_101_pseudo_1_1-30"))
        bam.write(make_read("ortholog_102_pseudo_1_1-30"))
        bam.write(make_read("ortholog_103_pseudo_1_1-30", secondary=True))
    bwa_runner.pysam.index(str(bam_path))

    _segments, event_support = bwa_runner.scan_bam(bam_path, "A" * 100)
    support = event_support[("snv", 20, 21, "A", "G")]
    rows = bwa_runner.make_bwa_event_rows(
        event_support,
        gene_id="1",
        target_meta={"genomic_begin": "100"},
        target_acc="NC_1",
        ortholog_meta_by_id={
            "101": {"tax_id": "1"},
            "102": {"tax_id": "2"},
        },
    )

    assert sorted(support) == ["101", "102", "103"]
    assert support["101"]["native_record_id"] == "ortholog_101_pseudo_1_1-30"
    assert support["101"]["is_primary"] is True
    assert support["103"]["is_primary"] is False
    assert len(rows) == 3
    assert len({row["event_id"] for row in rows}) == 3
    assert {row["ortholog_gene_id"]: row["qc_flags"] for row in rows} == {
        "101": "bwa_cigar_event",
        "102": "bwa_cigar_event",
        "103": "bwa_cigar_event,non_primary",
    }
