from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pytest


PROJECT_DIR = Path(__file__).resolve().parents[1]
PREPARE_SCRIPT = PROJECT_DIR / "bin" / "prepare_alignment_tasks.py"

from bin.alignment_task_io import load_task_context, materialize_task_fastas
from bin.prepare_alignment_tasks import iter_ortholog_groups


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_alignment_tasks_capture_target_and_ortholog_readiness(tmp_path: Path) -> None:
    genes = tmp_path / "genes.tsv.gz"
    gene_fields = [
        "gene_id",
        "symbol",
        "genomic_accession",
        "chromosome",
        "begin",
        "end",
        "orientation",
        "sequence_orientation",
        "sequence_length",
        "sequence_sha256",
    ]
    write_tsv_gz(
        genes,
        gene_fields,
        [
            {
                "gene_id": gene_id,
                "symbol": f"G{gene_id}",
                "genomic_accession": "NC_000001.11",
                "chromosome": "1",
                "begin": str(begin),
                "end": str(begin + 3),
                "orientation": "plus",
                "sequence_orientation": "plus",
                "sequence_length": "4",
                "sequence_sha256": "checksum",
            }
            for gene_id, begin in (("1", 100), ("2", 200))
        ],
    )

    orthologs = tmp_path / "orthologs.tsv.gz"
    ortholog_fields = [
        "query_gene_id",
        "ortholog_gene_id",
        "tax_id",
        "taxname",
        "symbol",
        "gene_type",
        "accession",
        "chromosome",
        "begin",
        "end",
        "range_text",
        "orientation",
        "source_complement",
        "sequence_length",
        "sequence_sha256",
    ]
    write_tsv_gz(
        orthologs,
        ortholog_fields,
        [
            {
                "query_gene_id": "1",
                "ortholog_gene_id": "101",
                "tax_id": "10090",
                "taxname": "Mus musculus",
                "sequence_length": "4",
            },
            {
                "query_gene_id": "1",
                "ortholog_gene_id": "102",
                "tax_id": "10116",
                "taxname": "Rattus norvegicus",
                "sequence_length": "6",
            },
            {
                "query_gene_id": "2",
                "ortholog_gene_id": "102",
                "tax_id": "10090",
                "taxname": "Mus musculus",
                "sequence_length": "4",
            },
        ],
    )

    sequences = tmp_path / "sequences"
    for directory, gene_ids in (("targets", ("1", "2")), ("orthologs", ("1",))):
        output_dir = sequences / directory
        output_dir.mkdir(parents=True)
        for gene_id in gene_ids:
            header = (
                f"target_{gene_id}"
                if directory == "targets"
                else f"query_{gene_id}|ortholog_gene_10{gene_id}|tax_id=10090"
            )
            with gzip.open(output_dir / f"{gene_id}.fa.gz", "wt") as handle:
                handle.write(f">{header}\nACGT\n")
                if directory == "orthologs" and gene_id == "1":
                    handle.write(">query_1|ortholog_gene_102|tax_id=10116\nACGTAC\n")

    outdir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(PREPARE_SCRIPT),
            "--genes-tsv",
            str(genes),
            "--orthologs-tsv",
            str(orthologs),
            "--outdir",
            str(outdir),
            "--sequences-dir",
            str(sequences),
            "--partition-size",
            "10",
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    task_rows = {
        row["gene_id"]: row
        for row in read_tsv_gz(outdir / "alignment_tasks.tsv.gz")
    }
    assert {gene_id: row["status"] for gene_id, row in task_rows.items()} == {
        "1": "ready",
        "2": "missing_ortholog_fasta",
    }
    assert {
        gene_id: (row["target_ready"], row["ortholog_ready"])
        for gene_id, row in task_rows.items()
    } == {
        "1": ("true", "true"),
        "2": ("true", "false"),
    }
    assert {
        gene_id: row["ortholog_sequence_bp"] for gene_id, row in task_rows.items()
    } == {"1": "10", "2": "4"}
    assert (outdir / "tasks" / "task_2" / "task.json").exists()

    task_dir = outdir / "tasks" / "task_1"
    manifest, target, ortholog_metadata = load_task_context(task_dir)
    assert manifest == {
        "gene_id": "1",
        "target": {
            "sequence_id": "target_1",
            "genomic_accession": "NC_000001.11",
            "genomic_begin": "100",
            "sequence_length": "4",
        },
    }
    assert ortholog_metadata == [
        {
            "sequence_id": "ortholog_101",
            "ortholog_gene_id": "101",
            "tax_id": "10090",
            "taxname": "Mus musculus",
            "sequence_length": "4",
        },
        {
            "sequence_id": "ortholog_102",
            "ortholog_gene_id": "102",
            "tax_id": "10116",
            "taxname": "Rattus norvegicus",
            "sequence_length": "6",
        },
    ]
    with (task_dir / "orthologs.metadata.tsv").open(newline="") as handle:
        assert next(csv.reader(handle, delimiter="\t")) == [
            "sequence_id",
            "ortholog_gene_id",
            "tax_id",
            "taxname",
            "sequence_length",
        ]
    with gzip.open(outdir / "alignment_tasks.tsv.gz", "rt", newline="") as handle:
        assert next(csv.reader(handle, delimiter="\t")) == [
            "gene_id",
            "partition_id",
            "ortholog_sequence_bp",
            "target_ready",
            "ortholog_ready",
            "status",
            "message",
        ]
    target_fasta, ortholog_fasta = materialize_task_fastas(
        sequences / "targets" / "1.fa.gz",
        sequences / "orthologs" / "1.fa.gz",
        manifest,
        ortholog_metadata,
        tmp_path / "materialized",
    )
    assert target_fasta.read_text() == ">target_1\nACGT\n"
    assert ortholog_fasta.read_text() == (
        ">ortholog_101\nACGT\n"
        ">ortholog_102\nACGTAC\n"
    )


def test_ortholog_groups_reject_noncontiguous_gene_rows(tmp_path: Path) -> None:
    path = tmp_path / "orthologs.tsv.gz"
    fields = [
        "query_gene_id",
        "ortholog_gene_id",
        "tax_id",
        "taxname",
        "sequence_length",
    ]
    write_tsv_gz(
        path,
        fields,
        [
            {"query_gene_id": "1", "ortholog_gene_id": "101"},
            {"query_gene_id": "2", "ortholog_gene_id": "201"},
            {"query_gene_id": "1", "ortholog_gene_id": "102"},
        ],
    )

    with pytest.raises(ValueError, match="not grouped by query_gene_id: 1"):
        list(iter_ortholog_groups(path))
