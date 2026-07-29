from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
PREPARE_SCRIPT = PROJECT_DIR / "bin" / "prepare_alignment_tasks.py"


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_partition_genes_contains_only_ready_alignment_tasks(tmp_path: Path) -> None:
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
                "query_gene_id": gene_id,
                "ortholog_gene_id": f"10{gene_id}",
                "tax_id": "10090",
            }
            for gene_id in ("1", "2")
        ],
    )

    target_features = tmp_path / "target_features.tsv.gz"
    write_tsv_gz(
        target_features,
        ["gene_id", "feature_type", "target_start0", "target_end0"],
        [
            {
                "gene_id": gene_id,
                "feature_type": "gene",
                "target_start0": "0",
                "target_end0": "4",
            }
            for gene_id in ("1", "2")
        ],
    )
    taxonomy = tmp_path / "taxonomy.tsv.gz"
    write_tsv_gz(taxonomy, ["tax_id", "preset_group", "minimap2_preset"], [])
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps({"orthologs_selected_grouped_by_query_gene_id": True}) + "\n"
    )

    sequences = tmp_path / "sequences"
    for directory, gene_ids in (("targets", ("1", "2")), ("orthologs", ("1",))):
        output_dir = sequences / directory
        output_dir.mkdir(parents=True)
        for gene_id in gene_ids:
            (output_dir / f"{gene_id}.fa.gz").write_text(">sequence\nACGT\n")

    outdir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(PREPARE_SCRIPT),
            "--genes-tsv",
            str(genes),
            "--orthologs-tsv",
            str(orthologs),
            "--fetch-manifest",
            str(manifest),
            "--target-features",
            str(target_features),
            "--taxonomy-presets",
            str(taxonomy),
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
    partition_rows = read_tsv_gz(
        outdir / "partition_genes" / "partition_000001.tsv.gz"
    )
    assert [row["gene_id"] for row in partition_rows] == ["1"]
    target_partition_rows = read_tsv_gz(
        outdir / "target_partition_genes" / "partition_000001.tsv.gz"
    )
    assert [row["gene_id"] for row in target_partition_rows] == ["1", "2"]
    assert (outdir / "tasks" / "task_2" / "task.json").exists()
