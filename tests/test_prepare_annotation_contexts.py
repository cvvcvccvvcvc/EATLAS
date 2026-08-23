from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_DIR / "bin" / "prepare_annotation_contexts.py"


def write_tsv(path: Path, rows: list[list[str]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "genomic_accession", "chromosome", "begin", "end"])
        writer.writerows(rows)


def write_partition(root: Path, partition_id: str, gene_ids: list[str]) -> None:
    partition = root / "partitions" / partition_id
    partition.mkdir(parents=True)
    (partition / "manifest.json").write_text(
        json.dumps(
            {
                "schema": "normalized_alignment_evidence_partition_v1",
                "partition_id": partition_id,
                "gene_count": len(gene_ids),
                "gene_ids": gene_ids,
            }
        )
        + "\n"
    )


def run_prepare(
    evidence: Path,
    genes: Path,
    targets: Path,
    outdir: Path,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--alignment-evidence-dir",
            str(evidence),
            "--genes-tsv",
            str(genes),
            "--target-sequences-dir",
            str(targets),
            "--outdir",
            str(outdir),
        ],
        cwd=PROJECT_DIR,
        text=True,
        capture_output=True,
        check=False,
    )


def test_prepares_exact_context_for_each_evidence_partition(tmp_path: Path) -> None:
    evidence = tmp_path / "evidence"
    write_partition(evidence, "partition_000001", ["2"])
    write_partition(evidence, "partition_000002", ["1", "3"])
    genes = tmp_path / "genes.tsv.gz"
    write_tsv(
        genes,
        [
            ["1", "NC_1", "1", "10", "20"],
            ["2", "NC_2", "2", "30", "40"],
            ["3", "NC_3", "3", "50", "60"],
        ],
    )
    targets = tmp_path / "targets"
    targets.mkdir()
    for gene_id in ["1", "2", "3"]:
        (targets / f"{gene_id}.fa.gz").write_bytes(f"target-{gene_id}".encode())

    completed = run_prepare(evidence, genes, targets, tmp_path / "contexts")

    assert completed.returncode == 0, completed.stderr
    first = tmp_path / "contexts" / "partition_000001"
    second = tmp_path / "contexts" / "partition_000002"
    with gzip.open(first / "genes.tsv.gz", "rt", newline="") as handle:
        assert [row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")] == ["2"]
    with gzip.open(second / "genes.tsv.gz", "rt", newline="") as handle:
        assert [row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")] == ["1", "3"]
    assert {path.name for path in (first / "targets").iterdir()} == {"2.fa.gz"}
    assert {path.name for path in (second / "targets").iterdir()} == {
        "1.fa.gz",
        "3.fa.gz",
    }


def test_rejects_noncanonical_partition_schema(tmp_path: Path) -> None:
    evidence = tmp_path / "evidence"
    write_partition(evidence, "partition_000001", ["1"])
    manifest_path = evidence / "partitions" / "partition_000001" / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    del manifest["schema"]
    manifest_path.write_text(json.dumps(manifest) + "\n")
    genes = tmp_path / "genes.tsv.gz"
    write_tsv(genes, [["1", "NC_1", "1", "10", "20"]])
    targets = tmp_path / "targets"
    targets.mkdir()
    (targets / "1.fa.gz").write_bytes(b"target")

    completed = run_prepare(evidence, genes, targets, tmp_path / "contexts")

    assert completed.returncode != 0
    assert "unsupported schema" in completed.stderr
