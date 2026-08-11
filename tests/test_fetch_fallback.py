from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest


PROJECT_DIR = Path(__file__).resolve().parents[1]
FETCH_SCRIPT = PROJECT_DIR / "bin" / "fetch_parse_chunk.py"
BUILD_SCRIPT = PROJECT_DIR / "bin" / "build_fetch_dataset.py"
sys.path.insert(0, str(PROJECT_DIR / "bin"))

import fetch_parse_chunk as fetch_chunk  # noqa: E402


def read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_fake_datasets(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import json
import sys
import zipfile
from pathlib import Path

if sys.argv[1:] == ["--version"]:
    print("datasets version: test")
    raise SystemExit(0)

args = sys.argv[1:]
ids_path = Path(args[args.index("--inputfile") + 1])
output_path = Path(args[args.index("--filename") + 1])
gene_ids = ids_path.read_text().split()

if gene_ids != ["1"]:
    print("simulated transient download failure", file=sys.stderr)
    raise SystemExit(1)

report = {
    "geneId": "1",
    "taxId": 9606,
    "taxname": "Homo sapiens",
    "symbol": "GENE1",
    "type": "protein-coding",
    "orientation": "plus",
    "annotations": [{
        "assemblyAccession": "GCF_000001405.40",
        "genomicLocations": [{
            "genomicAccessionVersion": "NC_000001.11",
            "genomicRange": {"begin": 100, "end": 103}
        }]
    }]
}

with zipfile.ZipFile(output_path, "w") as archive:
    archive.writestr(
        "ncbi_dataset/data/data_report.jsonl",
        json.dumps(report) + "\\n",
    )
    archive.writestr(
        "ncbi_dataset/data/gene.fna",
        ">NC_000001.11:100-103 [GeneID=1] "
        "[organism=Homo sapiens] [chromosome=1]\\nACGT\\n",
    )
"""
    )
    path.chmod(0o755)


def test_batch_failure_falls_back_to_singletons_and_builds_partial_dataset(
    tmp_path: Path,
) -> None:
    ids_file = tmp_path / "chunk_000001.ids.txt"
    ids_file.write_text("1\n2\n")
    fake_datasets = tmp_path / "datasets"
    write_fake_datasets(fake_datasets)
    chunk_dir = tmp_path / "fetch_chunk_000001"

    fetch = subprocess.run(
        [
            sys.executable,
            str(FETCH_SCRIPT),
            "--ids-file",
            str(ids_file),
            "--outdir",
            str(chunk_dir),
            "--datasets-bin",
            str(fake_datasets),
            "--download-retries",
            "1",
            "--download-retry-base-seconds",
            "0",
        ],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
    )

    assert fetch.returncode == 0, fetch.stderr
    chunk_manifest = json.loads((chunk_dir / "manifest.json").read_text())
    assert chunk_manifest["status"] == "partial"
    assert chunk_manifest["download_mode"] == "singleton_fallback"
    assert chunk_manifest["batch_download_attempts"] == 2
    assert chunk_manifest["singleton_download_attempts"] == 3
    assert chunk_manifest["target_gene_count"] == 1
    assert chunk_manifest["failure_count"] == 1
    assert [row["gene_id"] for row in read_tsv_gz(chunk_dir / "genes.tsv.gz")] == ["1"]
    assert read_tsv_gz(chunk_dir / "failures.tsv.gz") == [
        {
            "gene_id": "2",
            "failure_type": "ncbi_download_failed",
            "message": (
                "attempts=2; exit_code=1; "
                "error=simulated transient download failure"
            ),
        }
    ]

    ids_tsv = tmp_path / "input.ids.tsv"
    ids_tsv.write_text(
        "input_index\traw_token\tgene_id\taccepted\treason\tduplicate_of_input_index\n"
        "1\t1\t1\ttrue\t\t\n"
        "2\t2\t2\ttrue\t\t\n"
    )
    chunks_tsv = tmp_path / "chunks.tsv"
    chunks_tsv.write_text(
        "chunk_id\tchunk_file\tgene_count\tgene_ids\n"
        "chunk_000001\tchunk_000001.ids.txt\t2\t1,2\n"
    )
    gff3 = tmp_path / "genomic.gff"
    gff3.write_text(
        "##gff-version 3\n"
        "NC_000001.11\tRefSeq\texon\t101\t103\t.\t+\t.\t"
        "ID=exon1;Dbxref=GeneID:1\n"
    )
    fetch_dataset = tmp_path / "fetch"

    build = subprocess.run(
        [
            sys.executable,
            str(BUILD_SCRIPT),
            "--ids-tsv",
            str(ids_tsv),
            "--chunks-tsv",
            str(chunks_tsv),
            "--outdir",
            str(fetch_dataset),
            "--target-annotation-gff3",
            str(gff3),
            "--chunk-dir",
            str(chunk_dir),
        ],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
    )

    assert build.returncode == 0, build.stderr
    manifest = json.loads((fetch_dataset / "manifest.json").read_text())
    assert manifest["status"] == "partial"
    assert manifest["target_gene_count"] == 1
    assert manifest["failure_count"] == 1
    assert manifest["download_failed_gene_count"] == 1
    assert manifest["singleton_fallback_chunk_count"] == 1


def test_download_retries_apply_to_singleton_fallback(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ids_file = tmp_path / "gene.ids.txt"
    ids_file.write_text("1\n")
    zip_path = tmp_path / "ncbi_dataset.zip"
    calls = 0

    def fake_download(
        _datasets_bin: str,
        _ids_file: Path,
        output_path: Path,
    ) -> subprocess.CompletedProcess[str]:
        nonlocal calls
        calls += 1
        if calls == 1:
            return subprocess.CompletedProcess([], 1, "", "transient")
        with zipfile.ZipFile(output_path, "w") as archive:
            archive.writestr("ncbi_dataset/data/data_report.jsonl", "{}\n")
            archive.writestr("ncbi_dataset/data/gene.fna", ">x\nA\n")
        return subprocess.CompletedProcess([], 0, "", "")

    monkeypatch.setattr(fetch_chunk, "run_datasets_download", fake_download)
    monkeypatch.setattr(fetch_chunk.time, "sleep", lambda _seconds: None)
    attempts, wait_seconds = fetch_chunk.download_package(
        "datasets",
        ids_file,
        zip_path,
        retries=1,
        retry_base_seconds=0,
    )

    assert attempts == 2
    assert wait_seconds == 0
    assert calls == 2
