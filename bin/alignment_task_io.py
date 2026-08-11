"""Shared helpers for per-gene alignment task inputs."""

from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path


ORTHOLOG_GENE_RE = re.compile(r"(?:^|\|)ortholog_gene_(\d+)(?:\||$)|^ortholog_(\d+)(?:\s|$)")
FASTA_WIDTH = 80
TASK_FIELDS = {"gene_id", "target"}
TARGET_FIELDS = {"sequence_id", "genomic_accession", "genomic_begin", "sequence_length"}
ORTHOLOG_FIELDS = {
    "sequence_id",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "sequence_length",
}


def read_tsv(path: Path, required_fields: set[str]) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required_fields - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Task table {path} missing required columns: "
                + ", ".join(sorted(missing))
            )
        return [dict(row) for row in reader]


def iter_fasta(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header = None
    seq_parts: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            elif header is not None:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts)


def write_fasta_record(handle, sequence_id: str, seq: str) -> None:
    handle.write(f">{sequence_id}\n")
    for index in range(0, len(seq), FASTA_WIDTH):
        handle.write(seq[index : index + FASTA_WIDTH] + "\n")


def parse_ortholog_gene_id(header: str) -> str:
    match = ORTHOLOG_GENE_RE.search(header)
    if not match:
        return ""
    return next(group for group in match.groups() if group)


def load_task_context(task_dir: Path) -> tuple[dict[str, object], dict[str, str], list[dict[str, str]]]:
    """Load the metadata-only task manifest used by alignment strategies."""

    manifest = json.loads((task_dir / "task.json").read_text())
    missing = TASK_FIELDS - set(manifest)
    if missing:
        raise ValueError(
            f"Task manifest {task_dir / 'task.json'} missing fields: "
            + ", ".join(sorted(missing))
        )
    if not manifest["gene_id"]:
        raise ValueError(f"Task manifest {task_dir / 'task.json'} has an empty gene_id")
    if not isinstance(manifest["target"], dict):
        raise ValueError(f"Task manifest {task_dir / 'task.json'} target must be an object")
    target_meta = dict(manifest["target"])
    missing_target = TARGET_FIELDS - set(target_meta)
    if missing_target:
        raise ValueError(
            f"Task manifest {task_dir / 'task.json'} target missing fields: "
            + ", ".join(sorted(missing_target))
        )
    ortholog_meta = read_tsv(task_dir / "orthologs.metadata.tsv", ORTHOLOG_FIELDS)
    ortholog_ids = [row["ortholog_gene_id"] for row in ortholog_meta]
    sequence_ids = [row["sequence_id"] for row in ortholog_meta]
    if len(ortholog_ids) != len(set(ortholog_ids)):
        raise ValueError(f"Task {task_dir} contains duplicate ortholog_gene_id values")
    if len(sequence_ids) != len(set(sequence_ids)):
        raise ValueError(f"Task {task_dir} contains duplicate sequence_id values")
    return manifest, target_meta, ortholog_meta


def materialize_task_fastas(
    source_target_fasta: Path,
    source_ortholog_fasta: Path,
    manifest: dict[str, object],
    ortholog_meta: list[dict[str, str]],
    work_dir: Path,
) -> tuple[Path, Path]:
    """Write normalized uncompressed FASTA inputs for one aligner process."""

    work_dir.mkdir(parents=True, exist_ok=True)
    target_fasta = work_dir / "target.fa"
    ortholog_fasta = work_dir / "orthologs.fa"

    target_records = list(iter_fasta(source_target_fasta))
    if len(target_records) != 1:
        raise ValueError(f"Expected one target FASTA record in {source_target_fasta}, found {len(target_records)}")
    target_meta = manifest["target"]
    target_id = str(target_meta["sequence_id"])
    expected_target_length = int(target_meta["sequence_length"])
    if expected_target_length != len(target_records[0][1]):
        raise ValueError(
            f"Target length mismatch in {source_target_fasta}: "
            f"metadata={expected_target_length}, fasta={len(target_records[0][1])}"
        )
    with target_fasta.open("w") as handle:
        write_fasta_record(handle, target_id, target_records[0][1])

    expected_by_ortholog = {row["ortholog_gene_id"]: row for row in ortholog_meta}
    seen: set[str] = set()
    with ortholog_fasta.open("w") as handle:
        for header, seq in iter_fasta(source_ortholog_fasta):
            ortholog_gene_id = parse_ortholog_gene_id(header)
            row = expected_by_ortholog.get(ortholog_gene_id)
            if row is None:
                continue
            sequence_id = row["sequence_id"]
            expected_length = int(row["sequence_length"])
            if expected_length != len(seq):
                raise ValueError(
                    f"Ortholog {ortholog_gene_id} length mismatch in {source_ortholog_fasta}: "
                    f"metadata={expected_length}, fasta={len(seq)}"
                )
            write_fasta_record(handle, sequence_id, seq)
            seen.add(ortholog_gene_id)

    missing = sorted(set(expected_by_ortholog) - seen, key=lambda value: int(value) if value.isdigit() else value)
    if missing:
        preview = ", ".join(missing[:10])
        suffix = "..." if len(missing) > 10 else ""
        raise ValueError(f"Source ortholog FASTA is missing {len(missing)} selected records: {preview}{suffix}")

    return target_fasta, ortholog_fasta
