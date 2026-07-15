"""Shared helpers for per-gene alignment task inputs."""

from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path


ORTHOLOG_GENE_RE = re.compile(r"(?:^|\|)ortholog_gene_(\d+)(?:\||$)|^ortholog_(\d+)(?:\s|$)")
FASTA_WIDTH = 80


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


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
    target_meta = dict(manifest["target"])
    ortholog_metadata_path = task_dir / str(manifest.get("ortholog_metadata", "orthologs.metadata.tsv"))
    return manifest, target_meta, read_tsv(ortholog_metadata_path)


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
    target_id = str(manifest.get("target_id") or f"target_{manifest['gene_id']}")
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
            sequence_id = row.get("sequence_id") or f"ortholog_{ortholog_gene_id}"
            expected_length = row.get("sequence_length")
            if expected_length and int(expected_length) != len(seq):
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
