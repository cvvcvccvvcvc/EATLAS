#!/usr/bin/env python3
"""Create per-gene alignment task directories from fetch-stage outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from itertools import groupby
from pathlib import Path
from typing import Iterable


TSV_NULL = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes-tsv", required=True, type=Path)
    parser.add_argument("--orthologs-tsv", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sequences-dir", required=True, type=Path)
    parser.add_argument("--partition-size", required=True, type=int)
    return parser.parse_args()


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, TSV_NULL) for field in fields})
            count += 1
    return count


def write_tsv_gz(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, TSV_NULL) for field in fields})
            count += 1
    return count


def read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def gene_id_from_fasta_path(path: Path) -> str:
    name = path.name
    if name.endswith(".fa.gz"):
        return name[:-6]
    if name.endswith(".fasta.gz"):
        return name[:-9]
    if name.endswith(".fa"):
        return name[:-3]
    if name.endswith(".fasta"):
        return name[:-6]
    return path.stem


def fasta_paths_by_gene(directory: Path) -> dict[str, Path]:
    if not directory.is_dir():
        raise NotADirectoryError(f"Sequence directory does not exist: {directory}")
    paths = {}
    for path in directory.iterdir():
        if path.is_file() and path.name.endswith((".fa.gz", ".fasta.gz", ".fa", ".fasta")):
            gene_id = gene_id_from_fasta_path(path)
            if gene_id in paths:
                raise ValueError(f"Multiple FASTA files found for gene {gene_id} in {directory}")
            paths[gene_id] = path
    return paths


def iter_ortholog_groups(path: Path) -> Iterable[tuple[str, list[dict[str, str]]]]:
    required = {
        "query_gene_id",
        "ortholog_gene_id",
        "tax_id",
        "taxname",
        "sequence_length",
    }
    seen: set[str] = set()
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Ortholog table {path} missing required columns: "
                + ", ".join(sorted(missing))
            )
        for gene_id, group in groupby(reader, key=lambda row: row["query_gene_id"]):
            if not gene_id:
                raise ValueError(f"Ortholog table {path} contains an empty query_gene_id")
            if gene_id in seen:
                raise ValueError(
                    f"Ortholog table {path} is not grouped by query_gene_id: {gene_id}"
                )
            seen.add(gene_id)
            yield gene_id, list(group)


def target_metadata(gene: dict[str, str], target_id: str) -> dict[str, object]:
    return {
        "sequence_id": target_id,
        "genomic_accession": gene.get("genomic_accession", ""),
        "genomic_begin": gene.get("begin", ""),
        "sequence_length": gene.get("sequence_length", ""),
    }


def chromosome_sort_key(value: str) -> tuple[int, str]:
    chromosome = str(value or "").removeprefix("chr")
    if chromosome.isdigit():
        return int(chromosome), chromosome
    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(chromosome, 10**6), chromosome


def partition_ids(genes: dict[str, dict[str, str]], partition_size: int) -> dict[str, str]:
    if partition_size <= 0:
        raise ValueError("--partition-size must be a positive integer")

    def genomic_key(gene_id: str) -> tuple[tuple[int, str], int, int, tuple[int, str]]:
        gene = genes[gene_id]
        begin = int(gene.get("begin") or 0)
        end = int(gene.get("end") or 0)
        gene_key = (0, gene_id.zfill(20)) if gene_id.isdigit() else (1, gene_id)
        return chromosome_sort_key(gene.get("chromosome", "")), min(begin, end), max(begin, end), gene_key

    ordered_gene_ids = sorted(genes, key=genomic_key)
    return {
        gene_id: f"partition_{index // partition_size + 1:06d}"
        for index, gene_id in enumerate(ordered_gene_ids)
    }


def ortholog_metadata_row(source_meta: dict[str, str]) -> dict[str, object]:
    ortholog_gene_id = source_meta["ortholog_gene_id"]
    return {
        "sequence_id": f"ortholog_{ortholog_gene_id}",
        "ortholog_gene_id": ortholog_gene_id,
        "tax_id": source_meta.get("tax_id", ""),
        "taxname": source_meta.get("taxname", ""),
        "sequence_length": source_meta.get("sequence_length", ""),
    }


def prepare_gene_task(
    tasks_dir: Path,
    gene_id: str,
    gene: dict[str, str],
    source_orthologs: list[dict[str, str]],
    target_path: Path | None,
    ortholog_path: Path | None,
    partition_id: str,
) -> dict[str, object]:
    ortholog_ids = [row["ortholog_gene_id"] for row in source_orthologs]
    if any(not value for value in ortholog_ids):
        raise ValueError(f"Ortholog metadata for gene {gene_id} contains an empty ortholog_gene_id")
    if len(ortholog_ids) != len(set(ortholog_ids)):
        raise ValueError(
            f"Ortholog metadata for gene {gene_id} contains duplicate ortholog_gene_id values"
        )
    ortholog_sequence_bp = 0
    for row in source_orthologs:
        ortholog_id = row["ortholog_gene_id"]
        if not row.get("tax_id"):
            raise ValueError(f"Ortholog {ortholog_id} for gene {gene_id} has no tax_id")
        try:
            sequence_length = int(row["sequence_length"])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Ortholog {ortholog_id} for gene {gene_id} has invalid sequence_length"
            ) from exc
        if sequence_length <= 0:
            raise ValueError(
                f"Ortholog {ortholog_id} for gene {gene_id} has invalid sequence_length"
            )
        ortholog_sequence_bp += sequence_length
    ortholog_meta_by_id = {row["ortholog_gene_id"]: row for row in source_orthologs}
    target_ready = target_path is not None
    ortholog_ready = target_ready and ortholog_path is not None and bool(ortholog_meta_by_id)
    task_row = {
        "gene_id": gene_id,
        "partition_id": partition_id,
        "ortholog_sequence_bp": ortholog_sequence_bp,
        "target_ready": str(target_ready).lower(),
        "ortholog_ready": str(ortholog_ready).lower(),
        "status": "ready",
        "message": "",
    }
    if target_path is None:
        task_row.update({"status": "missing_target_fasta", "message": "No target FASTA for gene"})
        return task_row
    if ortholog_path is None:
        task_row.update({"status": "missing_ortholog_fasta", "message": "No ortholog FASTA for gene"})
    elif not ortholog_meta_by_id:
        task_row.update({"status": "no_ortholog_metadata", "message": "No ortholog metadata rows for gene"})

    task_dir = tasks_dir / f"task_{gene_id}"
    task_dir.mkdir(parents=True, exist_ok=True)
    target_id = f"target_{gene_id}"
    ortholog_meta_rows = [ortholog_metadata_row(source_meta) for source_meta in source_orthologs]
    ortholog_fields = [
        "sequence_id",
        "ortholog_gene_id",
        "tax_id",
        "taxname",
        "sequence_length",
    ]
    write_tsv(task_dir / "orthologs.metadata.tsv", ortholog_fields, ortholog_meta_rows)
    manifest = {
        "gene_id": gene_id,
        "target": target_metadata(gene, target_id),
    }
    (task_dir / "task.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return task_row


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    tasks_dir = args.outdir / "tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)

    genes = {row["gene_id"]: row for row in read_tsv_gz(args.genes_tsv)}
    gene_partitions = partition_ids(genes, args.partition_size)
    target_fastas = fasta_paths_by_gene(args.sequences_dir / "targets")
    ortholog_fastas = fasta_paths_by_gene(args.sequences_dir / "orthologs")
    for label, observed in (
        ("target FASTA", set(target_fastas)),
        ("ortholog FASTA", set(ortholog_fastas)),
    ):
        unexpected = observed - set(genes)
        if unexpected:
            raise ValueError(
                f"{label} inputs contain unknown gene_id values: "
                + ", ".join(sorted(unexpected))
            )

    task_rows: list[dict[str, object]] = []
    processed_gene_ids: set[str] = set()
    for gene_id, source_orthologs in iter_ortholog_groups(args.orthologs_tsv):
        if gene_id not in genes:
            raise ValueError(f"Ortholog table contains unknown query_gene_id: {gene_id}")
        task_rows.append(
            prepare_gene_task(
                tasks_dir,
                gene_id,
                genes[gene_id],
                source_orthologs,
                target_fastas.get(gene_id),
                ortholog_fastas.get(gene_id),
                gene_partitions[gene_id],
            )
        )
        processed_gene_ids.add(gene_id)

    for gene_id in set(genes) - processed_gene_ids:
        task_rows.append(
            prepare_gene_task(
                tasks_dir,
                gene_id,
                genes[gene_id],
                [],
                target_fastas.get(gene_id),
                ortholog_fastas.get(gene_id),
                gene_partitions[gene_id],
            )
        )

    task_rows.sort(
        key=lambda row: (
            (0, int(row["gene_id"]))
            if str(row["gene_id"]).isdigit()
            else (1, str(row["gene_id"]))
        )
    )
    task_fields = [
        "gene_id",
        "partition_id",
        "ortholog_sequence_bp",
        "target_ready",
        "ortholog_ready",
        "status",
        "message",
    ]
    write_tsv_gz(args.outdir / "alignment_tasks.tsv.gz", task_fields, task_rows)


if __name__ == "__main__":
    main()
