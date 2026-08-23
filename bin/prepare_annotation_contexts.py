#!/usr/bin/env python3
"""Materialize the target context required by each alignment evidence partition."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-evidence-dir", required=True, type=Path)
    parser.add_argument("--genes-tsv", required=True, type=Path)
    parser.add_argument("--target-sequences-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def load_genes(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(f"Target genes table not found: {path}")
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if "gene_id" not in fields:
            raise ValueError(f"Target genes table missing gene_id: {path}")
        genes: dict[str, dict[str, str]] = {}
        for row in reader:
            gene_id = str(row.get("gene_id") or "")
            if not gene_id:
                raise ValueError(f"Target genes table contains an empty gene_id: {path}")
            if gene_id in genes:
                raise ValueError(f"Target genes table contains duplicate gene_id={gene_id}: {path}")
            genes[gene_id] = row
    return fields, genes


def partition_gene_ids(partition: Path) -> list[str]:
    manifest_path = partition / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Alignment partition missing manifest.json: {partition}")
    manifest = json.loads(manifest_path.read_text())
    partition_id = str(manifest.get("partition_id") or "")
    if partition_id != partition.name:
        raise ValueError(
            f"Alignment partition identity mismatch: directory={partition.name!r}, "
            f"manifest={partition_id!r}"
        )
    if manifest.get("schema") != "normalized_alignment_evidence_partition_v2":
        raise ValueError(f"Alignment partition has unsupported schema: {partition}")
    gene_ids = manifest.get("gene_ids")
    if (
        not isinstance(gene_ids, list)
        or not gene_ids
        or not all(isinstance(gene_id, str) and gene_id for gene_id in gene_ids)
        or len(set(gene_ids)) != len(gene_ids)
    ):
        raise ValueError(f"Alignment partition has invalid gene_ids: {partition}")
    if manifest.get("gene_count") != len(gene_ids):
        raise ValueError(f"Alignment partition gene_count does not match gene_ids: {partition}")
    return gene_ids


def write_genes(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    partitions_root = args.alignment_evidence_dir / "partitions"
    if not partitions_root.is_dir():
        raise NotADirectoryError(
            f"Alignment evidence partitions directory not found: {partitions_root}"
        )
    if not args.target_sequences_dir.is_dir():
        raise NotADirectoryError(
            f"Target sequences directory not found: {args.target_sequences_dir}"
        )
    partitions = sorted(path for path in partitions_root.iterdir() if path.is_dir())
    if not partitions:
        raise ValueError(f"No alignment evidence partitions found in {partitions_root}")

    fields, genes = load_genes(args.genes_tsv)
    seen_gene_ids: set[str] = set()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for partition in partitions:
        gene_ids = partition_gene_ids(partition)
        duplicate_gene_ids = seen_gene_ids.intersection(gene_ids)
        if duplicate_gene_ids:
            raise ValueError(
                "Genes occur in multiple alignment evidence partitions: "
                + ", ".join(sorted(duplicate_gene_ids))
            )
        seen_gene_ids.update(gene_ids)

        missing_genes = [gene_id for gene_id in gene_ids if gene_id not in genes]
        if missing_genes:
            raise ValueError(
                f"Alignment partition {partition.name} references genes absent from genes.tsv.gz: "
                + ", ".join(missing_genes)
            )
        context_dir = args.outdir / partition.name
        targets_dir = context_dir / "targets"
        targets_dir.mkdir(parents=True, exist_ok=False)
        write_genes(
            context_dir / "genes.tsv.gz",
            fields,
            [genes[gene_id] for gene_id in gene_ids],
        )
        for gene_id in gene_ids:
            source = args.target_sequences_dir / f"{gene_id}.fa.gz"
            if not source.exists():
                raise FileNotFoundError(
                    f"Target FASTA not found for gene {gene_id}: {source}"
                )
            shutil.copy2(source, targets_dir / source.name)


if __name__ == "__main__":
    main()
