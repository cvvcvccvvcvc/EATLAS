#!/usr/bin/env python3
"""Create per-gene alignment task directories from fetch-stage outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable


TSV_NULL = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes-tsv", required=True, type=Path)
    parser.add_argument("--orthologs-tsv", required=True, type=Path)
    parser.add_argument("--taxonomy-presets", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--target-fasta", action="append", required=True, type=Path)
    parser.add_argument("--ortholog-fasta", action="append", required=True, type=Path)
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


def load_taxonomy(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv_gz(path)
    return {row.get("tax_id", ""): row for row in rows if row.get("tax_id")}


def target_metadata(gene_id: str, gene: dict[str, str], target_id: str) -> dict[str, object]:
    return {
        "gene_id": gene_id,
        "sequence_id": target_id,
        "genomic_accession": gene.get("genomic_accession", ""),
        "genomic_begin": gene.get("begin", ""),
        "genomic_end": gene.get("end", ""),
        "orientation": gene.get("orientation", ""),
        "sequence_orientation": gene.get("sequence_orientation", ""),
        "sequence_length": gene.get("sequence_length", ""),
        "sequence_sha256": gene.get("sequence_sha256", ""),
    }


def reconstructed_ortholog_header(row: dict[str, str]) -> str:
    taxname = row.get("taxname", "").replace(" ", "_")
    return (
        f"query_{row.get('query_gene_id', '')}|ortholog_gene_{row.get('ortholog_gene_id', '')}|"
        f"symbol={row.get('symbol', '')}|tax_id={row.get('tax_id', '')}|taxname={taxname}|"
        f"accession={row.get('accession', '')}|range={row.get('range_text', '')}|"
        f"orientation={row.get('orientation', '')}"
    )


def ortholog_metadata_row(gene_id: str, source_meta: dict[str, str], taxonomy: dict[str, dict[str, str]]) -> dict[str, object]:
    ortholog_gene_id = source_meta["ortholog_gene_id"]
    tax = taxonomy.get(source_meta.get("tax_id", ""), {})
    return {
        "gene_id": gene_id,
        "sequence_id": f"ortholog_{ortholog_gene_id}",
        "ortholog_gene_id": ortholog_gene_id,
        "tax_id": source_meta.get("tax_id", ""),
        "taxname": source_meta.get("taxname", ""),
        "symbol": source_meta.get("symbol", ""),
        "gene_type": source_meta.get("gene_type", ""),
        "accession": source_meta.get("accession", ""),
        "chromosome": source_meta.get("chromosome", ""),
        "begin": source_meta.get("begin", ""),
        "end": source_meta.get("end", ""),
        "orientation": source_meta.get("orientation", ""),
        "source_complement": source_meta.get("source_complement", ""),
        "sequence_length": source_meta.get("sequence_length", ""),
        "sequence_sha256": source_meta.get("sequence_sha256", ""),
        "taxonomy_group": tax.get("preset_group", "other_or_unknown"),
        "minimap2_preset": tax.get("minimap2_preset", "asm20"),
        "original_header": reconstructed_ortholog_header(source_meta),
    }


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    tasks_dir = args.outdir / "tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)

    genes = {row["gene_id"]: row for row in read_tsv_gz(args.genes_tsv)}
    ortholog_rows = read_tsv_gz(args.orthologs_tsv)
    orthologs_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in ortholog_rows:
        orthologs_by_gene[row["query_gene_id"]].append(row)

    taxonomy = load_taxonomy(args.taxonomy_presets)
    target_fastas = {gene_id_from_fasta_path(path): path for path in args.target_fasta}
    ortholog_fastas = {gene_id_from_fasta_path(path): path for path in args.ortholog_fasta}

    task_rows: list[dict[str, object]] = []
    for gene_id in sorted(genes, key=lambda value: int(value) if value.isdigit() else value):
        gene = genes[gene_id]
        target_path = target_fastas.get(gene_id)
        ortholog_path = ortholog_fastas.get(gene_id)
        ortholog_meta_by_id = {
            row["ortholog_gene_id"]: row for row in orthologs_by_gene.get(gene_id, [])
        }

        task_row = {
            "gene_id": gene_id,
            "symbol": gene.get("symbol", ""),
            "target_fasta": str(target_path or ""),
            "ortholog_fasta": str(ortholog_path or ""),
            "ortholog_count": len(ortholog_meta_by_id),
            "target_length": gene.get("sequence_length", ""),
            "status": "ready",
            "message": "",
        }
        if target_path is None:
            task_row.update({"status": "missing_target_fasta", "message": "No target FASTA for gene"})
            task_rows.append(task_row)
            continue
        if ortholog_path is None:
            task_row.update({"status": "missing_ortholog_fasta", "message": "No ortholog FASTA for gene"})
            task_rows.append(task_row)
            continue
        if not ortholog_meta_by_id:
            task_row.update({"status": "no_ortholog_metadata", "message": "No ortholog metadata rows for gene"})
            task_rows.append(task_row)
            continue

        task_dir = tasks_dir / f"task_{gene_id}"
        task_dir.mkdir(parents=True, exist_ok=True)
        target_id = f"target_{gene_id}"
        target_meta = target_metadata(gene_id, gene, target_id)
        ortholog_meta_rows = [
            ortholog_metadata_row(gene_id, source_meta, taxonomy)
            for source_meta in orthologs_by_gene.get(gene_id, [])
        ]

        ortholog_fields = [
            "gene_id",
            "sequence_id",
            "ortholog_gene_id",
            "tax_id",
            "taxname",
            "symbol",
            "gene_type",
            "accession",
            "chromosome",
            "begin",
            "end",
            "orientation",
            "source_complement",
            "sequence_length",
            "sequence_sha256",
            "taxonomy_group",
            "minimap2_preset",
            "original_header",
        ]
        write_tsv(task_dir / "orthologs.metadata.tsv", ortholog_fields, ortholog_meta_rows)
        manifest = {
            "gene_id": gene_id,
            "symbol": gene.get("symbol", ""),
            "target_id": target_id,
            "target_length": gene.get("sequence_length", ""),
            "target": target_meta,
            "ortholog_count": len(ortholog_meta_rows),
            "ortholog_metadata": "orthologs.metadata.tsv",
        }
        (task_dir / "task.json").write_text(json.dumps(manifest, indent=2) + "\n")
        task_row.update(
            {
                "target_fasta": str(Path("sequences") / "targets" / target_path.name),
                "ortholog_fasta": str(Path("sequences") / "orthologs" / ortholog_path.name),
                "ortholog_count": len(ortholog_meta_rows),
                "target_length": gene.get("sequence_length", ""),
            }
        )
        task_rows.append(task_row)

    task_fields = [
        "gene_id",
        "symbol",
        "target_fasta",
        "ortholog_fasta",
        "ortholog_count",
        "target_length",
        "status",
        "message",
    ]
    write_tsv_gz(args.outdir / "alignment_tasks.tsv.gz", task_fields, task_rows)


if __name__ == "__main__":
    main()
