#!/usr/bin/env python3
"""Create per-gene alignment task directories from fetch-stage outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from itertools import groupby
from pathlib import Path
from typing import Iterable


TSV_NULL = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes-tsv", required=True, type=Path)
    parser.add_argument("--orthologs-tsv", required=True, type=Path)
    parser.add_argument("--fetch-manifest", required=True, type=Path)
    parser.add_argument("--target-features", required=True, type=Path)
    parser.add_argument("--taxonomy-presets", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sequences-dir", required=True, type=Path)
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


def iter_tsv_gz(path: Path) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


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


def partition_target_features(path: Path, output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    partitions: dict[str, Path] = {}
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if "gene_id" not in fields:
            raise ValueError(f"Target features table missing gene_id column: {path}")
        for gene_id, rows in groupby(reader, key=lambda row: row["gene_id"]):
            if gene_id in partitions:
                raise ValueError(f"Target features are not grouped by gene_id in {path}: {gene_id}")
            partition = output_dir / f"{gene_id}.tsv.gz"
            write_tsv_gz(partition, fields, rows)
            partitions[gene_id] = partition
    return partitions


def iter_ortholog_groups(
    path: Path,
    grouped_by_gene: bool,
) -> Iterable[tuple[str, list[dict[str, str]]]]:
    rows = iter_tsv_gz(path)
    if grouped_by_gene:
        for gene_id, group in groupby(rows, key=lambda row: row["query_gene_id"]):
            yield gene_id, list(group)
        return

    print(
        f"Compatibility mode: loading legacy ungrouped ortholog table into memory: {path}",
        flush=True,
    )
    orthologs_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        orthologs_by_gene[row["query_gene_id"]].append(row)
    yield from orthologs_by_gene.items()


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


def prepare_gene_task(
    tasks_dir: Path,
    gene_id: str,
    gene: dict[str, str],
    source_orthologs: list[dict[str, str]],
    taxonomy: dict[str, dict[str, str]],
    target_path: Path | None,
    ortholog_path: Path | None,
    target_features_path: Path | None,
) -> dict[str, object]:
    ortholog_meta_by_id = {row["ortholog_gene_id"]: row for row in source_orthologs}
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
        return task_row
    if ortholog_path is None:
        task_row.update({"status": "missing_ortholog_fasta", "message": "No ortholog FASTA for gene"})
        return task_row
    if not ortholog_meta_by_id:
        task_row.update({"status": "no_ortholog_metadata", "message": "No ortholog metadata rows for gene"})
        return task_row
    if target_features_path is None:
        task_row.update({"status": "missing_target_features", "message": "No target features for gene"})
        return task_row

    task_dir = tasks_dir / f"task_{gene_id}"
    task_dir.mkdir(parents=True, exist_ok=True)
    target_features_path.replace(task_dir / "target_features.tsv.gz")
    target_id = f"target_{gene_id}"
    ortholog_meta_rows = [
        ortholog_metadata_row(gene_id, source_meta, taxonomy)
        for source_meta in source_orthologs
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
        "target": target_metadata(gene_id, gene, target_id),
        "ortholog_count": len(ortholog_meta_rows),
        "ortholog_metadata": "orthologs.metadata.tsv",
        "target_features": "target_features.tsv.gz",
    }
    (task_dir / "task.json").write_text(json.dumps(manifest, indent=2) + "\n")
    task_row.update(
        {
            "target_fasta": str(Path("sequences") / "targets" / target_path.name),
            "ortholog_fasta": str(Path("sequences") / "orthologs" / ortholog_path.name),
            "ortholog_count": len(ortholog_meta_rows),
        }
    )
    return task_row


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    tasks_dir = args.outdir / "tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)

    genes = {row["gene_id"]: row for row in read_tsv_gz(args.genes_tsv)}
    fetch_manifest = json.loads(args.fetch_manifest.read_text())
    grouped_orthologs = fetch_manifest.get("orthologs_selected_grouped_by_query_gene_id") is True
    taxonomy = load_taxonomy(args.taxonomy_presets)
    target_fastas = fasta_paths_by_gene(args.sequences_dir / "targets")
    ortholog_fastas = fasta_paths_by_gene(args.sequences_dir / "orthologs")
    target_features = partition_target_features(args.target_features, args.outdir / "target_feature_parts")

    task_rows: list[dict[str, object]] = []
    processed_gene_ids: set[str] = set()
    for gene_id, source_orthologs in iter_ortholog_groups(args.orthologs_tsv, grouped_orthologs):
        if gene_id not in genes:
            continue
        task_rows.append(
            prepare_gene_task(
                tasks_dir,
                gene_id,
                genes[gene_id],
                source_orthologs,
                taxonomy,
                target_fastas.get(gene_id),
                ortholog_fastas.get(gene_id),
                target_features.get(gene_id),
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
                taxonomy,
                target_fastas.get(gene_id),
                ortholog_fastas.get(gene_id),
                target_features.get(gene_id),
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
