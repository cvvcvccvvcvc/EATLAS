#!/usr/bin/env python3
"""Group per-gene alignment tasks by Ensembl Compara MAF source chunk."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from bin.ensembl_compara_maf import (
    STRATEGY_NAME,
    refseq_to_ensembl_seq_region,
    select_candidate_chunks,
)


CHUNK_TASK_FIELDS = [
    "chunk_id",
    "seq_region",
    "chunk_order",
    "source",
    "gene_count",
]

GENE_FIELDS = [
    "gene_id",
    "human_src",
    "genomic_accession",
    "target_origin1",
    "target_end1",
    "target_length",
]

FAILURE_GENE_FIELDS = [
    "gene_id",
    "failure_type",
    "message",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes-tsv", required=True, type=Path)
    parser.add_argument("--maf-manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--strategy", default=STRATEGY_NAME)
    parser.add_argument("--candidate-neighbors", type=int, default=1)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return len(rows)


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return len(rows)


def chunk_id_for(row: dict[str, str]) -> str:
    seq_region = row.get("seq_region") or "unknown"
    chunk_order = row.get("chunk_order") or ""
    if chunk_order:
        return f"chunk_{seq_region}_{int(chunk_order):06d}"
    digest = hashlib.sha1((row.get("source") or "").encode()).hexdigest()[:12]
    return f"chunk_{seq_region}_{digest}"


def target_bounds(target_meta: dict[str, str]) -> tuple[int, int]:
    values = [int(target_meta["begin"]), int(target_meta["end"])]
    return min(values), max(values)


def build_gene_row(target_meta: dict[str, str]) -> tuple[dict[str, object] | None, dict[str, object] | None]:
    gene_id = str(target_meta["gene_id"])
    genomic_accession = target_meta.get("genomic_accession", "")
    seq_region = refseq_to_ensembl_seq_region(genomic_accession)
    if not seq_region:
        return None, {
            "gene_id": gene_id,
            "failure_type": "unknown_ensembl_seq_region",
            "message": f"Could not map genomic_accession={genomic_accession!r} to an Ensembl seq_region",
        }

    target_origin1, target_end1 = target_bounds(target_meta)
    target_length = int(target_meta.get("sequence_length") or target_end1 - target_origin1 + 1)
    return {
        "gene_id": gene_id,
        "human_src": f"homo_sapiens.{seq_region}",
        "seq_region": seq_region,
        "genomic_accession": genomic_accession,
        "target_origin1": target_origin1,
        "target_end1": target_end1,
        "target_length": target_length,
    }, None


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    chunk_root = args.outdir / "maf_chunk_tasks"
    chunk_root.mkdir(parents=True, exist_ok=True)

    manifest_rows = read_tsv_gz(args.maf_manifest)
    grouped: dict[str, dict[str, object]] = {}
    failures: list[dict[str, object]] = []
    gene_count = 0

    for target_meta in read_tsv_gz(args.genes_tsv):
        gene_row, failure = build_gene_row(target_meta)
        if failure:
            failures.append(failure)
            continue
        assert gene_row is not None
        gene_count += 1
        candidates = select_candidate_chunks(
            manifest_rows,
            str(gene_row["seq_region"]),
            int(gene_row["target_origin1"]),
            int(gene_row["target_end1"]),
            args.candidate_neighbors,
        )
        if not candidates:
            failures.append(
                {
                    "gene_id": gene_row["gene_id"],
                    "failure_type": "no_candidate_maf_chunks",
                    "message": (
                        f"No MAF chunks in manifest overlap {gene_row['human_src']}:"
                        f"{gene_row['target_origin1']}-{gene_row['target_end1']}"
                    ),
                }
            )
            continue
        for candidate in candidates:
            chunk_id = chunk_id_for(candidate)
            group = grouped.setdefault(
                chunk_id,
                {
                    "chunk_id": chunk_id,
                    "seq_region": candidate.get("seq_region", ""),
                    "chunk_order": candidate.get("chunk_order", ""),
                    "source": candidate["source"],
                    "genes": [],
                },
            )
            group["genes"].append(gene_row)

    chunk_rows: list[dict[str, object]] = []
    for chunk_id, group in sorted(grouped.items()):
        chunk_dir = chunk_root / chunk_id
        genes = sorted(group["genes"], key=lambda row: (int(row["target_origin1"]), str(row["gene_id"])))
        write_tsv(chunk_dir / "genes.tsv", GENE_FIELDS, genes)
        chunk_manifest = {
            "task_type": "maf_chunk",
            "chunk_id": chunk_id,
            "seq_region": group["seq_region"],
            "chunk_order": group["chunk_order"],
            "source": group["source"],
            "strategy": args.strategy,
            "gene_count": len(genes),
            "genes_tsv": "genes.tsv",
        }
        (chunk_dir / "chunk.json").write_text(json.dumps(chunk_manifest, indent=2, sort_keys=True) + "\n")
        chunk_rows.append(
            {
                "chunk_id": chunk_id,
                "seq_region": group["seq_region"],
                "chunk_order": group["chunk_order"],
                "source": group["source"],
                "gene_count": len(genes),
            }
        )

    if failures:
        failure_dir = chunk_root / "chunk_failures"
        write_tsv(failure_dir / "genes.tsv", FAILURE_GENE_FIELDS, failures)
        failure_manifest = {
            "task_type": "failures",
            "chunk_id": "chunk_failures",
            "strategy": args.strategy,
            "gene_count": len(failures),
            "genes_tsv": "genes.tsv",
        }
        (failure_dir / "chunk.json").write_text(json.dumps(failure_manifest, indent=2, sort_keys=True) + "\n")
        chunk_rows.append(
            {
                "chunk_id": "chunk_failures",
                "seq_region": "",
                "chunk_order": "",
                "source": "",
                "gene_count": len(failures),
            }
        )

    write_tsv_gz(args.outdir / "maf_chunk_tasks.tsv.gz", CHUNK_TASK_FIELDS, chunk_rows)
    run_manifest = {
        "created_at": utc_now(),
        "strategy": args.strategy,
        "candidate_neighbors": args.candidate_neighbors,
        "input_gene_count": gene_count,
        "chunk_task_count": len(chunk_rows),
        "maf_chunk_task_count": len(grouped),
        "failed_gene_count": len(failures),
    }
    (args.outdir / "maf_chunk_tasks_manifest.json").write_text(
        json.dumps(run_manifest, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
