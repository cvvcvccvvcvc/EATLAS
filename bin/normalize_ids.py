#!/usr/bin/env python3
"""Normalize Entrez Gene IDs and split them into deterministic chunks."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


TOKEN_RE = re.compile(r"[,\s]+")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids-file", required=True, type=Path)
    parser.add_argument("--chunk-size", required=True, type=int)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def iter_tokens(path: Path):
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            raw_line = line.strip()
            if not raw_line or raw_line.startswith("#"):
                continue
            for raw in TOKEN_RE.split(raw_line):
                token = raw.strip()
                if token:
                    yield line_number, token


def normalize_gene_id(raw: str) -> int:
    if not raw.isdigit():
        raise ValueError(f"Invalid Entrez Gene ID {raw!r}: expected a positive integer")
    gene_id = int(raw)
    if gene_id <= 0:
        raise ValueError(f"Invalid Entrez Gene ID {raw!r}: expected a positive integer")
    return gene_id


def main() -> None:
    args = parse_args()
    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be a positive integer")

    args.outdir.mkdir(parents=True, exist_ok=True)
    chunks_dir = args.outdir / "chunks"
    chunks_dir.mkdir(parents=True, exist_ok=True)

    seen: dict[int, int] = {}
    accepted: list[int] = []
    rows: list[dict[str, object]] = []

    input_position = 0
    for line_number, raw in iter_tokens(args.ids_file):
        input_position += 1
        gene_id = normalize_gene_id(raw)
        if gene_id in seen:
            rows.append(
                {
                    "input_position": input_position,
                    "line_number": line_number,
                    "raw_value": raw,
                    "gene_id": gene_id,
                    "accepted": "false",
                    "accepted_index": "",
                    "duplicate_of_index": seen[gene_id],
                }
            )
            continue

        accepted_index = len(accepted) + 1
        seen[gene_id] = accepted_index
        accepted.append(gene_id)
        rows.append(
            {
                "input_position": input_position,
                "line_number": line_number,
                "raw_value": raw,
                "gene_id": gene_id,
                "accepted": "true",
                "accepted_index": accepted_index,
                "duplicate_of_index": "",
            }
        )

    if not accepted:
        raise ValueError(f"No Entrez Gene IDs found in {args.ids_file}")

    input_tsv = args.outdir / "input.ids.tsv"
    input_fields = [
        "input_position",
        "line_number",
        "raw_value",
        "gene_id",
        "accepted",
        "accepted_index",
        "duplicate_of_index",
    ]
    with input_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=input_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    chunks_tsv = args.outdir / "chunks.tsv"
    chunk_fields = ["chunk_id", "chunk_file", "gene_count", "first_gene_id", "last_gene_id"]
    with chunks_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=chunk_fields, delimiter="\t")
        writer.writeheader()
        for chunk_index, offset in enumerate(range(0, len(accepted), args.chunk_size), start=1):
            chunk = accepted[offset : offset + args.chunk_size]
            chunk_id = f"chunk_{chunk_index:06d}"
            chunk_path = chunks_dir / f"{chunk_id}.ids.txt"
            chunk_path.write_text("".join(f"{gene_id}\n" for gene_id in chunk))
            writer.writerow(
                {
                    "chunk_id": chunk_id,
                    "chunk_file": str(chunk_path.relative_to(args.outdir)),
                    "gene_count": len(chunk),
                    "first_gene_id": chunk[0],
                    "last_gene_id": chunk[-1],
                }
            )


if __name__ == "__main__":
    main()
