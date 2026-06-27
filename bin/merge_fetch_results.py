#!/usr/bin/env python3
"""Merge normalized fetch-stage chunk outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids-tsv", required=True, type=Path)
    parser.add_argument("--chunks-tsv", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--target-assembly-accession", required=True)
    parser.add_argument("--target-assembly-name", required=True)
    parser.add_argument("--target-tax-id", required=True)
    parser.add_argument("--chunk-dir", action="append", required=True, type=Path)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def count_tsv_gz_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with gzip.open(path, "rt") as handle:
        next(handle, None)
        return sum(1 for _ in handle)


def merge_tsv_gz(inputs: list[Path], output: Path) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    wrote_header = False
    count = 0
    with gzip.open(output, "wt", newline="") as out:
        writer = None
        for path in inputs:
            if not path.exists():
                continue
            with gzip.open(path, "rt", newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader, None)
                if header is None:
                    continue
                if not wrote_header:
                    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
                    writer.writerow(header)
                    wrote_header = True
                elif writer is None:
                    raise RuntimeError("Internal error: writer not initialized")
                for row in reader:
                    writer.writerow(row)
                    count += 1
    if not wrote_header:
        with gzip.open(output, "wt", newline="") as out:
            out.write("")
    return count


def copy_file_once(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        raise FileExistsError(f"Refusing to overwrite duplicate output: {dst}")
    shutil.copy2(src, dst)


def copy_or_keep(src: Path, dst: Path) -> None:
    if src.resolve() == dst.resolve():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def copy_sequences(chunk_dirs: list[Path], outdir: Path) -> tuple[int, int]:
    target_count = 0
    ortholog_count = 0
    for chunk_dir in chunk_dirs:
        targets_dir = chunk_dir / "sequences" / "targets"
        if targets_dir.exists():
            for src in sorted(targets_dir.glob("*.fa.gz")):
                copy_file_once(src, outdir / "sequences" / "targets" / src.name)
                target_count += 1

        orthologs_dir = chunk_dir / "sequences" / "orthologs"
        if orthologs_dir.exists():
            for src in sorted(orthologs_dir.glob("*.fa.gz")):
                copy_file_once(src, outdir / "sequences" / "orthologs" / src.name)
                ortholog_count += 1
    return target_count, ortholog_count


def read_input_counts(ids_tsv: Path) -> tuple[int, int]:
    total = 0
    accepted = 0
    with ids_tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total += 1
            if row.get("accepted") == "true":
                accepted += 1
    return total, accepted


def load_chunk_manifests(chunk_dirs: list[Path]) -> list[dict]:
    manifests = []
    for chunk_dir in chunk_dirs:
        path = chunk_dir / "manifest.json"
        if path.exists():
            manifests.append(json.loads(path.read_text()))
    return manifests


def main() -> None:
    args = parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    copy_or_keep(args.ids_tsv, outdir / "input.ids.tsv")
    copy_or_keep(args.chunks_tsv, outdir / "chunks.tsv")

    table_inputs = {
        "genes.tsv.gz": [chunk / "genes.tsv.gz" for chunk in args.chunk_dir],
        "orthologs.selected.tsv.gz": [chunk / "orthologs.selected.tsv.gz" for chunk in args.chunk_dir],
        "orthologs.candidates.tsv.gz": [chunk / "orthologs.candidates.tsv.gz" for chunk in args.chunk_dir],
        "failures.tsv.gz": [chunk / "failures.tsv.gz" for chunk in args.chunk_dir],
    }
    table_counts = {
        name: merge_tsv_gz(paths, outdir / name) for name, paths in table_inputs.items()
    }
    target_files, ortholog_files = copy_sequences(args.chunk_dir, outdir)

    input_total, input_unique = read_input_counts(args.ids_tsv)
    chunk_manifests = load_chunk_manifests(args.chunk_dir)
    datasets_versions = sorted(
        {manifest.get("datasets_version", "") for manifest in chunk_manifests if manifest.get("datasets_version")}
    )

    manifest = {
        "created_at": utc_now(),
        "stage": "fetch",
        "input_record_count": input_total,
        "unique_gene_count": input_unique,
        "chunk_count": len(args.chunk_dir),
        "target_gene_count": table_counts["genes.tsv.gz"],
        "selected_ortholog_count": table_counts["orthologs.selected.tsv.gz"],
        "candidate_record_count": table_counts["orthologs.candidates.tsv.gz"],
        "failure_count": table_counts["failures.tsv.gz"],
        "target_sequence_files": target_files,
        "ortholog_sequence_files": ortholog_files,
        "target_assembly_accession": args.target_assembly_accession,
        "target_assembly_name": args.target_assembly_name,
        "target_tax_id": args.target_tax_id,
        "ortholog_scope": "all",
        "datasets_versions": datasets_versions,
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
