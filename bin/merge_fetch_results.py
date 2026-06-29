#!/usr/bin/env python3
"""Merge normalized fetch-stage chunk outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids-tsv", required=True, type=Path)
    parser.add_argument("--chunks-tsv", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--target-assembly-accession", required=True)
    parser.add_argument("--target-assembly-name", required=True)
    parser.add_argument("--target-tax-id", required=True)
    parser.add_argument("--datasets-bin", required=True)
    parser.add_argument("--target-annotation-gff3", type=Path)
    parser.add_argument("--chunk-dir", action="append", required=True, type=Path)
    return parser.parse_args()


TARGET_FEATURE_FIELDS = [
    "gene_id",
    "feature_type",
    "feature_id",
    "genomic_accession",
    "genomic_start1",
    "genomic_end1",
    "target_start0",
    "target_end0",
    "length_bp",
    "strand",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def count_tsv_gz_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with gzip.open(path, "rt") as handle:
        next(handle, None)
        return sum(1 for _ in handle)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
            count += 1
    return count


def resolve_datasets_bin(raw: str) -> str:
    expanded = Path(raw).expanduser()
    if expanded.is_file():
        return str(expanded)
    found = shutil.which(raw)
    if found:
        return found
    raise FileNotFoundError(f"NCBI Datasets CLI not found: {raw!r}")


def download_genome_gff3(datasets_bin: str, assembly_accession: str, workdir: Path) -> tuple[Path, dict[str, object]]:
    zip_path = workdir / "target_annotation.zip"
    extract_dir = workdir / "target_annotation_unpacked"
    cmd = [
        datasets_bin,
        "download",
        "genome",
        "accession",
        assembly_accession,
        "--include",
        "gff3",
        "--filename",
        str(zip_path),
        "--no-progressbar",
    ]
    api_key = os.environ.get("NCBI_API_KEY") or os.environ.get("ENTREZ_API_KEY")
    if api_key:
        cmd.extend(["--api-key", api_key])

    result = None
    for attempt in range(1, 4):
        if zip_path.exists():
            zip_path.unlink()
        result = subprocess.run(cmd, text=True, capture_output=True)
        if result.returncode == 0:
            break
        if attempt < 3:
            message = (result.stderr or result.stdout or "").strip()
            print(
                f"datasets genome GFF3 download attempt {attempt}/3 failed; retrying: {message}",
                file=sys.stderr,
            )
            time.sleep(5 * attempt)

    if result is None or result.returncode != 0:
        if result and result.stdout:
            print(result.stdout, file=sys.stderr)
        if result and result.stderr:
            print(result.stderr, file=sys.stderr)
        exit_code = result.returncode if result else "unknown"
        raise RuntimeError(f"datasets genome GFF3 download failed with exit code {exit_code}")

    with zipfile.ZipFile(zip_path) as archive:
        archive.extractall(extract_dir)

    candidates = sorted(
        path
        for path in extract_dir.rglob("*")
        if path.is_file() and path.name.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz"))
    )
    if not candidates:
        raise FileNotFoundError(f"No GFF3 file found in downloaded annotation package: {zip_path}")

    return candidates[0], {
        "target_annotation_source": "ncbi_datasets_genome_gff3",
        "target_annotation_package_sha256": sha256_file(zip_path),
    }


def parse_gff3_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in raw.split(";"):
        if not part:
            continue
        key, sep, value = part.partition("=")
        if not sep:
            continue
        attrs[unquote(key)] = unquote(value)
    return attrs


def split_attr_values(value: str) -> list[str]:
    return [item for item in value.split(",") if item]


def gene_ids_from_attrs(attrs: dict[str, str]) -> set[str]:
    gene_ids: set[str] = set()
    for key in ("Dbxref", "db_xref"):
        for value in split_attr_values(attrs.get(key, "")):
            if value.startswith("GeneID:"):
                gene_ids.add(value.split(":", 1)[1])
    for key in ("gene_id", "GeneID"):
        value = attrs.get(key, "")
        if value.isdigit():
            gene_ids.add(value)
    return gene_ids


def read_genes(path: Path) -> dict[str, dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def to_int(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Expected integer value, got {value!r}") from error


def strand_from_gene(row: dict[str, str]) -> str:
    raw = (row.get("orientation") or "").lower()
    if raw in {"minus", "-", "reverse"}:
        return "-"
    if raw in {"plus", "+", "forward"}:
        return "+"
    return ""


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def subtract_intervals(
    bases: list[tuple[int, int]],
    masks: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    masks = merge_intervals(masks)
    result: list[tuple[int, int]] = []
    for base_start, base_end in merge_intervals(bases):
        cursor = base_start
        for mask_start, mask_end in masks:
            if mask_end <= cursor:
                continue
            if mask_start >= base_end:
                break
            if mask_start > cursor:
                result.append((cursor, min(mask_start, base_end)))
            cursor = max(cursor, mask_end)
            if cursor >= base_end:
                break
        if cursor < base_end:
            result.append((cursor, base_end))
    return result


def map_gff3_feature_ids(gff3_path: Path, target_gene_ids: set[str]) -> dict[str, str]:
    feature_to_gene: dict[str, str] = {}
    with open_maybe_gzip(gff3_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = parse_gff3_attributes(fields[8])
            direct_gene_ids = gene_ids_from_attrs(attrs) & target_gene_ids
            parent_gene_ids = {
                feature_to_gene[parent_id]
                for parent_id in split_attr_values(attrs.get("Parent", ""))
                if parent_id in feature_to_gene
            }
            resolved_gene_ids = direct_gene_ids or parent_gene_ids
            if len(resolved_gene_ids) != 1:
                continue
            gene_id = next(iter(resolved_gene_ids))
            for feature_id in split_attr_values(attrs.get("ID", "")):
                feature_to_gene[feature_id] = gene_id
    return feature_to_gene


def collect_gff3_intervals(
    gff3_path: Path,
    genes: dict[str, dict[str, str]],
) -> dict[str, dict[str, list[tuple[int, int]]]]:
    target_gene_ids = set(genes)
    feature_to_gene = map_gff3_feature_ids(gff3_path, target_gene_ids)
    intervals: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    wanted_types = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR"}

    with open_maybe_gzip(gff3_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _source, feature_type, start_text, end_text, _score, _strand, _phase, raw_attrs = fields
            if feature_type not in wanted_types:
                continue
            attrs = parse_gff3_attributes(raw_attrs)
            direct_gene_ids = gene_ids_from_attrs(attrs) & target_gene_ids
            parent_gene_ids = {
                feature_to_gene[parent_id]
                for parent_id in split_attr_values(attrs.get("Parent", ""))
                if parent_id in feature_to_gene
            }
            resolved_gene_ids = direct_gene_ids or parent_gene_ids
            if len(resolved_gene_ids) != 1:
                continue

            gene_id = next(iter(resolved_gene_ids))
            gene = genes[gene_id]
            if seqid != gene.get("genomic_accession"):
                continue

            gene_begin = min(to_int(gene["begin"]), to_int(gene["end"]))
            gene_length = to_int(gene["sequence_length"])
            start1 = to_int(start_text)
            end1 = to_int(end_text)
            start0 = max(0, min(start1, end1) - gene_begin)
            end0 = min(gene_length, max(start1, end1) - gene_begin + 1)
            if end0 <= start0:
                continue

            normalized_type = "utr" if feature_type.endswith("_UTR") else feature_type.lower()
            intervals[gene_id][normalized_type].append((start0, end0))

    return intervals


def feature_rows_for_gene(
    gene_id: str,
    gene: dict[str, str],
    raw_intervals: dict[str, list[tuple[int, int]]],
) -> list[dict[str, object]]:
    gene_begin = min(to_int(gene["begin"]), to_int(gene["end"]))
    gene_length = to_int(gene["sequence_length"])
    genomic_accession = gene.get("genomic_accession", "")
    strand = strand_from_gene(gene)

    exons = merge_intervals(raw_intervals.get("exon", []))
    cds = merge_intervals(raw_intervals.get("cds", []))
    utr = merge_intervals(raw_intervals.get("utr", []))
    if not utr and cds:
        utr = subtract_intervals(exons, cds)
    introns = subtract_intervals([(0, gene_length)], exons) if exons else []

    intervals_by_type = {
        "gene": [(0, gene_length)],
        "exon": exons,
        "cds": cds,
        "utr": utr,
        "intron": introns,
    }

    rows: list[dict[str, object]] = []
    for feature_type, intervals in intervals_by_type.items():
        for index, (start0, end0) in enumerate(intervals, start=1):
            rows.append(
                {
                    "gene_id": gene_id,
                    "feature_type": feature_type,
                    "feature_id": f"{feature_type}:{index:03d}",
                    "genomic_accession": genomic_accession,
                    "genomic_start1": gene_begin + start0,
                    "genomic_end1": gene_begin + end0 - 1,
                    "target_start0": start0,
                    "target_end0": end0,
                    "length_bp": end0 - start0,
                    "strand": strand,
                }
            )
    return rows


def build_target_features(genes_tsv: Path, gff3_path: Path, output: Path) -> tuple[int, dict[str, object]]:
    genes = read_genes(genes_tsv)
    if not genes:
        return write_tsv_gz(output, TARGET_FEATURE_FIELDS, []), {
            "target_feature_gene_count": 0,
            "target_genes_without_exon_features": 0,
        }

    raw_intervals = collect_gff3_intervals(gff3_path, genes)
    rows: list[dict[str, object]] = []
    genes_without_exons = 0
    for gene_id in sorted(genes, key=lambda value: int(value) if value.isdigit() else value):
        gene_intervals = raw_intervals.get(gene_id, {})
        if not gene_intervals.get("exon"):
            genes_without_exons += 1
        rows.extend(feature_rows_for_gene(gene_id, genes[gene_id], gene_intervals))

    if rows and all(row["feature_type"] == "gene" for row in rows):
        raise ValueError(f"No exon/CDS/UTR features matched target genes in {gff3_path}")

    count = write_tsv_gz(output, TARGET_FEATURE_FIELDS, rows)
    return count, {
        "target_feature_gene_count": len(genes),
        "target_genes_without_exon_features": genes_without_exons,
    }


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
    if args.target_annotation_gff3:
        gff3_path = args.target_annotation_gff3.expanduser()
        if not gff3_path.exists():
            raise FileNotFoundError(f"Target annotation GFF3 does not exist: {gff3_path}")
        annotation_manifest = {
            "target_annotation_source": "user_gff3",
            "target_annotation_gff3": str(gff3_path),
            "target_annotation_gff3_sha256": sha256_file(gff3_path),
        }
    else:
        datasets_bin = resolve_datasets_bin(args.datasets_bin)
        gff3_path, annotation_manifest = download_genome_gff3(
            datasets_bin,
            args.target_assembly_accession,
            Path("."),
        )
        annotation_manifest["target_annotation_datasets_bin"] = datasets_bin

    target_feature_count, feature_manifest = build_target_features(
        outdir / "genes.tsv.gz",
        gff3_path,
        outdir / "target_features.tsv.gz",
    )

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
        "target_feature_count": target_feature_count,
        "target_assembly_accession": args.target_assembly_accession,
        "target_assembly_name": args.target_assembly_name,
        "target_tax_id": args.target_tax_id,
        "ortholog_scope": "all",
        "datasets_versions": datasets_versions,
        **annotation_manifest,
        **feature_manifest,
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
