#!/usr/bin/env python3
"""Build the normalized fetch-stage dataset from per-chunk outputs."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import shutil
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids-tsv", required=True, type=Path)
    parser.add_argument("--chunks-tsv", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--target-annotation-gff3", required=True, type=Path)
    parser.add_argument("--chunk-dir", action="append", default=[], type=Path)
    parser.add_argument(
        "--chunk-root",
        type=Path,
        help="Directory containing staged per-chunk result directories.",
    )
    return parser.parse_args()


def resolve_chunk_dirs(explicit: list[Path], root: Path | None) -> list[Path]:
    chunk_dirs = list(explicit)
    if root:
        if not root.is_dir():
            raise NotADirectoryError(f"Chunk root is not a directory: {root}")
        chunk_dirs.extend(path for path in root.iterdir() if path.is_dir())
    resolved = sorted(set(chunk_dirs), key=lambda path: path.name)
    if not resolved:
        raise ValueError("Provide at least one --chunk-dir or a non-empty --chunk-root")
    return resolved


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

CHUNK_METRIC_FIELDS = [
    "chunk_id",
    "status",
    "download_mode",
    "batch_download_attempts",
    "singleton_download_attempts",
    "requested_gene_count",
    "target_gene_count",
    "selected_ortholog_count",
    "candidate_record_count",
    "failure_count",
    "gene_fna_uncompressed_bytes",
    "data_report_uncompressed_bytes",
    "ncbi_api_key_configured",
    "ncbi_contact_email_configured",
    "request_stagger_seconds",
    "request_stagger_wait_seconds",
    "timing_total_seconds",
    "timing_download_package_seconds",
    "timing_extract_package_seconds",
    "timing_load_report_seconds",
    "timing_scan_fasta_seconds",
    "timing_select_records_seconds",
    "timing_write_sequences_seconds",
    "timing_write_tables_seconds",
    "timing_package_sha256_seconds",
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


@dataclass(frozen=True)
class TargetInterval:
    gene_id: str
    start1: int
    end1: int


@dataclass(frozen=True)
class TargetIntervalIndex:
    intervals_by_accession: dict[str, list[TargetInterval]]
    starts_by_accession: dict[str, list[int]]

    @classmethod
    def from_genes(cls, genes: dict[str, dict[str, str]]) -> "TargetIntervalIndex":
        intervals_by_accession: dict[str, list[TargetInterval]] = defaultdict(list)
        for gene_id, gene in genes.items():
            accession = gene.get("genomic_accession", "")
            if not accession:
                continue
            begin = to_int(gene["begin"])
            end = to_int(gene["end"])
            start1 = min(begin, end)
            end1 = max(begin, end)
            intervals_by_accession[accession].append(TargetInterval(gene_id, start1, end1))

        starts_by_accession: dict[str, list[int]] = {}
        for accession, intervals in intervals_by_accession.items():
            intervals.sort(key=lambda item: (item.start1, item.end1, item.gene_id))
            starts_by_accession[accession] = [item.start1 for item in intervals]
        return cls(dict(intervals_by_accession), starts_by_accession)

    def overlapping_gene_ids(self, accession: str, start1: int, end1: int) -> set[str]:
        intervals = self.intervals_by_accession.get(accession)
        if not intervals:
            return set()
        starts = self.starts_by_accession[accession]
        limit = bisect.bisect_right(starts, end1)
        return {
            item.gene_id
            for item in intervals[:limit]
            if item.end1 >= start1
        }


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


def collect_gff3_intervals(
    gff3_path: Path,
    genes: dict[str, dict[str, str]],
) -> dict[str, dict[str, list[tuple[int, int]]]]:
    target_gene_ids = set(genes)
    interval_index = TargetIntervalIndex.from_genes(genes)
    feature_to_gene: dict[str, str] = {}
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
            if seqid not in interval_index.intervals_by_accession:
                continue
            start1 = to_int(start_text)
            end1 = to_int(end_text)
            feature_start1 = min(start1, end1)
            feature_end1 = max(start1, end1)
            overlapping_gene_ids = interval_index.overlapping_gene_ids(seqid, feature_start1, feature_end1)
            if not overlapping_gene_ids:
                continue

            attrs = parse_gff3_attributes(raw_attrs)
            direct_gene_ids = gene_ids_from_attrs(attrs) & target_gene_ids & overlapping_gene_ids
            parent_gene_ids = {
                feature_to_gene[parent_id]
                for parent_id in split_attr_values(attrs.get("Parent", ""))
                if parent_id in feature_to_gene and feature_to_gene[parent_id] in overlapping_gene_ids
            } & target_gene_ids
            resolved_gene_ids = direct_gene_ids or parent_gene_ids
            if len(resolved_gene_ids) != 1:
                continue

            gene_id = next(iter(resolved_gene_ids))
            for feature_id in split_attr_values(attrs.get("ID", "")):
                feature_to_gene[feature_id] = gene_id

            if feature_type not in wanted_types:
                continue

            gene = genes[gene_id]

            gene_begin = min(to_int(gene["begin"]), to_int(gene["end"]))
            gene_length = to_int(gene["sequence_length"])
            start0 = max(0, feature_start1 - gene_begin)
            end0 = min(gene_length, feature_end1 - gene_begin + 1)
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
    expected_header: list[str] | None = None
    count = 0
    with gzip.open(output, "wt", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        for path in inputs:
            if not path.exists():
                raise FileNotFoundError(f"Missing chunk table: {path}")
            with gzip.open(path, "rt", newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader, None)
                if not header or any(not field for field in header):
                    raise ValueError(f"Chunk table has no header: {path}")
                if expected_header is None:
                    expected_header = header
                    writer.writerow(header)
                elif header != expected_header:
                    raise ValueError(
                        f"Chunk table header mismatch in {path}: "
                        f"expected={expected_header}, observed={header}"
                    )
                for row_number, row in enumerate(reader, start=2):
                    if len(row) != len(expected_header):
                        raise ValueError(
                            f"Chunk table row width mismatch in {path} at line {row_number}: "
                            f"expected={len(expected_header)}, observed={len(row)}"
                        )
                    writer.writerow(row)
                    count += 1
    return count


def validate_fetch_counts(target_gene_count: int, selected_ortholog_count: int) -> None:
    if target_gene_count == 0:
        raise ValueError("Fetch produced no target genes; the pipeline cannot continue")
    if selected_ortholog_count == 0:
        raise ValueError("Fetch produced no selected orthologs; alignment cannot continue")


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


def read_accepted_gene_ids(ids_tsv: Path) -> set[str]:
    with ids_tsv.open(newline="") as handle:
        return {
            row["gene_id"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("accepted") == "true"
        }


def read_expected_chunk_ids(chunks_tsv: Path) -> set[str]:
    with chunks_tsv.open(newline="") as handle:
        return {
            row["chunk_id"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("chunk_id")
        }


def read_tsv_gz_column(path: Path, column: str) -> list[str]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if column not in (reader.fieldnames or []):
            raise ValueError(f"Missing {column!r} column in {path}")
        return [row[column] for row in reader if row.get(column)]


def load_chunk_manifests(chunk_dirs: list[Path]) -> list[dict]:
    manifests = []
    for chunk_dir in chunk_dirs:
        path = chunk_dir / "manifest.json"
        if not path.exists():
            raise FileNotFoundError(f"Missing chunk manifest: {path}")
        manifests.append(json.loads(path.read_text()))
    return manifests


def validate_chunk_manifests(expected_ids: set[str], manifests: list[dict]) -> None:
    observed = [str(manifest.get("chunk_id") or "") for manifest in manifests]
    if any(not chunk_id for chunk_id in observed):
        raise ValueError("Every chunk manifest must contain a non-empty chunk_id")
    if len(observed) != len(set(observed)):
        raise ValueError("Duplicate chunk_id values found in chunk manifests")
    if set(observed) != expected_ids:
        missing = sorted(expected_ids - set(observed))
        unexpected = sorted(set(observed) - expected_ids)
        raise ValueError(
            f"Chunk manifest mismatch: missing={missing or '-'}, "
            f"unexpected={unexpected or '-'}"
        )


def consistent_manifest_value(manifests: list[dict], field: str) -> str:
    values = {str(manifest.get(field) or "") for manifest in manifests}
    if "" in values or len(values) != 1:
        raise ValueError(
            f"Chunk manifests must agree on a non-empty {field}: "
            f"values={sorted(values)}"
        )
    return next(iter(values))


def validate_gene_outcomes(
    accepted_ids: set[str],
    gene_ids: list[str],
    failure_ids: list[str],
) -> None:
    if len(gene_ids) != len(set(gene_ids)):
        raise ValueError("Duplicate gene_id values found in genes.tsv.gz")
    if len(failure_ids) != len(set(failure_ids)):
        raise ValueError("Duplicate gene_id values found in failures.tsv.gz")

    successes = set(gene_ids)
    failures = set(failure_ids)
    overlap = successes & failures
    unknown = (successes | failures) - accepted_ids
    missing = accepted_ids - successes - failures
    if overlap or unknown or missing:
        raise ValueError(
            "Fetch gene outcome mismatch: "
            f"both_success_and_failure={sorted(overlap) or '-'}, "
            f"unknown={sorted(unknown) or '-'}, "
            f"missing={sorted(missing) or '-'}"
        )


def chunk_metric_rows(chunk_manifests: list[dict]) -> list[dict[str, object]]:
    rows = []
    for manifest in chunk_manifests:
        timings = manifest.get("step_timings_seconds") or {}
        rows.append(
            {
                "chunk_id": manifest.get("chunk_id", ""),
                "status": manifest.get("status", ""),
                "download_mode": manifest.get("download_mode", ""),
                "batch_download_attempts": manifest.get("batch_download_attempts", ""),
                "singleton_download_attempts": manifest.get("singleton_download_attempts", ""),
                "requested_gene_count": manifest.get("requested_gene_count", ""),
                "target_gene_count": manifest.get("target_gene_count", ""),
                "selected_ortholog_count": manifest.get("selected_ortholog_count", ""),
                "candidate_record_count": manifest.get("candidate_record_count", ""),
                "failure_count": manifest.get("failure_count", ""),
                "gene_fna_uncompressed_bytes": manifest.get("gene_fna_uncompressed_bytes", ""),
                "data_report_uncompressed_bytes": manifest.get("data_report_uncompressed_bytes", ""),
                "ncbi_api_key_configured": manifest.get("ncbi_api_key_configured", ""),
                "ncbi_contact_email_configured": manifest.get("ncbi_contact_email_configured", ""),
                "request_stagger_seconds": manifest.get("request_stagger_seconds", ""),
                "request_stagger_wait_seconds": manifest.get("request_stagger_wait_seconds", ""),
                "timing_total_seconds": timings.get("total_seconds", ""),
                "timing_download_package_seconds": timings.get("download_package", ""),
                "timing_extract_package_seconds": timings.get("extract_package", ""),
                "timing_load_report_seconds": timings.get("load_report", ""),
                "timing_scan_fasta_seconds": timings.get("scan_fasta", ""),
                "timing_select_records_seconds": timings.get("select_records", ""),
                "timing_write_sequences_seconds": timings.get("write_sequences", ""),
                "timing_write_tables_seconds": timings.get("write_tables", ""),
                "timing_package_sha256_seconds": timings.get("package_sha256", ""),
            }
        )
    return sorted(rows, key=lambda row: str(row.get("chunk_id", "")))


def main() -> None:
    args = parse_args()
    chunk_dirs = resolve_chunk_dirs(args.chunk_dir, args.chunk_root)
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    copy_or_keep(args.ids_tsv, outdir / "input.ids.tsv")
    copy_or_keep(args.chunks_tsv, outdir / "chunks.tsv")
    accepted_gene_ids = read_accepted_gene_ids(args.ids_tsv)
    expected_chunk_ids = read_expected_chunk_ids(args.chunks_tsv)
    chunk_manifests = load_chunk_manifests(chunk_dirs)
    validate_chunk_manifests(expected_chunk_ids, chunk_manifests)
    target_metadata = {
        field: consistent_manifest_value(chunk_manifests, field)
        for field in (
            "target_assembly_accession",
            "target_assembly_name",
            "target_tax_id",
        )
    }

    table_inputs = {
        "genes.tsv.gz": [chunk / "genes.tsv.gz" for chunk in chunk_dirs],
        "orthologs.selected.tsv.gz": [chunk / "orthologs.selected.tsv.gz" for chunk in chunk_dirs],
        "orthologs.candidates.tsv.gz": [chunk / "orthologs.candidates.tsv.gz" for chunk in chunk_dirs],
        "failures.tsv.gz": [chunk / "failures.tsv.gz" for chunk in chunk_dirs],
    }
    table_counts = {
        name: merge_tsv_gz(paths, outdir / name) for name, paths in table_inputs.items()
    }
    validate_fetch_counts(
        table_counts["genes.tsv.gz"],
        table_counts["orthologs.selected.tsv.gz"],
    )
    target_files, ortholog_files = copy_sequences(chunk_dirs, outdir)
    gene_ids = read_tsv_gz_column(outdir / "genes.tsv.gz", "gene_id")
    failure_ids = read_tsv_gz_column(outdir / "failures.tsv.gz", "gene_id")
    validate_gene_outcomes(accepted_gene_ids, gene_ids, failure_ids)
    if target_files != len(gene_ids):
        raise ValueError(
            f"Target FASTA count mismatch: files={target_files}, genes={len(gene_ids)}"
        )

    with gzip.open(outdir / "failures.tsv.gz", "rt", newline="") as handle:
        failure_rows = list(csv.DictReader(handle, delimiter="\t"))
    download_failed_gene_count = sum(
        row.get("failure_type") in {"ncbi_download_failed", "ncbi_package_invalid"}
        for row in failure_rows
    )
    singleton_fallback_chunk_count = sum(
        manifest.get("download_mode") == "singleton_fallback"
        for manifest in chunk_manifests
    )

    gff3_path = args.target_annotation_gff3.expanduser()
    if not gff3_path.exists():
        raise FileNotFoundError(f"Target annotation GFF3 does not exist: {gff3_path}")
    annotation_manifest = {
        "target_annotation_source": "user_gff3",
        "target_annotation_gff3": str(gff3_path),
        "target_annotation_gff3_sha256": sha256_file(gff3_path),
    }

    target_feature_count, feature_manifest = build_target_features(
        outdir / "genes.tsv.gz",
        gff3_path,
        outdir / "target_features.tsv.gz",
    )

    input_total, input_unique = read_input_counts(args.ids_tsv)
    chunk_metric_count = write_tsv_gz(
        outdir / "chunk_metrics.tsv.gz",
        CHUNK_METRIC_FIELDS,
        chunk_metric_rows(chunk_manifests),
    )
    datasets_versions = sorted(
        {manifest.get("datasets_version", "") for manifest in chunk_manifests if manifest.get("datasets_version")}
    )

    manifest = {
        "created_at": utc_now(),
        "stage": "fetch",
        "status": "partial" if failure_ids else "complete",
        "input_record_count": input_total,
        "unique_gene_count": input_unique,
        "chunk_count": len(chunk_dirs),
        "chunk_metric_count": chunk_metric_count,
        "target_gene_count": table_counts["genes.tsv.gz"],
        "selected_ortholog_count": table_counts["orthologs.selected.tsv.gz"],
        "candidate_record_count": table_counts["orthologs.candidates.tsv.gz"],
        "failure_count": table_counts["failures.tsv.gz"],
        "download_failed_gene_count": download_failed_gene_count,
        "singleton_fallback_chunk_count": singleton_fallback_chunk_count,
        "target_sequence_files": target_files,
        "ortholog_sequence_files": ortholog_files,
        "target_feature_count": target_feature_count,
        **target_metadata,
        "ortholog_scope": "all",
        "datasets_versions": datasets_versions,
        **annotation_manifest,
        **feature_manifest,
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
