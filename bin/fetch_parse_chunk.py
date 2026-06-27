#!/usr/bin/env python3
"""Fetch one NCBI Datasets gene package and normalize target/ortholog records."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


GENE_ID_RE = re.compile(r"\[GeneID=(\d+)\]")
ORGANISM_RE = re.compile(r"\[organism=([^\]]+)\]")
CHROMOSOME_RE = re.compile(r"\[chromosome=([^\]]+)\]")
HEADER_LOC_RE = re.compile(r"^(?P<acc>[^:\s]+):(?P<complement>c?)(?P<a>\d+)-(?P<b>\d+)")

TSV_NULL = ""
FASTA_WIDTH = 80
TARGET_SEQUENCE_ORIENTATION = "plus"


@dataclass
class FastaMeta:
    record_index: int
    header: str
    accession: str
    range_text: str
    begin: int | None
    end: int | None
    is_complement: bool
    gene_id: str
    organism: str
    chromosome: str
    length: int
    query_gene_id: str
    tax_id: str
    taxname: str
    symbol: str
    gene_type: str
    orientation: str
    annotation_index: int | None
    selected: bool = False
    selection_role: str = "candidate"
    reject_reason: str = ""
    priority: tuple = field(default_factory=tuple)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids-file", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--datasets-bin", required=True)
    parser.add_argument("--target-assembly-accession", required=True)
    parser.add_argument("--target-assembly-name", required=True)
    parser.add_argument("--target-tax-id", required=True)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_ids(path: Path) -> list[str]:
    ids = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not ids:
        raise ValueError(f"No IDs in chunk file: {path}")
    return ids


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text_gz(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    return gzip.open(path, "wt", newline="")


def write_tsv_gz(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> int:
    count = 0
    with open_text_gz(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, TSV_NULL) for field in fields})
            count += 1
    return count


def resolve_datasets_bin(raw: str) -> str:
    if raw:
        expanded = Path(raw).expanduser()
        if expanded.is_file():
            return str(expanded)
        found = shutil.which(raw)
        if found:
            return found
    raise FileNotFoundError(f"NCBI Datasets CLI not found: {raw!r}")


def datasets_version(datasets_bin: str) -> str:
    result = subprocess.run(
        [datasets_bin, "--version"],
        check=True,
        text=True,
        capture_output=True,
    )
    return (result.stdout or result.stderr).strip()


def download_package(datasets_bin: str, ids_file: Path, zip_path: Path) -> None:
    cmd = [
        datasets_bin,
        "download",
        "gene",
        "gene-id",
        "--inputfile",
        str(ids_file),
        "--ortholog",
        "all",
        "--include",
        "gene",
        "--filename",
        str(zip_path),
        "--no-progressbar",
    ]
    api_key = os.environ.get("NCBI_API_KEY") or os.environ.get("ENTREZ_API_KEY")
    if api_key:
        cmd.extend(["--api-key", api_key])

    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout, file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"datasets download failed with exit code {result.returncode}")


def extract_package(zip_path: Path, extract_dir: Path) -> tuple[Path, Path, Path | None]:
    with zipfile.ZipFile(zip_path) as archive:
        archive.extractall(extract_dir)

    data_dir = extract_dir / "ncbi_dataset" / "data"
    report_path = data_dir / "data_report.jsonl"
    gene_fna = data_dir / "gene.fna"
    catalog_path = data_dir / "dataset_catalog.json"
    if not report_path.exists():
        raise FileNotFoundError(f"Missing NCBI data report: {report_path}")
    if not gene_fna.exists():
        raise FileNotFoundError(f"Missing NCBI gene FASTA: {gene_fna}")
    return report_path, gene_fna, catalog_path if catalog_path.exists() else None


def load_report(report_path: Path, requested_ids: set[str]) -> tuple[dict[str, dict], dict[str, str]]:
    report_by_gene_id: dict[str, dict] = {}
    gene_to_query: dict[str, str] = {}

    with report_path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            gene_id = str(row.get("geneId") or "").strip()
            if not gene_id:
                continue

            report_by_gene_id[gene_id] = row
            if gene_id in requested_ids:
                gene_to_query[gene_id] = gene_id
                continue

            for group in row.get("geneGroups") or []:
                group_id = str(group.get("id") or "").strip()
                if group_id in requested_ids:
                    gene_to_query[gene_id] = group_id
                    break

    return report_by_gene_id, gene_to_query


def selected_target_location(row: dict, assembly_accession: str) -> tuple[dict, dict] | tuple[None, None]:
    for annotation in row.get("annotations") or []:
        if annotation.get("assemblyAccession") != assembly_accession:
            continue
        locations = annotation.get("genomicLocations") or []
        if locations:
            return annotation, locations[0]
    return None, None


def build_annotation_index(row: dict) -> dict[tuple[str, int, int], int]:
    index: dict[tuple[str, int, int], int] = {}
    for annotation_index, annotation in enumerate(row.get("annotations") or [], start=1):
        for location in annotation.get("genomicLocations") or []:
            genomic_range = location.get("genomicRange") or {}
            accession = location.get("genomicAccessionVersion") or ""
            begin = to_int(genomic_range.get("begin"))
            end = to_int(genomic_range.get("end"))
            if accession and begin is not None and end is not None:
                index[(accession, min(begin, end), max(begin, end))] = annotation_index
    return index


def to_int(value) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def parse_header(header: str) -> dict[str, object]:
    first = header.split()[0]
    loc_match = HEADER_LOC_RE.search(first)
    accession = ""
    range_text = ""
    begin = None
    end = None
    is_complement = False
    if loc_match:
        accession = loc_match.group("acc")
        a = int(loc_match.group("a"))
        b = int(loc_match.group("b"))
        begin = min(a, b)
        end = max(a, b)
        is_complement = bool(loc_match.group("complement"))
        range_text = f"{loc_match.group('complement')}{a}-{b}"

    gene_match = GENE_ID_RE.search(header)
    organism_match = ORGANISM_RE.search(header)
    chromosome_match = CHROMOSOME_RE.search(header)
    return {
        "accession": accession,
        "range_text": range_text,
        "begin": begin,
        "end": end,
        "is_complement": is_complement,
        "gene_id": gene_match.group(1) if gene_match else "",
        "organism": organism_match.group(1) if organism_match else "",
        "chromosome": chromosome_match.group(1) if chromosome_match else "",
    }


def iter_fasta_lengths(path: Path):
    header = None
    length = 0
    record_index = 0
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield record_index, header, length
                record_index += 1
                header = line[1:]
                length = 0
            elif header is not None:
                length += len(line.strip())
        if header is not None:
            yield record_index, header, length


def accession_rank(accession: str) -> int:
    if accession.startswith("NC_"):
        return 0
    if accession.startswith("NW_"):
        return 1
    return 2


def build_fasta_metadata(
    gene_fna: Path,
    report_by_gene_id: dict[str, dict],
    gene_to_query: dict[str, str],
) -> list[FastaMeta]:
    annotation_indexes = {
        gene_id: build_annotation_index(row) for gene_id, row in report_by_gene_id.items()
    }
    records: list[FastaMeta] = []

    for record_index, header, length in iter_fasta_lengths(gene_fna):
        parsed = parse_header(header)
        gene_id = str(parsed["gene_id"])
        report = report_by_gene_id.get(gene_id, {})
        query_gene_id = gene_to_query.get(gene_id, "")
        accession = str(parsed["accession"])
        begin = parsed["begin"]
        end = parsed["end"]
        ann_index = None
        if gene_id in annotation_indexes and begin is not None and end is not None:
            ann_index = annotation_indexes[gene_id].get((accession, begin, end))

        priority = (
            ann_index if ann_index is not None else 10**9,
            accession_rank(accession),
            -length,
            accession,
            str(parsed["range_text"]),
            record_index,
        )
        records.append(
            FastaMeta(
                record_index=record_index,
                header=header,
                accession=accession,
                range_text=str(parsed["range_text"]),
                begin=begin,
                end=end,
                is_complement=bool(parsed["is_complement"]),
                gene_id=gene_id,
                organism=str(parsed["organism"]),
                chromosome=str(parsed["chromosome"]),
                length=length,
                query_gene_id=query_gene_id,
                tax_id=str(report.get("taxId") or ""),
                taxname=str(report.get("taxname") or ""),
                symbol=str(report.get("symbol") or ""),
                gene_type=str(report.get("type") or ""),
                orientation=str(report.get("orientation") or ""),
                annotation_index=ann_index,
                priority=priority,
            )
        )
    return records


def matches_location(record: FastaMeta, location: dict) -> bool:
    genomic_range = location.get("genomicRange") or {}
    begin = to_int(genomic_range.get("begin"))
    end = to_int(genomic_range.get("end"))
    accession = location.get("genomicAccessionVersion") or ""
    return (
        record.accession == accession
        and begin is not None
        and end is not None
        and record.begin == min(begin, end)
        and record.end == max(begin, end)
    )


def select_records(
    records: list[FastaMeta],
    requested_ids: list[str],
    report_by_gene_id: dict[str, dict],
    target_assembly_accession: str,
    target_tax_id: str,
) -> list[dict[str, object]]:
    failures: list[dict[str, object]] = []
    target_selected: set[int] = set()

    records_by_gene_id: defaultdict[str, list[FastaMeta]] = defaultdict(list)
    for record in records:
        if record.gene_id:
            records_by_gene_id[record.gene_id].append(record)

    for query_gene_id in requested_ids:
        row = report_by_gene_id.get(query_gene_id)
        if not row:
            failures.append(
                failure_row(query_gene_id, "target_report_missing", "No NCBI gene report row for requested ID")
            )
            continue
        if str(row.get("taxId") or "") != target_tax_id:
            failures.append(
                failure_row(query_gene_id, "target_not_human", f"Requested gene taxId={row.get('taxId')}")
            )
            continue

        annotation, location = selected_target_location(row, target_assembly_accession)
        if not annotation or not location:
            failures.append(
                failure_row(
                    query_gene_id,
                    "target_assembly_missing",
                    f"No genomic location on {target_assembly_accession}",
                )
            )
            continue

        candidates = [record for record in records_by_gene_id.get(query_gene_id, []) if matches_location(record, location)]
        if not candidates:
            failures.append(
                failure_row(
                    query_gene_id,
                    "target_sequence_missing",
                    f"No gene.fna record matched {location.get('genomicAccessionVersion')}",
                )
            )
            continue

        selected = sorted(candidates, key=lambda record: record.priority)[0]
        selected.selected = True
        selected.selection_role = "target"
        selected.reject_reason = ""
        target_selected.add(selected.record_index)
        for record in candidates:
            if record.record_index != selected.record_index:
                record.reject_reason = "duplicate_target_lower_priority"

    ortholog_groups: defaultdict[tuple[str, str], list[FastaMeta]] = defaultdict(list)
    for record in records:
        if record.record_index in target_selected:
            continue
        if not record.query_gene_id:
            record.reject_reason = "no_query_mapping"
            continue
        if record.tax_id == target_tax_id:
            record.reject_reason = "human_excluded"
            continue
        if not record.gene_id:
            record.reject_reason = "missing_gene_id"
            continue
        ortholog_groups[(record.query_gene_id, record.gene_id)].append(record)

    for _, candidates in ortholog_groups.items():
        selected = sorted(candidates, key=lambda record: record.priority)[0]
        selected.selected = True
        selected.selection_role = "ortholog"
        selected.reject_reason = ""
        for record in candidates:
            if record.record_index != selected.record_index:
                record.reject_reason = "duplicate_gene_lower_priority"

    return failures


def failure_row(gene_id: str, failure_type: str, message: str) -> dict[str, object]:
    return {
        "gene_id": gene_id,
        "failure_type": failure_type,
        "message": message,
    }


def iter_fasta_sequences(path: Path):
    header = None
    seq_parts: list[str] = []
    record_index = 0
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield record_index, header, "".join(seq_parts)
                record_index += 1
                header = line[1:]
                seq_parts = []
            elif header is not None:
                seq_parts.append(line.strip())
        if header is not None:
            yield record_index, header, "".join(seq_parts)


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")
    return seq.translate(table)[::-1]


def wrap_fasta(seq: str) -> str:
    return "\n".join(seq[i : i + FASTA_WIDTH] for i in range(0, len(seq), FASTA_WIDTH)) + "\n"


def fasta_header_target(record: FastaMeta, row: dict, assembly_accession: str, assembly_name: str) -> str:
    return (
        f"gene_{record.gene_id}|symbol={record.symbol}|tax_id={record.tax_id}|"
        f"assembly={assembly_accession}|assembly_name={assembly_name}|"
        f"accession={record.accession}|range={record.begin}-{record.end}|"
        f"orientation={record.orientation}|sequence_orientation={TARGET_SEQUENCE_ORIENTATION}"
    )


def fasta_header_ortholog(record: FastaMeta) -> str:
    safe_taxname = record.taxname.replace(" ", "_")
    return (
        f"query_{record.query_gene_id}|ortholog_gene_{record.gene_id}|symbol={record.symbol}|"
        f"tax_id={record.tax_id}|taxname={safe_taxname}|accession={record.accession}|"
        f"range={record.range_text}|orientation={record.orientation}"
    )


def write_sequences(
    gene_fna: Path,
    records: list[FastaMeta],
    report_by_gene_id: dict[str, dict],
    outdir: Path,
    target_assembly_accession: str,
    target_assembly_name: str,
) -> tuple[dict[int, str], dict[int, str]]:
    target_by_index = {record.record_index: record for record in records if record.selection_role == "target" and record.selected}
    ortholog_by_index = {
        record.record_index: record for record in records if record.selection_role == "ortholog" and record.selected
    }

    target_checksums: dict[int, str] = {}
    ortholog_checksums: dict[int, str] = {}
    target_handles: dict[str, gzip.GzipFile] = {}
    ortholog_handles: dict[str, gzip.GzipFile] = {}

    try:
        for record_index, _header, seq in iter_fasta_sequences(gene_fna):
            if record_index in target_by_index:
                record = target_by_index[record_index]
                output_seq = reverse_complement(seq) if record.is_complement else seq
                target_checksums[record_index] = sha256_text(output_seq.upper())
                handle = target_handles.get(record.gene_id)
                if handle is None:
                    path = outdir / "sequences" / "targets" / f"{record.gene_id}.fa.gz"
                    path.parent.mkdir(parents=True, exist_ok=True)
                    handle = gzip.open(path, "wt")
                    target_handles[record.gene_id] = handle
                row = report_by_gene_id.get(record.gene_id, {})
                handle.write(f">{fasta_header_target(record, row, target_assembly_accession, target_assembly_name)}\n")
                handle.write(wrap_fasta(output_seq))

            if record_index in ortholog_by_index:
                record = ortholog_by_index[record_index]
                ortholog_checksums[record_index] = sha256_text(seq.upper())
                handle = ortholog_handles.get(record.query_gene_id)
                if handle is None:
                    path = outdir / "sequences" / "orthologs" / f"{record.query_gene_id}.fa.gz"
                    path.parent.mkdir(parents=True, exist_ok=True)
                    handle = gzip.open(path, "wt")
                    ortholog_handles[record.query_gene_id] = handle
                handle.write(f">{fasta_header_ortholog(record)}\n")
                handle.write(wrap_fasta(seq))
    finally:
        for handle in target_handles.values():
            handle.close()
        for handle in ortholog_handles.values():
            handle.close()

    return target_checksums, ortholog_checksums


def genes_rows(
    requested_ids: list[str],
    report_by_gene_id: dict[str, dict],
    records: list[FastaMeta],
    target_checksums: dict[int, str],
    target_assembly_accession: str,
    target_assembly_name: str,
) -> Iterable[dict[str, object]]:
    selected_by_gene = {
        record.gene_id: record for record in records if record.selection_role == "target" and record.selected
    }
    for gene_id in requested_ids:
        record = selected_by_gene.get(gene_id)
        if record is None:
            continue
        row = report_by_gene_id.get(gene_id, {})
        annotation, location = selected_target_location(row, target_assembly_accession)
        genomic_range = (location or {}).get("genomicRange") or {}
        yield {
            "gene_id": gene_id,
            "symbol": row.get("symbol") or "",
            "description": row.get("description") or "",
            "tax_id": row.get("taxId") or "",
            "taxname": row.get("taxname") or "",
            "gene_type": row.get("type") or "",
            "assembly_accession": target_assembly_accession,
            "assembly_name": target_assembly_name,
            "annotation_name": (annotation or {}).get("annotationName") or "",
            "annotation_release_date": (annotation or {}).get("annotationReleaseDate") or "",
            "genomic_accession": record.accession,
            "chromosome": record.chromosome,
            "begin": genomic_range.get("begin") or record.begin or "",
            "end": genomic_range.get("end") or record.end or "",
            "orientation": record.orientation,
            "sequence_orientation": TARGET_SEQUENCE_ORIENTATION,
            "sequence_length": record.length,
            "sequence_sha256": target_checksums.get(record.record_index, ""),
            "ensembl_gene_ids": ",".join(row.get("ensemblGeneIds") or row.get("ensembl_gene_ids") or []),
        }


def candidate_rows(records: list[FastaMeta], target_tax_id: str) -> Iterable[dict[str, object]]:
    for record in records:
        if not record.query_gene_id:
            continue
        if record.tax_id == target_tax_id:
            continue
        yield {
            "query_gene_id": record.query_gene_id,
            "ortholog_gene_id": record.gene_id,
            "record_index": record.record_index,
            "selected": "true" if record.selected else "false",
            "selection_role": record.selection_role if record.selected else "candidate",
            "reject_reason": record.reject_reason,
            "tax_id": record.tax_id,
            "taxname": record.taxname,
            "symbol": record.symbol,
            "gene_type": record.gene_type,
            "organism": record.organism,
            "accession": record.accession,
            "accession_kind": record.accession[:2] if record.accession else "",
            "chromosome": record.chromosome,
            "begin": record.begin or "",
            "end": record.end or "",
            "range_text": record.range_text,
            "orientation": record.orientation,
            "source_complement": "true" if record.is_complement else "false",
            "sequence_length": record.length,
            "annotation_index": record.annotation_index or "",
        }


def selected_ortholog_rows(records: list[FastaMeta], checksums: dict[int, str]) -> Iterable[dict[str, object]]:
    selected = [record for record in records if record.selection_role == "ortholog" and record.selected]
    for record in selected:
        yield {
            "query_gene_id": record.query_gene_id,
            "ortholog_gene_id": record.gene_id,
            "tax_id": record.tax_id,
            "taxname": record.taxname,
            "symbol": record.symbol,
            "gene_type": record.gene_type,
            "accession": record.accession,
            "chromosome": record.chromosome,
            "begin": record.begin or "",
            "end": record.end or "",
            "range_text": record.range_text,
            "orientation": record.orientation,
            "source_complement": "true" if record.is_complement else "false",
            "sequence_length": record.length,
            "sequence_sha256": checksums.get(record.record_index, ""),
        }


def main() -> None:
    args = parse_args()
    requested_ids = read_ids(args.ids_file)
    requested_set = set(requested_ids)
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    datasets_bin = resolve_datasets_bin(args.datasets_bin)
    version = datasets_version(datasets_bin)

    zip_path = Path("ncbi_dataset.zip")
    extract_dir = Path("ncbi_dataset_unpacked")
    download_package(datasets_bin, args.ids_file, zip_path)
    report_path, gene_fna, catalog_path = extract_package(zip_path, extract_dir)

    report_by_gene_id, gene_to_query = load_report(report_path, requested_set)
    fasta_records = build_fasta_metadata(gene_fna, report_by_gene_id, gene_to_query)
    failures = select_records(
        fasta_records,
        requested_ids,
        report_by_gene_id,
        args.target_assembly_accession,
        args.target_tax_id,
    )
    target_checksums, ortholog_checksums = write_sequences(
        gene_fna,
        fasta_records,
        report_by_gene_id,
        outdir,
        args.target_assembly_accession,
        args.target_assembly_name,
    )

    gene_fields = [
        "gene_id",
        "symbol",
        "description",
        "tax_id",
        "taxname",
        "gene_type",
        "assembly_accession",
        "assembly_name",
        "annotation_name",
        "annotation_release_date",
        "genomic_accession",
        "chromosome",
        "begin",
        "end",
        "orientation",
        "sequence_orientation",
        "sequence_length",
        "sequence_sha256",
        "ensembl_gene_ids",
    ]
    selected_fields = [
        "query_gene_id",
        "ortholog_gene_id",
        "tax_id",
        "taxname",
        "symbol",
        "gene_type",
        "accession",
        "chromosome",
        "begin",
        "end",
        "range_text",
        "orientation",
        "source_complement",
        "sequence_length",
        "sequence_sha256",
    ]
    candidate_fields = [
        "query_gene_id",
        "ortholog_gene_id",
        "record_index",
        "selected",
        "selection_role",
        "reject_reason",
        "tax_id",
        "taxname",
        "symbol",
        "gene_type",
        "organism",
        "accession",
        "accession_kind",
        "chromosome",
        "begin",
        "end",
        "range_text",
        "orientation",
        "source_complement",
        "sequence_length",
        "annotation_index",
    ]
    failure_fields = ["gene_id", "failure_type", "message"]

    genes_count = write_tsv_gz(
        outdir / "genes.tsv.gz",
        gene_fields,
        genes_rows(
            requested_ids,
            report_by_gene_id,
            fasta_records,
            target_checksums,
            args.target_assembly_accession,
            args.target_assembly_name,
        ),
    )
    selected_count = write_tsv_gz(
        outdir / "orthologs.selected.tsv.gz",
        selected_fields,
        selected_ortholog_rows(fasta_records, ortholog_checksums),
    )
    candidate_count = write_tsv_gz(
        outdir / "orthologs.candidates.tsv.gz",
        candidate_fields,
        candidate_rows(fasta_records, args.target_tax_id),
    )
    failure_count = write_tsv_gz(outdir / "failures.tsv.gz", failure_fields, failures)

    catalog = {}
    if catalog_path and catalog_path.exists():
        catalog = json.loads(catalog_path.read_text())

    manifest = {
        "created_at": utc_now(),
        "requested_gene_count": len(requested_ids),
        "target_gene_count": genes_count,
        "selected_ortholog_count": selected_count,
        "candidate_record_count": candidate_count,
        "failure_count": failure_count,
        "target_assembly_accession": args.target_assembly_accession,
        "target_assembly_name": args.target_assembly_name,
        "ortholog_scope": "all",
        "datasets_bin": datasets_bin,
        "datasets_version": version,
        "package_sha256": sha256_file(zip_path),
        "dataset_catalog": catalog,
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
