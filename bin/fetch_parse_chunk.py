#!/usr/bin/env python3
"""Fetch one NCBI Datasets gene package and normalize target/ortholog records."""

from __future__ import annotations

import argparse
import contextlib
import csv
import fcntl
import gzip
import hashlib
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
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
SEQUENCE_GZIP_COMPRESSLEVEL = 3
TARGET_ASSEMBLY_ACCESSION = "GCF_000001405.40"
TARGET_ASSEMBLY_NAME = "GRCh38.p14"
TARGET_TAX_ID = "9606"


@dataclass
class FastaMeta:
    record_index: int
    header: str
    record_start_offset: int
    record_end_offset: int
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


@dataclass
class NormalizedPackage:
    genes: list[dict[str, object]]
    orthologs_selected: list[dict[str, object]]
    orthologs_candidates: list[dict[str, object]]
    failures: list[dict[str, object]]
    catalog: dict
    package_sha256: str


class DownloadError(RuntimeError):
    def __init__(
        self,
        attempts: int,
        returncode: int,
        message: str,
        stagger_wait_seconds: float,
    ) -> None:
        super().__init__(
            f"datasets download failed after {attempts} attempts "
            f"with exit code {returncode}: {message}"
        )
        self.attempts = attempts
        self.returncode = returncode
        self.message = message
        self.stagger_wait_seconds = stagger_wait_seconds


class PackageValidationError(RuntimeError):
    pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids-file", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--datasets-bin", default="datasets")
    parser.add_argument("--request-stagger-seconds", type=float, default=5.0)
    parser.add_argument("--request-throttle-dir", type=Path)
    parser.add_argument("--download-retries", type=int, default=4)
    parser.add_argument("--download-retry-base-seconds", type=float, default=30.0)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@contextlib.contextmanager
def time_step(timings: dict[str, float], name: str):
    start = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        timings[name] = round(timings.get(name, 0.0) + elapsed, 3)


def env_configured(*names: str) -> bool:
    return any(bool(os.environ.get(name)) for name in names)


def throttle_ncbi_request(throttle_dir: Path | None, stagger_seconds: float) -> float:
    if throttle_dir is None or stagger_seconds <= 0:
        return 0.0

    throttle_dir.mkdir(parents=True, exist_ok=True)
    lock_path = throttle_dir / "request_schedule.lock"
    state_path = throttle_dir / "next_request_epoch.txt"

    with lock_path.open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        now = time.time()
        next_allowed = 0.0
        if state_path.exists():
            try:
                next_allowed = float(state_path.read_text().strip())
            except ValueError:
                next_allowed = 0.0
        start_at = max(now, next_allowed)
        state_path.write_text(f"{start_at + stagger_seconds:.6f}\n")
        fcntl.flock(lock, fcntl.LOCK_UN)

    wait_seconds = max(0.0, start_at - time.time())
    if wait_seconds > 0:
        time.sleep(wait_seconds)
    return round(wait_seconds, 3)


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


def run_datasets_download(datasets_bin: str, ids_file: Path, zip_path: Path) -> subprocess.CompletedProcess:
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

    return subprocess.run(cmd, text=True, capture_output=True)


def concise_download_error(result: subprocess.CompletedProcess) -> str:
    text = result.stderr or result.stdout or "datasets download returned no error message"
    return " ".join(text.split())[-500:]


def valid_package_archive(path: Path) -> bool:
    if not path.is_file() or not zipfile.is_zipfile(path):
        return False
    try:
        with zipfile.ZipFile(path) as archive:
            names = set(archive.namelist())
    except (OSError, zipfile.BadZipFile):
        return False
    return {
        "ncbi_dataset/data/data_report.jsonl",
        "ncbi_dataset/data/gene.fna",
    }.issubset(names)


def download_package(
    datasets_bin: str,
    ids_file: Path,
    zip_path: Path,
    retries: int,
    retry_base_seconds: float,
    throttle_dir: Path | None = None,
    stagger_seconds: float = 0.0,
) -> tuple[int, float]:
    attempts = max(1, retries + 1)
    last_result: subprocess.CompletedProcess | None = None
    last_error_message = ""
    stagger_wait_seconds = 0.0
    for attempt in range(1, attempts + 1):
        if zip_path.exists():
            zip_path.unlink()
        stagger_wait_seconds += throttle_ncbi_request(throttle_dir, stagger_seconds)
        result = run_datasets_download(datasets_bin, ids_file, zip_path)
        last_result = result
        if result.returncode == 0 and valid_package_archive(zip_path):
            return attempt, round(stagger_wait_seconds, 3)
        if result.returncode == 0:
            last_error_message = "datasets reported success but did not produce a valid ZIP package"
            print(last_error_message, file=sys.stderr)
        else:
            last_error_message = concise_download_error(result)

        if result.stdout:
            print(result.stdout, file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        if attempt >= attempts:
            break

        wait_seconds = retry_base_seconds * (2 ** (attempt - 1))
        wait_seconds += random.uniform(0, retry_base_seconds * 0.25)
        print(
            f"datasets download failed with exit code {result.returncode}; "
            f"retrying attempt {attempt + 1}/{attempts} after {wait_seconds:.1f}s",
            file=sys.stderr,
            flush=True,
        )
        time.sleep(wait_seconds)

    result = last_result
    if result is None:
        raise RuntimeError("datasets download did not run")
    raise DownloadError(
        attempts=attempts,
        returncode=result.returncode,
        message=last_error_message or concise_download_error(result),
        stagger_wait_seconds=round(stagger_wait_seconds, 3),
    )


def extract_package(zip_path: Path, extract_dir: Path) -> tuple[Path, Path, Path | None]:
    try:
        with zipfile.ZipFile(zip_path) as archive:
            archive.extractall(extract_dir)
    except (OSError, zipfile.BadZipFile) as error:
        raise PackageValidationError(f"Invalid NCBI package: {error}") from error

    data_dir = extract_dir / "ncbi_dataset" / "data"
    report_path = data_dir / "data_report.jsonl"
    gene_fna = data_dir / "gene.fna"
    catalog_path = data_dir / "dataset_catalog.json"
    if not report_path.exists():
        raise PackageValidationError(f"Missing NCBI data report: {report_path}")
    if not gene_fna.exists():
        raise PackageValidationError(f"Missing NCBI gene FASTA: {gene_fna}")
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


def iter_fasta_slices(path: Path):
    header = None
    length = 0
    record_index = 0
    record_start_offset = 0
    with path.open("rb") as handle:
        while True:
            line_start = handle.tell()
            line = handle.readline()
            if not line:
                break
            text = line.decode("utf-8").rstrip("\r\n")
            if text.startswith(">"):
                if header is not None:
                    yield record_index, header, length, record_start_offset, line_start
                record_index += 1
                header = text[1:]
                length = 0
                record_start_offset = line_start
            elif header is not None:
                length += len(text.strip())
        if header is not None:
            yield record_index, header, length, record_start_offset, handle.tell()


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

    for record_index, header, length, record_start_offset, record_end_offset in iter_fasta_slices(gene_fna):
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
                record_start_offset=record_start_offset,
                record_end_offset=record_end_offset,
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


def read_fasta_sequence_slice(path: Path, start_offset: int, end_offset: int) -> str:
    header = None
    seq_parts: list[str] = []
    with path.open("rb") as handle:
        handle.seek(start_offset)
        data = handle.read(end_offset - start_offset).decode("utf-8").splitlines()
        for line in data:
            if line.startswith(">"):
                header = line[1:]
            elif header is not None:
                seq_parts.append(line.strip())
    if header is None:
        raise ValueError(f"No FASTA header found at byte offset {start_offset} in {path}")
    return "".join(seq_parts)


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
        for record_index in sorted(set(target_by_index) | set(ortholog_by_index)):
            record = target_by_index.get(record_index) or ortholog_by_index[record_index]
            seq = read_fasta_sequence_slice(gene_fna, record.record_start_offset, record.record_end_offset)
            if record_index in target_by_index:
                record = target_by_index[record_index]
                output_seq = reverse_complement(seq) if record.is_complement else seq
                target_checksums[record_index] = sha256_text(output_seq.upper())
                handle = target_handles.get(record.gene_id)
                if handle is None:
                    path = outdir / "sequences" / "targets" / f"{record.gene_id}.fa.gz"
                    path.parent.mkdir(parents=True, exist_ok=True)
                    handle = gzip.open(path, "wt", compresslevel=SEQUENCE_GZIP_COMPRESSLEVEL)
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
                    handle = gzip.open(path, "wt", compresslevel=SEQUENCE_GZIP_COMPRESSLEVEL)
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
    selected = sorted(
        (record for record in records if record.selection_role == "ortholog" and record.selected),
        key=lambda record: (
            (0, int(record.query_gene_id))
            if record.query_gene_id.isdigit()
            else (1, record.query_gene_id)
        ),
    )
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


def normalize_package(
    zip_path: Path,
    extract_dir: Path,
    requested_ids: list[str],
    outdir: Path,
    target_assembly_accession: str,
    target_assembly_name: str,
    target_tax_id: str,
    timings: dict[str, float],
) -> NormalizedPackage:
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    with time_step(timings, "extract_package"):
        report_path, gene_fna, catalog_path = extract_package(zip_path, extract_dir)
    with time_step(timings, "load_report"):
        report_by_gene_id, gene_to_query = load_report(report_path, set(requested_ids))
    with time_step(timings, "scan_fasta"):
        fasta_records = build_fasta_metadata(gene_fna, report_by_gene_id, gene_to_query)
    with time_step(timings, "select_records"):
        failures = select_records(
            fasta_records,
            requested_ids,
            report_by_gene_id,
            target_assembly_accession,
            target_tax_id,
        )
    with time_step(timings, "write_sequences"):
        target_checksums, ortholog_checksums = write_sequences(
            gene_fna,
            fasta_records,
            report_by_gene_id,
            outdir,
            target_assembly_accession,
            target_assembly_name,
        )

    catalog = {}
    if catalog_path and catalog_path.exists():
        with time_step(timings, "load_dataset_catalog"):
            catalog = json.loads(catalog_path.read_text())
    with time_step(timings, "package_sha256"):
        package_sha256 = sha256_file(zip_path)

    return NormalizedPackage(
        genes=list(
            genes_rows(
                requested_ids,
                report_by_gene_id,
                fasta_records,
                target_checksums,
                target_assembly_accession,
                target_assembly_name,
            )
        ),
        orthologs_selected=list(selected_ortholog_rows(fasta_records, ortholog_checksums)),
        orthologs_candidates=list(candidate_rows(fasta_records, target_tax_id)),
        failures=failures,
        catalog=catalog,
        package_sha256=package_sha256,
    )


def main() -> None:
    total_start = time.perf_counter()
    timings: dict[str, float] = {}
    args = parse_args()
    with time_step(timings, "read_ids"):
        requested_ids = read_ids(args.ids_file)
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    with time_step(timings, "resolve_datasets_bin"):
        datasets_bin = resolve_datasets_bin(args.datasets_bin)
    with time_step(timings, "datasets_version"):
        version = datasets_version(datasets_bin)

    zip_path = Path("ncbi_dataset.zip")
    extract_dir = Path("ncbi_dataset_unpacked")
    packages: list[NormalizedPackage] = []
    failures: list[dict[str, object]] = []
    download_mode = "batch"
    batch_download_attempts = 0
    singleton_download_attempts = 0
    ncbi_request_wait_seconds = 0.0

    try:
        with time_step(timings, "download_package"):
            batch_download_attempts, wait_seconds = download_package(
                datasets_bin,
                args.ids_file,
                zip_path,
                retries=args.download_retries,
                retry_base_seconds=args.download_retry_base_seconds,
                throttle_dir=args.request_throttle_dir,
                stagger_seconds=args.request_stagger_seconds,
            )
        ncbi_request_wait_seconds += wait_seconds
        packages.append(
            normalize_package(
                zip_path,
                extract_dir,
                requested_ids,
                outdir,
                TARGET_ASSEMBLY_ACCESSION,
                TARGET_ASSEMBLY_NAME,
                TARGET_TAX_ID,
                timings,
            )
        )
    except (DownloadError, PackageValidationError) as batch_error:
        download_mode = "singleton_fallback"
        if isinstance(batch_error, DownloadError):
            batch_download_attempts = batch_error.attempts
            ncbi_request_wait_seconds += batch_error.stagger_wait_seconds
        print(
            f"Batch download failed for {len(requested_ids)} genes; "
            f"falling back to singleton downloads: {batch_error}",
            file=sys.stderr,
            flush=True,
        )
        zip_path.unlink(missing_ok=True)
        shutil.rmtree(extract_dir, ignore_errors=True)

        for gene_id in requested_ids:
            with tempfile.TemporaryDirectory(
                prefix=f"ncbi_gene_{gene_id}_",
                dir=Path.cwd(),
            ) as tmp:
                tmp_dir = Path(tmp)
                singleton_ids = tmp_dir / "gene.ids.txt"
                singleton_zip = tmp_dir / "ncbi_dataset.zip"
                singleton_extract = tmp_dir / "unpacked"
                singleton_ids.write_text(f"{gene_id}\n")
                try:
                    with time_step(timings, "download_package"):
                        attempts, wait_seconds = download_package(
                            datasets_bin,
                            singleton_ids,
                            singleton_zip,
                            retries=args.download_retries,
                            retry_base_seconds=args.download_retry_base_seconds,
                            throttle_dir=args.request_throttle_dir,
                            stagger_seconds=args.request_stagger_seconds,
                        )
                    singleton_download_attempts += attempts
                    ncbi_request_wait_seconds += wait_seconds
                    packages.append(
                        normalize_package(
                            singleton_zip,
                            singleton_extract,
                            [gene_id],
                            outdir,
                            TARGET_ASSEMBLY_ACCESSION,
                            TARGET_ASSEMBLY_NAME,
                            TARGET_TAX_ID,
                            timings,
                        )
                    )
                except DownloadError as error:
                    singleton_download_attempts += error.attempts
                    ncbi_request_wait_seconds += error.stagger_wait_seconds
                    failures.append(
                        failure_row(
                            gene_id,
                            "ncbi_download_failed",
                            (
                                f"attempts={error.attempts}; exit_code={error.returncode}; "
                                f"error={error.message}"
                            ),
                        )
                    )
                except PackageValidationError as error:
                    failures.append(
                        failure_row(
                            gene_id,
                            "ncbi_package_invalid",
                            str(error),
                        )
                    )

    for package in packages:
        failures.extend(package.failures)

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

    with time_step(timings, "write_tables"):
        genes_count = write_tsv_gz(
            outdir / "genes.tsv.gz",
            gene_fields,
            (row for package in packages for row in package.genes),
        )
        selected_count = write_tsv_gz(
            outdir / "orthologs.selected.tsv.gz",
            selected_fields,
            (row for package in packages for row in package.orthologs_selected),
        )
        candidate_count = write_tsv_gz(
            outdir / "orthologs.candidates.tsv.gz",
            candidate_fields,
            (row for package in packages for row in package.orthologs_candidates),
        )
        failure_count = write_tsv_gz(outdir / "failures.tsv.gz", failure_fields, failures)

    def catalog_size(catalog: dict, file_path: str) -> int:
        files = catalog.get("genes", {}).get("files", []) if catalog else []
        for item in files:
            if item.get("filePath") == file_path:
                return int(item.get("uncompressedLengthBytes") or 0)
        return 0

    gene_fna_uncompressed_bytes = sum(
        catalog_size(package.catalog, "gene.fna") for package in packages
    )
    report_uncompressed_bytes = sum(
        catalog_size(package.catalog, "data_report.jsonl") for package in packages
    )
    catalog = packages[0].catalog if download_mode == "batch" and packages else {}
    package_sha256 = packages[0].package_sha256 if download_mode == "batch" and packages else ""
    manifest = {
        "created_at": utc_now(),
        "chunk_id": outdir.name.removeprefix("fetch_"),
        "status": "partial" if failures else "complete",
        "download_mode": download_mode,
        "batch_download_attempts": batch_download_attempts,
        "singleton_download_attempts": singleton_download_attempts,
        "requested_gene_count": len(requested_ids),
        "target_gene_count": genes_count,
        "selected_ortholog_count": selected_count,
        "candidate_record_count": candidate_count,
        "failure_count": failure_count,
        "gene_fna_uncompressed_bytes": gene_fna_uncompressed_bytes,
        "data_report_uncompressed_bytes": report_uncompressed_bytes,
        "target_assembly_accession": TARGET_ASSEMBLY_ACCESSION,
        "target_assembly_name": TARGET_ASSEMBLY_NAME,
        "target_tax_id": TARGET_TAX_ID,
        "ortholog_scope": "all",
        "datasets_bin": datasets_bin,
        "datasets_version": version,
        "ncbi_api_key_configured": env_configured("NCBI_API_KEY", "ENTREZ_API_KEY"),
        "ncbi_contact_email_configured": env_configured("NCBI_EMAIL", "ENTREZ_EMAIL"),
        "request_stagger_seconds": args.request_stagger_seconds,
        "request_stagger_wait_seconds": round(ncbi_request_wait_seconds, 3),
        "download_retries": args.download_retries,
        "download_retry_base_seconds": args.download_retry_base_seconds,
        "sequence_gzip_compresslevel": SEQUENCE_GZIP_COMPRESSLEVEL,
        "package_sha256": package_sha256,
        "dataset_catalog": catalog,
    }
    manifest_path = outdir / "manifest.json"
    manifest["step_timings_seconds"] = timings
    manifest_text = json.dumps(manifest, indent=2) + "\n"
    write_start = time.perf_counter()
    manifest_path.write_text(manifest_text)
    timings["write_manifest"] = round(time.perf_counter() - write_start, 3)
    timings["total_seconds"] = round(time.perf_counter() - total_start, 3)
    manifest["step_timings_seconds"] = timings
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
