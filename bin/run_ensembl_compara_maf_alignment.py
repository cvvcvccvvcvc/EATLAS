#!/usr/bin/env python3
"""Normalize precomputed Ensembl Compara MAF alignments as a GAPH Stage 2 strategy."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO

from alignment_task_io import load_task_context
from feature_coverage import summarize_feature_coverage


SEGMENT_FIELDS = [
    "gene_id",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "strategy",
    "tool",
    "preset",
    "sequence_id",
    "target_id",
    "query_id",
    "target_start0",
    "target_end0",
    "query_start0",
    "query_end0",
    "strand",
    "matches",
    "block_length",
    "identity",
    "mapq",
    "is_primary",
    "divergence",
    "gap_compressed_divergence",
    "native_record_id",
    "qc_flags",
]

EVENT_FIELDS = [
    "gene_id",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "strategy",
    "tool",
    "preset",
    "event_id",
    "event_type",
    "target_start0",
    "target_end0",
    "genomic_accession",
    "genomic_start1",
    "genomic_end1",
    "ref",
    "alt",
    "query_id",
    "strand",
    "native_record_id",
    "qc_flags",
]

SUMMARY_FIELDS = [
    "gene_id",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "strategy",
    "tool",
    "preset",
    "status",
    "target_length",
    "query_length",
    "segment_count",
    "primary_segment_count",
    "secondary_segment_count",
    "aligned_target_bp",
    "aligned_query_bp",
    "target_coverage",
    "query_coverage",
    "best_identity",
    "mean_identity",
    "event_count",
    "qc_flags",
]

FAILURE_FIELDS = ["gene_id", "ortholog_gene_id", "strategy", "tool", "failure_type", "message"]

DNA_BASES = {"A", "C", "G", "T"}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
TOOL_NAME = "ensembl_compara_maf"
OUTPUT_GZIP_COMPRESSLEVEL = 3


@dataclass(frozen=True)
class MafSequence:
    src: str
    start0: int
    size: int
    strand: str
    src_size: int
    text: str

    def species_and_region(self) -> tuple[str, str]:
        if "." not in self.src:
            return self.src, ""
        species, seq_region = self.src.split(".", 1)
        return species, seq_region

    def forward_interval0(self) -> tuple[int, int]:
        if self.strand == "+":
            return self.start0, self.start0 + self.size
        if self.strand == "-":
            return self.src_size - (self.start0 + self.size), self.src_size - self.start0
        raise ValueError(f"Unsupported MAF strand for {self.src}: {self.strand}")

    def rest_strand(self) -> int:
        return 1 if self.strand == "+" else -1


@dataclass(frozen=True)
class AlignmentRow:
    species: str
    seq_region: str
    start1: int
    end1: int
    strand: int
    seq: str
    description: str

    def query_id(self) -> str:
        return f"{self.species}:{self.seq_region}:{self.start1}:{self.end1}:{self.strand}"


class TsvGzWriter:
    def __init__(self, path: Path, fields: list[str]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = gzip.open(path, "wt", newline="", compresslevel=OUTPUT_GZIP_COMPRESSLEVEL)
        self.writer = csv.DictWriter(self.handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        self.writer.writeheader()
        self.fields = fields
        self.count = 0

    def write(self, row: dict[str, object]) -> None:
        self.writer.writerow({field: row.get(field, "") for field in self.fields})
        self.count += 1

    def close(self) -> None:
        self.handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-dir", required=True, type=Path)
    parser.add_argument("--maf-manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--strategy", default="precomputed_ensembl_92_mammals_epo_extended")
    parser.add_argument("--release", default="116")
    parser.add_argument("--species-set", default="92_mammals.epo_extended")
    parser.add_argument("--method", default="EPO_EXTENDED")
    parser.add_argument("--target-features", type=Path)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--retry-base-seconds", type=float, default=2.0)
    parser.add_argument("--retry-max-seconds", type=float, default=30.0)
    parser.add_argument("--candidate-neighbors", type=int, default=1)
    return parser.parse_args()


def read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv_gz(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="", compresslevel=OUTPUT_GZIP_COMPRESSLEVEL) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
            count += 1
    return count


def interval_union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def refseq_to_ensembl_seq_region(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if text.startswith("chr"):
        text = text[3:]
    if text in {"X", "Y", "MT"} or text.isdigit():
        return str(int(text)) if text.isdigit() else text
    if text == "M":
        return "MT"
    if text.startswith("NC_"):
        base = text.split(".", 1)[0]
        try:
            number = int(base.split("_", 1)[1])
        except (IndexError, ValueError):
            return ""
        if 1 <= number <= 22:
            return str(number)
        if number == 23:
            return "X"
        if number == 24:
            return "Y"
        if number == 12920:
            return "MT"
    return ""


def truthy(value: str) -> bool:
    return str(value).lower() in {"1", "true", "yes", "y"}


def is_remote_source(source: str) -> bool:
    return source.startswith(("http://", "https://"))


def maf_source_name(source: str) -> str:
    if is_remote_source(source):
        return Path(urllib.parse.urlparse(source).path).name
    return Path(source).name


def open_maf_text(source: str, timeout: float) -> TextIO:
    if is_remote_source(source):
        request = urllib.request.Request(source, headers={"User-Agent": "gaph-ensembl-compara-maf/0.1"})
        response = urllib.request.urlopen(request, timeout=timeout)
        return gzip.open(response, "rt")
    return gzip.open(source, "rt")


def retryable_maf_error(exc: Exception) -> bool:
    return isinstance(
        exc,
        (
            EOFError,
            gzip.BadGzipFile,
            urllib.error.URLError,
            TimeoutError,
            OSError,
        ),
    )


def retry_sleep_seconds(args: argparse.Namespace, attempt: int) -> float:
    base = max(float(args.retry_base_seconds), 0.0)
    cap = max(float(args.retry_max_seconds), 0.0)
    delay = min(cap, base * (2 ** max(attempt - 1, 0))) if cap else base
    jitter = random.uniform(0.0, min(1.5, delay * 0.25)) if delay > 0 else 0.0
    return delay + jitter


def source_read_failure(
    args: argparse.Namespace,
    gene_id: str,
    source: str,
    attempts: int,
    completed_block_count: int,
    used_block_count: int,
    exc: Exception | None,
) -> dict[str, object]:
    error_text = f"{type(exc).__name__}: {exc}" if exc else "unknown error"
    return {
        "gene_id": gene_id,
        "ortholog_gene_id": "",
        "strategy": args.strategy,
        "tool": TOOL_NAME,
        "failure_type": "maf_source_read_failed",
        "message": (
            f"{source} failed after {attempts} attempts; "
            f"committed_blocks={completed_block_count}; "
            f"used_blocks={used_block_count}; "
            f"last_error={error_text}"
        ),
    }


def iter_maf_blocks(handle: TextIO) -> Iterable[list[MafSequence]]:
    block: list[MafSequence] = []
    for line in handle:
        line = line.rstrip("\n")
        if not line:
            if block:
                yield block
                block = []
            continue
        if line.startswith("s "):
            fields = line.split()
            if len(fields) < 7:
                continue
            block.append(
                MafSequence(
                    src=fields[1],
                    start0=int(fields[2]),
                    size=int(fields[3]),
                    strand=fields[4],
                    src_size=int(fields[5]),
                    text=fields[6],
                )
            )
    if block:
        yield block


def overlaps(start1: int, end1: int, query_start1: int, query_end1: int) -> bool:
    return end1 >= query_start1 and start1 <= query_end1


def reverse_complement_alignment(text: str) -> str:
    return text.translate(COMPLEMENT)[::-1]


def to_alignment_row(row: MafSequence, flip_orientation: bool) -> AlignmentRow:
    species, seq_region = row.species_and_region()
    start0, end0 = row.forward_interval0()
    strand = row.rest_strand()
    seq = row.text
    if flip_orientation:
        strand *= -1
        seq = reverse_complement_alignment(seq)
    return AlignmentRow(
        species=species,
        seq_region=seq_region,
        start1=start0 + 1,
        end1=end0,
        strand=strand,
        seq=seq,
        description=row.src,
    )


def is_ancestral(row: AlignmentRow) -> bool:
    return (
        "[" in row.species
        or "]" in row.species
        or "-" in row.species
        or row.species.startswith("ancestral")
        or row.species == "ancestral_sequences"
    )


def row_first_pos(row: AlignmentRow) -> int:
    return row.start1 if row.strand >= 0 else row.end1


def advance_pos(row: AlignmentRow, pos: int) -> int:
    return pos + 1 if row.strand >= 0 else pos - 1


def empty_summary(args: argparse.Namespace, row: AlignmentRow, target_length: int) -> dict[str, object]:
    query_length = sum(1 for char in row.seq if char != "-")
    return {
        "gene_id": "",
        "ortholog_gene_id": row.species,
        "tax_id": "",
        "taxname": row.species,
        "strategy": args.strategy,
        "tool": TOOL_NAME,
        "preset": f"{args.method}:{args.species_set}",
        "status": "not_run",
        "target_length": target_length,
        "query_length": query_length,
        "segment_count": 0,
        "primary_segment_count": 0,
        "secondary_segment_count": 0,
        "target_intervals": [],
        "query_intervals": [],
        "identities": [],
        "best_identity": 0.0,
        "event_count": 0,
        "qc_flags": set(),
    }


def finalize_summary(row: dict[str, object]) -> dict[str, object]:
    target_length = int(row.get("target_length") or 0)
    query_length = int(row.get("query_length") or 0)
    aligned_target = interval_union_length(row.pop("target_intervals"))
    aligned_query = interval_union_length(row.pop("query_intervals"))
    identities = row.pop("identities")
    flags = row.pop("qc_flags")
    if row["status"] == "not_run":
        row["status"] = "no_alignment"
        flags.add("no_alignment")
    mean_identity = sum(identities) / len(identities) if identities else 0.0
    row.update(
        {
            "aligned_target_bp": aligned_target,
            "aligned_query_bp": aligned_query,
            "target_coverage": f"{aligned_target / target_length if target_length else 0.0:.6f}",
            "query_coverage": f"{aligned_query / query_length if query_length else 0.0:.6f}",
            "best_identity": f"{float(row['best_identity']):.6f}",
            "mean_identity": f"{mean_identity:.6f}",
            "qc_flags": ",".join(sorted(flags)),
        }
    )
    return row


def append_segment(
    writer: TsvGzWriter,
    summary: dict[str, object],
    args: argparse.Namespace,
    gene_id: str,
    human_row: AlignmentRow,
    query_row: AlignmentRow,
    native_record_id: str,
    target_start0: int,
    target_end0: int,
    query_positions: list[int],
    matches: int,
    block_length: int,
    qc_flags: set[str],
) -> None:
    if target_end0 <= target_start0 or block_length <= 0:
        return
    query_start0 = min(query_positions) - 1 if query_positions else ""
    query_end0 = max(query_positions) if query_positions else ""
    identity = matches / block_length if block_length else 0.0
    query_id = query_row.query_id()
    writer.write(
        {
            "gene_id": gene_id,
            "ortholog_gene_id": query_row.species,
            "tax_id": "",
            "taxname": query_row.species,
            "strategy": args.strategy,
            "tool": TOOL_NAME,
            "preset": f"{args.method}:{args.species_set}",
            "sequence_id": query_id,
            "target_id": human_row.seq_region,
            "query_id": query_id,
            "target_start0": target_start0,
            "target_end0": target_end0,
            "query_start0": query_start0,
            "query_end0": query_end0,
            "strand": query_row.strand,
            "matches": matches,
            "block_length": block_length,
            "identity": f"{identity:.6f}",
            "mapq": "",
            "is_primary": "true",
            "divergence": "",
            "gap_compressed_divergence": "",
            "native_record_id": native_record_id,
            "qc_flags": ",".join(sorted(qc_flags)),
        }
    )
    summary["status"] = "aligned"
    summary["segment_count"] += 1
    summary["primary_segment_count"] += 1
    summary["target_intervals"].append((target_start0, target_end0))
    if query_positions:
        summary["query_intervals"].append((min(query_positions) - 1, max(query_positions)))
    summary["identities"].append(identity)
    summary["best_identity"] = max(float(summary["best_identity"]), identity)


def append_event(
    writer: TsvGzWriter,
    event_id: int,
    args: argparse.Namespace,
    gene_id: str,
    genomic_accession: str,
    target_origin1: int,
    query_row: AlignmentRow,
    event_type: str,
    target_start0: int,
    target_end0: int,
    ref: str,
    alt: str,
    native_record_id: str,
    qc_flags: set[str],
) -> None:
    if target_end0 > target_start0:
        genomic_start1 = target_origin1 + target_start0
        genomic_end1 = target_origin1 + target_end0 - 1
    else:
        genomic_start1 = target_origin1 + target_start0
        genomic_end1 = genomic_start1
    writer.write(
        {
            "gene_id": gene_id,
            "ortholog_gene_id": query_row.species,
            "tax_id": "",
            "taxname": query_row.species,
            "strategy": args.strategy,
            "tool": TOOL_NAME,
            "preset": f"{args.method}:{args.species_set}",
            "event_id": event_id,
            "event_type": event_type,
            "target_start0": target_start0,
            "target_end0": target_end0,
            "genomic_accession": genomic_accession,
            "genomic_start1": genomic_start1,
            "genomic_end1": genomic_end1,
            "ref": ref,
            "alt": alt,
            "query_id": query_row.query_id(),
            "strand": query_row.strand,
            "native_record_id": native_record_id,
            "qc_flags": ",".join(sorted(qc_flags)),
        }
    )


def emit_pending_indel(
    pending: dict[str, object] | None,
    event_writer: TsvGzWriter,
    event_id: int,
    args: argparse.Namespace,
    gene_id: str,
    genomic_accession: str,
    target_origin1: int,
    query_row: AlignmentRow,
    native_record_id: str,
    summary: dict[str, object],
) -> int:
    if not pending:
        return event_id
    append_event(
        event_writer,
        event_id,
        args,
        gene_id,
        genomic_accession,
        target_origin1,
        query_row,
        str(pending["event_type"]),
        int(pending["target_start0"]),
        int(pending["target_end0"]),
        str(pending.get("ref") or ""),
        str(pending.get("alt") or ""),
        native_record_id,
        set(pending.get("qc_flags") or set()),
    )
    summary["status"] = "aligned"
    summary["event_count"] += 1
    return event_id + 1


def convert_pair(
    args: argparse.Namespace,
    gene_id: str,
    genomic_accession: str,
    target_origin1: int,
    target_end1: int,
    human_row: AlignmentRow,
    query_row: AlignmentRow,
    native_record_id: str,
    summary: dict[str, object],
    event_id: int,
    segment_writer: TsvGzWriter,
    event_writer: TsvGzWriter,
) -> int:
    human_seq = human_row.seq.upper()
    query_seq = query_row.seq.upper()
    if len(human_seq) != len(query_seq):
        raise ValueError(f"MSA row length mismatch for {query_row.species}")

    target_length = target_end1 - target_origin1 + 1
    human_pos = row_first_pos(human_row)
    query_pos = row_first_pos(query_row)
    next_target0 = human_pos - target_origin1
    active_segment: dict[str, object] | None = None
    pending_indel: dict[str, object] | None = None

    def close_segment() -> None:
        nonlocal active_segment
        if not active_segment:
            return
        append_segment(
            segment_writer,
            summary,
            args,
            gene_id,
            human_row,
            query_row,
            native_record_id,
            int(active_segment["target_start0"]),
            int(active_segment["target_end0"]),
            list(active_segment["query_positions"]),
            int(active_segment["matches"]),
            int(active_segment["block_length"]),
            set(active_segment["qc_flags"]),
        )
        active_segment = None

    def close_indel() -> None:
        nonlocal pending_indel, event_id
        event_id = emit_pending_indel(
            pending_indel,
            event_writer,
            event_id,
            args,
            gene_id,
            genomic_accession,
            target_origin1,
            query_row,
            native_record_id,
            summary,
        )
        pending_indel = None

    for human_base, query_base in zip(human_seq, query_seq):
        human_has_base = human_base != "-"
        query_has_base = query_base != "-"
        target0 = human_pos - target_origin1 if human_has_base else next_target0
        in_target = 0 <= target0 < target_length if human_has_base else 0 <= target0 <= target_length
        current_query_pos = query_pos if query_has_base else None

        if not in_target:
            close_segment()
            close_indel()
        elif human_has_base and query_has_base:
            close_indel()
            qc_flags: set[str] = set()
            if human_base not in DNA_BASES or query_base not in DNA_BASES:
                qc_flags.add("ambiguous_base")
            if active_segment is None:
                active_segment = {
                    "target_start0": target0,
                    "target_end0": target0,
                    "query_positions": [],
                    "matches": 0,
                    "block_length": 0,
                    "qc_flags": set(),
                }
            active_segment["target_end0"] = target0 + 1
            active_segment["query_positions"].append(current_query_pos)
            active_segment["block_length"] += 1
            active_segment["qc_flags"].update(qc_flags)
            if human_base == query_base:
                active_segment["matches"] += 1
            elif not qc_flags:
                append_event(
                    event_writer,
                    event_id,
                    args,
                    gene_id,
                    genomic_accession,
                    target_origin1,
                    query_row,
                    "snv",
                    target0,
                    target0 + 1,
                    human_base,
                    query_base,
                    native_record_id,
                    set(),
                )
                event_id += 1
                summary["event_count"] += 1

        elif human_has_base and not query_has_base:
            close_segment()
            if pending_indel and pending_indel["event_type"] == "del" and pending_indel["target_end0"] == target0:
                pending_indel["target_end0"] = target0 + 1
                pending_indel["ref"] += human_base
            else:
                close_indel()
                pending_indel = {
                    "event_type": "del",
                    "target_start0": target0,
                    "target_end0": target0 + 1,
                    "ref": human_base,
                    "alt": "",
                    "qc_flags": set(),
                }

        elif not human_has_base and query_has_base:
            close_segment()
            if pending_indel and pending_indel["event_type"] == "ins" and pending_indel["target_start0"] == target0:
                pending_indel["alt"] += query_base
            else:
                close_indel()
                pending_indel = {
                    "event_type": "ins",
                    "target_start0": target0,
                    "target_end0": target0,
                    "ref": "",
                    "alt": query_base,
                    "qc_flags": set(),
                }
        else:
            close_segment()
            close_indel()

        if human_has_base:
            human_pos = advance_pos(human_row, human_pos)
            next_target0 = human_pos - target_origin1
        if query_has_base:
            query_pos = advance_pos(query_row, query_pos)

    close_segment()
    close_indel()
    return event_id


def select_candidate_chunks(
    manifest_rows: list[dict[str, str]],
    seq_region: str,
    start1: int,
    end1: int,
    neighbors: int,
) -> list[dict[str, str]]:
    region_rows = [
        row
        for row in manifest_rows
        if row.get("seq_region") == seq_region and row.get("range_start1") and row.get("range_end1")
    ]
    region_rows.sort(key=lambda row: int(row.get("chunk_order") or 0))
    selected_indices = set()
    for index, row in enumerate(region_rows):
        range_start = int(row["range_start1"])
        range_end = int(row["range_end1"])
        if overlaps(range_start, range_end, start1, end1):
            for neighbor_index in range(max(0, index - neighbors), min(len(region_rows), index + neighbors + 1)):
                selected_indices.add(neighbor_index)
    return [row for index, row in enumerate(region_rows) if index in selected_indices]


def scan_source(
    source: str,
    human_src: str,
    args: argparse.Namespace,
    gene_id: str,
    genomic_accession: str,
    target_origin1: int,
    target_end1: int,
    summaries: dict[str, dict[str, object]],
    segment_writer: TsvGzWriter,
    event_writer: TsvGzWriter,
    event_id: int,
) -> tuple[int, int, int, dict[str, object] | None]:
    source_name = maf_source_name(source)
    completed_block_count = 0
    used_block_count = 0
    row_count = 0
    last_error: Exception | None = None
    attempts = max(int(args.retries), 1)
    for attempt in range(1, attempts + 1):
        current_block_count = 0
        attempt_error: Exception | None = None
        try:
            handle = open_maf_text(source, args.timeout)
        except Exception as exc:
            if not retryable_maf_error(exc):
                raise
            attempt_error = exc
        else:
            with handle:
                block_iter = iter_maf_blocks(handle)
                while True:
                    try:
                        block = next(block_iter)
                    except StopIteration:
                        return event_id, used_block_count, row_count, None
                    except Exception as exc:
                        if not retryable_maf_error(exc):
                            raise
                        attempt_error = exc
                        break

                    current_block_count += 1
                    if current_block_count <= completed_block_count:
                        continue
                    human_rows = [row for row in block if row.src == human_src]
                    if not human_rows:
                        completed_block_count = current_block_count
                        continue
                    human_maf = human_rows[0]
                    human_start0, human_end0 = human_maf.forward_interval0()
                    if not overlaps(human_start0 + 1, human_end0, target_origin1, target_end1):
                        completed_block_count = current_block_count
                        continue
                    flip_orientation = human_maf.strand == "-"
                    human_row = to_alignment_row(human_maf, flip_orientation)
                    block_row_count = 0
                    for query_index, maf_row in enumerate(block, start=1):
                        if maf_row.src == human_src:
                            continue
                        query_row = to_alignment_row(maf_row, flip_orientation)
                        if is_ancestral(query_row):
                            continue
                        summary = summaries.setdefault(
                            query_row.species,
                            empty_summary(args, query_row, target_end1 - target_origin1 + 1),
                        )
                        summary["gene_id"] = gene_id
                        native_record_id = f"{source_name}:block{current_block_count}:row{query_index}"
                        event_id = convert_pair(
                            args,
                            gene_id,
                            genomic_accession,
                            target_origin1,
                            target_end1,
                            human_row,
                            query_row,
                            native_record_id,
                            summary,
                            event_id,
                            segment_writer,
                            event_writer,
                        )
                        block_row_count += 1
                    used_block_count += 1
                    row_count += block_row_count
                    completed_block_count = current_block_count

        if attempt_error is None:
            continue
        last_error = attempt_error
        message = (
            f"MAF source read attempt {attempt}/{attempts} failed for {source_name} "
            f"after {completed_block_count} committed blocks: "
            f"{type(attempt_error).__name__}: {attempt_error}"
        )
        print(message, file=sys.stderr)
        if attempt < attempts:
            time.sleep(retry_sleep_seconds(args, attempt))
            continue
        return (
            event_id,
            used_block_count,
            row_count,
            source_read_failure(args, gene_id, source, attempts, completed_block_count, used_block_count, attempt_error),
        )
    return (
        event_id,
        used_block_count,
        row_count,
        source_read_failure(args, gene_id, source, attempts, completed_block_count, used_block_count, last_error),
    )


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    task, target_meta, _ortholog_meta = load_task_context(args.task_dir)
    gene_id = str(task["gene_id"])
    genomic_accession = target_meta.get("genomic_accession", "")
    seq_region = refseq_to_ensembl_seq_region(genomic_accession)
    start_values = [int(target_meta["genomic_begin"]), int(target_meta["genomic_end"])]
    target_origin1 = min(start_values)
    target_end1 = max(start_values)
    target_length = int(target_meta.get("sequence_length") or target_end1 - target_origin1 + 1)
    human_src = f"homo_sapiens.{seq_region}" if seq_region else ""

    segment_writer = TsvGzWriter(args.outdir / "alignment_segments.tsv.gz", SEGMENT_FIELDS)
    event_writer = TsvGzWriter(args.outdir / "alignment_events.tsv.gz", EVENT_FIELDS)
    failures: list[dict[str, object]] = []
    summaries: dict[str, dict[str, object]] = {}
    source_count = 0
    source_failure_count = 0
    used_block_count = 0
    alignment_row_count = 0
    event_id = 1

    try:
        if not seq_region:
            raise ValueError(f"Could not map genomic_accession={genomic_accession!r} to an Ensembl seq_region")
        manifest_rows = read_tsv_gz(args.maf_manifest)
        candidates = select_candidate_chunks(
            manifest_rows,
            seq_region,
            target_origin1,
            target_end1,
            args.candidate_neighbors,
        )
        if not candidates:
            failures.append(
                {
                    "gene_id": gene_id,
                    "ortholog_gene_id": "",
                    "strategy": args.strategy,
                    "tool": TOOL_NAME,
                    "failure_type": "no_candidate_maf_chunks",
                    "message": f"No MAF chunks in manifest overlap {human_src}:{target_origin1}-{target_end1}",
                }
            )
        for candidate in candidates:
            source = candidate["source"]
            source_count += 1
            event_id, block_delta, row_delta, source_failure = scan_source(
                source,
                human_src,
                args,
                gene_id,
                genomic_accession,
                target_origin1,
                target_end1,
                summaries,
                segment_writer,
                event_writer,
                event_id,
            )
            used_block_count += block_delta
            alignment_row_count += row_delta
            if source_failure:
                failures.append(source_failure)
                source_failure_count += 1
    except Exception as exc:
        failures.append(
            {
                "gene_id": gene_id,
                "ortholog_gene_id": "",
                "strategy": args.strategy,
                "tool": TOOL_NAME,
                "failure_type": "ensembl_compara_maf_failed",
                "message": str(exc),
            }
        )
        raise
    finally:
        segment_writer.close()
        event_writer.close()

    if source_failure_count:
        for summary in summaries.values():
            summary["qc_flags"].add("maf_source_read_failed")

    summary_rows = [finalize_summary(row) for row in summaries.values()]
    write_tsv_gz(args.outdir / "ortholog_alignment_summary.tsv.gz", SUMMARY_FIELDS, summary_rows)
    write_tsv_gz(args.outdir / "failures.tsv.gz", FAILURE_FIELDS, failures)
    feature_coverage_count = None
    if args.target_features:
        feature_coverage_count = summarize_feature_coverage(
            args.target_features,
            args.outdir / "ortholog_alignment_summary.tsv.gz",
            args.outdir / "alignment_segments.tsv.gz",
            args.outdir / "feature_coverage.tsv.gz",
        )
    manifest = {
        "gene_id": gene_id,
        "strategy": args.strategy,
        "tool": TOOL_NAME,
        "release": args.release,
        "species_set": args.species_set,
        "method": args.method,
        "output_gzip_compresslevel": OUTPUT_GZIP_COMPRESSLEVEL,
        "human_src": human_src,
        "genomic_accession": genomic_accession,
        "target_start1": target_origin1,
        "target_end1": target_end1,
        "target_length": target_length,
        "candidate_source_count": source_count,
        "source_failure_count": source_failure_count,
        "used_block_count": used_block_count,
        "alignment_row_count": alignment_row_count,
        "summary_count": len(summary_rows),
        "segment_count": segment_writer.count,
        "event_count": event_writer.count,
        "feature_coverage_count": feature_coverage_count,
        "failure_count": len(failures),
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if failures:
        print(f"{TOOL_NAME} completed with failures for gene_id={gene_id}: {failures}", file=sys.stderr)


if __name__ == "__main__":
    main()
