#!/usr/bin/env python3
"""Run minimap2 for one gene task and normalize alignment evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import tempfile
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from bin.alignment_table_schema import (
    EVENT_FIELDS,
    FAILURE_FIELDS,
    SEGMENT_FIELDS,
    SUMMARY_FIELDS,
)
from bin.alignment_task_io import (
    iter_fasta,
    load_task_context,
    materialize_task_fastas,
    write_fasta_record,
)
from bin.alignment_runtime import minimap2_software


CS_OP_RE = re.compile(r"(:\d+|=[A-Za-z]+|\*[A-Za-z][A-Za-z]|[+\-][A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2})")
TSV_NULL = ""


@dataclass(frozen=True)
class QuerySlice:
    source_sequence_id: str
    source_start0: int
    source_end0: int
    source_length: int
    read_index: int
    is_pseudoread: bool


@dataclass(frozen=True)
class PseudoreadGeneration:
    total_reads: int
    query_slices: dict[str, QuerySlice]


@dataclass(frozen=True)
class BackboneSelection:
    accepted_record_ids: frozenset[str]
    input_alignment_count: int
    after_strand_count: int
    retained_alignment_count: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-dir", required=True, type=Path)
    parser.add_argument("--source-target-fasta", required=True, type=Path)
    parser.add_argument("--source-ortholog-fasta", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--strategy", required=True)
    parser.add_argument("--preset", choices=["asm10", "asm20", "map-ont"], required=True)
    parser.add_argument("--pseudoread-len", default=0, type=int)
    parser.add_argument("--pseudoread-step", default=0, type=int)
    parser.add_argument("--minimap2-bin", default="minimap2")
    parser.add_argument("--threads", default=1, type=int)
    return parser.parse_args()


def validate_query_mode(preset: str, pseudoread_len: int, pseudoread_step: int) -> None:
    if (pseudoread_len == 0) != (pseudoread_step == 0):
        raise ValueError("--pseudoread-len and --pseudoread-step must be set together")
    if pseudoread_len < 0 or pseudoread_step < 0:
        raise ValueError("Pseudoread length and step must be non-negative")
    if pseudoread_len and pseudoread_step > pseudoread_len:
        raise ValueError("Pseudoread step must not exceed pseudoread length")
    if preset == "map-ont" and not pseudoread_len:
        raise ValueError("The map-ont strategy requires pseudoread generation parameters")
    if pseudoread_len and preset != "map-ont":
        raise ValueError("Pseudoread generation is supported only for the map-ont strategy")


def pseudoread_starts(seq_len: int, read_len: int, step: int) -> list[int]:
    """Return gap-free deterministic starts, including the final full window."""
    if seq_len <= 0:
        raise ValueError("Pseudoread source sequence must be non-empty")
    if read_len <= 0 or step <= 0:
        raise ValueError("Pseudoread length and step must be positive")
    if step > read_len:
        raise ValueError("Pseudoread step must not exceed pseudoread length")
    if seq_len <= read_len:
        return [0]

    final_start = seq_len - read_len
    starts = list(range(0, final_start + 1, step))
    if starts[-1] != final_start:
        starts.append(final_start)
    return starts


def full_query_slices(meta_by_sequence: dict[str, dict[str, str]]) -> dict[str, QuerySlice]:
    return {
        sequence_id: QuerySlice(
            source_sequence_id=sequence_id,
            source_start0=0,
            source_end0=int(meta["sequence_length"]),
            source_length=int(meta["sequence_length"]),
            read_index=1,
            is_pseudoread=False,
        )
        for sequence_id, meta in meta_by_sequence.items()
    }


def generate_long_pseudoreads(
    orthologs_fasta: Path,
    out_fasta: Path,
    meta_by_sequence: dict[str, dict[str, str]],
    read_len: int,
    step: int,
) -> PseudoreadGeneration:
    query_slices: dict[str, QuerySlice] = {}
    seen_sources: set[str] = set()

    with out_fasta.open("w") as out:
        for header, sequence in iter_fasta(orthologs_fasta):
            source_id = header.split()[0]
            meta = meta_by_sequence.get(source_id)
            if meta is None:
                raise ValueError(f"Pseudoread source is missing metadata: {source_id}")
            source_length = int(meta["sequence_length"])
            if len(sequence) != source_length:
                raise ValueError(
                    f"Pseudoread source length mismatch for {source_id}: "
                    f"metadata={source_length}, fasta={len(sequence)}"
                )
            seen_sources.add(source_id)
            for read_index, start0 in enumerate(
                pseudoread_starts(source_length, read_len, step),
                start=1,
            ):
                end0 = min(start0 + read_len, source_length)
                query_name = f"{source_id}_long_{read_index}_{start0}-{end0}"
                write_fasta_record(out, query_name, sequence[start0:end0])
                query_slices[query_name] = QuerySlice(
                    source_sequence_id=source_id,
                    source_start0=start0,
                    source_end0=end0,
                    source_length=source_length,
                    read_index=read_index,
                    is_pseudoread=True,
                )

    missing = sorted(set(meta_by_sequence) - seen_sources)
    if missing:
        raise ValueError(
            "Pseudoread input is missing source sequence(s): " + ", ".join(missing[:10])
        )
    return PseudoreadGeneration(len(query_slices), query_slices)


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


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def is_primary(tags: dict[str, str]) -> bool:
    return tags.get("tp", "P") in {"P", "I"}


def interval_union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def genomic_coords(target_meta: dict[str, str], start0: int, end0: int) -> tuple[str, str]:
    begin_text = target_meta.get("genomic_begin") or ""
    if not begin_text:
        return "", ""
    begin = int(begin_text)
    genomic_start = begin + start0
    if end0 > start0:
        genomic_end = begin + end0 - 1
    else:
        genomic_end = genomic_start
    return str(genomic_start), str(genomic_end)


def paf_record_digest(line: str) -> str:
    return hashlib.sha256(line.encode("utf-8")).hexdigest()[:32]


def stable_paf_record_id(record_digest: str, occurrence: int) -> str:
    return f"paf:{record_digest}:{occurrence}"


def iter_paf_records(path: Path):
    record_occurrences: dict[str, int] = defaultdict(int)
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                preview = line[:200] + ("..." if len(line) > 200 else "")
                raise ValueError(
                    f"Malformed PAF record in {path} at line {line_number}: "
                    f"expected at least 12 tab-separated fields, observed "
                    f"{len(fields)}; record={preview!r}"
                )
            native_record_id = ""
            if fields[5] != "*" and fields[4] != "*":
                record_digest = paf_record_digest(line)
                record_occurrences[record_digest] += 1
                native_record_id = stable_paf_record_id(
                    record_digest,
                    record_occurrences[record_digest],
                )
            yield fields, native_record_id, line_number


def require_query_slice(
    query_slices: dict[str, QuerySlice],
    query_name: str,
    paf_path: Path,
    line_number: int,
) -> QuerySlice:
    query_slice = query_slices.get(query_name)
    if query_slice is None:
        raise ValueError(
            f"PAF record in {paf_path} at line {line_number} references "
            f"unknown query ID {query_name!r}"
        )
    return query_slice


def longest_increasing_indices(values: list[int]) -> set[int]:
    if not values:
        return set()

    tails_values: list[int] = []
    tails_indices: list[int] = []
    previous = [-1] * len(values)
    for index, value in enumerate(values):
        position = bisect_left(tails_values, value)
        if position == len(tails_values):
            tails_values.append(value)
            tails_indices.append(index)
        else:
            tails_values[position] = value
            tails_indices[position] = index
        if position > 0:
            previous[index] = tails_indices[position - 1]

    selected: set[int] = set()
    current = tails_indices[-1]
    while current != -1:
        selected.add(current)
        current = previous[current]
    return selected


def select_pseudoread_backbone(
    paf_path: Path,
    query_slices: dict[str, QuerySlice],
) -> BackboneSelection:
    records_by_source: dict[str, list[dict[str, object]]] = defaultdict(list)
    for fields, native_record_id, line_number in iter_paf_records(paf_path):
        query_slice = require_query_slice(
            query_slices,
            fields[0],
            paf_path,
            line_number,
        )
        if not native_record_id:
            continue
        records_by_source[query_slice.source_sequence_id].append(
            {
                "native_record_id": native_record_id,
                "read_index": query_slice.read_index,
                "target_start": int(fields[7]),
                "is_reverse": fields[4] == "-",
            }
        )

    accepted: set[str] = set()
    input_count = 0
    after_strand_count = 0
    for source_id in sorted(records_by_source):
        records = records_by_source[source_id]
        input_count += len(records)
        forward_count = sum(not bool(record["is_reverse"]) for record in records)
        reverse_count = len(records) - forward_count
        dominant_reverse = reverse_count > forward_count
        strand_records = [
            record
            for record in records
            if bool(record["is_reverse"]) == dominant_reverse
        ]
        strand_records.sort(
            key=lambda record: (
                int(record["target_start"]),
                str(record["native_record_id"]),
            )
        )
        after_strand_count += len(strand_records)
        source_order = [int(record["read_index"]) for record in strand_records]
        if dominant_reverse:
            source_order = [-value for value in source_order]
        backbone = longest_increasing_indices(source_order)
        accepted.update(
            str(record["native_record_id"])
            for index, record in enumerate(strand_records)
            if index in backbone
        )

    return BackboneSelection(
        accepted_record_ids=frozenset(accepted),
        input_alignment_count=input_count,
        after_strand_count=after_strand_count,
        retained_alignment_count=len(accepted),
    )


def cs_events(
    cs: str,
    target_start0: int,
    record: dict[str, object],
    target_meta: dict[str, str],
    event_id_prefix: str,
) -> list[dict[str, object]]:
    events: list[dict[str, object]] = []
    target_pos = target_start0
    event_index = 1

    for match in CS_OP_RE.finditer(cs):
        op = match.group(0)
        if op.startswith(":"):
            target_pos += int(op[1:])
        elif op.startswith("="):
            target_pos += len(op) - 1
        elif op.startswith("*"):
            ref = op[1].upper()
            alt = op[2].upper()
            genomic_start, genomic_end = genomic_coords(target_meta, target_pos, target_pos + 1)
            events.append(
                {
                    **record,
                    "event_id": f"{event_id_prefix}:{event_index}",
                    "event_type": "snv",
                    "target_start0": target_pos,
                    "target_end0": target_pos + 1,
                    "genomic_accession": target_meta.get("genomic_accession", ""),
                    "genomic_start1": genomic_start,
                    "genomic_end1": genomic_end,
                    "ref": ref,
                    "alt": alt,
                }
            )
            event_index += 1
            target_pos += 1
        elif op.startswith("+"):
            alt = op[1:].upper()
            genomic_start, genomic_end = genomic_coords(target_meta, target_pos, target_pos)
            events.append(
                {
                    **record,
                    "event_id": f"{event_id_prefix}:{event_index}",
                    "event_type": "ins",
                    "target_start0": target_pos,
                    "target_end0": target_pos,
                    "genomic_accession": target_meta.get("genomic_accession", ""),
                    "genomic_start1": genomic_start,
                    "genomic_end1": genomic_end,
                    "ref": "",
                    "alt": alt,
                }
            )
            event_index += 1
        elif op.startswith("-"):
            ref = op[1:].upper()
            genomic_start, genomic_end = genomic_coords(target_meta, target_pos, target_pos + len(ref))
            events.append(
                {
                    **record,
                    "event_id": f"{event_id_prefix}:{event_index}",
                    "event_type": "del",
                    "target_start0": target_pos,
                    "target_end0": target_pos + len(ref),
                    "genomic_accession": target_meta.get("genomic_accession", ""),
                    "genomic_start1": genomic_start,
                    "genomic_end1": genomic_end,
                    "ref": ref,
                    "alt": "",
                }
            )
            event_index += 1
            target_pos += len(ref)
        elif op.startswith("~"):
            intron_len = int(re.search(r"\d+", op).group(0))
            target_pos += intron_len

    return events


def run_minimap2(
    minimap2_bin: str,
    preset: str,
    target_fa: Path,
    query_fa: Path,
    paf_path: Path,
    threads: int,
) -> str:
    cmd = [
        minimap2_bin,
        "-t",
        str(threads),
        "-x",
        preset,
        "-c",
        "--cs=long",
        "--paf-no-hit",
        str(target_fa),
        str(query_fa),
    ]
    with paf_path.open("w") as handle:
        result = subprocess.run(cmd, text=True, stdout=handle, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(f"minimap2 failed for preset={preset}: {result.stderr.strip()}")
    return " ".join(cmd)


def parse_paf(
    paf_path: Path,
    gene_id: str,
    strategy: str,
    preset: str,
    target_meta: dict[str, str],
    meta_by_sequence: dict[str, dict[str, str]],
    summaries: dict[str, dict[str, object]],
    query_slices: dict[str, QuerySlice] | None = None,
    accepted_record_ids: frozenset[str] | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    segments: list[dict[str, object]] = []
    event_by_key: dict[tuple[object, ...], dict[str, object]] = {}
    query_slices = query_slices or full_query_slices(meta_by_sequence)

    for fields, native_record_id, line_number in iter_paf_records(paf_path):
        qname = fields[0]
        query_slice = require_query_slice(
            query_slices,
            qname,
            paf_path,
            line_number,
        )
        source_id = query_slice.source_sequence_id
        meta = meta_by_sequence.get(source_id)
        if meta is None:
            raise ValueError(
                f"PAF query {qname!r} in {paf_path} at line {line_number} "
                f"references source sequence {source_id!r} without metadata"
            )
        summary = summaries[source_id]

        if fields[5] == "*" or fields[4] == "*":
            if query_slice.is_pseudoread:
                summary["qc_flags"].add("has_unmapped_pseudoread")
            else:
                summary["status"] = "no_alignment"
                summary["qc_flags"].add("no_alignment")
            continue
        if accepted_record_ids is not None and native_record_id not in accepted_record_ids:
            continue

        qstart = query_slice.source_start0 + int(fields[2])
        qend = query_slice.source_start0 + int(fields[3])
        strand = fields[4]
        target_id = fields[5]
        target_length = int(fields[6])
        target_start = int(fields[7])
        target_end = int(fields[8])
        matches = int(fields[9])
        block_length = int(fields[10])
        mapq = int(fields[11])
        tags = parse_tags(fields)
        native_alignment_type = tags.get("tp", "")
        primary = is_primary(tags)
        identity = matches / block_length if block_length else 0.0
        flags = []
        if query_slice.is_pseudoread:
            flags.append("filtered_pseudoread")

        segment = {
            "gene_id": gene_id,
            "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
            "tax_id": meta.get("tax_id", ""),
            "taxname": meta.get("taxname", ""),
            "strategy": strategy,
            "tool": "minimap2",
            "preset": preset,
            "sequence_id": source_id,
            "target_id": target_id,
            "query_id": source_id,
            "target_start0": target_start,
            "target_end0": target_end,
            "query_start0": qstart,
            "query_end0": qend,
            "strand": strand,
            "matches": matches,
            "block_length": block_length,
            "identity": f"{identity:.6f}",
            "mapq": mapq,
            "is_primary": str(primary).lower(),
            "divergence": tags.get("dv", ""),
            "gap_compressed_divergence": tags.get("de", ""),
            "native_record_id": native_record_id,
            "qc_flags": ",".join(flags),
        }
        segments.append(segment)

        summary["status"] = "aligned"
        summary["segment_count"] += 1
        summary["target_length"] = target_length
        summary["query_length"] = query_slice.source_length
        summary["identities"].append(identity)
        summary["best_identity"] = max(summary["best_identity"], identity)
        if query_slice.is_pseudoread:
            summary["qc_flags"].add("pseudoread_backbone")
        if primary:
            summary["primary_segment_count"] += 1
            summary["target_intervals"].append((target_start, target_end))
            summary["query_intervals"].append((qstart, qend))
        else:
            summary["secondary_segment_count"] += 1
            summary["qc_flags"].add("has_secondary")

        cs = tags.get("cs")
        if cs:
            event_record = {
                "gene_id": gene_id,
                "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
                "tax_id": meta.get("tax_id", ""),
                "taxname": meta.get("taxname", ""),
                "strategy": strategy,
                "tool": "minimap2",
                "preset": preset,
                "query_id": source_id,
                "strand": strand,
                "mapq": mapq,
                "native_alignment_type": native_alignment_type,
                "native_record_id": native_record_id,
                "qc_flags": ",".join(flags),
            }
            new_events = cs_events(
                cs,
                target_start,
                event_record,
                target_meta,
                f"{strategy}:{native_record_id}",
            )
            for event in new_events:
                event["_is_primary"] = primary
                key = (
                    event["query_id"],
                    event["event_type"],
                    event["target_start0"],
                    event["target_end0"],
                    event["ref"],
                    event["alt"],
                )
                current = event_by_key.get(key)
                candidate_rank = (not primary, str(event["native_record_id"]))
                current_rank = (
                    not bool(current["_is_primary"]),
                    str(current["native_record_id"]),
                ) if current is not None else None
                if current_rank is None or candidate_rank < current_rank:
                    event_by_key[key] = event

    events = sorted(
        event_by_key.values(),
        key=lambda row: (
            str(row["query_id"]),
            int(row["target_start0"]),
            str(row["event_type"]),
            str(row["ref"]),
            str(row["alt"]),
        ),
    )
    for event in events:
        event.pop("_is_primary", None)
        summaries[str(event["query_id"])]["event_count"] += 1

    return segments, events


def empty_summary(gene_id: str, strategy: str, preset: str, meta: dict[str, str], target_length: int) -> dict[str, object]:
    return {
        "gene_id": gene_id,
        "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
        "tax_id": meta.get("tax_id", ""),
        "taxname": meta.get("taxname", ""),
        "strategy": strategy,
        "tool": "minimap2",
        "preset": preset,
        "status": "not_run",
        "target_length": target_length,
        "query_length": int(meta.get("sequence_length") or 0),
        "segment_count": 0,
        "primary_segment_count": 0,
        "secondary_segment_count": 0,
        "target_intervals": [],
        "query_intervals": [],
        "best_identity": 0.0,
        "identities": [],
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
            "target_coverage": f"{(aligned_target / target_length) if target_length else 0.0:.6f}",
            "query_coverage": f"{(aligned_query / query_length) if query_length else 0.0:.6f}",
            "best_identity": f"{float(row['best_identity']):.6f}",
            "mean_identity": f"{mean_identity:.6f}",
            "qc_flags": ",".join(sorted(flags)),
        }
    )
    return row


def main() -> None:
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")
    validate_query_mode(args.preset, args.pseudoread_len, args.pseudoread_step)
    software = minimap2_software(args.minimap2_bin)
    args.outdir.mkdir(parents=True, exist_ok=True)

    task, target_meta, ortholog_meta = load_task_context(args.task_dir)
    gene_id = task["gene_id"]
    meta_by_sequence = {row["sequence_id"]: row for row in ortholog_meta}
    target_length = int(target_meta["sequence_length"])

    summaries = {
        sequence_id: empty_summary(
            gene_id,
            args.strategy,
            args.preset,
            meta,
            target_length,
        )
        for sequence_id, meta in meta_by_sequence.items()
    }

    all_segments: list[dict[str, object]] = []
    all_events: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    commands: list[str] = []

    with tempfile.TemporaryDirectory(prefix=f"{args.strategy}_", dir=args.outdir) as tmp_name:
        work_dir = Path(tmp_name)
        target_fasta, orthologs_fasta = materialize_task_fastas(
            args.source_target_fasta,
            args.source_ortholog_fasta,
            task,
            ortholog_meta,
            work_dir,
        )
        query_fasta = orthologs_fasta
        query_slices = full_query_slices(meta_by_sequence)
        pseudoreads: PseudoreadGeneration | None = None
        backbone: BackboneSelection | None = None
        if args.pseudoread_len:
            query_fasta = work_dir / "long_pseudoreads.fa"
            pseudoreads = generate_long_pseudoreads(
                orthologs_fasta,
                query_fasta,
                meta_by_sequence,
                args.pseudoread_len,
                args.pseudoread_step,
            )
            query_slices = pseudoreads.query_slices

        paf_path = work_dir / f"{args.strategy}.{args.preset}.paf"
        try:
            command = run_minimap2(
                args.minimap2_bin,
                args.preset,
                target_fasta,
                query_fasta,
                paf_path,
                args.threads,
            )
            commands.append(command)
            accepted_record_ids: frozenset[str] | None = None
            if pseudoreads is not None:
                backbone = select_pseudoread_backbone(paf_path, query_slices)
                accepted_record_ids = backbone.accepted_record_ids
            segments, events = parse_paf(
                paf_path,
                gene_id,
                args.strategy,
                args.preset,
                target_meta,
                meta_by_sequence,
                summaries,
                query_slices,
                accepted_record_ids,
            )
            all_segments.extend(segments)
            all_events.extend(events)
        except Exception as exc:
            failures.append(
                {
                    "gene_id": gene_id,
                    "ortholog_gene_id": "",
                    "strategy": args.strategy,
                    "tool": "minimap2",
                    "failure_type": "minimap2_failed",
                    "message": str(exc),
                }
            )
            raise

    summary_rows = [finalize_summary(row) for row in summaries.values()]
    write_tsv_gz(args.outdir / "alignment_segments.tsv.gz", SEGMENT_FIELDS, all_segments)
    write_tsv_gz(args.outdir / "alignment_events.tsv.gz", EVENT_FIELDS, all_events)
    write_tsv_gz(args.outdir / "ortholog_alignment_summary.tsv.gz", SUMMARY_FIELDS, summary_rows)
    write_tsv_gz(args.outdir / "failures.tsv.gz", FAILURE_FIELDS, failures)
    strategy_parameters: dict[str, object] = {
        "preset": args.preset,
        "mapq_policy": "aligner_reported_unfiltered",
        "non_primary_policy": "retained_native_type",
        "software": software,
    }
    if pseudoreads is not None:
        strategy_parameters.update(
            {
                "pseudoread_len": args.pseudoread_len,
                "pseudoread_step": args.pseudoread_step,
                "sequencing_error_model": "none",
                "backbone_policy": "dominant_strand_lis",
            }
        )
    manifest = {
        "gene_ids": [gene_id],
        "strategies": [args.strategy],
        "strategy_parameters": {args.strategy: strategy_parameters},
        "tool": "minimap2",
        "commands": commands,
        "ortholog_alignment_summary_count": len(summary_rows),
        "alignment_segment_count": len(all_segments),
        "alignment_event_mode": "raw",
        "raw_alignment_event_count": len(all_events),
        "alignment_event_count": len(all_events),
        "failure_count": len(failures),
        "ortholog_count": len(ortholog_meta),
    }
    if pseudoreads is not None and backbone is not None:
        manifest.update(
            {
                "pseudoread_count": pseudoreads.total_reads,
                "pseudoread_alignment_record_count": backbone.input_alignment_count,
                "pseudoread_after_strand_count": backbone.after_strand_count,
                "pseudoread_backbone_record_count": backbone.retained_alignment_count,
            }
        )
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
