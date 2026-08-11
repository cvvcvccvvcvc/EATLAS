#!/usr/bin/env python3
"""Run nucmer for one gene task and normalize comparator alignment evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable

import pysam

from alignment_task_io import load_task_context, materialize_task_fastas
from feature_coverage import summarize_feature_coverage_rows


TSV_NULL = ""

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-dir", required=True, type=Path)
    parser.add_argument("--source-target-fasta", required=True, type=Path)
    parser.add_argument("--source-ortholog-fasta", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--nucmer-bin", default="nucmer")
    parser.add_argument("--threads", default=1, type=int)
    parser.add_argument("--target-features", type=Path)
    parser.add_argument("--keep-native", default="false")
    return parser.parse_args()


def truthy(value: str) -> bool:
    return str(value).lower() in {"1", "true", "yes", "y"}


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


def run_command(cmd: list[str]) -> str:
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {(result.stderr or result.stdout or '').strip()}")
    return " ".join(cmd)


def empty_summary(gene_id: str, meta: dict[str, str], target_length: int) -> dict[str, object]:
    return {
        "gene_id": gene_id,
        "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
        "tax_id": meta.get("tax_id", ""),
        "taxname": meta.get("taxname", ""),
        "strategy": "nucmer",
        "tool": "nucmer",
        "preset": "default",
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
        "qc_flags": {"unfiltered_nucmer"},
    }


def read_first_fasta_sequence(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    parts: list[str] = []
    found_record = False
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if found_record:
                    break
                found_record = True
                continue
            if found_record:
                parts.append(line)
    if not parts:
        raise ValueError(f"Target FASTA has no sequence: {path}")
    return "".join(parts).upper()


def cigar_block_length(read: pysam.AlignedSegment) -> int:
    return sum(length for op, length in (read.cigartuples or []) if op in {0, 1, 2, 7, 8})


def query_interval(read: pysam.AlignedSegment) -> tuple[int, int, int]:
    cigartuples = read.cigartuples or []
    query_length = sum(length for op, length in cigartuples if op in {0, 1, 4, 5, 7, 8})
    if query_length <= 0:
        query_length = read.infer_query_length(always=True) or len(read.query_sequence or "")

    leading_clip = 0
    for op, length in cigartuples:
        if op not in {4, 5}:
            break
        leading_clip += length

    trailing_clip = 0
    for op, length in reversed(cigartuples):
        if op not in {4, 5}:
            break
        trailing_clip += length

    oriented_start = leading_clip
    oriented_end = query_length - trailing_clip
    if read.is_reverse:
        return query_length - oriented_end, query_length - oriented_start, query_length
    return oriented_start, oriented_end, query_length


def append_event(
    event_by_key: dict[tuple[object, ...], dict[str, object]],
    *,
    gene_id: str,
    target_meta: dict[str, str],
    meta: dict[str, str],
    query_id: str,
    strand: str,
    native_record_id: int,
    event_ordinal: int,
    event_type: str,
    target_start0: int,
    target_end0: int,
    ref: str,
    alt: str,
    is_primary: bool,
) -> bool:
    ref = ref.upper()
    alt = alt.upper()
    if any(base not in DNA_BASES for base in ref + alt):
        return False

    genomic_start, genomic_end = genomic_coords(target_meta, target_start0, target_end0)
    flags = {"unfiltered_nucmer"}
    if not is_primary:
        flags.add("non_primary")
    event = {
        "gene_id": gene_id,
        "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
        "tax_id": meta.get("tax_id", ""),
        "taxname": meta.get("taxname", ""),
        "strategy": "nucmer",
        "tool": "nucmer",
        "preset": "default",
        "event_id": f"{native_record_id}:{event_ordinal}",
        "event_type": event_type,
        "target_start0": target_start0,
        "target_end0": target_end0,
        "genomic_accession": target_meta.get("genomic_accession", ""),
        "genomic_start1": genomic_start,
        "genomic_end1": genomic_end,
        "ref": ref,
        "alt": alt,
        "query_id": query_id,
        "strand": strand,
        "native_record_id": native_record_id,
        "qc_flags": ",".join(sorted(flags)),
        "_is_primary": is_primary,
    }
    key = (
        query_id,
        event_type,
        target_start0,
        target_end0,
        ref,
        alt,
    )
    current = event_by_key.get(key)
    if current is None or (is_primary and not bool(current["_is_primary"])):
        event_by_key[key] = event
    return True


def parse_sam(
    sam_path: Path,
    gene_id: str,
    target_meta: dict[str, str],
    target_seq: str,
    meta_by_sequence: dict[str, dict[str, str]],
    summaries: dict[str, dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], int]:
    segments: list[dict[str, object]] = []
    event_by_key: dict[tuple[object, ...], dict[str, object]] = {}
    ambiguous_event_allele_count = 0

    with pysam.AlignmentFile(str(sam_path), "r") as sam:
        for native_record_id, read in enumerate(sam.fetch(until_eof=True), start=1):
            if read.is_unmapped or read.reference_start is None or read.reference_end is None:
                continue
            query_id = read.query_name
            meta = meta_by_sequence.get(query_id)
            if meta is None:
                continue

            query_start0, query_end0, query_length = query_interval(read)
            target_start0 = int(read.reference_start)
            target_end0 = int(read.reference_end)
            strand = "-" if read.is_reverse else "+"
            block_length = cigar_block_length(read)
            try:
                edit_distance = int(read.get_tag("NM"))
            except KeyError:
                edit_distance = 0
            matches = max(0, block_length - edit_distance)
            identity = matches / block_length if block_length else 0.0
            is_primary = not read.is_secondary and not read.is_supplementary
            flags = {"unfiltered_nucmer"}
            if not is_primary:
                flags.add("non_primary")

            segment = {
                "gene_id": gene_id,
                "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
                "tax_id": meta.get("tax_id", ""),
                "taxname": meta.get("taxname", ""),
                "strategy": "nucmer",
                "tool": "nucmer",
                "preset": "default",
                "sequence_id": query_id,
                "target_id": read.reference_name or f"target_{gene_id}",
                "query_id": query_id,
                "target_start0": target_start0,
                "target_end0": target_end0,
                "query_start0": query_start0,
                "query_end0": query_end0,
                "strand": strand,
                "matches": matches,
                "block_length": block_length,
                "identity": f"{identity:.6f}",
                "mapq": "",
                "is_primary": str(is_primary).lower(),
                "divergence": "",
                "gap_compressed_divergence": "",
                "native_record_id": native_record_id,
                "qc_flags": ",".join(sorted(flags)),
            }
            segments.append(segment)

            summary = summaries[query_id]
            summary["status"] = "aligned"
            summary["target_length"] = len(target_seq)
            summary["query_length"] = query_length or summary["query_length"]
            summary["segment_count"] += 1
            if is_primary:
                summary["primary_segment_count"] += 1
            else:
                summary["secondary_segment_count"] += 1
            summary["target_intervals"].append((target_start0, target_end0))
            summary["query_intervals"].append((query_start0, query_end0))
            summary["identities"].append(identity)
            summary["best_identity"] = max(summary["best_identity"], identity)

            query_seq = (read.query_sequence or "").upper()
            if not query_seq:
                raise ValueError(f"Nucmer SAM record has no query sequence: {query_id}")
            ref_pos = target_start0
            query_pos = 0
            event_ordinal = 0
            for op, length in read.cigartuples or []:
                if op in {0, 7, 8}:
                    if op == 7:
                        ref_pos += length
                        query_pos += length
                        continue
                    for offset in range(length):
                        target_index = ref_pos + offset
                        query_index = query_pos + offset
                        if target_index >= len(target_seq) or query_index >= len(query_seq):
                            raise ValueError(
                                f"Nucmer SAM alignment exceeds sequence bounds for {query_id}"
                            )
                        ref = target_seq[target_index]
                        alt = query_seq[query_index]
                        if ref == alt:
                            continue
                        event_ordinal += 1
                        if not append_event(
                            event_by_key,
                            gene_id=gene_id,
                            target_meta=target_meta,
                            meta=meta,
                            query_id=query_id,
                            strand=strand,
                            native_record_id=native_record_id,
                            event_ordinal=event_ordinal,
                            event_type="snv",
                            target_start0=target_index,
                            target_end0=target_index + 1,
                            ref=ref,
                            alt=alt,
                            is_primary=is_primary,
                        ):
                            summary["qc_flags"].add("ambiguous_event_allele")
                            ambiguous_event_allele_count += 1
                    ref_pos += length
                    query_pos += length
                elif op == 1:
                    alt = query_seq[query_pos : query_pos + length]
                    if len(alt) != length:
                        raise ValueError(f"Truncated Nucmer insertion sequence for {query_id}")
                    event_ordinal += 1
                    if not append_event(
                        event_by_key,
                        gene_id=gene_id,
                        target_meta=target_meta,
                        meta=meta,
                        query_id=query_id,
                        strand=strand,
                        native_record_id=native_record_id,
                        event_ordinal=event_ordinal,
                        event_type="ins",
                        target_start0=ref_pos,
                        target_end0=ref_pos,
                        ref="",
                        alt=alt,
                        is_primary=is_primary,
                    ):
                        summary["qc_flags"].add("ambiguous_event_allele")
                        ambiguous_event_allele_count += 1
                    query_pos += length
                elif op == 2:
                    ref = target_seq[ref_pos : ref_pos + length]
                    if len(ref) != length:
                        raise ValueError(f"Truncated Nucmer deletion sequence for {query_id}")
                    event_ordinal += 1
                    if not append_event(
                        event_by_key,
                        gene_id=gene_id,
                        target_meta=target_meta,
                        meta=meta,
                        query_id=query_id,
                        strand=strand,
                        native_record_id=native_record_id,
                        event_ordinal=event_ordinal,
                        event_type="del",
                        target_start0=ref_pos,
                        target_end0=ref_pos + length,
                        ref=ref,
                        alt="",
                        is_primary=is_primary,
                    ):
                        summary["qc_flags"].add("ambiguous_event_allele")
                        ambiguous_event_allele_count += 1
                    ref_pos += length
                elif op == 3:
                    ref_pos += length
                elif op == 4:
                    query_pos += length
                elif op in {5, 6}:
                    continue

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
    return segments, events, ambiguous_event_allele_count


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


def gzip_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    with src.open("rb") as inp, gzip.open(dst, "wb") as out:
        shutil.copyfileobj(inp, out)


def main() -> None:
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")
    args.outdir.mkdir(parents=True, exist_ok=True)
    keep_native = truthy(args.keep_native)

    task, target_meta, ortholog_meta = load_task_context(args.task_dir)
    gene_id = task["gene_id"]
    meta_by_sequence = {row["sequence_id"]: row for row in ortholog_meta}
    target_length = int(target_meta["sequence_length"])
    summaries = {
        sequence_id: empty_summary(gene_id, meta, target_length)
        for sequence_id, meta in meta_by_sequence.items()
    }

    failures: list[dict[str, object]] = []
    commands: list[str] = []
    ambiguous_event_allele_count = 0

    with tempfile.TemporaryDirectory(prefix="nucmer_", dir=args.outdir) as tmp_name:
        work_dir = Path(tmp_name)
        target_fasta, orthologs_fasta = materialize_task_fastas(
            args.source_target_fasta,
            args.source_ortholog_fasta,
            task,
            ortholog_meta,
            work_dir,
        )
        sam_path = work_dir / "nucmer.sam"

        try:
            commands.append(
                run_command(
                    [
                        args.nucmer_bin,
                        "--threads",
                        str(args.threads),
                        f"--sam-long={sam_path}",
                        str(target_fasta),
                        str(orthologs_fasta),
                    ]
                )
            )
            target_seq = read_first_fasta_sequence(target_fasta)
            segments, events, ambiguous_event_allele_count = parse_sam(
                sam_path,
                gene_id,
                target_meta,
                target_seq,
                meta_by_sequence,
                summaries,
            )
            if keep_native:
                gzip_copy(sam_path, args.outdir / "native" / f"{gene_id}.sam.gz")
        except Exception as exc:
            failures.append(
                {
                    "gene_id": gene_id,
                    "ortholog_gene_id": "",
                    "strategy": "nucmer",
                    "tool": "nucmer",
                    "failure_type": "nucmer_failed",
                    "message": str(exc),
                }
            )
            raise

    summary_rows = [finalize_summary(row) for row in summaries.values()]
    write_tsv_gz(args.outdir / "alignment_segments.tsv.gz", SEGMENT_FIELDS, segments)
    write_tsv_gz(args.outdir / "alignment_events.tsv.gz", EVENT_FIELDS, events)
    write_tsv_gz(args.outdir / "ortholog_alignment_summary.tsv.gz", SUMMARY_FIELDS, summary_rows)
    write_tsv_gz(args.outdir / "failures.tsv.gz", FAILURE_FIELDS, failures)
    feature_coverage_count = None
    if args.target_features:
        feature_coverage_count = summarize_feature_coverage_rows(
            args.target_features,
            summary_rows,
            segments,
            args.outdir / "feature_coverage.tsv.gz",
        )
    manifest = {
        "gene_id": gene_id,
        "strategy": "nucmer",
        "tool": "nucmer",
        "commands": commands,
        "segment_count": len(segments),
        "event_count": len(events),
        "ambiguous_event_allele_count": ambiguous_event_allele_count,
        "feature_coverage_count": feature_coverage_count,
        "ortholog_count": len(ortholog_meta),
        "keep_native": keep_native,
        "filtering": "no global one-to-one filtering; SAM/CIGAR records are evaluated per ortholog",
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
