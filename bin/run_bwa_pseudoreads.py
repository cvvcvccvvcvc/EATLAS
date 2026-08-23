#!/usr/bin/env python3
"""Run BWA pseudoread alignment and emit normalized Stage 2 evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pysam

from bin.alignment_table_schema import (
    EVENT_FIELDS,
    FAILURE_FIELDS,
    SEGMENT_FIELDS,
    SUMMARY_FIELDS,
)
from bin.alignment_task_io import load_task_context, materialize_task_fastas
from bin import bwa_pseudoread_filter

EventKey = tuple[str, int, int, str, str]
EventSupport = dict[EventKey, dict[str, dict[str, object]]]


@dataclass(frozen=True)
class PseudoreadGeneration:
    total_reads: int
    query_lengths: dict[str, int]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-dir", required=True, type=Path)
    parser.add_argument("--source-target-fasta", required=True, type=Path)
    parser.add_argument("--source-ortholog-fasta", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--strategy", required=True)
    parser.add_argument("--bwa-bin", default="bwa")
    parser.add_argument("--samtools-bin", default="samtools")
    parser.add_argument("--threads", default=2, type=int)
    parser.add_argument("--pseudoread-len", required=True, type=int)
    parser.add_argument("--pseudoread-step", required=True, type=int)
    parser.add_argument("--pseudoread-phred", required=True, type=int)
    return parser.parse_args()


def run_checked(cmd: list[str], *, stdout=None) -> None:
    print("Running: " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, stdout=stdout)


def run_bwa_mem_pipeline(
    bwa_bin: str,
    samtools_bin: str,
    target_fasta: Path,
    fastq: Path,
    sorted_bam: Path,
    threads: int,
) -> int:
    if threads < 2:
        raise ValueError("--threads must be at least 2 for concurrent BWA and samtools")

    bwa_threads = threads - 1
    bwa_cmd = [
        bwa_bin,
        "mem",
        "-t",
        str(bwa_threads),
        str(target_fasta),
        str(fastq),
    ]
    sort_cmd = [samtools_bin, "sort", "-o", str(sorted_bam), "-"]

    run_checked([bwa_bin, "index", str(target_fasta)])
    print(
        "Running: "
        + " ".join(bwa_cmd)
        + " | "
        + " ".join(sort_cmd),
        flush=True,
    )
    bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
    assert bwa_proc.stdout is not None
    sort_proc = subprocess.Popen(sort_cmd, stdin=bwa_proc.stdout)
    bwa_proc.stdout.close()

    sort_code = sort_proc.wait()
    bwa_code = bwa_proc.wait()
    failed = [
        (name, code)
        for name, code in [("bwa mem", bwa_code), ("samtools sort", sort_code)]
        if code != 0
    ]
    if failed:
        details = ", ".join(f"{name} exit {code}" for name, code in failed)
        raise subprocess.CalledProcessError(failed[0][1], details)
    run_checked([samtools_bin, "index", str(sorted_bam)])
    return bwa_threads


def iter_fasta(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header = None
    seq_parts: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


def write_tsv_gz(path: Path, headers: list[str], rows) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in headers})
            count += 1
    return count


def pseudoread_starts(
    seq_len: int,
    read_len: int,
    step: int,
    min_read_len: int = 20,
) -> list[int]:
    """Return sliding-window starts, including the sequence's final window."""
    if step <= 0:
        raise ValueError("Pseudoread step must be positive")
    if seq_len < min_read_len or read_len < min_read_len:
        return []
    if seq_len <= read_len:
        return [0]

    final_start = seq_len - read_len
    starts = list(range(0, final_start + 1, step))
    if starts[-1] != final_start:
        starts.append(final_start)
    return starts


def expected_pseudoreads(seq_len: int, read_len: int, step: int) -> int:
    """Count pseudo-reads produced for one sequence."""
    return len(pseudoread_starts(seq_len, read_len, step))


def generate_pseudoreads(
    orthologs_fa: Path,
    out_fastq: Path,
    read_len: int,
    step: int,
    phred: int,
) -> PseudoreadGeneration:
    if read_len < 20:
        raise ValueError("Pseudoread length must be at least 20")
    if not 0 <= phred <= 93:
        raise ValueError("Pseudoread PHRED must be between 0 and 93")
    phred_char = chr(phred + 33)
    total_reads = 0
    query_lengths: dict[str, int] = {}
    with out_fastq.open("w") as out:
        for header, seq in iter_fasta(orthologs_fa):
            ortholog_id = header.split()[0]
            query_lengths[ortholog_id.removeprefix("ortholog_")] = len(seq)
            read_index = 1
            for start in pseudoread_starts(len(seq), read_len, step):
                read_seq = seq[start : start + read_len]
                qual = phred_char * len(read_seq)
                end = start + len(read_seq)
                out.write(f"@{ortholog_id}_pseudo_{read_index}_{start + 1}-{end}\n{read_seq}\n+\n{qual}\n")
                read_index += 1
                total_reads += 1
    return PseudoreadGeneration(total_reads, query_lengths)


def read_ortholog_gene_id(read_name: str) -> str:
    prefix = read_name.split("_pseudo_", 1)[0]
    return prefix.removeprefix("ortholog_")


def genomic_coords(target_meta: dict[str, str], start0: int, end0: int) -> tuple[str, str]:
    begin_text = target_meta.get("genomic_begin") or ""
    if not begin_text:
        return "", ""
    begin = int(begin_text)
    genomic_start = begin + start0
    genomic_end = begin + end0 - 1 if end0 > start0 else genomic_start
    return str(genomic_start), str(genomic_end)


def interval_union_length(intervals: list[tuple[int, int]]) -> int:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def fmt_fraction(numerator: int | float, denominator: int | float) -> str:
    if denominator <= 0:
        return "0.000000"
    return f"{numerator / denominator:.6f}"


def read_target_sequence(target_fasta: Path) -> str:
    for _, seq in iter_fasta(target_fasta):
        return seq
    raise ValueError(f"Target FASTA has no records: {target_fasta}")


def sequence_length_from_read(read: pysam.AlignedSegment) -> int:
    return read.infer_query_length(always=True) or len(read.query_sequence or "")


def cigar_block_length(read: pysam.AlignedSegment) -> int:
    return sum(length for op, length in (read.cigartuples or []) if op in {0, 1, 2, 7, 8})


def identity_from_read(read: pysam.AlignedSegment) -> tuple[int, int, str]:
    block_length = cigar_block_length(read)
    if block_length <= 0:
        return 0, 0, "0.000000"
    try:
        nm = int(read.get_tag("NM"))
    except KeyError:
        nm = 0
    matches = max(0, block_length - nm)
    return matches, block_length, fmt_fraction(matches, block_length)


def sam_alignment_type(read: pysam.AlignedSegment) -> str:
    if read.is_secondary and read.is_supplementary:
        return "secondary_supplementary"
    if read.is_secondary:
        return "secondary"
    if read.is_supplementary:
        return "supplementary"
    return "primary"


def add_event_support(
    event_support: EventSupport,
    key: EventKey,
    ortholog_id: str,
    strand: str,
    native_record_id: str,
    mapq: int,
    native_alignment_type: str,
    is_primary: bool,
) -> None:
    support = {
        "ortholog_gene_id": ortholog_id,
        "strand": strand,
        "native_record_id": native_record_id,
        "mapq": mapq,
        "native_alignment_type": native_alignment_type,
        "is_primary": is_primary,
    }
    current = event_support[key].get(ortholog_id)
    if current is None or (
        not is_primary,
        str(support["native_record_id"]),
        str(support["strand"]),
    ) < (
        not bool(current["is_primary"]),
        str(current["native_record_id"]),
        str(current["strand"]),
    ):
        event_support[key][ortholog_id] = support


def scan_bam(
    bam_path: Path,
    target_seq: str,
) -> tuple[list[dict[str, object]], EventSupport]:
    segments: list[dict[str, object]] = []
    event_support: EventSupport = defaultdict(dict)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.reference_start is None or read.reference_end is None:
                continue
            ortholog_id = read_ortholog_gene_id(read.query_name)
            strand = "-" if read.is_reverse else "+"
            primary = not read.is_secondary and not read.is_supplementary
            native_alignment_type = sam_alignment_type(read)
            matches, block_length, identity = identity_from_read(read)
            segments.append(
                {
                    "ortholog_gene_id": ortholog_id,
                    "sequence_id": f"ortholog_{ortholog_id}",
                    "query_id": f"ortholog_{ortholog_id}",
                    "target_start0": read.reference_start,
                    "target_end0": read.reference_end,
                    "query_start0": read.query_alignment_start or 0,
                    "query_end0": read.query_alignment_end or sequence_length_from_read(read),
                    "strand": strand,
                    "matches": matches,
                    "block_length": block_length,
                    "identity": identity,
                    "mapq": read.mapping_quality,
                    "is_primary": str(primary).lower(),
                    "native_record_id": read.query_name,
                    "qc_flags": "filtered_pseudoread",
                }
            )

            q_seq = read.query_sequence or ""
            r_pos = read.reference_start
            q_pos = 0
            for op, length in read.cigartuples or []:
                if op in {0, 7, 8}:
                    if op == 7:
                        r_pos += length
                        q_pos += length
                        continue
                    for offset in range(length):
                        ref_index = r_pos + offset
                        query_index = q_pos + offset
                        if ref_index >= len(target_seq) or query_index >= len(q_seq):
                            continue
                        ref = target_seq[ref_index].upper()
                        alt = q_seq[query_index].upper()
                        if ref == alt:
                            continue
                        key = ("snv", ref_index, ref_index + 1, ref, alt)
                        add_event_support(
                            event_support,
                            key,
                            ortholog_id,
                            strand,
                            read.query_name,
                            read.mapping_quality,
                            native_alignment_type,
                            primary,
                        )
                    r_pos += length
                    q_pos += length
                elif op == 1:
                    alt = q_seq[q_pos : q_pos + length].upper()
                    key = ("ins", r_pos, r_pos, "", alt)
                    add_event_support(
                        event_support,
                        key,
                        ortholog_id,
                        strand,
                        read.query_name,
                        read.mapping_quality,
                        native_alignment_type,
                        primary,
                    )
                    q_pos += length
                elif op in {2, 3}:
                    if op == 2:
                        ref = target_seq[r_pos : r_pos + length].upper()
                        key = ("del", r_pos, r_pos + length, ref, "")
                        add_event_support(
                            event_support,
                            key,
                            ortholog_id,
                            strand,
                            read.query_name,
                            read.mapping_quality,
                            native_alignment_type,
                            primary,
                        )
                    r_pos += length
                elif op == 4:
                    q_pos += length
                elif op in {5, 6}:
                    continue
    return segments, event_support


def make_segment_rows(
    base_segments: list[dict[str, object]],
    gene_id: str,
    target_id: str,
    ortholog_meta_by_id: dict[str, dict[str, str]],
    strategy: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for segment in base_segments:
        ortholog_id = str(segment["ortholog_gene_id"])
        meta = ortholog_meta_by_id.get(ortholog_id, {})
        row = {
            "gene_id": gene_id,
            "tax_id": meta.get("tax_id", ""),
            "taxname": meta.get("taxname", ""),
            "strategy": strategy,
            "tool": "bwa",
            "preset": "pseudo",
            "target_id": target_id,
            "divergence": "",
            "gap_compressed_divergence": "",
        }
        row.update(segment)
        rows.append(row)
    return rows


def make_event_row(
    *,
    gene_id: str,
    target_meta: dict[str, str],
    target_acc: str,
    ortholog_meta_by_id: dict[str, dict[str, str]],
    strategy: str,
    event_type: str,
    start0: int,
    end0: int,
    ref: str,
    alt: str,
    support: dict[str, object],
    qc_flags: str,
) -> dict[str, object]:
    ortholog_id = str(support["ortholog_gene_id"])
    meta = ortholog_meta_by_id.get(ortholog_id, {})
    genomic_start, genomic_end = genomic_coords(target_meta, start0, end0)
    event_id = f"{strategy}:{ortholog_id}:{event_type}:{start0}:{end0}:{ref}>{alt}"
    return {
        "gene_id": gene_id,
        "ortholog_gene_id": ortholog_id,
        "tax_id": meta.get("tax_id", ""),
        "taxname": meta.get("taxname", ""),
        "strategy": strategy,
        "tool": "bwa",
        "preset": "pseudo",
        "event_id": event_id,
        "event_type": event_type,
        "target_start0": start0,
        "target_end0": end0,
        "genomic_accession": target_acc,
        "genomic_start1": genomic_start,
        "genomic_end1": genomic_end,
        "ref": ref,
        "alt": alt,
        "query_id": f"ortholog_{ortholog_id}",
        "strand": support.get("strand", ""),
        "mapq": support.get("mapq", ""),
        "native_alignment_type": support.get("native_alignment_type", ""),
        "native_record_id": support.get("native_record_id", ""),
        "qc_flags": qc_flags,
    }


def make_bwa_event_rows(
    event_support: EventSupport,
    gene_id: str,
    target_meta: dict[str, str],
    target_acc: str,
    ortholog_meta_by_id: dict[str, dict[str, str]],
    strategy: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (event_type, start0, end0, ref, alt), support_by_ortholog in sorted(event_support.items()):
        for ortholog_id in sorted(support_by_ortholog, key=int):
            support = support_by_ortholog[ortholog_id]
            flags = ["bwa_cigar_event"]
            rows.append(
                make_event_row(
                    gene_id=gene_id,
                    target_meta=target_meta,
                    target_acc=target_acc,
                    ortholog_meta_by_id=ortholog_meta_by_id,
                    strategy=strategy,
                    event_type=event_type,
                    start0=start0,
                    end0=end0,
                    ref=ref,
                    alt=alt,
                    support=support,
                    qc_flags=",".join(flags),
                )
            )
    return rows


def make_summary_rows(
    gene_id: str,
    target_length: int,
    ortholog_meta: list[dict[str, str]],
    query_lengths: dict[str, int],
    segment_rows: list[dict[str, object]],
    event_rows: list[dict[str, object]],
    strategy: str,
) -> list[dict[str, object]]:
    segments_by_key: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in segment_rows:
        segments_by_key[(str(row["strategy"]), str(row["ortholog_gene_id"]))].append(row)
    event_counts: dict[tuple[str, str], int] = defaultdict(int)
    for row in event_rows:
        event_counts[(str(row["strategy"]), str(row["ortholog_gene_id"]))] += 1

    summaries: list[dict[str, object]] = []
    for meta in ortholog_meta:
        ortholog_id = meta["ortholog_gene_id"]
        segments = segments_by_key.get((strategy, ortholog_id), [])
        target_intervals = [(int(row["target_start0"]), int(row["target_end0"])) for row in segments]
        query_intervals = [(int(row["query_start0"]), int(row["query_end0"])) for row in segments]
        identities = [float(row["identity"]) for row in segments if row.get("identity") not in {"", None}]
        primary_count = sum(1 for row in segments if str(row.get("is_primary", "")).lower() == "true")
        secondary_count = len(segments) - primary_count
        aligned_target_bp = interval_union_length(target_intervals)
        query_length = int(query_lengths.get(ortholog_id) or meta.get("sequence_length") or 0)
        aligned_query_bp = interval_union_length(query_intervals)
        event_count = event_counts.get((strategy, ortholog_id), 0)
        summaries.append(
            {
                "gene_id": gene_id,
                "ortholog_gene_id": ortholog_id,
                "tax_id": meta.get("tax_id", ""),
                "taxname": meta.get("taxname", ""),
                "strategy": strategy,
                "tool": "bwa",
                "preset": "pseudo",
                "status": "aligned" if segments else "no_alignment",
                "target_length": target_length,
                "query_length": query_length,
                "segment_count": len(segments),
                "primary_segment_count": primary_count,
                "secondary_segment_count": secondary_count,
                "aligned_target_bp": aligned_target_bp,
                "aligned_query_bp": aligned_query_bp,
                "target_coverage": fmt_fraction(aligned_target_bp, target_length),
                "query_coverage": fmt_fraction(aligned_query_bp, query_length),
                "best_identity": f"{max(identities):.6f}" if identities else "0.000000",
                "mean_identity": f"{sum(identities) / len(identities):.6f}" if identities else "0.000000",
                "event_count": event_count,
                "qc_flags": "pseudoread_bam_segments" if segments else "no_filtered_pseudoread_segments",
            }
        )
    return summaries


def main() -> None:
    args = parse_args()
    task_dir = args.task_dir
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    manifest, target_meta, ortholog_meta = load_task_context(task_dir)
    gene_id = manifest["gene_id"]
    target_id = target_meta["sequence_id"]
    ortholog_meta_by_id = {row["ortholog_gene_id"]: row for row in ortholog_meta}
    target_acc = target_meta.get("genomic_accession", "")

    with tempfile.TemporaryDirectory(prefix=f"{args.strategy}_", dir=outdir) as tmp_name:
        work_dir = Path(tmp_name)
        local_target_fasta, local_orthologs_fasta = materialize_task_fastas(
            args.source_target_fasta,
            args.source_ortholog_fasta,
            manifest,
            ortholog_meta,
            work_dir,
        )
        target_seq = read_target_sequence(local_target_fasta)

        pseudoreads_fastq = work_dir / "pseudo_reads.fastq"
        pseudoreads = generate_pseudoreads(
            local_orthologs_fasta,
            pseudoreads_fastq,
            read_len=args.pseudoread_len,
            step=args.pseudoread_step,
            phred=args.pseudoread_phred,
        )

        sorted_bam = work_dir / "aln.sorted.bam"
        run_bwa_mem_pipeline(
            args.bwa_bin,
            args.samtools_bin,
            local_target_fasta,
            pseudoreads_fastq,
            sorted_bam,
            args.threads,
        )

        bwa_pseudoread_filter.filter_bam_for_gene(work_dir)
        filtered_bam = work_dir / "aln.filtered.lis.bam"

        base_segments, event_support = scan_bam(filtered_bam, target_seq)
        segment_rows = make_segment_rows(
            base_segments,
            gene_id,
            target_id,
            ortholog_meta_by_id,
            args.strategy,
        )

        event_rows = make_bwa_event_rows(
            event_support,
            gene_id,
            target_meta,
            target_acc,
            ortholog_meta_by_id,
            args.strategy,
        )

        summary_rows = make_summary_rows(
            gene_id,
            len(target_seq),
            ortholog_meta,
            pseudoreads.query_lengths,
            segment_rows,
            event_rows,
            args.strategy,
        )

        write_tsv_gz(outdir / "alignment_segments.tsv.gz", SEGMENT_FIELDS, segment_rows)
        write_tsv_gz(outdir / "alignment_events.tsv.gz", EVENT_FIELDS, event_rows)
        write_tsv_gz(outdir / "ortholog_alignment_summary.tsv.gz", SUMMARY_FIELDS, summary_rows)
        write_tsv_gz(outdir / "failures.tsv.gz", FAILURE_FIELDS, [])
        manifest_out = {
            "gene_ids": [gene_id],
            "strategies": [args.strategy],
            "strategy_parameters": {
                args.strategy: {
                    "pseudoread_len": args.pseudoread_len,
                    "pseudoread_step": args.pseudoread_step,
                    "pseudoread_phred": args.pseudoread_phred,
                }
            },
            "tool": "bwa",
            "ortholog_alignment_summary_count": len(summary_rows),
            "alignment_segment_count": len(segment_rows),
            "alignment_event_mode": "raw",
            "raw_alignment_event_count": len(event_rows),
            "alignment_event_count": len(event_rows),
            "failure_count": 0,
            "ortholog_count": len(ortholog_meta),
            "pseudoread_count": pseudoreads.total_reads,
        }
        (outdir / "manifest.json").write_text(json.dumps(manifest_out, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
