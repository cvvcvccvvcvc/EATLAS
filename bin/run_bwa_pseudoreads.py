#!/usr/bin/env python3
"""Run BWA pseudoread alignment and emit normalized Stage 2 evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pysam

from alignment_task_io import load_task_context, materialize_task_fastas
import bam_filtering_v1


BWA_STRATEGY = "bwa_pseudoreads"
BWA_VARSCAN_STRATEGY = "bwa_pseudoreads_varscan"
BWA_STRATEGIES = {BWA_STRATEGY, BWA_VARSCAN_STRATEGY}

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-dir", required=True, type=Path)
    parser.add_argument("--source-target-fasta", required=True, type=Path)
    parser.add_argument("--source-ortholog-fasta", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--strategies", required=True, help="Comma-separated BWA strategy names to emit.")
    parser.add_argument("--bwa-bin", default="bwa")
    parser.add_argument("--samtools-bin", default="samtools")
    parser.add_argument("--varscan-jar", default="VarScan.v2.4.6.jar")
    parser.add_argument("--varscan-min-coverage", default=2, type=int)
    parser.add_argument("--varscan-min-reads2", default=1, type=int)
    parser.add_argument("--varscan-min-var-freq", default=0.01, type=float)
    parser.add_argument("--pseudoread-len", default=75, type=int)
    parser.add_argument("--pseudoread-step", default=35, type=int)
    parser.add_argument("--pseudoread-phred", default=30, type=int)
    parser.add_argument("--keep-native", default="false")
    return parser.parse_args()


def parse_strategies(raw: str) -> list[str]:
    selected = [item.strip() for item in raw.split(",") if item.strip()]
    unknown = sorted(set(selected) - BWA_STRATEGIES)
    if unknown:
        raise ValueError(f"Unknown BWA strategy value(s): {', '.join(unknown)}")
    if not selected:
        raise ValueError("At least one BWA strategy is required")
    return selected


def truthy(raw: str | bool) -> bool:
    if isinstance(raw, bool):
        return raw
    return raw.strip().lower() in {"1", "true", "yes", "y", "on"}


def run_checked(cmd: list[str], *, stdout=None) -> None:
    print("Running: " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, stdout=stdout)


def run_bwa_mem_pipeline(bwa_bin: str, samtools_bin: str, target_fasta: Path, fastq: Path, sorted_bam: Path) -> None:
    run_checked([bwa_bin, "index", str(target_fasta)])
    print(
        "Running: "
        + " ".join([bwa_bin, "mem", str(target_fasta), str(fastq)])
        + " | "
        + " ".join([samtools_bin, "view", "-bS", "-"])
        + " | "
        + " ".join([samtools_bin, "sort", "-o", str(sorted_bam), "-"]),
        flush=True,
    )
    bwa_proc = subprocess.Popen([bwa_bin, "mem", str(target_fasta), str(fastq)], stdout=subprocess.PIPE)
    assert bwa_proc.stdout is not None
    view_proc = subprocess.Popen([samtools_bin, "view", "-bS", "-"], stdin=bwa_proc.stdout, stdout=subprocess.PIPE)
    bwa_proc.stdout.close()
    assert view_proc.stdout is not None
    sort_proc = subprocess.Popen([samtools_bin, "sort", "-o", str(sorted_bam), "-"], stdin=view_proc.stdout)
    view_proc.stdout.close()

    sort_code = sort_proc.wait()
    view_code = view_proc.wait()
    bwa_code = bwa_proc.wait()
    failed = [
        (name, code)
        for name, code in [("bwa mem", bwa_code), ("samtools view", view_code), ("samtools sort", sort_code)]
        if code != 0
    ]
    if failed:
        details = ", ".join(f"{name} exit {code}" for name, code in failed)
        raise subprocess.CalledProcessError(failed[0][1], details)
    run_checked([samtools_bin, "index", str(sorted_bam)])


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


def generate_pseudoreads(
    orthologs_fa: Path,
    out_fastq: Path,
    read_len: int,
    step: int,
    phred: int,
) -> int:
    phred_char = chr(phred + 33)
    total_reads = 0
    with out_fastq.open("w") as out:
        for header, seq in iter_fasta(orthologs_fa):
            ortholog_id = header.split()[0]
            read_index = 1
            for start in range(0, max(1, len(seq) - read_len + 1), step):
                read_seq = seq[start : start + read_len]
                if len(read_seq) < 20:
                    continue
                qual = phred_char * len(read_seq)
                end = start + len(read_seq)
                out.write(f"@{ortholog_id}_pseudo_{read_index}_{start + 1}-{end}\n{read_seq}\n+\n{qual}\n")
                read_index += 1
                total_reads += 1
    return total_reads


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


def query_lengths_by_ortholog(orthologs_fasta: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for header, seq in iter_fasta(orthologs_fasta):
        sequence_id = header.split()[0]
        lengths[sequence_id.removeprefix("ortholog_")] = len(seq)
    return lengths


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


def scan_bam(
    bam_path: Path,
    target_seq: str,
) -> tuple[list[dict[str, object]], dict[tuple[str, int, int, str, str], list[dict[str, object]]]]:
    segments: list[dict[str, object]] = []
    event_support: dict[tuple[str, int, int, str, str], list[dict[str, object]]] = defaultdict(list)
    seen_events: set[tuple[tuple[str, int, int, str, str], str, str]] = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.reference_start is None or read.reference_end is None:
                continue
            ortholog_id = read_ortholog_gene_id(read.query_name)
            strand = "-" if read.is_reverse else "+"
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
                    "is_primary": str(not read.is_secondary and not read.is_supplementary).lower(),
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
                        support_key = (key, ortholog_id, read.query_name)
                        if support_key not in seen_events:
                            event_support[key].append(
                                {"ortholog_gene_id": ortholog_id, "strand": strand, "native_record_id": read.query_name}
                            )
                            seen_events.add(support_key)
                    r_pos += length
                    q_pos += length
                elif op == 1:
                    alt = q_seq[q_pos : q_pos + length].upper()
                    key = ("ins", r_pos, r_pos, "", alt)
                    support_key = (key, ortholog_id, read.query_name)
                    if support_key not in seen_events:
                        event_support[key].append(
                            {"ortholog_gene_id": ortholog_id, "strand": strand, "native_record_id": read.query_name}
                        )
                        seen_events.add(support_key)
                    q_pos += length
                elif op in {2, 3}:
                    if op == 2:
                        ref = target_seq[r_pos : r_pos + length].upper()
                        key = ("del", r_pos, r_pos + length, ref, "")
                        support_key = (key, ortholog_id, read.query_name)
                        if support_key not in seen_events:
                            event_support[key].append(
                                {"ortholog_gene_id": ortholog_id, "strand": strand, "native_record_id": read.query_name}
                            )
                            seen_events.add(support_key)
                    r_pos += length
                elif op == 4:
                    q_pos += length
                elif op in {5, 6}:
                    continue
    return segments, event_support


def make_segment_rows(
    base_segments: list[dict[str, object]],
    selected_strategies: list[str],
    gene_id: str,
    target_id: str,
    ortholog_meta_by_id: dict[str, dict[str, str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for strategy in selected_strategies:
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
        "native_record_id": support.get("native_record_id", ""),
        "qc_flags": qc_flags,
    }


def make_bwa_event_rows(
    event_support: dict[tuple[str, int, int, str, str], list[dict[str, object]]],
    gene_id: str,
    target_meta: dict[str, str],
    target_acc: str,
    ortholog_meta_by_id: dict[str, dict[str, str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (event_type, start0, end0, ref, alt), support_rows in sorted(event_support.items()):
        for support in support_rows:
            rows.append(
                make_event_row(
                    gene_id=gene_id,
                    target_meta=target_meta,
                    target_acc=target_acc,
                    ortholog_meta_by_id=ortholog_meta_by_id,
                    strategy=BWA_STRATEGY,
                    event_type=event_type,
                    start0=start0,
                    end0=end0,
                    ref=ref,
                    alt=alt,
                    support=support,
                    qc_flags="bwa_cigar_event",
                )
            )
    return rows


def varscan_command(varscan_jar: str) -> list[str]:
    if varscan_jar == "varscan":
        return ["varscan"]
    jar = Path(varscan_jar)
    if jar.exists():
        return ["java", "-jar", str(jar)]
    if shutil.which("varscan"):
        return ["varscan"]
    return ["java", "-jar", varscan_jar]


def run_varscan(
    samtools_bin: str,
    varscan_jar: str,
    target_fasta: Path,
    filtered_bam: Path,
    work_dir: Path,
    min_coverage: int,
    min_reads2: int,
    min_var_freq: float,
) -> list[Path]:
    pileup = work_dir / "gene.mpileup"
    snps_vcf = work_dir / "gene_snps.vcf"
    indels_vcf = work_dir / "gene_indels.vcf"
    run_checked([samtools_bin, "faidx", str(target_fasta)])
    with pileup.open("w") as out:
        run_checked([samtools_bin, "mpileup", "-f", str(target_fasta), str(filtered_bam)], stdout=out)

    base_cmd = varscan_command(varscan_jar)
    common_args = [
        str(pileup),
        "--min-coverage",
        str(min_coverage),
        "--min-reads2",
        str(min_reads2),
        "--min-var-freq",
        str(min_var_freq),
        "--output-vcf",
        "1",
    ]
    with snps_vcf.open("w") as out:
        run_checked([*base_cmd, "mpileup2snp", *common_args], stdout=out)
    with indels_vcf.open("w") as out:
        run_checked([*base_cmd, "mpileup2indel", *common_args], stdout=out)
    return [snps_vcf, indels_vcf]


def normalize_vcf_variant(pos1: int, ref: str, alt: str) -> tuple[str, int, int, str, str]:
    pos0 = pos1 - 1
    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
        inserted = alt[1:]
        start0 = pos0 + 1
        return "ins", start0, start0, "", inserted
    if len(alt) == 1 and len(ref) > 1 and ref.startswith(alt):
        deleted = ref[1:]
        start0 = pos0 + 1
        return "del", start0, start0 + len(deleted), deleted, ""
    if len(ref) == 1 and len(alt) == 1:
        return "snv", pos0, pos0 + 1, ref, alt
    return "complex", pos0, pos0 + len(ref), ref, alt


def iter_varscan_vcf_events(vcf_paths: list[Path]):
    for vcf_path in vcf_paths:
        with vcf_path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                pos1 = int(parts[1])
                ref = parts[3].upper()
                for alt in parts[4].upper().split(","):
                    yield normalize_vcf_variant(pos1, ref, alt)


def make_varscan_event_rows(
    vcf_paths: list[Path],
    event_support: dict[tuple[str, int, int, str, str], list[dict[str, object]]],
    gene_id: str,
    target_meta: dict[str, str],
    target_acc: str,
    ortholog_meta_by_id: dict[str, dict[str, str]],
) -> tuple[list[dict[str, object]], int]:
    rows: list[dict[str, object]] = []
    unassigned = 0
    for key in iter_varscan_vcf_events(vcf_paths):
        support_rows = event_support.get(key, [])
        if not support_rows:
            unassigned += 1
            continue
        event_type, start0, end0, ref, alt = key
        for support in support_rows:
            rows.append(
                make_event_row(
                    gene_id=gene_id,
                    target_meta=target_meta,
                    target_acc=target_acc,
                    ortholog_meta_by_id=ortholog_meta_by_id,
                    strategy=BWA_VARSCAN_STRATEGY,
                    event_type=event_type,
                    start0=start0,
                    end0=end0,
                    ref=ref,
                    alt=alt,
                    support=support,
                    qc_flags="varscan_call_with_bam_read_support",
                )
            )
    return rows, unassigned


def make_summary_rows(
    selected_strategies: list[str],
    gene_id: str,
    target_length: int,
    ortholog_meta: list[dict[str, str]],
    query_lengths: dict[str, int],
    segment_rows: list[dict[str, object]],
    event_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    segments_by_key: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in segment_rows:
        segments_by_key[(str(row["strategy"]), str(row["ortholog_gene_id"]))].append(row)
    event_counts: dict[tuple[str, str], int] = defaultdict(int)
    for row in event_rows:
        event_counts[(str(row["strategy"]), str(row["ortholog_gene_id"]))] += 1

    summaries: list[dict[str, object]] = []
    for strategy in selected_strategies:
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


def keep_native_outputs(work_dir: Path, outdir: Path) -> None:
    native_dir = outdir / "native"
    native_dir.mkdir(parents=True, exist_ok=True)
    for name in [
        "pseudo_reads.fastq",
        "aln.sorted.bam",
        "aln.sorted.bam.bai",
        "aln.filtered.lis.bam",
        "aln.filtered.lis.bam.bai",
        "bam_filtering_stats.json",
        "bam_filtering_overall.json",
        "gene.mpileup",
        "gene_snps.vcf",
        "gene_indels.vcf",
    ]:
        src = work_dir / name
        if src.exists():
            shutil.copy2(src, native_dir / name)


def main() -> None:
    args = parse_args()
    selected_strategies = parse_strategies(args.strategies)
    keep_native = truthy(args.keep_native)
    task_dir = args.task_dir
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    manifest, target_meta, ortholog_meta = load_task_context(task_dir)
    gene_id = manifest["gene_id"]
    target_id = manifest.get("target_id", f"target_{gene_id}")
    ortholog_meta_by_id = {row["ortholog_gene_id"]: row for row in ortholog_meta}
    target_acc = target_meta.get("genomic_accession", "")

    with tempfile.TemporaryDirectory(prefix="bwa_pseudoreads_", dir=outdir) as tmp_name:
        work_dir = Path(tmp_name)
        local_target_fasta, local_orthologs_fasta = materialize_task_fastas(
            args.source_target_fasta,
            args.source_ortholog_fasta,
            manifest,
            ortholog_meta,
            work_dir,
        )
        target_seq = read_target_sequence(local_target_fasta)
        query_lengths = query_lengths_by_ortholog(local_orthologs_fasta)

        pseudoreads_fastq = work_dir / "pseudo_reads.fastq"
        pseudoread_count = generate_pseudoreads(
            local_orthologs_fasta,
            pseudoreads_fastq,
            read_len=args.pseudoread_len,
            step=args.pseudoread_step,
            phred=args.pseudoread_phred,
        )

        sorted_bam = work_dir / "aln.sorted.bam"
        run_bwa_mem_pipeline(args.bwa_bin, args.samtools_bin, local_target_fasta, pseudoreads_fastq, sorted_bam)

        filter_cfg = {
            "wrong_strand": True,
            "lis": True,
            "overlap": True,
            "min_mapped_pct_of_generated": 0.0,
            "max_pct_filtered": 100.0,
            "min_kept_pct_of_reference": 0.0,
        }
        bam_filtering_v1.filter_bam_for_gene(
            work_dir=work_dir,
            filtering_cfg=filter_cfg,
            read_len=args.pseudoread_len,
            step=args.pseudoread_step,
        )
        filtered_bam = work_dir / "aln.filtered.lis.bam"

        base_segments, event_support = scan_bam(filtered_bam, target_seq)
        segment_rows = make_segment_rows(
            base_segments,
            selected_strategies,
            gene_id,
            target_id,
            ortholog_meta_by_id,
        )

        event_rows: list[dict[str, object]] = []
        varscan_unassigned = 0
        if BWA_STRATEGY in selected_strategies:
            event_rows.extend(make_bwa_event_rows(event_support, gene_id, target_meta, target_acc, ortholog_meta_by_id))
        if BWA_VARSCAN_STRATEGY in selected_strategies:
            vcf_paths = run_varscan(
                args.samtools_bin,
                args.varscan_jar,
                local_target_fasta,
                filtered_bam,
                work_dir,
                min_coverage=args.varscan_min_coverage,
                min_reads2=args.varscan_min_reads2,
                min_var_freq=args.varscan_min_var_freq,
            )
            varscan_rows, varscan_unassigned = make_varscan_event_rows(
                vcf_paths,
                event_support,
                gene_id,
                target_meta,
                target_acc,
                ortholog_meta_by_id,
            )
            event_rows.extend(varscan_rows)

        summary_rows = make_summary_rows(
            selected_strategies,
            gene_id,
            len(target_seq),
            ortholog_meta,
            query_lengths,
            segment_rows,
            event_rows,
        )

        write_tsv_gz(outdir / "alignment_segments.tsv.gz", SEGMENT_FIELDS, segment_rows)
        write_tsv_gz(outdir / "alignment_events.tsv.gz", EVENT_FIELDS, event_rows)
        write_tsv_gz(outdir / "ortholog_alignment_summary.tsv.gz", SUMMARY_FIELDS, summary_rows)
        write_tsv_gz(outdir / "failures.tsv.gz", FAILURE_FIELDS, [])

        if keep_native:
            keep_native_outputs(work_dir, outdir)

        manifest_out = {
            "gene_id": gene_id,
            "strategy": "bwa_pseudoreads",
            "strategies": selected_strategies,
            "tool": "bwa",
            "segment_count": len(segment_rows),
            "event_count": len(event_rows),
            "ortholog_count": len(ortholog_meta),
            "pseudoread_count": pseudoread_count,
            "keep_native": keep_native,
            "varscan_unassigned_event_count": varscan_unassigned,
            "pseudoread_len": args.pseudoread_len,
            "pseudoread_step": args.pseudoread_step,
            "pseudoread_phred": args.pseudoread_phred,
            "varscan_min_coverage": args.varscan_min_coverage,
            "varscan_min_reads2": args.varscan_min_reads2,
            "varscan_min_var_freq": args.varscan_min_var_freq,
        }
        (outdir / "manifest.json").write_text(json.dumps(manifest_out, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
