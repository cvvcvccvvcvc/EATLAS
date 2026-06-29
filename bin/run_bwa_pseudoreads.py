#!/usr/bin/env python3
"""Run BWA pseudoreads alignment and variant extraction (pysam and varscan)."""
import argparse
import gzip
import json
import subprocess
import csv
from pathlib import Path
from collections import defaultdict
import pysam
import re

import bam_filtering_v1


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
    parser.add_argument("--bwa-bin", default="bwa", type=str)
    parser.add_argument("--samtools-bin", default="samtools", type=str)
    parser.add_argument("--varscan-jar", default="VarScan.v2.4.6.jar", type=str)
    return parser.parse_args()


def run_cmd(cmd: str):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def iter_fasta(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    header = None
    seq_parts = []
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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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



def generate_pseudoreads(orthologs_fa: Path, out_fastq: Path, read_len=75, step=35, phred=30):
    phred_char = chr(phred + 33)
    total_reads = 0
    with open(out_fastq, "w") as out:
        for header, seq in iter_fasta(orthologs_fa):
            ortholog_id = header.split()[0]
            n = len(seq)
            read_index = 1
            for start in range(0, max(1, n - read_len + 1), step):
                read_seq = seq[start : start + read_len]
                if len(read_seq) < 20:
                    continue
                qual = phred_char * len(read_seq)
                read_header = f"@{ortholog_id}_pseudo_{read_index}_{start+1}-{start+len(read_seq)}"
                out.write(f"{read_header}\n{read_seq}\n+\n{qual}\n")
                read_index += 1
                total_reads += 1
    return total_reads


def extract_pysam_variants(bam_path: Path, target_seq: str, target_acc: str, gene_id: str):
    """Extract variants from BAM using pysam."""
    variants = set() # (ortholog_id, ref_pos, ref_base, alt_base)
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    for read in bam.fetch():
        if read.is_unmapped:
            continue
        # read_name format: ortholog_<id>_pseudo_...
        # We need to extract the ortholog id
        match = re.match(r"^(ortholog_\d+)_pseudo_", read.query_name)
        if match:
            ortholog_id = match.group(1).replace("ortholog_", "")
        else:
            continue

        pairs = read.get_aligned_pairs(with_seq=True)
        for q_pos, r_pos, ref_base in pairs:
            if q_pos is not None and r_pos is not None:
                query_base = read.query_sequence[q_pos]
                if ref_base and query_base.upper() != ref_base.upper():
                    # SNV
                    variants.add((ortholog_id, r_pos, ref_base.upper(), query_base.upper(), 'snv'))
            elif q_pos is not None and r_pos is None:
                # Insertion (simplified)
                # To get exact ins sequence, we'd need to reconstruct it from consecutive q_pos
                # For basic counting/comparison, we just mark an insertion at the LAST r_pos
                pass # Simplified parsing might miss multi-base indels properly. Let's use pileup or just basic CIGAR.

    # To get proper variants matching standard gaph_v2 output, we should parse the CIGAR properly
    # Let's write a better CIGAR parser for events
    events = set()
    
    # We will rewind and parse CIGARs properly
    bam.reset()
    for read in bam.fetch():
        if read.is_unmapped:
            continue
        match = re.match(r"^(ortholog_\d+)_pseudo_", read.query_name)
        if not match:
            continue
        ortholog_id = match.group(1).replace("ortholog_", "")
        
        q_seq = read.query_sequence
        r_pos = read.reference_start
        q_pos = 0
        
        for op, length in read.cigartuples:
            if op == 0: # M
                for i in range(length):
                    rb = target_seq[r_pos + i].upper()
                    qb = q_seq[q_pos + i].upper()
                    if rb != qb:
                        events.add((ortholog_id, r_pos + i, r_pos + i + 1, rb, qb, 'snv'))
                r_pos += length
                q_pos += length
            elif op == 1: # I (Insertion into reference)
                qb = q_seq[q_pos : q_pos + length].upper()
                events.add((ortholog_id, r_pos, r_pos, "", qb, 'ins'))
                q_pos += length
            elif op == 2: # D (Deletion from reference)
                rb = target_seq[r_pos : r_pos + length].upper()
                events.add((ortholog_id, r_pos, r_pos + length, rb, "", 'del'))
                r_pos += length
            elif op == 4: # S
                q_pos += length
            elif op == 5: # H
                pass

    return events


def write_tsv_gz(path: Path, headers: list[str], rows: list[dict]):
    with gzip.open(path, "wt", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    task_dir = args.task_dir
    manifest_path = task_dir / "task.json"
    with open(manifest_path) as f:
        manifest = json.load(f)

    gene_id = manifest["gene_id"]
    target_fasta = task_dir / manifest["target_fasta"]
    orthologs_fasta = task_dir / manifest["ortholog_fasta"]

    target_meta = read_tsv(task_dir / "target.metadata.tsv")[0]
    ortholog_meta = read_tsv(task_dir / "orthologs.metadata.tsv")
    ortholog_meta_by_id = {row["ortholog_gene_id"]: row for row in ortholog_meta}

    # Target info
    target_acc = target_meta.get("genomic_accession", "unknown")
    target_seq = ""
    for _, seq in iter_fasta(target_fasta):
        target_seq = seq
        break
    query_lengths = {}
    for header, seq in iter_fasta(orthologs_fasta):
        sequence_id = header.split()[0]
        ortholog_id = sequence_id.replace("ortholog_", "")
        query_lengths[ortholog_id] = len(seq)

    def make_event_row(ortholog_id, ref_start, ref_end, ref, alt, event_type, strategy):
        meta = ortholog_meta_by_id.get(ortholog_id, {})
        genomic_start, genomic_end = genomic_coords(target_meta, ref_start, ref_end)
        return {
            "gene_id": gene_id,
            "ortholog_gene_id": ortholog_id,
            "tax_id": meta.get("tax_id", ""),
            "taxname": meta.get("taxname", ""),
            "strategy": strategy,
            "tool": "bwa",
            "preset": "pseudo",
            "event_id": f"{ortholog_id}_{ref_start}_{ref}_{alt}",
            "event_type": event_type,
            "target_start0": ref_start,
            "target_end0": ref_end,
            "genomic_accession": target_acc,
            "genomic_start1": genomic_start,
            "genomic_end1": genomic_end,
            "ref": ref,
            "alt": alt,
            "query_id": f"ortholog_{ortholog_id}",
            "strand": "+",  # BWA maps forward pseudo reads directly to genomic target strand
            "native_record_id": "",
            "qc_flags": ""
        }

    pysam_var_dict = defaultdict(list)
    event_rows = []

    # 1. Generate Pseudoreads
    pseudoreads_fastq = task_dir / "pseudo_reads.fastq"
    generate_pseudoreads(orthologs_fasta, pseudoreads_fastq)

    # 2. BWA Alignment
    run_cmd(f"{args.bwa_bin} index {target_fasta}")
    aln_sam = task_dir / "aln.sam"
    aln_bam = task_dir / "aln.bam"
    aln_sorted_bam = task_dir / "aln.sorted.bam"
    run_cmd(f"{args.bwa_bin} mem {target_fasta} {pseudoreads_fastq} > {aln_sam}")
    run_cmd(f"{args.samtools_bin} view -bS {aln_sam} > {aln_bam}")
    run_cmd(f"{args.samtools_bin} sort {aln_bam} -o {aln_sorted_bam}")
    run_cmd(f"{args.samtools_bin} index {aln_sorted_bam}")

    # 3. Filter BAM using bam_filtering_v1
    filter_cfg = {
        "wrong_strand": True,
        "lis": True,
        "overlap": True,
        "min_mapped_pct_of_generated": 0.0,
        "max_pct_filtered": 100.0,
        "min_kept_pct_of_reference": 0.0
    }
    bam_filtering_v1.filter_bam_for_gene(
        work_dir=task_dir,
        filtering_cfg=filter_cfg
    )
    filtered_bam = task_dir / "aln.filtered.lis.bam"

    # 4. Extract pysam events
    pysam_events = extract_pysam_variants(filtered_bam, target_seq, target_acc, gene_id)

    # We will format events for output
    event_rows = []
    
    # helper for rows
    pysam_var_dict = defaultdict(list)
    for (ortholog_id, ref_start, ref_end, ref, alt, event_type) in pysam_events:
        event_rows.append(make_event_row(ortholog_id, ref_start, ref_end, ref, alt, event_type, "bwa_pseudoreads_pysam"))
        # Store for varscan intersection
        pysam_var_dict[(ref_start, ref, alt)].append((ortholog_id, ref_end, event_type))

    # 5. VarScan Calling
    pileup = task_dir / "gene.mpileup"
    snps_vcf = task_dir / "gene_snps.vcf"
    indels_vcf = task_dir / "gene_indels.vcf"
    
    run_cmd(f"{args.samtools_bin} faidx {target_fasta}")
    run_cmd(f"{args.samtools_bin} mpileup -f {target_fasta} {filtered_bam} > {pileup}")
    
    has_varscan = False
    
    import shutil
    java_bin = shutil.which("java")
    if not java_bin and Path("/opt/homebrew/opt/openjdk/bin/java").exists():
        java_bin = "/opt/homebrew/opt/openjdk/bin/java"
    if not java_bin:
        java_bin = "java"

    varscan_cmd = "varscan" if args.varscan_jar == "varscan" else f"{java_bin} -jar {args.varscan_jar}"
    
    try:
        # First check if varscan is in path if it's set to "varscan"
        if args.varscan_jar != "varscan" and not Path(args.varscan_jar).exists():
            # If jar doesn't exist, try just 'varscan'
            varscan_cmd = "varscan"
            
        run_cmd(f"{varscan_cmd} mpileup2snp {pileup} --min-coverage 2 --min-reads2 1 --min-var-freq 0.01 --output-vcf 1 > {snps_vcf}")
        run_cmd(f"{varscan_cmd} mpileup2indel {pileup} --min-coverage 2 --min-reads2 1 --min-var-freq 0.01 --output-vcf 1 > {indels_vcf}")
        has_varscan = True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("VarScan failed or not found. Skipping varscan strategy.")

    if has_varscan:
        # Parse VarScan VCFs
        for vcf in [snps_vcf, indels_vcf]:
            with open(vcf) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 5:
                        continue
                    # 1-based pos
                    pos = int(parts[1]) - 1
                    ref = parts[3]
                    alt = parts[4]
                    
                    # VarScan output format mapping
                    # VarScan indels include the preceding reference base (e.g. ref=C alt=CA)
                    # Our minimap2 format (and pysam parser) uses empty string for insertion ref,
                    # and empty string for deletion alt.
                    v_ref, v_alt = ref, alt
                    v_pos = pos
                    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
                        # Insertion. Minimap2 format: target_start0=pos+1, target_end0=pos+1, ref="", alt=ins_seq
                        v_ref = ""
                        v_alt = alt[1:]
                        v_pos = pos + 1
                    elif len(alt) == 1 and len(ref) > 1 and ref.startswith(alt):
                        # Deletion. Minimap2 format: target_start0=pos+1, target_end0=pos+len, ref=del_seq, alt=""
                        v_ref = ref[1:]
                        v_alt = ""
                        v_pos = pos + 1
                    
                    # Look up in pysam dict
                    key = (v_pos, v_ref, v_alt)
                    if key in pysam_var_dict:
                        for (ortholog_id, ref_end, event_type) in pysam_var_dict[key]:
                            event_rows.append(make_event_row(ortholog_id, v_pos, ref_end, v_ref, v_alt, event_type, "bwa_pseudoreads_varscan"))

    # Write events
    write_tsv_gz(task_dir / "alignment_events.tsv.gz", EVENT_FIELDS, event_rows)
    
    # Write empty segments and summaries for now, so MERGE_ALIGNMENT doesn't crash
    write_tsv_gz(task_dir / "alignment_segments.tsv.gz", SEGMENT_FIELDS, [])
    
    write_tsv_gz(task_dir / "failures.tsv.gz", FAILURE_FIELDS, [])
    
    # Dummy summary for orthologs found
    summaries = []
    unique_orthologs = set(ortholog_id for (ortholog_id, *_) in pysam_events)
    for oid in unique_orthologs:
        for strat in ["bwa_pseudoreads_pysam", "bwa_pseudoreads_varscan"]:
            if not has_varscan and strat == "bwa_pseudoreads_varscan":
                continue
            summaries.append({
                "gene_id": gene_id,
                "ortholog_gene_id": oid,
                "tax_id": "",
                "taxname": "",
                "strategy": strat,
                "tool": "bwa",
                "preset": "pseudo",
                "status": "aligned",
                "target_length": len(target_seq),
                "query_length": query_lengths.get(oid, 0),
                "segment_count": 0,
                "primary_segment_count": 0,
                "secondary_segment_count": 0,
                "aligned_target_bp": 0,
                "aligned_query_bp": 0,
                "target_coverage": "0.000000",
                "query_coverage": "0.000000",
                "best_identity": "0.000000",
                "mean_identity": "0.000000",
                "event_count": sum(1 for e in event_rows if e["ortholog_gene_id"] == oid and e["strategy"] == strat),
                "qc_flags": "no_segments_bwa_pseudoreads"
            })
    write_tsv_gz(task_dir / "ortholog_alignment_summary.tsv.gz", SUMMARY_FIELDS, summaries)
    
    # Overwrite task.json with updated strategy
    manifest["strategy"] = "bwa_pseudoreads"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    manifest_out = {
        "gene_id": gene_id,
        "strategy": "bwa_pseudoreads",
        "strategies": ["bwa_pseudoreads_pysam"] + (["bwa_pseudoreads_varscan"] if has_varscan else []),
        "tool": "bwa",
        "segment_count": 0,
        "event_count": len(event_rows),
        "ortholog_count": len(ortholog_meta),
        "keep_native": False,
    }
    (task_dir / "manifest.json").write_text(json.dumps(manifest_out, indent=2) + "\n")

if __name__ == "__main__":
    main()
