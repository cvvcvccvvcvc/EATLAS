"""Production BAM filtering for the main pipeline (no plotting/CLI helpers).

Filtering order per homologue:
1) Keep dominant strand only
2) Keep monotonic LIS/LDS backbone in alignment order
3) Deduplicate overlapping alignments at the same position
4) Optional whole-homologue drops (mapped%, pct filtered, kept/reference%)
"""

from __future__ import annotations

import json
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pysam
import logging

class LoguruDummyLogger:
    def __init__(self):
        self.logger = logging.getLogger("bam_filtering")
    def info(self, msg, *args, **kwargs):
        if args or kwargs:
            self.logger.info(msg.format(*args, **kwargs))
        else:
            self.logger.info(msg)
    def warning(self, msg, *args, **kwargs):
        if args or kwargs:
            self.logger.warning(msg.format(*args, **kwargs))
        else:
            self.logger.warning(msg)
    def success(self, msg, *args, **kwargs):
        self.info(msg, *args, **kwargs)

logger = LoguruDummyLogger()



@dataclass
class FilterResult:
    """Output metadata for one gene-level BAM filtering run."""

    input_bam: Path
    output_bam: Path
    output_bai: Path
    per_homologue_stats_json: Path
    overall_stats_json: Path
    filtering_stats: Dict[str, Dict[str, Any]]
    overall: Dict[str, Any]


def _read_key(read: pysam.AlignedSegment) -> Tuple[Any, ...]:
    """Stable identity for one BAM alignment record."""
    return (
        read.query_name,
        read.flag,
        read.reference_id,
        read.reference_start,
        read.reference_end,
        read.cigarstring,
        read.mapping_quality,
    )


def parse_read_name(read_name: str) -> Tuple[str, int, int, int]:
    """
    Parse read name like:
    gene_100439533_pseudo_2384_83406-83480
    """
    parts = read_name.split("_pseudo_")
    if len(parts) != 2:
        raise ValueError(f"Unexpected read name format: {read_name}")

    homologue_id = parts[0]
    read_num_str, pos_range = parts[1].split("_", 1)
    start_str, end_str = pos_range.split("-")

    return homologue_id, int(read_num_str), int(start_str), int(end_str)


def _build_read_record(
    read: pysam.AlignedSegment,
    parsed: Tuple[str, int, int, int],
) -> Dict[str, Any]:
    """Convert one alignment row to compact dict used by filters."""
    return {
        "read_name": read.query_name,
        "read_key": _read_key(read),
        "homologue_id": parsed[0],
        "actual_read_num": parsed[1],
        "original_start": parsed[2],
        "original_end": parsed[3],
        "alignment_pos": read.reference_start,
        "is_reverse": read.is_reverse,
    }


def collect_homologue_data(
    bam_path: Path,
) -> Tuple[Dict[str, List[Dict[str, Any]]], Dict[str, Dict[str, Any]]]:
    """Collect per-homologue reads and strand stats in a single BAM pass."""
    stats = defaultdict(
        lambda: {
            "total_reads": 0,
            "forward_reads": 0,
            "reverse_reads": 0,
        }
    )
    reads_by_homologue: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            parsed = parse_read_name(read.query_name)
            homologue_id = parsed[0]

            stats[homologue_id]["total_reads"] += 1
            if read.is_reverse:
                stats[homologue_id]["reverse_reads"] += 1
            else:
                stats[homologue_id]["forward_reads"] += 1

            reads_by_homologue[homologue_id].append(_build_read_record(read, parsed))

    for homologue_id, s in stats.items():
        s["dominant_strand"] = "forward" if s["forward_reads"] >= s["reverse_reads"] else "reverse"
        reads_by_homologue[homologue_id].sort(key=lambda r: r["alignment_pos"])

    return dict(reads_by_homologue), dict(stats)


def _assign_observed_numbers(reads: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Assign observed numbers based on sorted alignment positions."""
    current_pos = None
    current_num = 0
    for read in reads:
        if read["alignment_pos"] != current_pos:
            current_num += 1
            current_pos = read["alignment_pos"]
        read["observed_num"] = current_num
    return reads


def _assign_aligned_read_numbers(reads: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Re-enumerate aligned reads by original pseudo-read order."""
    sorted_by_actual = sorted(reads, key=lambda r: r["actual_read_num"])
    current_actual = None
    current_aligned = 0
    actual_to_aligned: Dict[int, int] = {}

    for read in sorted_by_actual:
        if read["actual_read_num"] != current_actual:
            current_aligned += 1
            current_actual = read["actual_read_num"]
            actual_to_aligned[current_actual] = current_aligned

    for read in reads:
        read["aligned_read_num"] = actual_to_aligned[read["actual_read_num"]]

    return reads


def _lis_indices(values: List[int]) -> Set[int]:
    """Indices of one longest strictly increasing subsequence."""
    n = len(values)
    if n == 0:
        return set()

    tails_values: List[int] = []
    tails_indices: List[int] = []
    prev = [-1] * n

    for i, value in enumerate(values):
        pos = bisect_left(tails_values, value)
        if pos == len(tails_values):
            tails_values.append(value)
            tails_indices.append(i)
        else:
            tails_values[pos] = value
            tails_indices[pos] = i

        if pos > 0:
            prev[i] = tails_indices[pos - 1]

    out: Set[int] = set()
    cur = tails_indices[-1]
    while cur != -1:
        out.add(cur)
        cur = prev[cur]

    return out


def _backbone_indices(reads: List[Dict[str, Any]], dominant_strand: str) -> Set[int]:
    """
    Longest monotonic backbone in alignment order.

    - forward: LIS on aligned_read_num
    - reverse: LDS via LIS on negated aligned_read_num
    """
    sequence = [r["aligned_read_num"] for r in reads]
    if dominant_strand == "reverse":
        sequence = [-x for x in sequence]
    return _lis_indices(sequence)


def _filter_homologue_reads(
    reads: List[Dict[str, Any]],
    dominant_strand: str,
    apply_wrong_strand: bool = True,
    apply_lis: bool = True,
    apply_overlap: bool = True,
) -> Tuple[Set[Tuple[Any, ...]], Dict[str, Any]]:
    """Apply dominant strand + LIS/LDS + overlap filtering to one homologue."""
    if not reads:
        return set(), {}

    initial_count = len(reads)
    forward_initial = sum(1 for r in reads if not r["is_reverse"])
    reverse_initial = sum(1 for r in reads if r["is_reverse"])

    reads = _assign_observed_numbers(reads)
    reads = _assign_aligned_read_numbers(reads)

    if apply_wrong_strand:
        if dominant_strand == "forward":
            reads_after_strand = [r for r in reads if not r["is_reverse"]]
        else:
            reads_after_strand = [r for r in reads if r["is_reverse"]]
    else:
        reads_after_strand = list(reads)

    if apply_lis:
        keep_idx = _backbone_indices(reads_after_strand, dominant_strand)
        reads_after_order = [r for i, r in enumerate(reads_after_strand) if i in keep_idx]
    else:
        keep_idx = set(range(len(reads_after_strand)))
        reads_after_order = list(reads_after_strand)

    if apply_overlap:
        position_groups: Dict[int, List[Dict[str, Any]]] = defaultdict(list)
        for read in reads_after_order:
            position_groups[read["alignment_pos"]].append(read)

        reads_final: List[Dict[str, Any]] = []
        for reads_at_pos in position_groups.values():
            if len(reads_at_pos) == 1:
                reads_final.append(reads_at_pos[0])
            else:
                reads_final.append(min(reads_at_pos, key=lambda r: r["actual_read_num"]))
    else:
        reads_final = list(reads_after_order)

    keep_keys = {r["read_key"] for r in reads_final}

    stats: Dict[str, Any] = {
        "initial_count": initial_count,
        "forward_initial": forward_initial,
        "reverse_initial": reverse_initial,
        "dominant_strand": dominant_strand,
        "order_method": "lis" if apply_lis else "disabled",
        "wrong_strand_enabled": apply_wrong_strand,
        "lis_enabled": apply_lis,
        "overlap_enabled": apply_overlap,
        "after_strand_filter": len(reads_after_strand),
        "after_order_filter": len(reads_after_order),
        "after_overlap_filter": len(reads_final),
        "filtered_by_strand": initial_count - len(reads_after_strand),
        "filtered_by_order": len(reads_after_strand) - len(reads_after_order),
        "filtered_by_overlap": len(reads_after_order) - len(reads_final),
        "lis_backbone_size": len(keep_idx),
        "total_filtered": initial_count - len(reads_final),
        "pct_kept": 100.0 * len(reads_final) / initial_count if initial_count > 0 else 0.0,
    }

    return keep_keys, stats


def expected_pseudoreads(seq_len: int, read_len: int, step: int) -> int:
    """Count sliding-window pseudoreads for one sequence."""
    if seq_len < read_len:
        return 0
    return 1 + (seq_len - read_len) // step


def _generated_from_genes_fastq(path: Path, read_len: int, step: int) -> Dict[str, int]:
    """Generate per-homologue expected pseudo-read counts from source gene FASTQ."""
    counts: Dict[str, int] = {}
    with path.open() as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().rstrip("\n")
            plus = f.readline()
            qual = f.readline()
            if not plus or not qual:
                break
            rid = header[1:].strip().split()[0] if header.startswith("@") else header.strip().split()[0]
            counts[rid] = expected_pseudoreads(len(seq), read_len=read_len, step=step)
    return counts


def _generated_from_pseudoreads_fastq(path: Path) -> Dict[str, int]:
    """Generate per-homologue pseudo-read counts by counting pseudo read headers."""
    counts: Dict[str, int] = {}
    with path.open() as f:
        line_no = 0
        for line in f:
            if line_no % 4 == 0:
                rid = line[1:].strip().split()[0] if line.startswith("@") else line.strip().split()[0]
                hid = rid.split("_pseudo_")[0]
                counts[hid] = counts.get(hid, 0) + 1
            line_no += 1
    return counts


def load_generated_pseudoreads(
    work_dir: Path,
    read_len: int,
    step: int,
) -> Tuple[Dict[str, int], Optional[str]]:
    """
    Load generated pseudo-read counts from run directory.

    Priority:
    1) genes_from_coordinates.fastq (small, deterministic)
    2) pseudo_reads.fastq
    """
    genes_fastq = work_dir / "genes_from_coordinates.fastq"
    pseudo_fastq = work_dir / "pseudo_reads.fastq"

    if genes_fastq.exists():
        return _generated_from_genes_fastq(genes_fastq, read_len=read_len, step=step), str(genes_fastq)
    if pseudo_fastq.exists():
        return _generated_from_pseudoreads_fastq(pseudo_fastq), str(pseudo_fastq)
    return {}, None


def _first_fasta_seq_len(path: Path) -> Optional[int]:
    """Length of first FASTA sequence, or None if no valid sequence exists."""
    seq_len = 0
    in_seq = False
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if in_seq:
                    break
                in_seq = True
                continue
            if in_seq:
                seq_len += len(line)
    return seq_len if in_seq else None


def load_reference_pseudoreads(
    work_dir: Path,
    bam_path: Path,
    read_len: int,
    step: int,
) -> Tuple[Optional[int], Optional[int], Optional[str]]:
    """
    Load reference pseudo-read denominator N.

    Priority:
    1) gene_seq.fasta
    2) BAM header reference length
    """
    ref_fasta = work_dir / "gene_seq.fasta"
    if ref_fasta.exists():
        ref_len = _first_fasta_seq_len(ref_fasta)
        if ref_len is not None:
            return expected_pseudoreads(ref_len, read_len=read_len, step=step), ref_len, str(ref_fasta)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if bam.lengths:
            ref_len = int(bam.lengths[0])
            return expected_pseudoreads(ref_len, read_len=read_len, step=step), ref_len, "bam_header"

    return None, None, None


def filter_bam_for_gene(
    work_dir: Path,
    filtering_cfg: Dict[str, Any],
    read_len: int = 75,
    step: int = 35,
    verbose: bool = False,
) -> FilterResult:
    """
    Apply LIS-based BAM filtering in production pipeline.

    Reads:
      - <work_dir>/aln.sorted.bam
    Writes:
      - <work_dir>/aln.filtered.lis.bam (+ .bai)
      - <work_dir>/bam_filtering_stats.json
      - <work_dir>/bam_filtering_overall.json
    """
    work_dir = Path(work_dir)
    input_bam = work_dir / "aln.sorted.bam"
    output_bam = work_dir / "aln.filtered.lis.bam"
    per_homologue_stats_json = work_dir / "bam_filtering_stats.json"
    overall_stats_json = work_dir / "bam_filtering_overall.json"

    if not input_bam.exists():
        raise FileNotFoundError(f"Input BAM not found: {input_bam}")

    wrong_strand = bool(filtering_cfg.get("wrong_strand", True))
    lis = bool(filtering_cfg.get("lis", True))
    overlap = bool(filtering_cfg.get("overlap", True))
    min_mapped_pct = filtering_cfg.get("min_mapped_pct_of_generated")
    max_pct_filtered = filtering_cfg.get("max_pct_filtered")
    min_kept_pct_of_reference = filtering_cfg.get("min_kept_pct_of_reference")

    logger.info(
        "BAM filtering started: input={}, wrong_strand={}, lis={}, overlap={}, "
        "min_mapped_pct_of_generated={}, max_pct_filtered={}, min_kept_pct_of_reference={}, "
        "read_len={}, step={}",
        input_bam,
        wrong_strand,
        lis,
        overlap,
        min_mapped_pct,
        max_pct_filtered,
        min_kept_pct_of_reference,
        read_len,
        step,
    )

    generated_counts, generated_source = load_generated_pseudoreads(work_dir, read_len=read_len, step=step)
    if generated_source:
        logger.info(
            "Generated pseudo-read counts loaded: source={} homologues={}",
            generated_source,
            len(generated_counts),
        )
    else:
        logger.warning(
            "Generated pseudo-read source not found (genes_from_coordinates.fastq / pseudo_reads.fastq)"
        )

    reference_pseudoreads, reference_seq_len, reference_source = load_reference_pseudoreads(
        work_dir, input_bam, read_len=read_len, step=step
    )
    logger.info(
        "Reference pseudo-read denominator: source={} reference_len={} N={}",
        reference_source,
        reference_seq_len,
        reference_pseudoreads,
    )

    reads_by_homologue, homologue_stats = collect_homologue_data(input_bam)
    if verbose:
        logger.info("Collected homologue data: {} homologues", len(homologue_stats))

    all_keep_keys: Set[Tuple[Any, ...]] = set()
    filtering_stats: Dict[str, Dict[str, Any]] = {}
    missing_generated_for_mapped_filter: List[str] = []

    for homologue_id in sorted(homologue_stats.keys()):
        reads = reads_by_homologue.get(homologue_id, [])
        dominant_strand = homologue_stats[homologue_id]["dominant_strand"]
        keep_keys, stats = _filter_homologue_reads(
            reads,
            dominant_strand,
            apply_wrong_strand=wrong_strand,
            apply_lis=lis,
            apply_overlap=overlap,
        )

        generated = generated_counts.get(homologue_id)
        mapped_pct = None
        if generated is not None and generated > 0:
            mapped_pct = 100.0 * stats["initial_count"] / generated
        if min_mapped_pct is not None and generated is None:
            missing_generated_for_mapped_filter.append(homologue_id)

        kept_pct_of_reference = None
        if reference_pseudoreads and reference_pseudoreads > 0:
            kept_pct_of_reference = 100.0 * stats["after_overlap_filter"] / reference_pseudoreads

        stats["generated_pseudoreads"] = generated
        stats["mapped_pct_of_generated"] = mapped_pct
        stats["reference_pseudoreads"] = reference_pseudoreads
        stats["base_kept_pct_of_reference"] = kept_pct_of_reference
        stats["kept_pct_of_reference"] = kept_pct_of_reference

        stats["base_kept_count"] = stats["after_overlap_filter"]
        stats["base_total_filtered"] = stats["total_filtered"]
        stats["base_pct_filtered"] = (
            100.0 * stats["base_total_filtered"] / stats["initial_count"] if stats["initial_count"] > 0 else 0.0
        )

        drop_reason = None
        if (
            min_mapped_pct is not None
            and mapped_pct is not None
            and mapped_pct < float(min_mapped_pct)
        ):
            drop_reason = "low_mapped_pct_of_generated"

        if (
            drop_reason is None
            and max_pct_filtered is not None
            and stats["base_pct_filtered"] > float(max_pct_filtered)
        ):
            drop_reason = "high_pct_filtered"

        if (
            drop_reason is None
            and min_kept_pct_of_reference is not None
            and kept_pct_of_reference is not None
            and kept_pct_of_reference < float(min_kept_pct_of_reference)
        ):
            drop_reason = "low_kept_pct_of_reference"

        stats["dropped_homologue"] = False
        stats["drop_reason"] = None
        stats["filtered_by_homologue_drop"] = 0
        stats["filtered_by_low_mapped_pct"] = 0
        stats["filtered_by_high_pct_filtered"] = 0
        stats["filtered_by_low_kept_pct_of_reference"] = 0

        if drop_reason is not None:
            additionally_removed = len(keep_keys)
            keep_keys = set()
            stats["dropped_homologue"] = True
            stats["drop_reason"] = drop_reason
            stats["filtered_by_homologue_drop"] = additionally_removed

            if drop_reason == "low_mapped_pct_of_generated":
                stats["filtered_by_low_mapped_pct"] = additionally_removed
            elif drop_reason == "high_pct_filtered":
                stats["filtered_by_high_pct_filtered"] = additionally_removed
            elif drop_reason == "low_kept_pct_of_reference":
                stats["filtered_by_low_kept_pct_of_reference"] = additionally_removed

            stats["after_overlap_filter"] = 0
            stats["total_filtered"] = stats["initial_count"]
            stats["pct_kept"] = 0.0
            stats["kept_pct_of_reference"] = 0.0 if reference_pseudoreads else None

        stats["pct_filtered"] = 100.0 - stats["pct_kept"]

        all_keep_keys.update(keep_keys)
        filtering_stats[homologue_id] = stats

    if min_mapped_pct is not None and missing_generated_for_mapped_filter:
        missing_generated_for_mapped_filter = sorted(missing_generated_for_mapped_filter)
        sample = ", ".join(missing_generated_for_mapped_filter[:10])
        extra = len(missing_generated_for_mapped_filter) - 10
        if extra > 0:
            sample = f"{sample}, ... (+{extra} more)"
        logger.error(
            "Missing generated pseudo-read counts for {} homologues while "
            "min_mapped_pct_of_generated is enabled. Examples: {}",
            len(missing_generated_for_mapped_filter),
            sample,
        )
        raise RuntimeError(
            "Generated pseudo-read counts are required for all homologues when "
            "bam_filtering.min_mapped_pct_of_generated is enabled"
        )

    with pysam.AlignmentFile(str(input_bam), "rb") as inbam:
        with pysam.AlignmentFile(str(output_bam), "wb", header=inbam.header) as outbam:
            for read in inbam.fetch():
                if _read_key(read) in all_keep_keys:
                    outbam.write(read)

    pysam.index(str(output_bam))
    output_bai = Path(str(output_bam) + ".bai")

    total_initial = sum(s["initial_count"] for s in filtering_stats.values())
    total_kept = sum(s["after_overlap_filter"] for s in filtering_stats.values())
    total_filtered = total_initial - total_kept
    dropped_homologues = sum(1 for s in filtering_stats.values() if s.get("dropped_homologue"))

    overall: Dict[str, Any] = {
        "input_bam": str(input_bam),
        "output_bam": str(output_bam),
        "per_homologue_stats_json": str(per_homologue_stats_json),
        "order_method": "lis" if lis else "disabled",
        "pseudo_read_len": read_len,
        "pseudo_read_step": step,
        "wrong_strand": wrong_strand,
        "lis": lis,
        "overlap": overlap,
        "min_mapped_pct_of_generated": min_mapped_pct,
        "max_pct_filtered": max_pct_filtered,
        "min_kept_pct_of_reference": min_kept_pct_of_reference,
        "homologue_count": len(filtering_stats),
        "dropped_homologues": dropped_homologues,
        "total_initial_reads": total_initial,
        "total_kept_reads": total_kept,
        "total_filtered_reads": total_filtered,
        "pct_kept": (100.0 * total_kept / total_initial) if total_initial else 0.0,
        "pct_filtered": (100.0 * total_filtered / total_initial) if total_initial else 0.0,
        "generated_counts_source": generated_source,
        "homologues_with_generated_counts": len(generated_counts),
        "reference_length_source": reference_source,
        "reference_sequence_length": reference_seq_len,
        "reference_pseudoreads": reference_pseudoreads,
    }
    if generated_counts:
        total_generated_known = sum(generated_counts.values())
        overall["total_generated_pseudoreads_known"] = total_generated_known
        overall["mapped_pct_of_generated_total"] = (
            100.0 * total_initial / total_generated_known if total_generated_known else 0.0
        )

    with per_homologue_stats_json.open("w") as f:
        json.dump(filtering_stats, f, indent=2)
    with overall_stats_json.open("w") as f:
        json.dump(overall, f, indent=2)

    logger.info(
        "BAM filtering complete: output={} kept={}/{} ({:.2f}%) dropped_homologues={}",
        output_bam,
        total_kept,
        total_initial,
        overall["pct_kept"],
        dropped_homologues,
    )

    return FilterResult(
        input_bam=input_bam,
        output_bam=output_bam,
        output_bai=output_bai,
        per_homologue_stats_json=per_homologue_stats_json,
        overall_stats_json=overall_stats_json,
        filtering_stats=filtering_stats,
        overall=overall,
    )
