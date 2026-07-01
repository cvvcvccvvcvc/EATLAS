#!/usr/bin/env python3
"""Annotate TSV variants with ClinVar and gnomAD."""

import argparse
import bisect
import csv
import gzip
import logging
import sys
import concurrent.futures
from collections import Counter, defaultdict
from pathlib import Path

# Add bin to path so we can import fetch_gnomad_variants
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from fetch_gnomad_variants import fetch_region_variants_recursive, _select_af_metrics
except ImportError as e:
    print(f"Error importing fetch_gnomad_variants: {e}")
    sys.exit(1)

try:
    import pysam
except ImportError:
    print("pysam is required but not installed.")
    sys.exit(1)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--events-tsv", required=True, type=Path)
    parser.add_argument("--genes-tsv", required=False, type=Path)
    parser.add_argument("--target-sequences-dir", required=False, type=Path)
    parser.add_argument("--clinvar-vcf", required=False, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()

def _refseq_accession_to_gnomad_chrom(chr_acc: str) -> str | None:
    if not chr_acc: return None
    c = str(chr_acc).strip()
    if c.startswith("chr"): c = c[3:]
    if c in {"X", "Y", "MT", "M"}: return "MT" if c == "M" else c
    if c.isdigit():
        val = int(c)
        if val == 23: return "X"
        if val == 24: return "Y"
        return str(val)
    import re
    match = re.search(r"NC_0+(\d+)\.", c)
    if not match: return None
    num = int(match.group(1))
    if num == 23: return "X"
    if num == 24: return "Y"
    if num in {12920, 1807}: return "MT"
    return str(num)

def cluster_positions(positions: list[int], max_gap: int = 100000) -> list[tuple[int, int]]:
    if not positions: return []
    positions = sorted(positions)
    clusters = []
    start = positions[0]
    last = positions[0]
    for p in positions[1:]:
        if p - last > max_gap:
            clusters.append((start, last))
            start = p
        last = p
    clusters.append((start, last))
    return clusters


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def read_fasta_sequence(path: Path) -> str:
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            chunks.append(line.strip())
    return "".join(chunks).upper()


def normalize_chrom(value: str | None) -> str | None:
    if not value:
        return None
    c = str(value).strip()
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        return "MT"
    return c


def load_target_contexts(genes_tsv: Path | None, target_sequences_dir: Path | None) -> dict[str, dict]:
    if not genes_tsv and not target_sequences_dir:
        logger.warning("No target context provided; ClinVar/gnomAD lookup will use raw event keys.")
        return {}
    if not genes_tsv or not target_sequences_dir:
        raise ValueError("--genes-tsv and --target-sequences-dir must be provided together.")
    if not genes_tsv.exists():
        raise FileNotFoundError(f"Target genes table not found: {genes_tsv}")
    if not target_sequences_dir.exists():
        raise FileNotFoundError(f"Target sequences directory not found: {target_sequences_dir}")

    contexts = {}
    with open_text(genes_tsv) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", "genomic_accession", "chromosome", "begin", "end"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Target genes table missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            gene_id = row["gene_id"]
            fasta_path = target_sequences_dir / f"{gene_id}.fa.gz"
            if not fasta_path.exists():
                raise FileNotFoundError(f"Target FASTA not found for gene {gene_id}: {fasta_path}")
            contexts[gene_id] = {
                "gene_id": gene_id,
                "accession": row["genomic_accession"],
                "chrom": normalize_chrom(row["chromosome"]) or _refseq_accession_to_gnomad_chrom(row["genomic_accession"]),
                "begin": int(row["begin"]),
                "end": int(row["end"]),
                "fasta_path": fasta_path,
            }
    logger.info(f"Loaded target context for {len(contexts)} gene(s).")
    return contexts


def context_sequence(context: dict) -> str:
    seq = context.get("seq")
    if seq is None:
        seq = read_fasta_sequence(context["fasta_path"])
        context["seq"] = seq
    return seq


def build_context_index(contexts: dict[str, dict]) -> dict[str, tuple[list[dict], list[int]]]:
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for context in contexts.values():
        chrom = context.get("chrom")
        if chrom:
            by_chrom[chrom].append(context)

    index: dict[str, tuple[list[dict], list[int]]] = {}
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda row: (int(row["begin"]), int(row["end"]), row["gene_id"]))
        index[chrom] = (rows, [int(row["begin"]) for row in rows])
    return index


def normalize_vcf_key_for_context(
    context: dict,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> tuple[tuple[str, int, str, str] | None, str]:
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    if not ref or not alt:
        return None, "empty_vcf_allele"

    seq = context_sequence(context)
    pos0 = pos - int(context["begin"])
    if pos0 < 0 or pos0 + len(ref) > len(seq):
        return None, "out_of_target"
    if seq[pos0 : pos0 + len(ref)] != ref:
        return None, "ref_mismatch"

    if len(ref) != len(alt):
        while pos0 > 0 and ref[-1] == alt[-1]:
            prev = seq[pos0 - 1]
            ref = prev + ref[:-1]
            alt = prev + alt[:-1]
            pos0 -= 1

    return (chrom, int(context["begin"]) + pos0, ref, alt), "ok"


def event_vcf_key(row: dict, contexts: dict[str, dict]) -> tuple[tuple[str, int, str, str] | None, str]:
    gene_id = row.get("gene_id", "")
    context = contexts.get(gene_id)
    chrom = _refseq_accession_to_gnomad_chrom(row.get("genomic_accession", ""))
    if not chrom:
        return None, "unknown_chrom"

    try:
        raw_pos = int(row.get("genomic_start1", 0))
    except ValueError:
        return None, "bad_position"
    raw_key = (chrom, raw_pos, row.get("ref", "").upper(), row.get("alt", "").upper())

    if not context:
        return raw_key, "raw_no_context"

    try:
        start0 = int(row.get("target_start0", 0))
    except ValueError:
        return raw_key, "bad_target_position"

    seq = context_sequence(context)
    event_type = row.get("event_type", "")
    ref = row.get("ref", "").upper()
    alt = row.get("alt", "").upper()

    if event_type == "snv":
        if len(ref) != 1 or len(alt) != 1:
            return raw_key, "bad_snv_allele"
        vcf = (chrom, int(context["begin"]) + start0, ref, alt)
    elif event_type == "del":
        if not ref or alt:
            return raw_key, "bad_del_allele"
        if start0 <= 0:
            return raw_key, "missing_left_anchor"
        anchor = seq[start0 - 1]
        vcf = (chrom, int(context["begin"]) + start0 - 1, anchor + ref, anchor)
    elif event_type == "ins":
        if ref or not alt:
            return raw_key, "bad_ins_allele"
        if start0 <= 0:
            return raw_key, "missing_left_anchor"
        anchor = seq[start0 - 1]
        vcf = (chrom, int(context["begin"]) + start0 - 1, anchor, anchor + alt)
    else:
        return raw_key, "unsupported_event_type"

    normalized, status = normalize_vcf_key_for_context(context, *vcf)
    return (normalized or raw_key), status


def contexts_for_variant(
    context_index: dict[str, tuple[list[dict], list[int]]],
    chrom: str,
    pos: int,
) -> list[dict]:
    rows, starts = context_index.get(chrom, ([], []))
    limit = bisect.bisect_right(starts, pos)
    return [context for context in rows[:limit] if int(context["end"]) >= pos]


def add_annotation_cache_entry(
    cache: dict,
    key: tuple[str, int, str, str],
    value,
    contexts: dict[str, dict],
    context_index: dict[str, tuple[list[dict], list[int]]],
    status_counts: Counter,
) -> None:
    cache[key] = value
    chrom, pos, ref, alt = key
    matched_contexts = contexts_for_variant(context_index, chrom, pos)
    if not matched_contexts:
        status_counts["raw_no_context"] += 1
        return
    for context in matched_contexts:
        normalized, status = normalize_vcf_key_for_context(context, chrom, pos, ref, alt)
        status_counts[status] += 1
        if normalized:
            cache[normalized] = value


def fetch_gnomad_for_cluster(chrom: str, start: int, end: int) -> list[dict]:
    # pad by 100 bases
    return fetch_region_variants_recursive(chrom, max(1, start - 100), end + 100)


def format_info_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return "|".join(str(item) for item in value)
    return str(value)


def build_clinvar_cache(
    clinvar,
    accession_positions: dict[str, set[int]],
    contexts: dict[str, dict],
    context_index: dict[str, tuple[list[dict], list[int]]],
) -> dict[tuple[str, int, str, str], tuple[str, str, str]]:
    cache = {}
    if clinvar is None:
        return cache
    status_counts = Counter()
    for acc, positions in accession_positions.items():
        chrom = _refseq_accession_to_gnomad_chrom(acc)
        if not chrom:
            continue
        for start, end in cluster_positions(list(positions), max_gap=200000):
            try:
                for rec in clinvar.fetch(chrom, max(0, start - 1), end + 1):
                    rec_alts = rec.alts or ()
                    clin_sig = format_info_value(rec.info.get("CLNSIG"))
                    clin_rev = format_info_value(rec.info.get("CLNREVSTAT"))
                    clin_id = rec.id or ""
                    for alt in rec_alts:
                        add_annotation_cache_entry(
                            cache,
                            (chrom, rec.pos, rec.ref, alt),
                            (clin_sig, clin_rev, clin_id),
                            contexts,
                            context_index,
                            status_counts,
                        )
            except ValueError:
                logger.warning(f"ClinVar contig not found for chrom={chrom}; skipping region {start}-{end}.")
    if status_counts:
        logger.info(f"ClinVar key normalization status: {dict(status_counts)}")
    return cache


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = args.outdir / "alignment_events_annotated.tsv.gz"
    contexts = load_target_contexts(args.genes_tsv, args.target_sequences_dir)
    context_index = build_context_index(contexts)
    
    # 1. Read variants to find regions to query
    accession_positions = defaultdict(set)
    event_key_status_counts = Counter()
    with gzip.open(args.events_tsv, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames
        required = {"gene_id", "event_type", "target_start0", "genomic_accession", "genomic_start1", "ref", "alt"}
        missing = required - set(header or [])
        if missing:
            raise ValueError(f"Events table missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            acc = row["genomic_accession"]
            if acc:
                pos = int(row["genomic_start1"])
                accession_positions[acc].add(pos)
                lookup_key, status = event_vcf_key(row, contexts)
                event_key_status_counts[status] += 1
                if lookup_key:
                    accession_positions[acc].add(int(lookup_key[1]))
    logger.info(f"Event key normalization status: {dict(event_key_status_counts)}")
                
    # 2. Determine gnomAD clusters
    gnomad_tasks = []
    for acc, positions in accession_positions.items():
        chrom = _refseq_accession_to_gnomad_chrom(acc)
        if not chrom: continue
        clusters = cluster_positions(list(positions), max_gap=200000)
        for c_start, c_end in clusters:
            gnomad_tasks.append((chrom, c_start, c_end))
            
    logger.info(f"Will fetch {len(gnomad_tasks)} region(s) from gnomAD API.")
    
    # 3. Fetch gnomAD in parallel and cache
    gnomad_cache = {}
    gnomad_key_status_counts = Counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        future_to_task = {
            executor.submit(fetch_gnomad_for_cluster, chrom, start, end): (chrom, start, end)
            for chrom, start, end in gnomad_tasks
        }
        for future in concurrent.futures.as_completed(future_to_task):
            task = future_to_task[future]
            try:
                vars_list = future.result()
                for v in vars_list:
                    key = (v["chrom"], int(v["pos"]), v["ref"], v["alt"])
                    add_annotation_cache_entry(
                        gnomad_cache,
                        key,
                        v,
                        contexts,
                        context_index,
                        gnomad_key_status_counts,
                    )
            except Exception as exc:
                logger.error(f"gnomAD fetch failed for {task}: {exc}")

    logger.info(f"Cached {len(gnomad_cache)} gnomAD variants.")
    if gnomad_key_status_counts:
        logger.info(f"gnomAD key normalization status: {dict(gnomad_key_status_counts)}")

    # 4. Open ClinVar
    clinvar = None
    if args.clinvar_vcf:
        if not args.clinvar_vcf.exists():
            raise FileNotFoundError(f"ClinVar VCF not found: {args.clinvar_vcf}")
        clinvar = pysam.VariantFile(str(args.clinvar_vcf))
    else:
        logger.warning("No clinvar VCF provided. ClinVar annotation will be empty.")
    clinvar_cache = build_clinvar_cache(clinvar, accession_positions, contexts, context_index)
    logger.info(f"Cached {len(clinvar_cache)} ClinVar variants.")

    # 5. Annotate rows
    new_header = list(header) + [
        "clinvar_sig", "clinvar_revstat", "clinvar_id",
        "gnomad_af", "gnomad_af_source", "gnomad_csq"
    ]
    
    with gzip.open(out_tsv, "wt", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=new_header, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        with gzip.open(args.events_tsv, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                lookup_key, _ = event_vcf_key(row, contexts)

                clin_sig = ""
                clin_rev = ""
                clin_id = ""
                gnom_af = ""
                gnom_src = ""
                gnom_csq = ""

                # ClinVar Lookup
                if lookup_key:
                    clin_sig, clin_rev, clin_id = clinvar_cache.get(lookup_key, ("", "", ""))

                # gnomAD Lookup
                if lookup_key and lookup_key in gnomad_cache:
                    v = gnomad_cache[lookup_key]
                    af, af_src, _, _, _, _, _ = _select_af_metrics(v)
                    gnom_af = f"{af:.6g}" if af is not None else ""
                    gnom_src = af_src or ""
                    gnom_csq = v.get("consequence", "")

                new_row = dict(row)
                new_row["clinvar_sig"] = clin_sig
                new_row["clinvar_revstat"] = clin_rev
                new_row["clinvar_id"] = clin_id
                new_row["gnomad_af"] = gnom_af
                new_row["gnomad_af_source"] = gnom_src
                new_row["gnomad_csq"] = gnom_csq

                writer.writerow(new_row)
            
    logger.info(f"Saved annotated events to {out_tsv}")

if __name__ == "__main__":
    main()
