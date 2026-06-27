#!/usr/bin/env python3
"""Annotate TSV variants with ClinVar and gnomAD."""

import argparse
import csv
import gzip
import logging
import sys
import os
import concurrent.futures
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

def fetch_gnomad_for_cluster(chrom: str, start: int, end: int) -> list[dict]:
    # pad by 100 bases
    return fetch_region_variants_recursive(chrom, max(1, start - 100), end + 100)

def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = args.outdir / "alignment_events_annotated.tsv.gz"
    
    # 1. Read variants to find regions to query
    accession_positions = {}
    rows = []
    with gzip.open(args.events_tsv, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames
        for row in reader:
            rows.append(row)
            acc = row["genomic_accession"]
            if acc:
                pos = int(row["genomic_start1"])
                if acc not in accession_positions:
                    accession_positions[acc] = set()
                accession_positions[acc].add(pos)
                
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
                    gnomad_cache[key] = v
            except Exception as exc:
                logger.error(f"gnomAD fetch failed for {task}: {exc}")

    logger.info(f"Cached {len(gnomad_cache)} gnomAD variants.")

    # 4. Open ClinVar
    clinvar = None
    if args.clinvar_vcf and args.clinvar_vcf.exists():
        clinvar = pysam.VariantFile(str(args.clinvar_vcf))
    else:
        logger.warning("No clinvar VCF provided or found. ClinVar annotation will be empty.")

    # 5. Annotate rows
    new_header = list(header) + [
        "clinvar_sig", "clinvar_revstat", "clinvar_id",
        "gnomad_af", "gnomad_af_source", "gnomad_csq"
    ]
    
    with gzip.open(out_tsv, "wt", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=new_header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        
        for row in rows:
            acc = row.get("genomic_accession")
            pos = int(row.get("genomic_start1", 0))
            ref = row.get("ref", "")
            alt = row.get("alt", "")
            chrom = _refseq_accession_to_gnomad_chrom(acc)
            
            clin_sig = ""
            clin_rev = ""
            clin_id = ""
            gnom_af = ""
            gnom_src = ""
            gnom_csq = ""
            
            # ClinVar Lookup
            if clinvar and chrom:
                try:
                    for rec in clinvar.fetch(chrom, max(0, pos-1), pos+1):
                        if rec.pos == pos and rec.ref == ref and alt in rec.alts:
                            clin_sig = "|".join(rec.info.get("CLNSIG", []))
                            clin_rev = "|".join(rec.info.get("CLNREVSTAT", []))
                            clin_id = rec.id or ""
                            break
                except ValueError:
                    pass # contig not found or something
                    
            # gnomAD Lookup
            if chrom:
                key = (chrom, pos, ref, alt)
                if key in gnomad_cache:
                    v = gnomad_cache[key]
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
