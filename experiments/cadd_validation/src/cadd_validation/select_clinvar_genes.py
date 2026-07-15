#!/usr/bin/env python3
"""Select ClinVar-rich Entrez Gene IDs for a real GAPH validation pilot."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

from .build_variant_universe import label_from_clnsig, parse_info, review_stars
from .io import write_tsv


GENE_FIELDS = [
    "gene_id",
    "symbol",
    "pathogenic_count",
    "benign_count",
    "total_count",
]
VARIANT_FIELDS = [
    "variant_id",
    "label",
    "gene_id",
    "symbol",
    "genomic_accession",
    "genomic_start1",
    "ref",
    "alt",
    "clinvar_id",
    "clinvar_sig",
    "clinvar_revstat",
    "clinvar_stars",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinvar-vcf", required=True, type=Path)
    parser.add_argument("--out-genes", required=True, type=Path)
    parser.add_argument("--out-gene-summary", required=True, type=Path)
    parser.add_argument("--out-variants", type=Path)
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--min-review-stars", type=int, default=1)
    parser.add_argument("--min-per-label", type=int, default=3)
    parser.add_argument("--max-genes", type=int, default=10)
    parser.add_argument("--max-variants-per-label-per-gene", type=int, default=25)
    parser.add_argument("--snv-only", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def open_vcf(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_geneinfo(raw: str) -> list[tuple[str, str]]:
    genes = []
    for part in (raw or "").split("|"):
        if not part or ":" not in part:
            continue
        symbol, gene_id = part.rsplit(":", 1)
        gene_id = gene_id.strip()
        symbol = symbol.strip()
        if gene_id and gene_id != "-1":
            genes.append((gene_id, symbol))
    return genes


def variant_id(accession: str, pos1: str, ref: str, alt: str) -> str:
    return f"{accession}:{pos1}:{ref}>{alt}"


def iter_labeled_gene_variants(args: argparse.Namespace):
    skipped = Counter()
    with open_vcf(args.clinvar_vcf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                skipped["malformed_vcf_row"] += 1
                continue
            chrom, pos, clinvar_id, ref, alts, _qual, _filter, raw_info = fields[:8]
            info = parse_info(raw_info)
            stars = review_stars(info.get("CLNREVSTAT", ""))
            if stars < args.min_review_stars:
                skipped["low_review_status"] += 1
                continue
            label = label_from_clnsig(info.get("CLNSIG", ""))
            if label is None:
                skipped["unsupported_label"] += 1
                continue
            genes = parse_geneinfo(info.get("GENEINFO", ""))
            if not genes:
                skipped["missing_geneinfo"] += 1
                continue
            for alt in alts.split(","):
                if args.snv_only and (len(ref) != 1 or len(alt) != 1):
                    skipped["non_snv"] += 1
                    continue
                for gene_id, symbol in genes:
                    yield {
                        "variant_id": variant_id(chrom, pos, ref, alt),
                        "label": label,
                        "gene_id": gene_id,
                        "symbol": symbol,
                        "genomic_accession": chrom,
                        "genomic_start1": pos,
                        "ref": ref,
                        "alt": alt,
                        "clinvar_id": clinvar_id,
                        "clinvar_sig": info.get("CLNSIG", ""),
                        "clinvar_revstat": info.get("CLNREVSTAT", ""),
                        "clinvar_stars": stars,
                    }
    yield {"__skipped__": dict(skipped)}


def main() -> None:
    args = parse_args()
    counts: dict[str, Counter] = defaultdict(Counter)
    symbols: dict[str, Counter] = defaultdict(Counter)
    variants_by_gene_label: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    skipped: dict[str, int] = {}

    for row in iter_labeled_gene_variants(args):
        if "__skipped__" in row:
            skipped = row["__skipped__"]
            continue
        gene_id = str(row["gene_id"])
        label = str(row["label"])
        counts[gene_id][label] += 1
        symbols[gene_id][str(row.get("symbol", ""))] += 1
        bucket = variants_by_gene_label[(gene_id, label)]
        if len(bucket) < args.max_variants_per_label_per_gene:
            bucket.append(row)

    candidates = []
    for gene_id, counter in counts.items():
        pathogenic = counter.get("pathogenic", 0)
        benign = counter.get("benign", 0)
        if pathogenic < args.min_per_label or benign < args.min_per_label:
            continue
        symbol = symbols[gene_id].most_common(1)[0][0] if symbols[gene_id] else ""
        candidates.append(
            {
                "gene_id": gene_id,
                "symbol": symbol,
                "pathogenic_count": pathogenic,
                "benign_count": benign,
                "total_count": pathogenic + benign,
            }
        )

    candidates.sort(
        key=lambda row: (
            min(int(row["pathogenic_count"]), int(row["benign_count"])),
            int(row["total_count"]),
            int(row["pathogenic_count"]),
        ),
        reverse=True,
    )
    selected = candidates[: args.max_genes]
    selected_ids = {str(row["gene_id"]) for row in selected}

    args.out_genes.parent.mkdir(parents=True, exist_ok=True)
    args.out_genes.write_text("".join(f"{row['gene_id']}\n" for row in selected))
    write_tsv(args.out_gene_summary, selected, GENE_FIELDS)

    variant_rows = []
    if args.out_variants:
        for gene_id in selected_ids:
            for label in ("pathogenic", "benign"):
                variant_rows.extend(variants_by_gene_label.get((gene_id, label), []))
        variant_rows.sort(key=lambda row: (str(row["gene_id"]), str(row["label"]), str(row["variant_id"])))
        write_tsv(args.out_variants, variant_rows, VARIANT_FIELDS)

    summary = {
        "clinvar_vcf": str(args.clinvar_vcf),
        "selected_gene_count": len(selected),
        "selected_variant_count": len(variant_rows),
        "candidate_gene_count": len(candidates),
        "min_review_stars": args.min_review_stars,
        "min_per_label": args.min_per_label,
        "max_genes": args.max_genes,
        "skipped": skipped,
    }
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

