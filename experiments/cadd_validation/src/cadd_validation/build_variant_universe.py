#!/usr/bin/env python3
"""Build an independent labeled variant universe from ClinVar VCF."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from .fetch_cadd_scores import cadd_chrom
from .io import iter_tsv, write_tsv


LABEL_SPLIT_RE = re.compile(r"[/|,;]")
OUTPUT_FIELDS = [
    "variant_id",
    "label",
    "gene_id",
    "genomic_accession",
    "genomic_start1",
    "ref",
    "alt",
    "target_start0",
    "clinvar_id",
    "clinvar_sig",
    "clinvar_revstat",
    "clinvar_stars",
]


@dataclass(frozen=True)
class GeneInterval:
    gene_id: str
    genomic_accession: str
    genomic_start1: int
    genomic_end1: int

    @property
    def low(self) -> int:
        return min(self.genomic_start1, self.genomic_end1)

    @property
    def high(self) -> int:
        return max(self.genomic_start1, self.genomic_end1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinvar-vcf", required=True, type=Path)
    parser.add_argument("--target-features-tsv", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--min-review-stars", type=int, default=1)
    parser.add_argument("--snv-only", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def open_vcf(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def load_gene_intervals(path: Path) -> dict[str, list[GeneInterval]]:
    by_accession: dict[str, list[GeneInterval]] = defaultdict(list)
    for row in iter_tsv(path):
        if row.get("feature_type") != "gene":
            continue
        gene_id = row.get("gene_id", "")
        accession = row.get("genomic_accession", "")
        if not gene_id or not accession:
            continue
        try:
            start1 = int(row["genomic_start1"])
            end1 = int(row["genomic_end1"])
        except (KeyError, TypeError, ValueError):
            continue
        interval = GeneInterval(gene_id, accession, start1, end1)
        aliases = {accession}
        chrom = cadd_chrom(accession)
        if chrom:
            aliases.add(chrom)
            aliases.add(f"chr{chrom}")
        for alias in aliases:
            by_accession[alias].append(interval)
    for intervals in by_accession.values():
        intervals.sort(key=lambda item: (item.low, item.high, item.gene_id))
    return by_accession


def parse_info(raw: str) -> dict[str, str]:
    info = {}
    for part in raw.split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            info[key] = unquote(value)
        else:
            info[part] = "true"
    return info


def review_stars(revstat: str) -> int:
    text = revstat.lower()
    if "practice_guideline" in text:
        return 4
    if "reviewed_by_expert_panel" in text:
        return 3
    if "multiple_submitters" in text and "no_conflicts" in text:
        return 2
    if "criteria_provided" in text and "single_submitter" in text:
        return 1
    if "criteria_provided" in text:
        return 1
    return 0


def label_from_clnsig(clnsig: str) -> str | None:
    text = clnsig.lower()
    if not text or "uncertain" in text or "conflicting" in text or "not_provided" in text:
        return None
    tokens = {token.strip() for token in LABEL_SPLIT_RE.split(text) if token.strip()}
    pathogenic = bool({"pathogenic", "likely_pathogenic"} & tokens)
    benign = bool({"benign", "likely_benign"} & tokens)
    if pathogenic and not benign:
        return "pathogenic"
    if benign and not pathogenic:
        return "benign"
    return None


def overlapping_genes(intervals: dict[str, list[GeneInterval]], accession: str, pos1: int) -> list[GeneInterval]:
    return [item for item in intervals.get(accession, []) if item.low <= pos1 <= item.high]


def variant_id(accession: str, pos1: int, ref: str, alt: str) -> str:
    return f"{accession}:{pos1}:{ref}>{alt}"


def build_rows(args: argparse.Namespace) -> tuple[list[dict[str, object]], dict[str, object]]:
    intervals = load_gene_intervals(args.target_features_tsv)
    rows = []
    skipped = Counter()
    with open_vcf(args.clinvar_vcf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                skipped["malformed_vcf_row"] += 1
                continue
            chrom, pos_text, clinvar_id, ref, alts, _qual, _filter, raw_info = fields[:8]
            try:
                pos1 = int(pos_text)
            except ValueError:
                skipped["bad_position"] += 1
                continue
            info = parse_info(raw_info)
            clnsig = info.get("CLNSIG", "")
            revstat = info.get("CLNREVSTAT", "")
            stars = review_stars(revstat)
            if stars < args.min_review_stars:
                skipped["low_review_status"] += 1
                continue
            label = label_from_clnsig(clnsig)
            if label is None:
                skipped["unsupported_label"] += 1
                continue
            genes = overlapping_genes(intervals, chrom, pos1)
            if not genes:
                skipped["outside_target_genes"] += 1
                continue
            for alt in alts.split(","):
                if args.snv_only and (len(ref) != 1 or len(alt) != 1):
                    skipped["non_snv"] += 1
                    continue
                for gene in genes:
                    rows.append(
                        {
                            "variant_id": variant_id(chrom, pos1, ref, alt),
                            "label": label,
                            "gene_id": gene.gene_id,
                            "genomic_accession": chrom,
                            "genomic_start1": pos1,
                            "ref": ref,
                            "alt": alt,
                            "target_start0": pos1 - gene.low,
                            "clinvar_id": clinvar_id,
                            "clinvar_sig": clnsig,
                            "clinvar_revstat": revstat,
                            "clinvar_stars": stars,
                        }
                    )
    summary = {
        "clinvar_vcf": str(args.clinvar_vcf),
        "target_features_tsv": str(args.target_features_tsv),
        "row_count": len(rows),
        "skipped": dict(skipped),
        "min_review_stars": args.min_review_stars,
        "snv_only": args.snv_only,
    }
    return rows, summary


def main() -> None:
    args = parse_args()
    rows, summary = build_rows(args)
    if not rows:
        raise ValueError("No labeled ClinVar rows overlapped target genes")
    write_tsv(args.out_tsv, rows, OUTPUT_FIELDS)
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
