#!/usr/bin/env python3
"""Fetch a GERP conservation-score region from an Ensembl Compara bigWig."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path


DEFAULT_BIGWIG_URL = (
    "https://ftp.ensembl.org/pub/release-116/compara/conservation_scores/"
    "92_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bigwig", default=DEFAULT_BIGWIG_URL, help="Local path or http(s) URL to a GERP bigWig")
    parser.add_argument("--chrom", required=True, help="Ensembl chromosome/seq_region name, for example 4")
    parser.add_argument("--start1", required=True, type=int)
    parser.add_argument("--end1", required=True, type=int)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--score-source", default="ensembl_92_mammals_gerp")
    return parser.parse_args()


def require_pybigwig():
    try:
        import pyBigWig  # type: ignore
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "pyBigWig is required for bigWig region access. "
            "Install it in the active environment or set PYTHONPATH to a temporary target."
        ) from exc
    return pyBigWig


def finite_values(values: list[float | None]) -> list[float]:
    result: list[float] = []
    for value in values:
        if value is None:
            continue
        numeric = float(value)
        if math.isfinite(numeric):
            result.append(numeric)
    return result


def main() -> None:
    args = parse_args()
    if args.end1 < args.start1:
        raise ValueError("--end1 must be >= --start1")

    pyBigWig = require_pybigwig()
    start0 = args.start1 - 1
    end0 = args.end1
    bw = pyBigWig.open(args.bigwig)
    if not bw:
        raise RuntimeError(f"Could not open bigWig: {args.bigwig}")
    chroms = bw.chroms()
    if args.chrom not in chroms:
        bw.close()
        raise ValueError(f"Chromosome {args.chrom!r} not found in bigWig. Available example keys: {list(chroms)[:10]}")
    if end0 > chroms[args.chrom]:
        bw.close()
        raise ValueError(f"Requested end {args.end1} exceeds chromosome length {chroms[args.chrom]}")

    intervals = bw.intervals(args.chrom, start0, end0) or []
    values = finite_values(bw.values(args.chrom, start0, end0, numpy=False))
    bw.close()

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.out_tsv, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "genomic_accession",
                "genomic_start1",
                "genomic_end1",
                "score_source",
                "score",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for interval_start0, interval_end0, score in intervals:
            writer.writerow(
                {
                    "genomic_accession": args.chrom,
                    "genomic_start1": int(interval_start0) + 1,
                    "genomic_end1": int(interval_end0),
                    "score_source": args.score_source,
                    "score": f"{float(score):.6g}",
                }
            )

    summary = {
        "bigwig": args.bigwig,
        "chrom": args.chrom,
        "start1": args.start1,
        "end1": args.end1,
        "length_bp": args.end1 - args.start1 + 1,
        "interval_count": len(intervals),
        "value_count": len(values),
        "missing_count": (args.end1 - args.start1 + 1) - len(values),
        "min": min(values) if values else None,
        "max": max(values) if values else None,
        "mean": (sum(values) / len(values)) if values else None,
        "score_source": args.score_source,
        "out_tsv": str(args.out_tsv),
    }
    args.summary_json.parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
