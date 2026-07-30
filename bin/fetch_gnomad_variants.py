#!/usr/bin/env python3
"""Fetch gnomAD variants for one region and write a compact VCF."""

import argparse
import logging
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomics.gnomad import fetch_region_to_vcf


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--out-vcf", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    fetch_region_to_vcf(args.chrom, args.start, args.end, args.out_vcf)


if __name__ == "__main__":
    main()
