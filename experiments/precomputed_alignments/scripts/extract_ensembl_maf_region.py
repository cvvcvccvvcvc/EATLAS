#!/usr/bin/env python3
"""Extract Ensembl Compara MAF blocks overlapping a human region.

The output JSON intentionally matches the shape consumed by
rest_alignment_to_gaph_tables.py:

[
  {
    "alignments": [
      {
        "species": "homo_sapiens",
        "seq_region": "4",
        "start": 122612501,
        "end": 122612650,
        "strand": 1,
        "seq": "ACGT..."
      }
    ]
  }
]
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class MafSequence:
    src: str
    start0: int
    size: int
    strand: str
    src_size: int
    text: str

    def species_and_region(self) -> tuple[str, str]:
        if "." not in self.src:
            return self.src, ""
        species, seq_region = self.src.split(".", 1)
        return species, seq_region

    def forward_interval0(self) -> tuple[int, int]:
        if self.strand == "+":
            return self.start0, self.start0 + self.size
        if self.strand == "-":
            return self.src_size - (self.start0 + self.size), self.src_size - self.start0
        raise ValueError(f"Unsupported MAF strand for {self.src}: {self.strand}")

    def rest_strand(self) -> int:
        return 1 if self.strand == "+" else -1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", action="append", required=True, help="Local path or http(s) URL to .maf.gz")
    parser.add_argument("--human-src", default="homo_sapiens.4", help="MAF source name for the human chromosome")
    parser.add_argument("--start1", required=True, type=int, help="1-based inclusive human start")
    parser.add_argument("--end1", required=True, type=int, help="1-based inclusive human end")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-blocks", type=int, default=0, help="Stop after this many overlapping blocks; 0 means no limit")
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--progress-every", type=int, default=10000)
    return parser.parse_args()


def open_text_maf(source: str, timeout: float) -> TextIO:
    if source.startswith(("http://", "https://")):
        request = urllib.request.Request(source, headers={"User-Agent": "gaph-maf-region-extractor/0.1"})
        response = urllib.request.urlopen(request, timeout=timeout)
        return gzip.open(response, "rt")
    return gzip.open(source, "rt")


def iter_maf_blocks(handle: TextIO) -> Iterable[list[MafSequence]]:
    block: list[MafSequence] = []
    for line in handle:
        line = line.rstrip("\n")
        if not line:
            if block:
                yield block
                block = []
            continue
        if line.startswith("s "):
            fields = line.split()
            if len(fields) < 7:
                continue
            block.append(
                MafSequence(
                    src=fields[1],
                    start0=int(fields[2]),
                    size=int(fields[3]),
                    strand=fields[4],
                    src_size=int(fields[5]),
                    text=fields[6],
                )
            )
    if block:
        yield block


def overlaps(start0: int, end0: int, query_start0: int, query_end0: int) -> bool:
    return end0 > query_start0 and start0 < query_end0


def reverse_complement_alignment(text: str) -> str:
    return text.translate(COMPLEMENT)[::-1]


def to_rest_row(row: MafSequence, flip_orientation: bool) -> dict[str, object]:
    species, seq_region = row.species_and_region()
    start0, end0 = row.forward_interval0()
    strand = row.rest_strand()
    text = row.text
    if flip_orientation:
        strand *= -1
        text = reverse_complement_alignment(text)
    return {
        "species": species,
        "seq_region": seq_region,
        "start": start0 + 1,
        "end": end0,
        "strand": strand,
        "seq": text,
        "description": row.src,
    }


def extract_blocks(args: argparse.Namespace) -> list[dict[str, object]]:
    query_start0 = args.start1 - 1
    query_end0 = args.end1
    extracted: list[dict[str, object]] = []
    scanned_blocks = 0
    scanned_files = 0

    for source in args.maf:
        scanned_files += 1
        print(f"Scanning {source}", file=sys.stderr)
        with open_text_maf(source, args.timeout) as handle:
            for block in iter_maf_blocks(handle):
                scanned_blocks += 1
                if args.progress_every and scanned_blocks % args.progress_every == 0:
                    print(f"Scanned {scanned_blocks:,} blocks", file=sys.stderr)
                human_rows = [row for row in block if row.src == args.human_src]
                if not human_rows:
                    continue
                human = human_rows[0]
                human_start0, human_end0 = human.forward_interval0()
                if not overlaps(human_start0, human_end0, query_start0, query_end0):
                    continue
                flip_orientation = human.strand == "-"
                extracted.append(
                    {
                        "source_maf": source,
                        "human_src": args.human_src,
                        "requested_start1": args.start1,
                        "requested_end1": args.end1,
                        "human_block_start1": human_start0 + 1,
                        "human_block_end1": human_end0,
                        "alignments": [to_rest_row(row, flip_orientation) for row in block],
                    }
                )
                if args.max_blocks and len(extracted) >= args.max_blocks:
                    print(
                        f"Stopping after {len(extracted)} overlapping blocks; scanned {scanned_blocks:,} blocks",
                        file=sys.stderr,
                    )
                    return extracted
    print(
        f"Finished scanning {scanned_files} file(s), {scanned_blocks:,} blocks, extracted {len(extracted)} blocks",
        file=sys.stderr,
    )
    return extracted


def main() -> None:
    args = parse_args()
    if args.end1 < args.start1:
        raise ValueError("--end1 must be >= --start1")
    blocks = extract_blocks(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(blocks, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.output}")
    print(json.dumps({"block_count": len(blocks), "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
