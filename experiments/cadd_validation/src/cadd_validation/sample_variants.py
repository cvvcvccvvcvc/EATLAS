#!/usr/bin/env python3
"""Create a reproducible balanced sample from a variant TSV."""

from __future__ import annotations

import argparse
import random
from collections import defaultdict
from pathlib import Path

from .io import read_tsv, write_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants-tsv", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--group-columns", default="gene_id,label")
    parser.add_argument("--max-per-group", type=int, default=20)
    parser.add_argument("--random-seed", type=int, default=13)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_tsv(args.variants_tsv)
    if not rows:
        raise ValueError("Input variants table has no rows")
    group_columns = [col.strip() for col in args.group_columns.split(",") if col.strip()]
    if not group_columns:
        raise ValueError("--group-columns must contain at least one column")
    missing = [col for col in group_columns if col not in rows[0]]
    if missing:
        raise ValueError("Missing group column(s): " + ", ".join(missing))
    rng = random.Random(args.random_seed)
    groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = tuple(row.get(col, "") for col in group_columns)
        groups[key].append(row)
    sampled = []
    for key in sorted(groups):
        bucket = groups[key]
        bucket = sorted(bucket, key=lambda row: row.get("variant_id", ""))
        if len(bucket) > args.max_per_group:
            bucket = rng.sample(bucket, args.max_per_group)
            bucket.sort(key=lambda row: row.get("variant_id", ""))
        sampled.extend(bucket)
    sampled.sort(key=lambda row: tuple(row.get(col, "") for col in group_columns) + (row.get("variant_id", ""),))
    write_tsv(args.out_tsv, sampled, list(rows[0].keys()))


if __name__ == "__main__":
    main()

