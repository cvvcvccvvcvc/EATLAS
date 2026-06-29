#!/usr/bin/env python3
"""Join GAPH features with external baseline annotations and labels."""

from __future__ import annotations

import argparse
from pathlib import Path

from .io import read_tsv, write_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gaph-features-tsv", required=True, type=Path)
    parser.add_argument("--baseline-tsv", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--left-key", default="variant_id")
    parser.add_argument("--right-key", default="variant_id")
    parser.add_argument("--how", choices=["inner", "left"], default="inner")
    parser.add_argument(
        "--baseline-prefix",
        default="",
        help="Optional prefix applied to baseline columns that collide with GAPH columns.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    gaph_rows = read_tsv(args.gaph_features_tsv)
    baseline_rows = read_tsv(args.baseline_tsv)
    if not gaph_rows:
        raise ValueError("GAPH features table has no rows")
    if not baseline_rows:
        raise ValueError("Baseline table has no rows")
    baseline_index = {}
    for row in baseline_rows:
        key = row.get(args.right_key, "")
        if not key:
            continue
        if key in baseline_index:
            raise ValueError(f"Duplicate baseline key: {key}")
        baseline_index[key] = row

    gaph_fields = list(gaph_rows[0].keys())
    baseline_fields = [field for field in baseline_rows[0].keys() if field != args.right_key]
    output_fields = gaph_fields[:]
    for field in baseline_fields:
        out_field = field
        if out_field in output_fields:
            out_field = args.baseline_prefix + field if args.baseline_prefix else f"baseline_{field}"
        output_fields.append(out_field)

    out_rows: list[dict[str, object]] = []
    for row in gaph_rows:
        key = row.get(args.left_key, "")
        baseline = baseline_index.get(key)
        if baseline is None and args.how == "inner":
            continue
        out = dict(row)
        for field in baseline_fields:
            out_field = field
            if out_field in gaph_fields:
                out_field = args.baseline_prefix + field if args.baseline_prefix else f"baseline_{field}"
            out[out_field] = baseline.get(field, "") if baseline else ""
        out_rows.append(out)

    if not out_rows:
        raise ValueError("Join produced no rows")
    write_tsv(args.out_tsv, out_rows, output_fields)


if __name__ == "__main__":
    main()

