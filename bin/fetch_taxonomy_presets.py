#!/usr/bin/env python3
"""Build compact taxonomy-to-minimap2-preset metadata for ortholog taxa using offline dictionary."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path
from typing import Iterable


TSV_NULL = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthologs-tsv", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--taxonomy-classes", required=True, type=Path, help="Path to taxonomy_classes.json.gz")
    # Kept for compatibility but ignored
    parser.add_argument("--datasets-bin", help=argparse.SUPPRESS)
    parser.add_argument("--batch-size", help=argparse.SUPPRESS)
    parser.add_argument("--taxonomy-file", help=argparse.SUPPRESS)
    parser.add_argument("--fast-mock", action="store_true", help=argparse.SUPPRESS)
    return parser.parse_args()


def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def write_tsv_gz(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, TSV_NULL) for field in fields})
            count += 1
    return count


def read_unique_tax_ids(path: Path) -> list[str]:
    tax_ids: set[str] = set()
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tax_id = (row.get("tax_id") or "").strip()
            if tax_id:
                tax_ids.add(tax_id)
    return sorted(tax_ids, key=lambda value: int(value) if value.isdigit() else value)


def load_taxonomy_dict(path: Path) -> dict[str, str]:
    with gzip.open(path, "rt") as f:
        return json.load(f)


def taxonomy_row(tax_id: str, preset_group: str) -> dict[str, object]:
    preset = "asm10" if preset_group == "primates" else "asm20"
    
    return {
        "tax_id": tax_id,
        "scientific_name": "",
        "rank": "",
        "group_name": "",
        "class_id": "",
        "class_name": "",
        "order_id": "",
        "order_name": "",
        "family_id": "",
        "family_name": "",
        "genus_id": "",
        "genus_name": "",
        "species_id": "",
        "species_name": "",
        "is_hominidae": str(preset_group == "primates").lower(), # Approximated
        "is_primate": str(preset_group == "primates").lower(),
        "is_mammal": str(preset_group in ("primates", "other_mammals")).lower(),
        "is_vertebrate": str(preset_group in ("primates", "other_mammals", "other_vertebrates")).lower(),
        "preset_group": preset_group,
        "minimap2_preset": preset,
        "parent_tax_ids": "",
    }


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    fields = [
        "tax_id",
        "scientific_name",
        "rank",
        "group_name",
        "class_id",
        "class_name",
        "order_id",
        "order_name",
        "family_id",
        "family_name",
        "genus_id",
        "genus_name",
        "species_id",
        "species_name",
        "is_hominidae",
        "is_primate",
        "is_mammal",
        "is_vertebrate",
        "preset_group",
        "minimap2_preset",
        "parent_tax_ids",
    ]
    failure_fields = ["tax_id", "failure_type", "message"]

    tax_ids = read_unique_tax_ids(args.orthologs_tsv)
    
    taxonomy_dict = load_taxonomy_dict(args.taxonomy_classes)
    
    rows = []
    failures = []
    
    for tax_id in tax_ids:
        preset_group = taxonomy_dict.get(str(tax_id), "other_or_unknown")
        rows.append(taxonomy_row(tax_id, preset_group))

    write_tsv_gz(args.outdir / "taxonomy_presets.tsv.gz", fields, rows)
    write_tsv_gz(args.outdir / "taxonomy_failures.tsv.gz", failure_fields, failures)

if __name__ == "__main__":
    main()
