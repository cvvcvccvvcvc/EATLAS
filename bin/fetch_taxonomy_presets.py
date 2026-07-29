#!/usr/bin/env python3
"""Build compact taxonomy and minimap2-preset metadata for ortholog taxa."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
from pathlib import Path
from typing import Iterable

from taxonomic_evidence import (
    build_taxonomy_summary_rows,
    load_taxonomy_profiles,
    write_taxonomy_summary,
)


ANCESTOR_IDS = {
    "hominidae": "9604",
    "primates": "9443",
    "mammalia": "40674",
    "vertebrata": "7742",
}


TSV_NULL = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthologs-tsv", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--taxonomy-classes", required=True, type=Path, help="Path to taxonomy_classes.json.gz")
    parser.add_argument("--datasets-bin", default="datasets")
    return parser.parse_args()


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


def fetch_taxonomy_records(
    tax_ids: list[str],
    datasets_bin: str,
    outdir: Path,
) -> dict[str, dict]:
    ids_path = outdir / "taxonomy.ids.txt"
    jsonl_path = outdir / "taxonomy.jsonl"
    ids_path.write_text("\n".join(tax_ids) + "\n")
    with jsonl_path.open("w") as output:
        result = subprocess.run(
            [
                datasets_bin,
                "summary",
                "taxonomy",
                "taxon",
                "--inputfile",
                str(ids_path),
                "--as-json-lines",
            ],
            text=True,
            stdout=output,
            stderr=subprocess.PIPE,
        )
    if result.returncode != 0:
        detail = (result.stderr or "").strip()
        raise RuntimeError(f"datasets taxonomy summary failed: {detail}")

    records: dict[str, dict] = {}
    with jsonl_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            payload = json.loads(line)
            record = payload.get("taxonomy") or {}
            tax_id = str(record.get("taxId") or record.get("tax_id") or "")
            if not tax_id:
                raise ValueError(f"Taxonomy JSONL row {line_number} has no tax_id")
            if tax_id in records:
                raise ValueError(f"Duplicate taxonomy response for tax_id {tax_id}")
            records[tax_id] = record
    return records


def parent_ids(record: dict) -> list[str]:
    values = []
    for item in record.get("parents") or []:
        value = (
            item.get("taxId") or item.get("tax_id") or item.get("id")
            if isinstance(item, dict)
            else item
        )
        if value is not None:
            values.append(str(value))
    tax_id = str(record.get("taxId") or record.get("tax_id") or "")
    if tax_id and tax_id not in values:
        values.append(tax_id)
    return values


def classification_value(record: dict, rank: str, field: str) -> str:
    return str(((record.get("classification") or {}).get(rank) or {}).get(field) or "")


def taxonomy_row(tax_id: str, preset_group: str, record: dict | None) -> dict[str, object]:
    preset = "asm10" if preset_group == "primates" else "asm20"
    record = record or {}
    ancestors = set(parent_ids(record))
    return {
        "tax_id": tax_id,
        "scientific_name": (
            (record.get("currentScientificName") or {}).get("name")
            or (record.get("current_scientific_name") or {}).get("name")
            or record.get("organism_name", "")
        ),
        "rank": record.get("rank", ""),
        "group_name": record.get("groupName") or record.get("group_name") or preset_group,
        "class_id": classification_value(record, "class", "id"),
        "class_name": classification_value(record, "class", "name"),
        "order_id": classification_value(record, "order", "id"),
        "order_name": classification_value(record, "order", "name"),
        "family_id": classification_value(record, "family", "id"),
        "family_name": classification_value(record, "family", "name"),
        "genus_id": classification_value(record, "genus", "id"),
        "genus_name": classification_value(record, "genus", "name"),
        "species_id": classification_value(record, "species", "id"),
        "species_name": classification_value(record, "species", "name"),
        "is_hominidae": str(ANCESTOR_IDS["hominidae"] in ancestors).lower(),
        "is_primate": str(ANCESTOR_IDS["primates"] in ancestors).lower(),
        "is_mammal": str(ANCESTOR_IDS["mammalia"] in ancestors).lower(),
        "is_vertebrate": str(ANCESTOR_IDS["vertebrata"] in ancestors).lower(),
        "preset_group": preset_group,
        "minimap2_preset": preset,
        "parent_tax_ids": ",".join(parent_ids(record)),
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
    taxonomy_records = fetch_taxonomy_records(tax_ids, args.datasets_bin, args.outdir)

    rows = []
    failures = []
    for tax_id in tax_ids:
        preset_group = taxonomy_dict.get(str(tax_id), "other_or_unknown")
        record = taxonomy_records.get(tax_id)
        if record is None:
            failures.append(
                {
                    "tax_id": tax_id,
                    "failure_type": "taxonomy_summary_missing",
                    "message": "NCBI Datasets taxonomy summary returned no record",
                }
            )
        rows.append(taxonomy_row(tax_id, preset_group, record))

    write_tsv_gz(args.outdir / "taxonomy_presets.tsv.gz", fields, rows)
    write_tsv_gz(args.outdir / "taxonomy_failures.tsv.gz", failure_fields, failures)
    taxonomy_profiles_path = args.outdir / "taxonomy_presets.tsv.gz"
    taxonomy_profiles = load_taxonomy_profiles(taxonomy_profiles_path)
    with gzip.open(args.orthologs_tsv, "rt", newline="") as handle:
        summary_rows = build_taxonomy_summary_rows(
            csv.DictReader(handle, delimiter="\t"),
            taxonomy_profiles,
        )
        write_taxonomy_summary(args.outdir / "taxonomy_summary.tsv.gz", summary_rows)

if __name__ == "__main__":
    main()
