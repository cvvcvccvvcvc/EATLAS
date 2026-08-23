#!/usr/bin/env python3
"""Build canonical wide taxonomy metadata for ortholog taxa."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
from pathlib import Path
from typing import Iterable

from genomics.taxonomy import TAXONOMY_FIELDS

TSV_NULL = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthologs-tsv", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
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


def lineage_tax_ids(record: dict) -> list[str]:
    values: list[str] = []
    seen: set[str] = set()
    for item in record.get("parents") or []:
        value = (
            item.get("taxId") or item.get("tax_id") or item.get("id")
            if isinstance(item, dict)
            else item
        )
        if value is not None and str(value) not in seen:
            values.append(str(value))
            seen.add(str(value))
    tax_id = str(record.get("taxId") or record.get("tax_id") or "")
    if tax_id and tax_id not in seen:
        values.append(tax_id)
    return values


def classification_value(record: dict, rank: str, field: str) -> str:
    return str(((record.get("classification") or {}).get(rank) or {}).get(field) or "")


def taxonomy_row(tax_id: str, record: dict | None) -> dict[str, object]:
    taxonomy_status = "resolved" if record is not None else "not_returned"
    record = record or {}
    return {
        "tax_id": tax_id,
        "taxonomy_status": taxonomy_status,
        "scientific_name": (
            (record.get("currentScientificName") or {}).get("name")
            or (record.get("current_scientific_name") or {}).get("name")
            or record.get("organism_name", "")
        ),
        "rank": record.get("rank", ""),
        "group_name": record.get("groupName") or record.get("group_name") or "",
        "domain_id": classification_value(record, "domain", "id"),
        "domain_name": classification_value(record, "domain", "name"),
        "kingdom_id": classification_value(record, "kingdom", "id"),
        "kingdom_name": classification_value(record, "kingdom", "name"),
        "phylum_id": classification_value(record, "phylum", "id"),
        "phylum_name": classification_value(record, "phylum", "name"),
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
        "lineage_tax_ids": ",".join(lineage_tax_ids(record)),
    }


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    failure_fields = ["tax_id", "failure_type", "message"]

    tax_ids = read_unique_tax_ids(args.orthologs_tsv)
    taxonomy_records = fetch_taxonomy_records(tax_ids, args.datasets_bin, args.outdir)

    rows = []
    failures = []
    for tax_id in tax_ids:
        record = taxonomy_records.get(tax_id)
        if record is None:
            failures.append(
                {
                    "tax_id": tax_id,
                    "failure_type": "taxonomy_summary_missing",
                    "message": "NCBI Datasets taxonomy summary returned no record",
                }
            )
        rows.append(taxonomy_row(tax_id, record))

    write_tsv_gz(args.outdir / "taxonomy.tsv.gz", TAXONOMY_FIELDS, rows)
    write_tsv_gz(args.outdir / "taxonomy_failures.tsv.gz", failure_fields, failures)

if __name__ == "__main__":
    main()
