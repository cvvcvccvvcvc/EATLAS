"""Small TSV helpers used by the CADD validation utilities."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Iterable


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def iter_tsv(path: Path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield dict(row)


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
            count += 1
    return count

