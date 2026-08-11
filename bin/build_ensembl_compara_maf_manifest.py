#!/usr/bin/env python3
"""Build a small run-specific manifest for Ensembl Compara MAF chunks."""

from __future__ import annotations

import argparse
import csv
import gzip
import html.parser
import json
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, TextIO

from ensembl_compara_maf import BASE_URL, RELEASE, SPECIES_SET


FIELDS = [
    "release",
    "species_set",
    "base_url",
    "file_name",
    "source",
    "human_src",
    "seq_region",
    "chunk_order",
    "first_start1",
    "first_end1",
    "range_start1",
    "range_end1",
    "chrom_length",
]

FAILURE_FIELDS = ["seq_region", "file_name", "source", "failure_type", "message"]


@dataclass(frozen=True)
class ChunkProbe:
    file_name: str
    source: str
    human_src: str
    seq_region: str
    chunk_order: int
    first_start1: int
    first_end1: int
    chrom_length: int


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        attr_map = dict(attrs)
        href = attr_map.get("href")
        if href:
            self.links.append(href)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes-tsv", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--release", default=RELEASE)
    parser.add_argument("--species-set", default=SPECIES_SET)
    parser.add_argument(
        "--base-url",
        default=BASE_URL,
        help="HTTP(S) directory URL or local directory containing the selected MAF set.",
    )
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--retries", type=int, default=3)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv_gz(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
            count += 1
    return count


def refseq_to_ensembl_seq_region(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if text.startswith("chr"):
        text = text[3:]
    if text in {"X", "Y", "MT"} or text.isdigit():
        return str(int(text)) if text.isdigit() else text
    if text == "M":
        return "MT"
    if text.startswith("NC_"):
        base = text.split(".", 1)[0]
        try:
            number = int(base.split("_", 1)[1])
        except (IndexError, ValueError):
            return ""
        if 1 <= number <= 22:
            return str(number)
        if number == 23:
            return "X"
        if number == 24:
            return "Y"
        if number == 12920:
            return "MT"
    return ""


def needed_seq_regions(genes_tsv: Path) -> list[str]:
    regions = set()
    for row in read_tsv_gz(genes_tsv):
        region = refseq_to_ensembl_seq_region(row.get("genomic_accession", ""))
        if region:
            regions.add(region)
    return sorted(regions, key=lambda item: (not item.isdigit(), int(item) if item.isdigit() else item))


def source_join(base_url: str, file_name: str) -> str:
    if base_url.startswith(("http://", "https://")):
        return base_url.rstrip("/") + "/" + file_name
    return str(Path(base_url) / file_name)


def fetch_url_text(url: str, timeout: float, retries: int) -> str:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "gaph-ensembl-maf-manifest/0.1"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(2.0 * attempt)
                continue
            raise RuntimeError(f"Could not fetch {url}: {exc}") from exc
    raise RuntimeError(f"Could not fetch {url}: {last_error}")


def list_remote_files(base_url: str, timeout: float, retries: int) -> list[str]:
    parser = LinkParser()
    parser.feed(fetch_url_text(base_url.rstrip("/") + "/", timeout, retries))
    return [
        link
        for link in parser.links
        if link.endswith(".maf.gz") and not link.startswith(("../", "http://", "https://"))
    ]


def list_local_files(base_url: str) -> list[str]:
    return sorted(path.name for path in Path(base_url).glob("*.maf.gz"))


def list_maf_files(base_url: str, timeout: float, retries: int) -> list[str]:
    if base_url.startswith(("http://", "https://")):
        return list_remote_files(base_url, timeout, retries)
    return list_local_files(base_url)


def chunk_order(file_name: str, species_set: str, seq_region: str) -> int | None:
    escaped = re.escape(species_set)
    match = re.fullmatch(rf"{escaped}\.{re.escape(seq_region)}_(\d+)\.maf\.gz", file_name)
    if not match:
        return None
    return int(match.group(1))


def files_for_region(file_names: list[str], species_set: str, seq_region: str) -> list[tuple[int, str]]:
    rows: list[tuple[int, str]] = []
    for file_name in file_names:
        order = chunk_order(file_name, species_set, seq_region)
        if order is not None:
            rows.append((order, file_name))
    return sorted(rows)


def open_maf_text(source: str, timeout: float):
    if source.startswith(("http://", "https://")):
        request = urllib.request.Request(source, headers={"User-Agent": "gaph-ensembl-maf-manifest/0.1"})
        response = urllib.request.urlopen(request, timeout=timeout)
        return gzip.open(response, "rt")
    return gzip.open(source, "rt")


def forward_interval0(start0: int, size: int, strand: str, src_size: int) -> tuple[int, int]:
    if strand == "+":
        return start0, start0 + size
    if strand == "-":
        return src_size - (start0 + size), src_size - start0
    raise ValueError(f"Unsupported MAF strand: {strand}")


def probe_first_human_row(
    source: str,
    file_name: str,
    human_src: str,
    seq_region: str,
    order: int,
    timeout: float,
    retries: int,
) -> ChunkProbe:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with open_maf_text(source, timeout) as handle:
                for line in handle:
                    if not line.startswith("s "):
                        continue
                    fields = line.split()
                    if len(fields) < 7 or fields[1] != human_src:
                        continue
                    start0 = int(fields[2])
                    size = int(fields[3])
                    strand = fields[4]
                    src_size = int(fields[5])
                    first_start0, first_end0 = forward_interval0(start0, size, strand, src_size)
                    return ChunkProbe(
                        file_name=file_name,
                        source=source,
                        human_src=human_src,
                        seq_region=seq_region,
                        chunk_order=order,
                        first_start1=first_start0 + 1,
                        first_end1=first_end0,
                        chrom_length=src_size,
                    )
            raise RuntimeError(f"No {human_src} row found")
        except (urllib.error.URLError, TimeoutError, OSError, RuntimeError) as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(2.0 * attempt)
                continue
            raise RuntimeError(f"Could not probe {source}: {exc}") from exc
    raise RuntimeError(f"Could not probe {source}: {last_error}")


def manifest_rows_for_region(
    base_url: str,
    release: str,
    species_set: str,
    seq_region: str,
    file_names: list[str],
    timeout: float,
    retries: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    failures: list[dict[str, object]] = []
    probes: list[ChunkProbe] = []
    human_src = f"homo_sapiens.{seq_region}"
    for order, file_name in files_for_region(file_names, species_set, seq_region):
        source = source_join(base_url, file_name)
        try:
            probes.append(probe_first_human_row(source, file_name, human_src, seq_region, order, timeout, retries))
        except Exception as exc:
            failures.append(
                {
                    "seq_region": seq_region,
                    "file_name": file_name,
                    "source": source,
                    "failure_type": "chunk_probe_failed",
                    "message": str(exc),
                }
            )

    probes.sort(key=lambda row: (row.first_start1, row.chunk_order))
    rows: list[dict[str, object]] = []
    for index, probe in enumerate(probes):
        next_start = probes[index + 1].first_start1 if index + 1 < len(probes) else probe.chrom_length + 1
        rows.append(
            {
                "release": release,
                "species_set": species_set,
                "base_url": base_url,
                "file_name": probe.file_name,
                "source": probe.source,
                "human_src": probe.human_src,
                "seq_region": probe.seq_region,
                "chunk_order": probe.chunk_order,
                "first_start1": probe.first_start1,
                "first_end1": probe.first_end1,
                "range_start1": probe.first_start1,
                "range_end1": max(probe.first_start1, next_start - 1),
                "chrom_length": probe.chrom_length,
            }
        )
    if not probes:
        failures.append(
            {
                "seq_region": seq_region,
                "file_name": "",
                "source": base_url,
                "failure_type": "no_chunks_for_seq_region",
                "message": f"No MAF chunks were found for seq_region={seq_region}",
            }
        )
    return rows, failures


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    seq_regions = needed_seq_regions(args.genes_tsv)
    file_names = list_maf_files(args.base_url, args.timeout, args.retries)

    all_rows: list[dict[str, object]] = []
    all_failures: list[dict[str, object]] = []
    for seq_region in seq_regions:
        rows, failures = manifest_rows_for_region(
            args.base_url,
            args.release,
            args.species_set,
            seq_region,
            file_names,
            args.timeout,
            args.retries,
        )
        all_rows.extend(rows)
        all_failures.extend(failures)

    manifest_path = args.outdir / "ensembl_compara_maf_manifest.tsv.gz"
    failures_path = args.outdir / "ensembl_compara_maf_manifest_failures.tsv.gz"
    row_count = write_tsv_gz(manifest_path, FIELDS, all_rows)
    failure_count = write_tsv_gz(failures_path, FAILURE_FIELDS, all_failures)
    run_manifest = {
        "created_at": utc_now(),
        "release": args.release,
        "species_set": args.species_set,
        "base_url": args.base_url,
        "seq_regions": seq_regions,
        "chunk_count": row_count,
        "failure_count": failure_count,
        "manifest": str(manifest_path),
        "failures": str(failures_path),
    }
    (args.outdir / "ensembl_compara_maf_manifest.json").write_text(json.dumps(run_manifest, indent=2) + "\n")
    if failure_count:
        raise RuntimeError(
            f"Failed to build complete Ensembl Compara MAF manifest; "
            f"{failure_count} failure(s) written to {failures_path}"
        )


if __name__ == "__main__":
    main()
