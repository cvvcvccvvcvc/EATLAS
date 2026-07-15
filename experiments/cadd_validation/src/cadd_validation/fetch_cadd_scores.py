#!/usr/bin/env python3
"""Fetch CADD raw/PHRED scores for a variant TSV using the public CADD API."""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path

from .io import read_tsv, write_tsv


NCBI_CHROM_RE = re.compile(r"^NC_(\d{6})(?:\.\d+)?$")
OUTPUT_FIELDS = [
    "variant_id",
    "CADD_RAW",
    "CADD_PHRED",
    "cadd_chrom",
    "cadd_pos",
    "cadd_ref",
    "cadd_alt",
    "cadd_query",
    "cadd_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants-tsv", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--genome-build", default="GRCh38-v1.7")
    parser.add_argument("--api-template", default="https://cadd.gs.washington.edu/api/v1.0/{build}/{query}")
    parser.add_argument("--sleep-seconds", type=float, default=0.2)
    parser.add_argument("--timeout-seconds", type=float, default=30.0)
    parser.add_argument("--max-variants", type=int)
    parser.add_argument("--cache-jsonl", type=Path)
    return parser.parse_args()


def cadd_chrom(value: str) -> str:
    raw = str(value or "").strip()
    match = NCBI_CHROM_RE.match(raw)
    if match:
        number = int(match.group(1))
        if 1 <= number <= 22:
            return str(number)
        if number == 23:
            return "X"
        if number == 24:
            return "Y"
    if raw == "NC_012920.1":
        return "MT"
    if raw.lower().startswith("chr"):
        return raw[3:]
    return raw


def make_query(row: dict[str, str]) -> tuple[str, str, str, str]:
    chrom = cadd_chrom(row.get("genomic_accession") or row.get("chrom") or row.get("#CHROM") or "")
    pos = str(row.get("genomic_start1") or row.get("pos") or row.get("POS") or "").strip()
    ref = str(row.get("ref") or row.get("REF") or "").strip().upper()
    alt = str(row.get("alt") or row.get("ALT") or "").strip().upper()
    return chrom, pos, ref, alt


def variant_id(row: dict[str, str], chrom: str, pos: str, ref: str, alt: str) -> str:
    return row.get("variant_id") or f"{row.get('genomic_accession') or chrom}:{pos}:{ref}>{alt}"


def api_url(template: str, build: str, chrom: str, pos: str, ref: str, alt: str) -> tuple[str, str]:
    query = f"{chrom}:{pos}_{ref}_{alt}"
    return template.format(build=build, query=query), query


def normalize_response(raw: str) -> tuple[str, str, str]:
    text = raw.strip()
    if not text:
        return "", "", "empty_response"
    try:
        payload = json.loads(text)
        return score_from_json(payload)
    except json.JSONDecodeError:
        return score_from_tsv(text)


def score_from_json(payload) -> tuple[str, str, str]:
    records = payload if isinstance(payload, list) else [payload]
    for record in records:
        if not isinstance(record, dict):
            continue
        lowered = {str(key).lower(): value for key, value in record.items()}
        raw = first_value(lowered, ["rawscore", "raw_score", "raw", "cadd_raw"])
        phred = first_value(lowered, ["phred", "phredscore", "cadd_phred"])
        if raw not in (None, "") and phred not in (None, ""):
            return str(raw), str(phred), "ok"
    return "", "", "score_not_found"


def first_value(record: dict[str, object], keys: list[str]):
    for key in keys:
        if key in record:
            return record[key]
    return None


def score_from_tsv(text: str) -> tuple[str, str, str]:
    lines = [line for line in text.splitlines() if line and not line.startswith("#")]
    if not lines:
        return "", "", "empty_tsv_response"
    reader = csv.DictReader(lines, delimiter="\t")
    if reader.fieldnames:
        for row in reader:
            lowered = {str(key).lower(): value for key, value in row.items()}
            raw = first_value(lowered, ["rawscore", "raw_score", "raw", "cadd_raw"])
            phred = first_value(lowered, ["phred", "phredscore", "cadd_phred"])
            if raw not in (None, "") and phred not in (None, ""):
                return str(raw), str(phred), "ok"
    return "", "", "score_not_found"


def fetch_url(url: str, timeout: float) -> tuple[str, str]:
    request = urllib.request.Request(url, headers={"User-Agent": "gaph-v2-cadd-validation/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.read().decode("utf-8"), "http_ok"
    except urllib.error.HTTPError as exc:
        return exc.read().decode("utf-8", errors="replace"), f"http_{exc.code}"
    except urllib.error.URLError as exc:
        return str(exc), "url_error"


def main() -> None:
    args = parse_args()
    rows = read_tsv(args.variants_tsv)
    if args.max_variants is not None:
        rows = rows[: args.max_variants]

    out_rows = []
    cache_handle = None
    if args.cache_jsonl:
        args.cache_jsonl.parent.mkdir(parents=True, exist_ok=True)
        cache_handle = args.cache_jsonl.open("w")
    try:
        for index, row in enumerate(rows):
            chrom, pos, ref, alt = make_query(row)
            vid = variant_id(row, chrom, pos, ref, alt)
            if not chrom or not pos or not ref or not alt:
                out_rows.append(
                    {
                        "variant_id": vid,
                        "CADD_RAW": "",
                        "CADD_PHRED": "",
                        "cadd_chrom": chrom,
                        "cadd_pos": pos,
                        "cadd_ref": ref,
                        "cadd_alt": alt,
                        "cadd_query": "",
                        "cadd_status": "missing_variant_key",
                    }
                )
                continue
            url, query = api_url(args.api_template, args.genome_build, chrom, pos, ref, alt)
            body, http_status = fetch_url(url, args.timeout_seconds)
            raw_score, phred, parse_status = normalize_response(body) if http_status == "http_ok" else ("", "", http_status)
            status = parse_status if parse_status != "ok" else "ok"
            out_rows.append(
                {
                    "variant_id": vid,
                    "CADD_RAW": raw_score,
                    "CADD_PHRED": phred,
                    "cadd_chrom": chrom,
                    "cadd_pos": pos,
                    "cadd_ref": ref,
                    "cadd_alt": alt,
                    "cadd_query": query,
                    "cadd_status": status,
                }
            )
            if cache_handle:
                cache_handle.write(
                    json.dumps(
                        {
                            "variant_id": vid,
                            "url": url,
                            "http_status": http_status,
                            "parse_status": parse_status,
                            "body": body[:2000],
                        },
                        sort_keys=True,
                    )
                    + "\n"
                )
            if args.sleep_seconds and index + 1 < len(rows):
                time.sleep(args.sleep_seconds)
    finally:
        if cache_handle:
            cache_handle.close()
    write_tsv(args.out_tsv, out_rows, OUTPUT_FIELDS)


if __name__ == "__main__":
    main()
