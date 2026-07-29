#!/usr/bin/env python3
"""Retry failed gnomAD regions without replacing the original annotation output."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import json
import os
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from fetch_gnomad_variants import fetch_region_variants_recursive, select_af_metrics
from gnomad_cache import GnomadRegionCache

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analytics.core.variant_keys import (  # noqa: E402
    build_context_index,
    contexts_for_variant,
    load_target_contexts,
    normalize_chrom,
    normalize_vcf_key_for_context,
    parse_variant_key,
)


REQUIRED_FILES = (
    "variant_annotations.tsv.gz",
    "failures.tsv.gz",
    "manifest.json",
)
FAILURE_FIELDS = ["source", "scope", "chrom", "start", "end", "failure_type", "message"]
GNOMAD_COLUMNS = ["gnomad_af", "gnomad_af_source", "gnomad_csq"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--source-annotation-dir", type=Path)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument(
        "--gnomad-cache-dir",
        type=Path,
        default=os.environ.get("GAPH_GNOMAD_CACHE_DIR") or None,
        help="Optional shared directory for resumable gnomAD regional responses.",
    )
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def atomic_path(destination: Path, suffix: str) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(dir=destination.parent, suffix=suffix, delete=False)
    handle.close()
    return Path(handle.name)


def read_failures(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != FAILURE_FIELDS:
            raise ValueError(
                f"Unexpected failure table header in {path}: {reader.fieldnames}"
            )
        return list(reader)


def failed_gnomad_region(row: dict[str, str]) -> tuple[str, int, int] | None:
    if row.get("source") != "gnomad" or row.get("scope") != "region":
        return None
    try:
        start = int(row["start"])
        end = int(row["end"])
    except (KeyError, TypeError, ValueError):
        return None
    chrom = normalize_chrom(row.get("chrom"))
    if not chrom or start < 1 or end < start:
        return None
    return chrom, start, end


def fetch_failed_regions(
    regions: list[tuple[str, int, int]],
    workers: int,
    region_cache: GnomadRegionCache,
) -> tuple[dict[tuple[str, int, int], list[dict]], dict[tuple[str, int, int], Exception]]:
    successes: dict[tuple[str, int, int], list[dict]] = {}
    failures: dict[tuple[str, int, int], Exception] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                region_cache.fetch_region,
                chrom,
                max(1, start - 100),
                end + 100,
            ): (chrom, start, end)
            for chrom, start, end in regions
        }
        for future in concurrent.futures.as_completed(futures):
            region = futures[future]
            try:
                successes[region] = future.result()
            except Exception as exc:
                failures[region] = exc
    return successes, failures


def build_gnomad_cache(
    records_by_region: dict[tuple[str, int, int], list[dict]],
    contexts: dict[str, dict],
) -> tuple[dict[tuple[str, int, str, str], dict], Counter]:
    cache: dict[tuple[str, int, str, str], dict] = {}
    status_counts: Counter = Counter()
    context_index = build_context_index(contexts)
    for records in records_by_region.values():
        for record in records:
            key = (
                normalize_chrom(record.get("chrom")) or "",
                int(record.get("pos", 0)),
                str(record.get("ref", "")).upper(),
                str(record.get("alt", "")).upper(),
            )
            cache[key] = record
            chrom, pos, ref, alt = key
            matched_contexts = contexts_for_variant(context_index, chrom, pos)
            if not matched_contexts:
                status_counts["raw_no_context"] += 1
                continue
            for context in matched_contexts:
                normalized, status = normalize_vcf_key_for_context(
                    context,
                    chrom,
                    pos,
                    ref,
                    alt,
                )
                status_counts[status] += 1
                if normalized:
                    cache[normalized] = record
    return cache, status_counts


def rewrite_annotations(
    source: Path,
    destination: Path,
    cache: dict[tuple[str, int, str, str], dict],
) -> tuple[int, int, Counter]:
    temporary = atomic_path(destination, ".tsv.gz")
    row_count = 0
    updated_count = 0
    nonempty_counts: Counter = Counter()
    try:
        with gzip.open(source, "rt", newline="") as in_handle, gzip.open(
            temporary, "wt", newline=""
        ) as out_handle:
            reader = csv.DictReader(in_handle, delimiter="\t")
            header = reader.fieldnames or []
            required = {"variant_key", *GNOMAD_COLUMNS}
            missing = required - set(header)
            if missing:
                raise ValueError(
                    "Variant annotations missing required columns: "
                    + ", ".join(sorted(missing))
                )
            writer = csv.DictWriter(
                out_handle,
                fieldnames=header,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                row_count += 1
                key = parse_variant_key(row.get("variant_key", ""))
                record = cache.get(key) if key is not None else None
                if record is not None:
                    annotation = gnomad_annotation(record)
                    if any(row.get(column, "") != annotation[column] for column in GNOMAD_COLUMNS):
                        updated_count += 1
                    row.update(annotation)
                for column in header:
                    if column.startswith(("clinvar_", "gnomad_")) and row[column]:
                        nonempty_counts[column] += 1
                writer.writerow(row)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return row_count, updated_count, nonempty_counts


def gnomad_annotation(record: dict) -> dict[str, str]:
    af, source, *_rest = select_af_metrics(record)
    return {
        "gnomad_af": f"{af:.6g}" if af is not None else "",
        "gnomad_af_source": source or "",
        "gnomad_csq": str(record.get("consequence") or ""),
    }


def write_failures(path: Path, rows: list[dict[str, str]]) -> None:
    temporary = atomic_path(path, ".tsv.gz")
    try:
        with gzip.open(temporary, "wt", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=FAILURE_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_manifest(path: Path, manifest: dict) -> None:
    temporary = atomic_path(path, ".json")
    try:
        temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def complete_gnomad_annotation(
    run_dir: Path,
    source_annotation_dir: Path | None = None,
    outdir: Path | None = None,
    workers: int = 5,
    gnomad_cache_dir: Path | None = None,
) -> dict:
    if workers < 1:
        raise ValueError("--workers must be >= 1")
    run_dir = run_dir.expanduser().resolve()
    original_source = (source_annotation_dir or run_dir / "annotation").expanduser().resolve()
    outdir = (outdir or run_dir / "annotation_gnomad_complete").expanduser().resolve()

    existing = [path for name in REQUIRED_FILES if (path := outdir / name).exists()]
    if existing and len(existing) != len(REQUIRED_FILES):
        raise ValueError(f"Incomplete existing output directory: {outdir}")
    source = outdir if len(existing) == len(REQUIRED_FILES) else original_source
    for name in REQUIRED_FILES:
        path = source / name
        if not path.exists():
            raise FileNotFoundError(f"Missing annotation input: {path}")

    failure_rows = read_failures(source / "failures.tsv.gz")
    region_rows = [
        (row, region)
        for row in failure_rows
        if (region := failed_gnomad_region(row)) is not None
    ]
    regions = sorted({region for _row, region in region_rows})
    region_cache = GnomadRegionCache(
        gnomad_cache_dir,
        fetcher=fetch_region_variants_recursive,
    )
    successes, fetch_failures = fetch_failed_regions(regions, workers, region_cache)

    contexts = load_target_contexts(
        run_dir / "fetch" / "genes.tsv.gz",
        run_dir / "fetch" / "sequences" / "targets",
    )
    cache, key_status_counts = build_gnomad_cache(successes, contexts)

    retained_failures = []
    for row in failure_rows:
        region = failed_gnomad_region(row)
        if region is None or region not in successes:
            if region in fetch_failures:
                exc = fetch_failures[region]
                row = {
                    **row,
                    "failure_type": type(exc).__name__,
                    "message": str(exc),
                }
            retained_failures.append(row)

    outdir.mkdir(parents=True, exist_ok=True)
    row_count, updated_count, nonempty_counts = rewrite_annotations(
        source / "variant_annotations.tsv.gz",
        outdir / "variant_annotations.tsv.gz",
        cache,
    )
    write_failures(outdir / "failures.tsv.gz", retained_failures)

    source_manifest = json.loads((source / "manifest.json").read_text())
    recovered_failure_rows = len(region_rows) - sum(
        failed_gnomad_region(row) is not None for row in retained_failures
    )
    manifest = {
        **source_manifest,
        "created_at": utc_now(),
        "annotated_variant_context_count": row_count,
        "failure_count": len(retained_failures),
        "annotation_nonempty_counts": dict(nonempty_counts),
        "gnomad_region_success_count": int(
            source_manifest.get("gnomad_region_success_count") or 0
        )
        + recovered_failure_rows,
        "gnomad_region_failure_count": int(
            source_manifest.get("gnomad_region_failure_count") or 0
        )
        - recovered_failure_rows,
        "gnomad_raw_variant_count": int(source_manifest.get("gnomad_raw_variant_count") or 0)
        + sum(len(records) for records in successes.values()),
        "gnomad_cached_variant_count": int(
            source_manifest.get("gnomad_cached_variant_count") or 0
        )
        + len(cache),
        "gnomad_key_status_counts": dict(
            Counter(source_manifest.get("gnomad_key_status_counts") or {}) + key_status_counts
        ),
        "gnomad_completion": {
            "source_annotation_dir": str(original_source),
            "attempted_region_count": len(regions),
            "recovered_region_count": len(successes),
            "remaining_region_count": len(fetch_failures),
            "updated_variant_context_count": updated_count,
            "shared_cache": region_cache.snapshot(),
        },
    }
    write_manifest(outdir / "manifest.json", manifest)
    return manifest


def main() -> None:
    args = parse_args()
    manifest = complete_gnomad_annotation(
        run_dir=args.run_dir,
        source_annotation_dir=args.source_annotation_dir,
        outdir=args.outdir,
        workers=args.workers,
        gnomad_cache_dir=args.gnomad_cache_dir,
    )
    completion = manifest["gnomad_completion"]
    print(
        "gnomAD completion finished: "
        f"recovered={completion['recovered_region_count']}, "
        f"remaining={completion['remaining_region_count']}, "
        f"updated_rows={completion['updated_variant_context_count']}"
    )


if __name__ == "__main__":
    main()
