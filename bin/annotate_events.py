#!/usr/bin/env python3
"""Annotate TSV variants with ClinVar and gnomAD."""

import argparse
import csv
import json
import gzip
import logging
import os
import sys
import time
import concurrent.futures
from collections import Counter, defaultdict
from collections.abc import Sequence
from pathlib import Path

import pysam

from genomics.clinvar import review_stars as clinvar_review_stars
from genomics.gnomad import (
    GNOMAD_API_URL,
    GNOMAD_DATASET,
    fetch_region_variants_recursive,
    select_af_metrics,
)
from genomics.gnomad_cache import GnomadRegionCache
from genomics.variants import (
    CLINVAR_ANNOTATION_FIELDS,
    GNOMAD_ANNOTATION_FIELDS,
    add_context_normalized_record,
    build_context_index,
    event_vcf_key,
    load_target_contexts,
    refseq_accession_to_chrom,
    variant_aggregate_key,
    variant_key_text,
)


csv.field_size_limit(sys.maxsize)


logger = logging.getLogger(__name__)


def start_phase(name: str) -> float:
    logger.info("Starting phase %s", name)
    return time.perf_counter()


def finish_phase(timings: dict[str, float], name: str, started_at: float) -> None:
    elapsed = round(time.perf_counter() - started_at, 3)
    timings[name] = elapsed
    logger.info("Timing %s: %.3f seconds", name, elapsed)

VARIANT_ANNOTATION_FIELDS = [
    "variant_key",
    "gene_id",
    "event_type",
    "ref",
    "alt",
    "lookup_status",
    "strategies",
    *CLINVAR_ANNOTATION_FIELDS,
    *GNOMAD_ANNOTATION_FIELDS,
]

VARIANT_ANNOTATION_SHARD_SIZE = 250_000
VARIANT_ANNOTATION_SHARD_FORMAT = "tsv_gzip_v1"

EVENT_VARIANT_MAP_FIELDS = [
    "event_group_id",
    "variant_key",
    "normalization_status",
]

FAILURE_FIELDS = ["source", "scope", "chrom", "start", "end", "failure_type", "message"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-manifest", required=True, type=Path)
    parser.add_argument("--events-tsv", required=True, type=Path)
    parser.add_argument("--genes-tsv", required=True, type=Path)
    parser.add_argument("--target-sequences-dir", required=True, type=Path)
    parser.add_argument("--clinvar-vcf", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--partition-id", default="")
    parser.add_argument(
        "--gnomad-cache-dir",
        type=Path,
        default=os.environ.get("GAPH_GNOMAD_CACHE_DIR") or None,
        help="Optional shared directory for resumable gnomAD regional responses.",
    )
    return parser.parse_args()

def cluster_positions(positions: list[int], max_gap: int = 100000) -> list[tuple[int, int]]:
    if not positions: return []
    positions = sorted(positions)
    clusters = []
    start = positions[0]
    last = positions[0]
    for p in positions[1:]:
        if p - last > max_gap:
            clusters.append((start, last))
            start = p
        last = p
    clusters.append((start, last))
    return clusters


def manifest_count(manifest: dict, field: str) -> int:
    value = manifest.get(field)
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"Alignment manifest has invalid {field}: {value!r}")
    return value


def load_alignment_manifest(path: Path, partition_id: str) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Alignment manifest not found: {path}")
    manifest = json.loads(path.read_text())
    if manifest.get("stage") != "alignment":
        raise ValueError(f"Alignment manifest has invalid stage: {manifest.get('stage')!r}")
    if manifest.get("alignment_event_mode") != "compact_support":
        raise ValueError(
            "Annotation requires compact_support alignment events, observed "
            f"{manifest.get('alignment_event_mode')!r}"
        )
    observed_partition_id = str(manifest.get("partition_id") or "")
    if observed_partition_id != partition_id:
        raise ValueError(
            "Alignment manifest partition mismatch: "
            f"expected {partition_id!r}, observed {observed_partition_id!r}"
        )
    if not partition_id:
        raise ValueError("Annotation requires a partition-scoped alignment manifest")
    expected_schema = "normalized_alignment_evidence_partition_v2"
    if manifest.get("schema") != expected_schema:
        raise ValueError(
            "Alignment manifest schema mismatch: "
            f"expected {expected_schema!r}, observed {manifest.get('schema')!r}"
        )
    manifest_count(manifest, "alignment_event_count")
    return manifest


def write_tsv_gz(
    path: Path,
    fields: list[str],
    rows: list[dict],
    *,
    include_header: bool = True,
) -> int:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        if include_header:
            writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return len(rows)


def write_variant_annotation_shards(
    directory: Path,
    rows: list[dict],
    *,
    shard_size: int = VARIANT_ANNOTATION_SHARD_SIZE,
) -> dict[str, object]:
    """Write one ordered, headered shard set without duplicating the full table."""

    if shard_size < 1:
        raise ValueError("Variant annotation shard size must be >= 1")
    directory.mkdir(parents=True, exist_ok=True)
    if any(directory.iterdir()):
        raise ValueError(f"Variant annotation shard directory is not empty: {directory}")

    shards = []
    shard_count = max(1, (len(rows) + shard_size - 1) // shard_size)
    for index in range(shard_count):
        shard_id = f"shard_{index + 1:06d}"
        path = directory / f"{shard_id}.tsv.gz"
        shard_rows = rows[index * shard_size : (index + 1) * shard_size]
        write_tsv_gz(path, VARIANT_ANNOTATION_FIELDS, shard_rows)
        shards.append(
            {
                "shard_id": shard_id,
                "path": path.name,
                "row_count": len(shard_rows),
                "size_bytes": path.stat().st_size,
            }
        )
    return {
        "layout": "partitioned",
        "format": VARIANT_ANNOTATION_SHARD_FORMAT,
        "path": directory.name,
        "row_count": len(rows),
        "shard_size": shard_size,
        "shard_count": len(shards),
        "fields": VARIANT_ANNOTATION_FIELDS,
        "shards": shards,
    }


def int_or_default(value, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def event_variant_map_row(
    event_group_id: int,
    lookup_key: tuple[str, int, str, str] | None,
    normalization_status: str,
) -> dict[str, object]:
    return {
        "event_group_id": event_group_id,
        "variant_key": variant_key_text(lookup_key),
        "normalization_status": normalization_status,
    }


def path_metadata(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": stat.st_size,
        "mtime": int(stat.st_mtime),
    }


def failure_row(
    source: str,
    scope: str,
    chrom: str | None,
    start: int | str | None,
    end: int | str | None,
    failure_type: str,
    message: str,
) -> dict[str, object]:
    return {
        "source": source,
        "scope": scope,
        "chrom": chrom or "",
        "start": start if start is not None else "",
        "end": end if end is not None else "",
        "failure_type": failure_type,
        "message": message,
    }


def fetch_gnomad_for_cluster(
    region_cache: GnomadRegionCache,
    chrom: str,
    start: int,
    end: int,
) -> list[dict]:
    # pad by 100 bases
    return region_cache.fetch_region(chrom, max(1, start - 100), end + 100)


def format_info_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return "|".join(str(item) for item in value)
    return str(value)


def format_review_status_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ",".join(str(item) for item in value)
    return str(value)


def empty_annotation(columns: Sequence[str]) -> dict[str, str]:
    return {column: "" for column in columns}


def count_pipe_values(value: str) -> str:
    if not value:
        return "0"
    return str(len([item for item in value.split("|") if item]))


def format_float(value) -> str:
    return f"{value:.6g}" if value is not None else ""


def clinvar_annotation_from_record(rec) -> dict[str, str]:
    review_status = format_review_status_value(rec.info.get("CLNREVSTAT"))
    stars, stars_status = clinvar_review_stars(review_status)
    scv_accessions = format_info_value(rec.info.get("CLNSIGSCV"))
    return {
        "clinvar_sig": format_info_value(rec.info.get("CLNSIG")),
        "clinvar_revstat": review_status,
        "clinvar_review_stars": stars,
        "clinvar_review_stars_status": stars_status,
        "clinvar_id": rec.id or "",
        "clinvar_allele_id": format_info_value(rec.info.get("ALLELEID")),
        "clinvar_scv_count": count_pipe_values(scv_accessions),
        "clinvar_hgvs": format_info_value(rec.info.get("CLNHGVS")),
        "clinvar_disease": format_info_value(rec.info.get("CLNDN")),
        "clinvar_variant_type": format_info_value(rec.info.get("CLNVC")),
    }


def gnomad_annotation_from_variant(variant: dict) -> dict[str, str]:
    af, af_source, *_ = select_af_metrics(variant)
    return {
        "gnomad_af": format_float(af),
        "gnomad_af_source": af_source or "",
        "gnomad_csq": str(variant.get("consequence") or ""),
    }


def build_clinvar_cache(
    clinvar,
    accession_positions: dict[str, set[int]],
    contexts: dict[str, dict],
    context_index: dict[str, tuple[list[dict], list[int]]],
    failures: list[dict],
) -> tuple[dict[tuple[str, int, str, str], dict[str, str]], Counter]:
    cache = {}
    status_counts = Counter()
    for acc, positions in accession_positions.items():
        chrom = refseq_accession_to_chrom(acc)
        if not chrom:
            failures.append(
                failure_row(
                    "clinvar",
                    "accession",
                    "",
                    "",
                    "",
                    "unknown_chrom",
                    f"Could not map genomic accession to chromosome: {acc}",
                )
            )
            continue
        for start, end in cluster_positions(list(positions), max_gap=200000):
            try:
                for rec in clinvar.fetch(chrom, max(0, start - 1), end + 1):
                    rec_alts = rec.alts or ()
                    annotation = clinvar_annotation_from_record(rec)
                    for alt in rec_alts:
                        add_context_normalized_record(
                            cache,
                            (chrom, rec.pos, rec.ref, alt),
                            annotation,
                            contexts,
                            context_index,
                            status_counts,
                        )
            except ValueError as exc:
                message = f"ClinVar contig not found for chrom={chrom}; skipping region {start}-{end}."
                logger.warning(message)
                failures.append(
                    failure_row("clinvar", "region", chrom, start, end, "contig_not_found", str(exc))
                )
    if status_counts:
        logger.info(f"ClinVar key normalization status: {dict(status_counts)}")
    return cache, status_counts


def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    alignment_manifest = load_alignment_manifest(args.alignment_manifest, args.partition_id)
    timings_seconds: dict[str, float] = {}
    args.outdir.mkdir(parents=True, exist_ok=True)
    variant_annotations_dir = args.outdir / "variant_annotation_shards"
    event_variant_map_tsv = args.outdir / "event_variant_map.tsv.gz"
    failures_tsv = args.outdir / "failures.tsv.gz"
    manifest_json = args.outdir / "manifest.json"
    required_inputs = [
        args.genes_tsv,
        args.target_sequences_dir,
    ]
    for path in required_inputs:
        if not path.exists():
            raise FileNotFoundError(f"Required annotation input not found: {path}")
    if not args.clinvar_vcf.exists():
        raise FileNotFoundError(f"ClinVar VCF not found: {args.clinvar_vcf}")
    clinvar_tbi = Path(f"{args.clinvar_vcf}.tbi")
    if not clinvar_tbi.exists():
        raise FileNotFoundError(f"ClinVar VCF index not found: {clinvar_tbi}")

    failures: list[dict] = []
    phase_started = start_phase("load_target_context")
    contexts = load_target_contexts(args.genes_tsv, args.target_sequences_dir)
    logger.info("Loaded target context for %s gene(s).", len(contexts))
    context_index = build_context_index(contexts)
    finish_phase(timings_seconds, "load_target_context", phase_started)

    # 1. Read events once, collect lookup regions, and collapse repeated support rows.
    phase_started = start_phase("collapse_events")
    accession_positions = defaultdict(set)
    event_key_status_counts = Counter()
    variant_aggregates: dict[tuple, dict] = {}
    input_row_count = 0
    with gzip.open(event_variant_map_tsv, "wt", newline="") as map_handle:
        map_writer = csv.DictWriter(
            map_handle,
            fieldnames=EVENT_VARIANT_MAP_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        map_writer.writeheader()
        with gzip.open(args.events_tsv, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            required = {
                "event_group_id",
                "gene_id",
                "event_type",
                "target_start0",
                "genomic_accession",
                "genomic_start1",
                "ref",
                "alt",
                "strategy",
            }
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"Events table missing required columns: {', '.join(sorted(missing))}"
                )
            for row in reader:
                input_row_count += 1
                event_group_id = int_or_default(row.get("event_group_id"), 0)
                if event_group_id != input_row_count:
                    raise ValueError(
                        "Compact event_group_id values must be consecutive from 1; "
                        f"expected {input_row_count}, observed {event_group_id}"
                    )
                gene_id = str(row.get("gene_id") or "")
                if gene_id not in contexts:
                    raise ValueError(
                        f"Alignment event references gene {gene_id!r} outside the supplied target context"
                    )

                acc = row["genomic_accession"]
                lookup_key, status = event_vcf_key(row, contexts)
                event_key_status_counts[status] += 1
                map_writer.writerow(event_variant_map_row(event_group_id, lookup_key, status))
                if status == "non_concrete_allele":
                    continue
                raw_pos = int_or_default(row.get("genomic_start1"), -1)
                if acc and raw_pos > 0:
                    accession_positions[acc].add(raw_pos)
                if acc and lookup_key:
                    accession_positions[acc].add(int(lookup_key[1]))

                variant_key = variant_key_text(lookup_key)
                aggregate_key = variant_aggregate_key(row, variant_key)
                aggregate = variant_aggregates.get(aggregate_key)
                if aggregate is None:
                    aggregate = {
                        "variant_key": variant_key,
                        "gene_id": row.get("gene_id", ""),
                        "event_type": row.get("event_type", ""),
                        "target_start0": row.get("target_start0", ""),
                        "ref": row.get("ref", ""),
                        "alt": row.get("alt", ""),
                        "lookup_status": status,
                        "_lookup_key": lookup_key,
                        "_strategies": set(),
                    }
                    variant_aggregates[aggregate_key] = aggregate
                strategy = str(row.get("strategy") or "")
                if not strategy:
                    raise ValueError("Alignment event requires one strategy")
                aggregate["_strategies"].add(strategy)
    expected_event_count = manifest_count(alignment_manifest, "alignment_event_count")
    if input_row_count != expected_event_count:
        raise ValueError(
            "Alignment event row count does not match alignment manifest: "
            f"rows={input_row_count}, manifest={expected_event_count}"
        )
    logger.info(f"Event key normalization status: {dict(event_key_status_counts)}")
    logger.info(f"Collapsed {input_row_count} event row(s) to {len(variant_aggregates)} variant-context row(s).")
    finish_phase(timings_seconds, "collapse_events", phase_started)

    # 2. Determine gnomAD clusters
    phase_started = start_phase("gnomad_lookup")
    gnomad_tasks = []
    for acc, positions in accession_positions.items():
        chrom = refseq_accession_to_chrom(acc)
        if not chrom:
            continue
        clusters = cluster_positions(list(positions), max_gap=200000)
        for c_start, c_end in clusters:
            gnomad_tasks.append((chrom, c_start, c_end))

    logger.info(f"Will fetch {len(gnomad_tasks)} region(s) from gnomAD API.")

    # 3. Fetch gnomAD in parallel and cache
    gnomad_region_cache = GnomadRegionCache(
        args.gnomad_cache_dir,
        fetcher=fetch_region_variants_recursive,
    )
    gnomad_cache = {}
    gnomad_key_status_counts = Counter()
    gnomad_region_success_count = 0
    gnomad_raw_variant_count = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        future_to_task = {
            executor.submit(fetch_gnomad_for_cluster, gnomad_region_cache, chrom, start, end): (
                chrom,
                start,
                end,
            )
            for chrom, start, end in gnomad_tasks
        }
        for future in concurrent.futures.as_completed(future_to_task):
            task = future_to_task[future]
            try:
                vars_list = future.result()
                gnomad_region_success_count += 1
                gnomad_raw_variant_count += len(vars_list)
                for v in vars_list:
                    key = (v["chrom"], int(v["pos"]), v["ref"], v["alt"])
                    add_context_normalized_record(
                        gnomad_cache,
                        key,
                        v,
                        contexts,
                        context_index,
                        gnomad_key_status_counts,
                    )
            except Exception as exc:
                logger.error(f"gnomAD fetch failed for {task}: {exc}")
                chrom, start, end = task
                failures.append(
                    failure_row(
                        "gnomad",
                        "region",
                        chrom,
                        start,
                        end,
                        type(exc).__name__,
                        str(exc),
                    )
                )

    logger.info(f"Cached {len(gnomad_cache)} gnomAD variants.")
    if gnomad_key_status_counts:
        logger.info(f"gnomAD key normalization status: {dict(gnomad_key_status_counts)}")
    finish_phase(timings_seconds, "gnomad_lookup", phase_started)

    # 4. Open ClinVar
    phase_started = start_phase("clinvar_lookup")
    clinvar = pysam.VariantFile(str(args.clinvar_vcf))
    clinvar_cache, clinvar_key_status_counts = build_clinvar_cache(
        clinvar,
        accession_positions,
        contexts,
        context_index,
        failures,
    )
    logger.info(f"Cached {len(clinvar_cache)} ClinVar variants.")
    finish_phase(timings_seconds, "clinvar_lookup", phase_started)

    # 5. Annotate unique variant-context rows.
    phase_started = start_phase("write_variant_annotations")
    variant_rows: list[dict] = []
    for aggregate in variant_aggregates.values():
        lookup_key = aggregate["_lookup_key"]
        clinvar_annotation = empty_annotation(CLINVAR_ANNOTATION_FIELDS)
        gnomad_annotation = empty_annotation(GNOMAD_ANNOTATION_FIELDS)

        if lookup_key:
            clinvar_annotation = clinvar_cache.get(lookup_key, clinvar_annotation)
        if lookup_key and lookup_key in gnomad_cache:
            gnomad_annotation = gnomad_annotation_from_variant(gnomad_cache[lookup_key])

        row = {field: aggregate.get(field, "") for field in VARIANT_ANNOTATION_FIELDS}
        row["strategies"] = ",".join(sorted(aggregate["_strategies"]))
        row.update(clinvar_annotation)
        row.update(gnomad_annotation)
        variant_rows.append(row)

    variant_rows.sort(
        key=lambda row: (
            int_or_default(row.get("gene_id"), 10**18),
            row.get("variant_key", ""),
            row.get("event_type", ""),
        )
    )
    variant_annotations = write_variant_annotation_shards(
        variant_annotations_dir,
        variant_rows,
    )
    output_row_count = int(variant_annotations["row_count"])
    finish_phase(timings_seconds, "write_variant_annotations", phase_started)

    phase_started = start_phase("write_failures")
    failure_count = write_tsv_gz(failures_tsv, FAILURE_FIELDS, failures)
    finish_phase(timings_seconds, "write_failures", phase_started)
    manifest = {
        "stage": "annotation",
        "schema": "normalized_annotation_evidence_partition_v2",
        "partition_id": args.partition_id,
        "event_row_count": input_row_count,
        "event_variant_map_count": input_row_count,
        "variant_context_count": len(variant_aggregates),
        "annotated_variant_context_count": output_row_count,
        "variant_annotations": variant_annotations,
        "clinvar_vcf": path_metadata(args.clinvar_vcf),
        "clinvar_tbi": path_metadata(clinvar_tbi),
        "clinvar_cached_variant_count": len(clinvar_cache),
        "gnomad_api_url": GNOMAD_API_URL,
        "gnomad_dataset": GNOMAD_DATASET,
        "gnomad_region_count": len(gnomad_tasks),
        "gnomad_region_success_count": gnomad_region_success_count,
        "gnomad_region_failure_count": len(gnomad_tasks) - gnomad_region_success_count,
        "gnomad_raw_variant_count": gnomad_raw_variant_count,
        "gnomad_cached_variant_count": len(gnomad_cache),
        "gnomad_shared_cache": gnomad_region_cache.snapshot(),
        "failure_count": failure_count,
        "event_key_status_counts": dict(event_key_status_counts),
        "gnomad_key_status_counts": dict(gnomad_key_status_counts),
        "clinvar_key_status_counts": dict(clinvar_key_status_counts),
        "timings_seconds": timings_seconds,
    }
    manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    logger.info(
        "Saved %s variant annotation row(s) in %s shard(s) to %s",
        output_row_count,
        variant_annotations["shard_count"],
        variant_annotations_dir,
    )
    logger.info(f"Saved event-to-variant lineage to {event_variant_map_tsv}")
    logger.info(f"Saved annotation failures to {failures_tsv}")
    logger.info(f"Saved annotation manifest to {manifest_json}")

if __name__ == "__main__":
    main()
