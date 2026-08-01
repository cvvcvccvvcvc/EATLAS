#!/usr/bin/env python3
"""Annotate TSV variants with ClinVar and gnomAD."""

import argparse
import csv
import json
import gzip
import logging
import os
import sys
import concurrent.futures
from collections.abc import Iterable
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pysam

if __package__ in {None, ""}:
    runtime_root = Path.cwd()
    if not (runtime_root / "genomics").is_dir():
        runtime_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(runtime_root))

from feature_coverage import load_snv_site_depth, site_aligned_ortholog_counts
from genomics.clinvar import review_stars as clinvar_review_stars
from genomics.gnomad import GNOMAD_API_URL, fetch_region_variants_recursive, select_af_metrics
from genomics.gnomad_cache import GnomadRegionCache
from genomics.variants import (
    add_context_normalized_record,
    build_context_index,
    event_vcf_key,
    load_target_contexts,
    normalize_chrom,
    refseq_accession_to_chrom,
    variant_key_text,
)
from ortholog_evidence_summary import write_ortholog_evidence_summary

logger = logging.getLogger(__name__)

CLINVAR_COLUMNS = [
    "clinvar_sig",
    "clinvar_revstat",
    "clinvar_review_stars",
    "clinvar_review_stars_status",
    "clinvar_id",
    "clinvar_allele_id",
    "clinvar_scv_count",
    "clinvar_hgvs",
    "clinvar_disease",
    "clinvar_variant_type",
]

GNOMAD_COLUMNS = [
    "gnomad_af",
    "gnomad_af_source",
    "gnomad_csq",
]

ANNOTATION_COLUMNS = CLINVAR_COLUMNS + GNOMAD_COLUMNS

VARIANT_ANNOTATION_FIELDS = [
    "variant_key",
    "gene_id",
    "event_type",
    "ref",
    "alt",
    "lookup_status",
    "support_row_count",
    "support_ortholog_count",
    "strategies",
    *ANNOTATION_COLUMNS,
]

VARIANT_STRATEGY_SUPPORT_FIELDS = [
    "variant_key",
    "gene_id",
    "strategy",
    "alt_support_row_count",
    "alt_support_ortholog_count",
    "site_aligned_ortholog_count",
]

VARIANT_ORTHOLOG_SUPPORT_FIELDS = [
    "variant_key",
    "gene_id",
    "strategy",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "support_row_count",
]

FAILURE_FIELDS = ["source", "scope", "chrom", "start", "end", "failure_type", "message"]
GNOMAD_DATASET = "gnomad_r4"

@dataclass(slots=True)
class StrategySupport:
    row_count: int = 0
    orthologs: set[str] = field(default_factory=set)
    ortholog_count_hint: int = 0


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--events-tsv", required=True, type=Path)
    parser.add_argument("--event-ortholog-support-tsv", type=Path)
    depth_input = parser.add_mutually_exclusive_group(required=True)
    depth_input.add_argument("--segments-tsv", type=Path)
    depth_input.add_argument("--snv-site-depth-tsv", type=Path)
    parser.add_argument("--snv-taxonomic-depth-tsv", type=Path)
    parser.add_argument("--snv-alt-taxonomic-support-tsv", type=Path)
    parser.add_argument("--genes-tsv", required=False, type=Path)
    parser.add_argument("--target-sequences-dir", required=False, type=Path)
    parser.add_argument("--target-features-dir", type=Path)
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


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict]) -> int:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return len(rows)


def split_values(value: str | None) -> set[str]:
    if not value:
        return set()
    return {item for item in str(value).replace("|", ",").split(",") if item}


def int_or_default(value, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def add_strategy_support(aggregate: dict, row: dict[str, str]) -> None:
    singular_strategies = split_values(row.get("strategy"))
    strategies = singular_strategies | split_values(row.get("strategies"))
    if not strategies:
        raise ValueError("Event row does not identify an alignment strategy")
    if len(strategies) > 1 and not singular_strategies:
        raise ValueError(
            "Event row contains support aggregated across multiple strategies; "
            "per-strategy support counts cannot be recovered"
        )

    support_row_count = int_or_default(row.get("support_row_count"), 1)
    ortholog_gene_id = row.get("ortholog_gene_id", "")
    ortholog_count_hint = int_or_default(row.get("support_ortholog_count"), 0)
    support_by_strategy = aggregate["_support_by_strategy"]
    for strategy in strategies:
        support = support_by_strategy.get(strategy)
        if support is None:
            support = StrategySupport()
            support_by_strategy[strategy] = support
        support.row_count += support_row_count
        if ortholog_gene_id:
            support.orthologs.add(ortholog_gene_id)
        else:
            support.ortholog_count_hint += ortholog_count_hint


def add_ortholog_support(aggregate: dict, row: dict[str, str]) -> None:
    strategy = str(row.get("strategy") or "")
    ortholog_gene_id = str(row.get("ortholog_gene_id") or "")
    if not strategy or not ortholog_gene_id:
        raise ValueError("Ortholog support row requires strategy and ortholog_gene_id")
    support_row_count = int_or_default(row.get("support_row_count"), 1)
    if support_row_count < 1:
        raise ValueError("Ortholog support_row_count must be positive")

    support_by_ortholog = aggregate["_ortholog_support"]
    key = (strategy, ortholog_gene_id)
    support = support_by_ortholog.get(key)
    if support is None:
        support_by_ortholog[key] = {
            "strategy": strategy,
            "ortholog_gene_id": ortholog_gene_id,
            "tax_id": str(row.get("tax_id") or ""),
            "taxname": str(row.get("taxname") or ""),
            "support_row_count": support_row_count,
        }
        return

    for field in ("tax_id", "taxname"):
        observed = str(row.get(field) or "")
        current = str(support.get(field) or "")
        if current and observed and current != observed:
            raise ValueError(
                f"Conflicting {field} for strategy={strategy}, "
                f"ortholog_gene_id={ortholog_gene_id}: {current!r} != {observed!r}"
            )
        if not current and observed:
            support[field] = observed
    support["support_row_count"] += support_row_count


def build_variant_ortholog_support(
    aggregates: Iterable[dict],
) -> tuple[list[dict[str, object]], int]:
    rows: list[dict[str, object]] = []
    missing_key_count = 0
    for aggregate in aggregates:
        variant_key = aggregate.get("variant_key", "")
        ortholog_support = aggregate["_ortholog_support"]
        if not variant_key:
            missing_key_count += len(ortholog_support)
            continue
        for support in ortholog_support.values():
            rows.append(
                {
                    "variant_key": variant_key,
                    "gene_id": aggregate.get("gene_id", ""),
                    **support,
                }
            )
    rows.sort(
        key=lambda row: (
            int_or_default(row.get("gene_id"), 10**18),
            row["variant_key"],
            row["strategy"],
            row["ortholog_gene_id"],
        )
    )
    return rows, missing_key_count


def validate_ortholog_support_totals(
    strategy_rows: Iterable[dict[str, object]],
    ortholog_rows: Iterable[dict[str, object]],
) -> None:
    observed: dict[tuple[str, str, str], tuple[int, int]] = {}
    for row in ortholog_rows:
        key = (str(row["variant_key"]), str(row["gene_id"]), str(row["strategy"]))
        ortholog_count, row_count = observed.get(key, (0, 0))
        observed[key] = (
            ortholog_count + 1,
            row_count + int(row["support_row_count"]),
        )
    for row in strategy_rows:
        key = (str(row["variant_key"]), str(row["gene_id"]), str(row["strategy"]))
        actual_orthologs, actual_rows = observed.get(key, (0, 0))
        expected_orthologs = int(row["alt_support_ortholog_count"])
        expected_rows = int(row["alt_support_row_count"])
        if (actual_orthologs, actual_rows) != (expected_orthologs, expected_rows):
            raise ValueError(
                "Variant ortholog support does not match strategy totals for "
                f"{key}: orthologs={actual_orthologs}/{expected_orthologs}, "
                f"rows={actual_rows}/{expected_rows}"
            )


def build_variant_strategy_support(
    aggregates: Iterable[dict],
    site_depths: dict[tuple[str, str, int], int] | None = None,
) -> tuple[list[dict[str, object]], int]:
    site_depths = site_depths or {}
    rows: list[dict[str, object]] = []
    missing_key_count = 0
    for aggregate in aggregates:
        variant_key = aggregate.get("variant_key", "")
        support_by_strategy = aggregate["_support_by_strategy"]
        if not variant_key:
            missing_key_count += len(support_by_strategy)
            continue
        for strategy, support in support_by_strategy.items():
            alt_support_count = max(
                len(support.orthologs),
                support.ortholog_count_hint,
            )
            site_depth: int | str = ""
            if aggregate.get("event_type") == "snv":
                depth_key = (
                    str(aggregate.get("gene_id") or ""),
                    strategy,
                    int(aggregate.get("target_start0") or 0),
                )
                if depth_key not in site_depths:
                    raise ValueError(f"Missing site ortholog depth for SNV {depth_key}")
                site_depth = site_depths[depth_key]
                if alt_support_count > site_depth:
                    raise ValueError(
                        "ALT-support ortholog count exceeds site-aligned ortholog count for "
                        f"{depth_key}: {alt_support_count} > {site_depth}"
                    )
            rows.append(
                {
                    "variant_key": variant_key,
                    "gene_id": aggregate.get("gene_id", ""),
                    "strategy": strategy,
                    "alt_support_row_count": support.row_count,
                    "alt_support_ortholog_count": alt_support_count,
                    "site_aligned_ortholog_count": site_depth,
                }
            )
    rows.sort(
        key=lambda row: (
            int_or_default(row.get("gene_id"), 10**18),
            row["variant_key"],
            row["strategy"],
        )
    )
    return rows, missing_key_count


def iter_variant_strategy_snv_sites(
    aggregates: Iterable[dict],
) -> Iterable[dict[str, object]]:
    for aggregate in aggregates:
        if aggregate.get("event_type") != "snv" or not aggregate.get("variant_key"):
            continue
        for strategy in aggregate["_support_by_strategy"]:
            yield {
                "gene_id": aggregate.get("gene_id", ""),
                "strategy": strategy,
                "target_start0": aggregate.get("target_start0", ""),
            }


def variant_aggregate_key(row: dict[str, str], variant_key: str) -> tuple:
    gene_id = row.get("gene_id", "")
    if variant_key:
        return "canonical", gene_id, variant_key
    return (
        "raw",
        gene_id,
        row.get("event_type", ""),
        row.get("target_start0", ""),
        row.get("target_end0", ""),
        row.get("genomic_accession", ""),
        row.get("genomic_start1", ""),
        row.get("genomic_end1", ""),
        row.get("ref", ""),
        row.get("alt", ""),
    )


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


def empty_annotation(columns: list[str]) -> dict[str, str]:
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


def build_gnomad_statuses(
    aggregates: Iterable[dict],
    gnomad_cache: dict[tuple[str, int, str, str], dict],
    failures: Iterable[dict],
) -> dict[tuple[str, int, str, str], str]:
    failed_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for failure in failures:
        if failure.get("source") != "gnomad" or failure.get("scope") != "region":
            continue
        chrom = normalize_chrom(str(failure.get("chrom") or ""))
        try:
            start = int(failure.get("start") or 0)
            end = int(failure.get("end") or 0)
        except (TypeError, ValueError):
            continue
        if chrom and start > 0 and end >= start:
            failed_by_chrom[chrom].append((start, end))
    for chrom in failed_by_chrom:
        failed_by_chrom[chrom].sort()

    statuses: dict[tuple[str, int, str, str], str] = {}
    for aggregate in aggregates:
        if aggregate.get("event_type") != "snv":
            continue
        target_key = (
            str(aggregate.get("gene_id") or ""),
            int(aggregate.get("target_start0") or 0),
            str(aggregate.get("ref") or "").upper(),
            str(aggregate.get("alt") or "").upper(),
        )
        lookup_key = aggregate.get("_lookup_key")
        found = False
        if lookup_key in gnomad_cache:
            found = bool(gnomad_annotation_from_variant(gnomad_cache[lookup_key])["gnomad_af"])
        if found:
            status = "found"
        elif aggregate.get("lookup_status") != "ok" or lookup_key is None:
            status = "lookup_failed"
        else:
            chrom, position, _ref, _alt = lookup_key
            status = "not_found"
            for start, end in failed_by_chrom.get(normalize_chrom(chrom) or "", []):
                if start <= position <= end:
                    status = "lookup_failed"
                    break
                if start > position:
                    break
        previous = statuses.setdefault(target_key, status)
        if previous != status:
            raise ValueError(f"Conflicting gnomAD status for target SNV {target_key}")
    return statuses


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
    args.outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = args.outdir / "variant_annotations.tsv.gz"
    support_tsv = args.outdir / "variant_strategy_support.tsv.gz"
    ortholog_support_tsv = args.outdir / "variant_ortholog_support.tsv.gz"
    ortholog_evidence_tsv = args.outdir / "ortholog_evidence_summary.tsv.gz"
    failures_tsv = args.outdir / "failures.tsv.gz"
    manifest_json = args.outdir / "manifest.json"
    if args.segments_tsv is not None and not args.segments_tsv.exists():
        raise FileNotFoundError(f"Alignment segments TSV not found: {args.segments_tsv}")
    if args.snv_site_depth_tsv is not None and not args.snv_site_depth_tsv.exists():
        raise FileNotFoundError(f"SNV site-depth TSV not found: {args.snv_site_depth_tsv}")
    if args.event_ortholog_support_tsv is not None and not args.event_ortholog_support_tsv.exists():
        raise FileNotFoundError(
            f"Event ortholog support TSV not found: {args.event_ortholog_support_tsv}"
        )
    taxonomic_inputs = [
        args.snv_taxonomic_depth_tsv,
        args.snv_alt_taxonomic_support_tsv,
        args.target_features_dir,
    ]
    if any(path is not None for path in taxonomic_inputs) and not all(
        path is not None for path in taxonomic_inputs
    ):
        raise ValueError(
            "Taxonomic ortholog evidence requires site depth, ALT support, and target features"
        )
    for path in taxonomic_inputs:
        if path is not None and not path.exists():
            raise FileNotFoundError(f"Taxonomic ortholog evidence input not found: {path}")
    if not args.clinvar_vcf.exists():
        raise FileNotFoundError(f"ClinVar VCF not found: {args.clinvar_vcf}")
    clinvar_tbi = Path(f"{args.clinvar_vcf}.tbi")
    if not clinvar_tbi.exists():
        raise FileNotFoundError(f"ClinVar VCF index not found: {clinvar_tbi}")

    failures: list[dict] = []
    if bool(args.genes_tsv) != bool(args.target_sequences_dir):
        raise ValueError("--genes-tsv and --target-sequences-dir must be provided together.")
    if args.genes_tsv and args.target_sequences_dir:
        contexts = load_target_contexts(args.genes_tsv, args.target_sequences_dir)
        logger.info("Loaded target context for %s gene(s).", len(contexts))
    else:
        contexts = {}
        logger.warning("No target context provided; ClinVar/gnomAD lookup will use raw event keys.")
    context_index = build_context_index(contexts)

    # 1. Read events once, collect lookup regions, and collapse repeated support rows.
    accession_positions = defaultdict(set)
    event_key_status_counts = Counter()
    unique_lookup_status_counts = Counter()
    variant_aggregates: dict[tuple, dict] = {}
    input_row_count = 0
    with gzip.open(args.events_tsv, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames
        required = {"gene_id", "event_type", "target_start0", "genomic_accession", "genomic_start1", "ref", "alt"}
        missing = required - set(header or [])
        if missing:
            raise ValueError(f"Events table missing required columns: {', '.join(sorted(missing))}")
        if not {"strategy", "strategies"} & set(header or []):
            raise ValueError("Events table must include strategy or strategies")
        events_have_ortholog_identity = "ortholog_gene_id" in set(header or [])
        if not events_have_ortholog_identity and args.event_ortholog_support_tsv is None:
            raise ValueError(
                "Compact events require --event-ortholog-support-tsv to publish exact supporters"
            )
        for row in reader:
            input_row_count += 1
            acc = row["genomic_accession"]
            lookup_key, status = event_vcf_key(row, contexts)
            event_key_status_counts[status] += 1
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
                    "support_row_count": 0,
                    "_lookup_key": lookup_key,
                    "_support_by_strategy": {},
                    "_ortholog_support": {},
                }
                variant_aggregates[aggregate_key] = aggregate
                unique_lookup_status_counts[status] += 1

            aggregate["support_row_count"] += int_or_default(row.get("support_row_count"), 1)
            add_strategy_support(aggregate, row)
            if events_have_ortholog_identity and args.event_ortholog_support_tsv is None:
                add_ortholog_support(aggregate, row)
    logger.info(f"Event key normalization status: {dict(event_key_status_counts)}")
    logger.info(f"Collapsed {input_row_count} event row(s) to {len(variant_aggregates)} variant-context row(s).")

    if args.event_ortholog_support_tsv is not None:
        with open_text(args.event_ortholog_support_tsv) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "gene_id",
                "event_type",
                "target_start0",
                "genomic_accession",
                "genomic_start1",
                "ref",
                "alt",
                "strategy",
                "ortholog_gene_id",
                "tax_id",
                "taxname",
            }
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    "Event ortholog support table missing required columns: "
                    + ", ".join(sorted(missing))
                )
            for row in reader:
                lookup_key, status = event_vcf_key(row, contexts)
                if status == "non_concrete_allele":
                    continue
                variant_key = variant_key_text(lookup_key)
                aggregate_key = variant_aggregate_key(row, variant_key)
                aggregate = variant_aggregates.get(aggregate_key)
                if aggregate is None:
                    raise ValueError(
                        "Event ortholog support row has no matching aggregate event: "
                        f"gene_id={row.get('gene_id', '')}, strategy={row.get('strategy', '')}, "
                        f"variant_key={variant_key}"
                    )
                add_ortholog_support(aggregate, row)

    if args.snv_site_depth_tsv is not None:
        site_depths = load_snv_site_depth(args.snv_site_depth_tsv)
    else:
        site_depths = site_aligned_ortholog_counts(
            args.segments_tsv,
            iter_variant_strategy_snv_sites(variant_aggregates.values()),
            args.outdir,
        )
    logger.info(f"Calculated site-aligned ortholog depth for {len(site_depths)} variant-strategy SNV(s).")

    # 2. Determine gnomAD clusters
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

    # 4. Open ClinVar
    clinvar = pysam.VariantFile(str(args.clinvar_vcf))
    clinvar_cache, clinvar_key_status_counts = build_clinvar_cache(
        clinvar,
        accession_positions,
        contexts,
        context_index,
        failures,
    )
    logger.info(f"Cached {len(clinvar_cache)} ClinVar variants.")

    # 5. Annotate unique variant-context rows.
    annotation_value_counts = Counter()
    variant_rows: list[dict] = []
    for aggregate in variant_aggregates.values():
        lookup_key = aggregate["_lookup_key"]
        clinvar_annotation = empty_annotation(CLINVAR_COLUMNS)
        gnomad_annotation = empty_annotation(GNOMAD_COLUMNS)

        if lookup_key:
            clinvar_annotation = clinvar_cache.get(lookup_key, clinvar_annotation)
        if lookup_key and lookup_key in gnomad_cache:
            gnomad_annotation = gnomad_annotation_from_variant(gnomad_cache[lookup_key])

        row = {field: aggregate.get(field, "") for field in VARIANT_ANNOTATION_FIELDS}
        support_by_strategy = aggregate["_support_by_strategy"]
        orthologs = set()
        ortholog_count_hint = 0
        for support in support_by_strategy.values():
            orthologs.update(support.orthologs)
            ortholog_count_hint += support.ortholog_count_hint
        row["support_ortholog_count"] = max(
            len(orthologs),
            ortholog_count_hint,
        )
        row["strategies"] = ",".join(sorted(support_by_strategy))
        row.update(clinvar_annotation)
        row.update(gnomad_annotation)
        for column in ANNOTATION_COLUMNS:
            if row[column]:
                annotation_value_counts[column] += 1
        variant_rows.append(row)

    variant_rows.sort(
        key=lambda row: (
            int_or_default(row.get("gene_id"), 10**18),
            row.get("variant_key", ""),
            row.get("event_type", ""),
        )
    )
    output_row_count = write_tsv_gz(out_tsv, VARIANT_ANNOTATION_FIELDS, variant_rows)
    strategy_support_rows, strategy_support_missing_key_count = build_variant_strategy_support(
        variant_aggregates.values(),
        site_depths,
    )
    strategy_support_count = write_tsv_gz(
        support_tsv,
        VARIANT_STRATEGY_SUPPORT_FIELDS,
        strategy_support_rows,
    )
    ortholog_support_rows, ortholog_support_missing_key_count = build_variant_ortholog_support(
        variant_aggregates.values()
    )
    validate_ortholog_support_totals(strategy_support_rows, ortholog_support_rows)
    ortholog_support_count = write_tsv_gz(
        ortholog_support_tsv,
        VARIANT_ORTHOLOG_SUPPORT_FIELDS,
        ortholog_support_rows,
    )
    ortholog_evidence_summary_count = 0
    if args.snv_taxonomic_depth_tsv is not None:
        target_feature_paths = sorted(args.target_features_dir.glob("*.tsv.gz"))
        if not target_feature_paths:
            raise ValueError(f"No target feature tables found in {args.target_features_dir}")
        ortholog_evidence_summary_count = write_ortholog_evidence_summary(
            args.snv_taxonomic_depth_tsv,
            args.snv_alt_taxonomic_support_tsv,
            target_feature_paths,
            build_gnomad_statuses(variant_aggregates.values(), gnomad_cache, failures),
            ortholog_evidence_tsv,
        )

    failure_count = write_tsv_gz(failures_tsv, FAILURE_FIELDS, failures)
    manifest = {
        "output_mode": "unique_variant_context",
        "partition_id": args.partition_id,
        "event_row_count": input_row_count,
        "excluded_non_concrete_event_count": event_key_status_counts["non_concrete_allele"],
        "variant_context_count": len(variant_aggregates),
        "annotated_variant_context_count": output_row_count,
        "variant_strategy_support_count": strategy_support_count,
        "variant_strategy_support_missing_key_count": strategy_support_missing_key_count,
        "variant_ortholog_support_count": ortholog_support_count,
        "variant_ortholog_support_missing_key_count": ortholog_support_missing_key_count,
        "variant_strategy_site_depth_count": len(site_depths),
        "ortholog_evidence_summary_count": ortholog_evidence_summary_count,
        "target_context_count": len(contexts),
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
        "annotation_nonempty_counts": dict(annotation_value_counts),
        "event_key_status_counts": dict(event_key_status_counts),
        "unique_lookup_status_counts": dict(unique_lookup_status_counts),
        "gnomad_key_status_counts": dict(gnomad_key_status_counts),
        "clinvar_key_status_counts": dict(clinvar_key_status_counts),
    }
    manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    logger.info(f"Saved variant annotations to {out_tsv}")
    logger.info(f"Saved variant-strategy support to {support_tsv}")
    logger.info(f"Saved variant-ortholog support to {ortholog_support_tsv}")
    logger.info(f"Saved annotation failures to {failures_tsv}")
    logger.info(f"Saved annotation manifest to {manifest_json}")

if __name__ == "__main__":
    main()
