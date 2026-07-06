#!/usr/bin/env python3
"""Annotate TSV variants with ClinVar and gnomAD."""

import argparse
import bisect
import csv
import json
import gzip
import logging
import sys
import concurrent.futures
from collections import Counter, defaultdict
from pathlib import Path

# Add bin to path so we can import fetch_gnomad_variants
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from fetch_gnomad_variants import GNOMAD_API_URL, fetch_region_variants_recursive, _select_af_metrics
except ImportError as e:
    print(f"Error importing fetch_gnomad_variants: {e}")
    sys.exit(1)

try:
    import pysam
except ImportError:
    print("pysam is required but not installed.")
    sys.exit(1)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

CLINVAR_COLUMNS = [
    "clinvar_sig",
    "clinvar_revstat",
    "clinvar_review_stars",
    "clinvar_review_stars_status",
    "clinvar_id",
    "clinvar_allele_id",
    "clinvar_sig_conflict",
    "clinvar_scv_count",
    "clinvar_scv_accessions",
    "clinvar_hgvs",
    "clinvar_geneinfo",
    "clinvar_disease",
    "clinvar_disease_db",
    "clinvar_variant_type",
    "clinvar_variant_type_so",
    "clinvar_origin",
    "clinvar_rs",
]

GNOMAD_COLUMNS = [
    "gnomad_af",
    "gnomad_af_source",
    "gnomad_csq",
    "gnomad_variant_id",
    "gnomad_af_exome",
    "gnomad_af_genome",
    "gnomad_af_joint",
    "gnomad_ac_joint",
    "gnomad_an_joint",
    "gnomad_hgvsc",
    "gnomad_hgvsp",
]

ANNOTATION_COLUMNS = CLINVAR_COLUMNS + GNOMAD_COLUMNS

VARIANT_ANNOTATION_FIELDS = [
    "variant_key",
    "gene_id",
    "event_type",
    "target_start0",
    "target_end0",
    "genomic_accession",
    "genomic_start1",
    "genomic_end1",
    "ref",
    "alt",
    "lookup_chrom",
    "lookup_pos",
    "lookup_ref",
    "lookup_alt",
    "lookup_status",
    "support_row_count",
    "support_ortholog_count",
    "support_strategy_count",
    "strategies",
    "tools",
    "presets",
    "tax_id_count",
    "taxname_count",
    *ANNOTATION_COLUMNS,
]

FAILURE_FIELDS = ["source", "scope", "chrom", "start", "end", "failure_type", "message"]
GNOMAD_DATASET = "gnomad_r4"

CLINVAR_REVIEW_STARS = {
    "practice_guideline": "4",
    "reviewed_by_expert_panel": "3",
    "criteria_provided,_multiple_submitters,_no_conflicts": "2",
    "criteria_provided,_multiple_submitters": "2",
    "criteria_provided,_conflicting_classifications": "1",
    "criteria_provided,_conflicting_interpretations": "1",
    "criteria_provided,_single_submitter": "1",
    "no_assertion_criteria_provided": "0",
    "no_assertion_provided": "0",
    "no_classification_provided": "0",
    "no_classification_for_the_individual_variant": "0",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--events-tsv", required=True, type=Path)
    parser.add_argument("--genes-tsv", required=False, type=Path)
    parser.add_argument("--target-sequences-dir", required=False, type=Path)
    parser.add_argument("--clinvar-vcf", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()

def _refseq_accession_to_gnomad_chrom(chr_acc: str) -> str | None:
    if not chr_acc: return None
    c = str(chr_acc).strip()
    if c.startswith("chr"): c = c[3:]
    if c in {"X", "Y", "MT", "M"}: return "MT" if c == "M" else c
    if c.isdigit():
        val = int(c)
        if val == 23: return "X"
        if val == 24: return "Y"
        return str(val)
    import re
    match = re.search(r"NC_0+(\d+)\.", c)
    if not match: return None
    num = int(match.group(1))
    if num == 23: return "X"
    if num == 24: return "Y"
    if num in {12920, 1807}: return "MT"
    return str(num)

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


def lookup_key_text(key: tuple[str, int, str, str] | None) -> str:
    if not key:
        return ""
    chrom, pos, ref, alt = key
    return f"{chrom}:{pos}:{ref}>{alt}"


def int_or_default(value, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


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


def read_fasta_sequence(path: Path) -> str:
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            chunks.append(line.strip())
    return "".join(chunks).upper()


def normalize_chrom(value: str | None) -> str | None:
    if not value:
        return None
    c = str(value).strip()
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        return "MT"
    return c


def load_target_contexts(genes_tsv: Path | None, target_sequences_dir: Path | None) -> dict[str, dict]:
    if not genes_tsv and not target_sequences_dir:
        logger.warning("No target context provided; ClinVar/gnomAD lookup will use raw event keys.")
        return {}
    if not genes_tsv or not target_sequences_dir:
        raise ValueError("--genes-tsv and --target-sequences-dir must be provided together.")
    if not genes_tsv.exists():
        raise FileNotFoundError(f"Target genes table not found: {genes_tsv}")
    if not target_sequences_dir.exists():
        raise FileNotFoundError(f"Target sequences directory not found: {target_sequences_dir}")

    contexts = {}
    with open_text(genes_tsv) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", "genomic_accession", "chromosome", "begin", "end"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Target genes table missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            gene_id = row["gene_id"]
            fasta_path = target_sequences_dir / f"{gene_id}.fa.gz"
            if not fasta_path.exists():
                raise FileNotFoundError(f"Target FASTA not found for gene {gene_id}: {fasta_path}")
            contexts[gene_id] = {
                "gene_id": gene_id,
                "accession": row["genomic_accession"],
                "chrom": normalize_chrom(row["chromosome"]) or _refseq_accession_to_gnomad_chrom(row["genomic_accession"]),
                "begin": int(row["begin"]),
                "end": int(row["end"]),
                "fasta_path": fasta_path,
            }
    logger.info(f"Loaded target context for {len(contexts)} gene(s).")
    return contexts


def context_sequence(context: dict) -> str:
    seq = context.get("seq")
    if seq is None:
        seq = read_fasta_sequence(context["fasta_path"])
        context["seq"] = seq
    return seq


def build_context_index(contexts: dict[str, dict]) -> dict[str, tuple[list[dict], list[int]]]:
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for context in contexts.values():
        chrom = context.get("chrom")
        if chrom:
            by_chrom[chrom].append(context)

    index: dict[str, tuple[list[dict], list[int]]] = {}
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda row: (int(row["begin"]), int(row["end"]), row["gene_id"]))
        index[chrom] = (rows, [int(row["begin"]) for row in rows])
    return index


def normalize_vcf_key_for_context(
    context: dict,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> tuple[tuple[str, int, str, str] | None, str]:
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    if not ref or not alt:
        return None, "empty_vcf_allele"

    seq = context_sequence(context)
    pos0 = pos - int(context["begin"])
    if pos0 < 0 or pos0 + len(ref) > len(seq):
        return None, "out_of_target"
    if seq[pos0 : pos0 + len(ref)] != ref:
        return None, "ref_mismatch"

    if len(ref) != len(alt):
        while pos0 > 0 and ref[-1] == alt[-1]:
            prev = seq[pos0 - 1]
            ref = prev + ref[:-1]
            alt = prev + alt[:-1]
            pos0 -= 1

    return (chrom, int(context["begin"]) + pos0, ref, alt), "ok"


def event_vcf_key(row: dict, contexts: dict[str, dict]) -> tuple[tuple[str, int, str, str] | None, str]:
    gene_id = row.get("gene_id", "")
    context = contexts.get(gene_id)
    chrom = _refseq_accession_to_gnomad_chrom(row.get("genomic_accession", ""))
    if not chrom:
        return None, "unknown_chrom"

    try:
        raw_pos = int(row.get("genomic_start1", 0))
    except ValueError:
        return None, "bad_position"
    raw_key = (chrom, raw_pos, row.get("ref", "").upper(), row.get("alt", "").upper())

    if not context:
        return raw_key, "raw_no_context"

    try:
        start0 = int(row.get("target_start0", 0))
    except ValueError:
        return raw_key, "bad_target_position"

    seq = context_sequence(context)
    event_type = row.get("event_type", "")
    ref = row.get("ref", "").upper()
    alt = row.get("alt", "").upper()

    if event_type == "snv":
        if len(ref) != 1 or len(alt) != 1:
            return raw_key, "bad_snv_allele"
        vcf = (chrom, int(context["begin"]) + start0, ref, alt)
    elif event_type == "del":
        if not ref or alt:
            return raw_key, "bad_del_allele"
        if start0 <= 0:
            return raw_key, "missing_left_anchor"
        anchor = seq[start0 - 1]
        vcf = (chrom, int(context["begin"]) + start0 - 1, anchor + ref, anchor)
    elif event_type == "ins":
        if ref or not alt:
            return raw_key, "bad_ins_allele"
        if start0 <= 0:
            return raw_key, "missing_left_anchor"
        anchor = seq[start0 - 1]
        vcf = (chrom, int(context["begin"]) + start0 - 1, anchor, anchor + alt)
    else:
        return raw_key, "unsupported_event_type"

    normalized, status = normalize_vcf_key_for_context(context, *vcf)
    return (normalized or raw_key), status


def contexts_for_variant(
    context_index: dict[str, tuple[list[dict], list[int]]],
    chrom: str,
    pos: int,
) -> list[dict]:
    rows, starts = context_index.get(chrom, ([], []))
    limit = bisect.bisect_right(starts, pos)
    return [context for context in rows[:limit] if int(context["end"]) >= pos]


def add_annotation_cache_entry(
    cache: dict,
    key: tuple[str, int, str, str],
    value,
    contexts: dict[str, dict],
    context_index: dict[str, tuple[list[dict], list[int]]],
    status_counts: Counter,
) -> None:
    cache[key] = value
    chrom, pos, ref, alt = key
    matched_contexts = contexts_for_variant(context_index, chrom, pos)
    if not matched_contexts:
        status_counts["raw_no_context"] += 1
        return
    for context in matched_contexts:
        normalized, status = normalize_vcf_key_for_context(context, chrom, pos, ref, alt)
        status_counts[status] += 1
        if normalized:
            cache[normalized] = value


def fetch_gnomad_for_cluster(chrom: str, start: int, end: int) -> list[dict]:
    # pad by 100 bases
    return fetch_region_variants_recursive(chrom, max(1, start - 100), end + 100)


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


def normalize_review_status(value: str) -> str:
    return str(value or "").strip().lower().replace(" ", "_")


def clinvar_review_stars(review_status: str) -> tuple[str, str]:
    if not review_status:
        return "", "missing"
    statuses = [normalize_review_status(item) for item in review_status.split("|") if item]
    if not statuses:
        return "", "missing"

    stars = []
    for status in statuses:
        star_value = CLINVAR_REVIEW_STARS.get(status)
        if star_value is None:
            return "", f"unmapped:{status}"
        stars.append(star_value)

    unique_stars = sorted(set(stars))
    if len(unique_stars) != 1:
        return "", "ambiguous_multiple_review_statuses"
    return unique_stars[0], "mapped"


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
        "clinvar_sig_conflict": format_info_value(rec.info.get("CLNSIGCONF")),
        "clinvar_scv_count": count_pipe_values(scv_accessions),
        "clinvar_scv_accessions": scv_accessions,
        "clinvar_hgvs": format_info_value(rec.info.get("CLNHGVS")),
        "clinvar_geneinfo": format_info_value(rec.info.get("GENEINFO")),
        "clinvar_disease": format_info_value(rec.info.get("CLNDN")),
        "clinvar_disease_db": format_info_value(rec.info.get("CLNDISDB")),
        "clinvar_variant_type": format_info_value(rec.info.get("CLNVC")),
        "clinvar_variant_type_so": format_info_value(rec.info.get("CLNVCSO")),
        "clinvar_origin": format_info_value(rec.info.get("ORIGIN")),
        "clinvar_rs": format_info_value(rec.info.get("RS")),
    }


def gnomad_annotation_from_variant(variant: dict) -> dict[str, str]:
    af, af_source, af_exome, af_genome, af_joint, an_joint, ac_joint = _select_af_metrics(variant)
    return {
        "gnomad_af": format_float(af),
        "gnomad_af_source": af_source or "",
        "gnomad_csq": str(variant.get("consequence") or ""),
        "gnomad_variant_id": str(variant.get("variant_id") or ""),
        "gnomad_af_exome": format_float(af_exome),
        "gnomad_af_genome": format_float(af_genome),
        "gnomad_af_joint": format_float(af_joint),
        "gnomad_ac_joint": str(ac_joint) if ac_joint is not None else "",
        "gnomad_an_joint": str(an_joint) if an_joint is not None else "",
        "gnomad_hgvsc": str(variant.get("hgvsc") or ""),
        "gnomad_hgvsp": str(variant.get("hgvsp") or ""),
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
        chrom = _refseq_accession_to_gnomad_chrom(acc)
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
                        add_annotation_cache_entry(
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
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = args.outdir / "variant_annotations.tsv.gz"
    failures_tsv = args.outdir / "failures.tsv.gz"
    manifest_json = args.outdir / "manifest.json"
    if not args.clinvar_vcf.exists():
        raise FileNotFoundError(f"ClinVar VCF not found: {args.clinvar_vcf}")
    clinvar_tbi = Path(f"{args.clinvar_vcf}.tbi")
    if not clinvar_tbi.exists():
        raise FileNotFoundError(f"ClinVar VCF index not found: {clinvar_tbi}")

    failures: list[dict] = []
    contexts = load_target_contexts(args.genes_tsv, args.target_sequences_dir)
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
        for row in reader:
            input_row_count += 1
            acc = row["genomic_accession"]
            raw_pos = int_or_default(row.get("genomic_start1"), -1)
            if acc and raw_pos > 0:
                accession_positions[acc].add(raw_pos)
            lookup_key, status = event_vcf_key(row, contexts)
            event_key_status_counts[status] += 1
            if acc and lookup_key:
                accession_positions[acc].add(int(lookup_key[1]))

            lookup_chrom = lookup_key[0] if lookup_key else ""
            lookup_pos = lookup_key[1] if lookup_key else ""
            lookup_ref = lookup_key[2] if lookup_key else ""
            lookup_alt = lookup_key[3] if lookup_key else ""
            variant_key = lookup_key_text(lookup_key)
            aggregate_key = (
                row.get("gene_id", ""),
                row.get("event_type", ""),
                row.get("target_start0", ""),
                row.get("target_end0", ""),
                row.get("genomic_accession", ""),
                row.get("genomic_start1", ""),
                row.get("genomic_end1", ""),
                row.get("ref", ""),
                row.get("alt", ""),
                variant_key,
            )
            aggregate = variant_aggregates.get(aggregate_key)
            if aggregate is None:
                aggregate = {
                    "variant_key": variant_key,
                    "gene_id": row.get("gene_id", ""),
                    "event_type": row.get("event_type", ""),
                    "target_start0": row.get("target_start0", ""),
                    "target_end0": row.get("target_end0", ""),
                    "genomic_accession": row.get("genomic_accession", ""),
                    "genomic_start1": row.get("genomic_start1", ""),
                    "genomic_end1": row.get("genomic_end1", ""),
                    "ref": row.get("ref", ""),
                    "alt": row.get("alt", ""),
                    "lookup_chrom": lookup_chrom,
                    "lookup_pos": lookup_pos,
                    "lookup_ref": lookup_ref,
                    "lookup_alt": lookup_alt,
                    "lookup_status": status,
                    "support_row_count": 0,
                    "_lookup_key": lookup_key,
                    "_orthologs": set(),
                    "_ortholog_count_hint": 0,
                    "_strategies": set(),
                    "_strategy_count_hint": 0,
                    "_tools": set(),
                    "_presets": set(),
                    "_tax_ids": set(),
                    "_tax_id_count_hint": 0,
                    "_taxnames": set(),
                    "_taxname_count_hint": 0,
                }
                variant_aggregates[aggregate_key] = aggregate
                unique_lookup_status_counts[status] += 1

            aggregate["support_row_count"] += int_or_default(row.get("support_row_count"), 1)
            if row.get("ortholog_gene_id"):
                aggregate["_orthologs"].add(row["ortholog_gene_id"])
            else:
                aggregate["_ortholog_count_hint"] += int_or_default(row.get("support_ortholog_count"), 0)

            aggregate["_strategies"].update(split_values(row.get("strategy")))
            aggregate["_strategies"].update(split_values(row.get("strategies")))
            aggregate["_strategy_count_hint"] += int_or_default(row.get("support_strategy_count"), 0)

            aggregate["_tools"].update(split_values(row.get("tool")))
            aggregate["_tools"].update(split_values(row.get("tools")))
            aggregate["_presets"].update(split_values(row.get("preset")))
            aggregate["_presets"].update(split_values(row.get("presets")))

            if row.get("tax_id"):
                aggregate["_tax_ids"].add(row["tax_id"])
            else:
                aggregate["_tax_id_count_hint"] += int_or_default(row.get("tax_id_count"), 0)
            if row.get("taxname"):
                aggregate["_taxnames"].add(row["taxname"])
            else:
                aggregate["_taxname_count_hint"] += int_or_default(row.get("taxname_count"), 0)
    logger.info(f"Event key normalization status: {dict(event_key_status_counts)}")
    logger.info(f"Collapsed {input_row_count} event row(s) to {len(variant_aggregates)} variant-context row(s).")

    # 2. Determine gnomAD clusters
    gnomad_tasks = []
    for acc, positions in accession_positions.items():
        chrom = _refseq_accession_to_gnomad_chrom(acc)
        if not chrom:
            continue
        clusters = cluster_positions(list(positions), max_gap=200000)
        for c_start, c_end in clusters:
            gnomad_tasks.append((chrom, c_start, c_end))

    logger.info(f"Will fetch {len(gnomad_tasks)} region(s) from gnomAD API.")

    # 3. Fetch gnomAD in parallel and cache
    gnomad_cache = {}
    gnomad_key_status_counts = Counter()
    gnomad_region_success_count = 0
    gnomad_raw_variant_count = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        future_to_task = {
            executor.submit(fetch_gnomad_for_cluster, chrom, start, end): (chrom, start, end)
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
                    add_annotation_cache_entry(
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
        row["support_ortholog_count"] = max(
            len(aggregate["_orthologs"]),
            int_or_default(aggregate["_ortholog_count_hint"]),
        )
        row["support_strategy_count"] = max(
            len(aggregate["_strategies"]),
            int_or_default(aggregate["_strategy_count_hint"]),
        )
        row["strategies"] = ",".join(sorted(aggregate["_strategies"]))
        row["tools"] = ",".join(sorted(aggregate["_tools"]))
        row["presets"] = ",".join(sorted(aggregate["_presets"]))
        row["tax_id_count"] = max(len(aggregate["_tax_ids"]), int_or_default(aggregate["_tax_id_count_hint"]))
        row["taxname_count"] = max(len(aggregate["_taxnames"]), int_or_default(aggregate["_taxname_count_hint"]))
        row.update(clinvar_annotation)
        row.update(gnomad_annotation)
        for column in ANNOTATION_COLUMNS:
            if row[column]:
                annotation_value_counts[column] += 1
        variant_rows.append(row)

    variant_rows.sort(
        key=lambda row: (
            int_or_default(row.get("gene_id"), 10**18),
            int_or_default(row.get("target_start0"), 10**18),
            row.get("event_type", ""),
            row.get("variant_key", ""),
        )
    )
    output_row_count = write_tsv_gz(out_tsv, VARIANT_ANNOTATION_FIELDS, variant_rows)

    failure_count = write_tsv_gz(failures_tsv, FAILURE_FIELDS, failures)
    manifest = {
        "output_mode": "unique_variant_context",
        "event_row_count": input_row_count,
        "variant_context_count": len(variant_aggregates),
        "annotated_variant_context_count": output_row_count,
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
        "failure_count": failure_count,
        "annotation_nonempty_counts": dict(annotation_value_counts),
        "event_key_status_counts": dict(event_key_status_counts),
        "unique_lookup_status_counts": dict(unique_lookup_status_counts),
        "gnomad_key_status_counts": dict(gnomad_key_status_counts),
        "clinvar_key_status_counts": dict(clinvar_key_status_counts),
    }
    manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    logger.info(f"Saved variant annotations to {out_tsv}")
    logger.info(f"Saved annotation failures to {failures_tsv}")
    logger.info(f"Saved annotation manifest to {manifest_json}")

if __name__ == "__main__":
    main()
