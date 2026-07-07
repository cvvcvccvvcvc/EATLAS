"""ClinVar ALT-observed validation for completed GAPH runs."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .stats import enrichment_result
from .variant_keys import (
    build_context_index,
    contexts_for_variant,
    load_target_contexts,
    normalize_chrom,
    normalize_vcf_key_for_context,
    variant_key_text,
    variant_type,
)


UNIVERSE_FIELDS = [
    "variant_key",
    "variant_type",
    "chrom",
    "pos",
    "ref",
    "alt",
    "label_class",
    "clinvar_ids",
    "clinvar_sigs",
    "gene_ids",
]
CACHE_VERSION = 2
VALIDATION_TYPES = ["snv", "indel"]


@dataclass(frozen=True)
class ClinvarValidation:
    universe_path: Path
    manifest_path: Path
    universe: pd.DataFrame
    strategy_results: pd.DataFrame
    manifest: dict
    observed_by_strategy_type: dict[tuple[str, str], set[str]]


def build_validation(
    *,
    run_dir: Path,
    variant_annotations_tsv: Path,
    genes_tsv: Path,
    target_sequences_dir: Path,
    clinvar_vcf: Path,
    strategies: list[str],
) -> ClinvarValidation:
    analytics_dir = run_dir / "analytics"
    analytics_dir.mkdir(parents=True, exist_ok=True)
    universe_path = analytics_dir / "clinvar_universe.snv_indel.tsv.gz"
    manifest_path = analytics_dir / "clinvar_universe.snv_indel.manifest.json"
    regions_path = analytics_dir / "clinvar_target_regions.bed"

    manifest = build_or_load_clinvar_universe(
        genes_tsv=genes_tsv,
        target_sequences_dir=target_sequences_dir,
        clinvar_vcf=clinvar_vcf,
        universe_path=universe_path,
        manifest_path=manifest_path,
        regions_path=regions_path,
    )
    universe = pd.read_csv(universe_path, sep="\t", compression="gzip", keep_default_na=False)
    observed_by_strategy_type = collect_observed_keys_by_strategy_type(
        universe=universe,
        variant_annotations_tsv=variant_annotations_tsv,
        strategies=strategies,
    )
    strategy_results = compute_strategy_results(
        universe=universe,
        strategies=strategies,
        observed_by_strategy_type=observed_by_strategy_type,
    )
    return ClinvarValidation(
        universe_path,
        manifest_path,
        universe,
        strategy_results,
        manifest,
        observed_by_strategy_type,
    )


def build_or_load_clinvar_universe(
    *,
    genes_tsv: Path,
    target_sequences_dir: Path,
    clinvar_vcf: Path,
    universe_path: Path,
    manifest_path: Path,
    regions_path: Path,
) -> dict:
    expected_inputs = {
        "genes_tsv": path_metadata(genes_tsv),
        "target_sequences_dir": directory_metadata(target_sequences_dir),
        "clinvar_vcf": path_metadata(clinvar_vcf),
        "clinvar_tbi": path_metadata(Path(f"{clinvar_vcf}.tbi")),
        "mode": "snv_indel",
        "cache_version": CACHE_VERSION,
    }
    if universe_path.exists() and manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("inputs") == expected_inputs:
            return manifest

    genes = read_genes(genes_tsv)
    contexts = load_target_contexts(genes_tsv, target_sequences_dir)
    context_index = build_context_index(contexts)
    intervals = merged_intervals(genes)
    write_regions_bed(regions_path, intervals)
    rows, counts = query_clinvar_variant_universe(clinvar_vcf, regions_path, context_index)
    write_universe(universe_path, rows)
    row_counts = Counter((row["variant_type"], row["label_class"]) for row in rows)
    manifest = {
        "inputs": expected_inputs,
        "target_gene_count": len(genes),
        "target_region_count": len(intervals),
        "raw_allele_count": counts["raw_allele_count"],
        "raw_snv_allele_count": counts["raw_snv_allele_count"],
        "raw_indel_allele_count": counts["raw_indel_allele_count"],
        "usable_allele_count": len(rows),
        "usable_snv_allele_count": sum(1 for row in rows if row["variant_type"] == "snv"),
        "usable_indel_allele_count": sum(1 for row in rows if row["variant_type"] == "indel"),
        "benign_count": sum(1 for row in rows if row["label_class"] == "benign"),
        "pathogenic_count": sum(1 for row in rows if row["label_class"] == "pathogenic"),
        "benign_snv_count": row_counts[("snv", "benign")],
        "pathogenic_snv_count": row_counts[("snv", "pathogenic")],
        "benign_indel_count": row_counts[("indel", "benign")],
        "pathogenic_indel_count": row_counts[("indel", "pathogenic")],
        "excluded_vus_count": counts["excluded_vus_count"],
        "excluded_missing_count": counts["excluded_missing_count"],
        "excluded_other_count": counts["excluded_other_count"],
        "excluded_vus_snv_count": counts["excluded_vus_snv_count"],
        "excluded_missing_snv_count": counts["excluded_missing_snv_count"],
        "excluded_other_snv_count": counts["excluded_other_snv_count"],
        "excluded_vus_indel_count": counts["excluded_vus_indel_count"],
        "excluded_missing_indel_count": counts["excluded_missing_indel_count"],
        "excluded_other_indel_count": counts["excluded_other_indel_count"],
        "excluded_normalization_snv_count": counts["excluded_normalization_snv_count"],
        "excluded_normalization_indel_count": counts["excluded_normalization_indel_count"],
        "excluded_complex_allele_count": counts["excluded_complex_allele_count"],
        "excluded_unsupported_allele_count": counts["excluded_unsupported_allele_count"],
        "ambiguous_mixed_label_count": counts["ambiguous_mixed_label_count"],
        "ambiguous_mixed_label_snv_count": counts["ambiguous_mixed_label_snv_count"],
        "ambiguous_mixed_label_indel_count": counts["ambiguous_mixed_label_indel_count"],
        "duplicate_usable_key_count": counts["duplicate_usable_key_count"],
        "regions_bed": str(regions_path),
        "universe_tsv": str(universe_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def read_genes(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", "chromosome", "begin", "end"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Genes table missing required columns: {', '.join(sorted(missing))}")
        return [row for row in reader]


def merged_intervals(genes: list[dict[str, str]]) -> list[tuple[str, int, int]]:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in genes:
        chrom = normalize_chrom(row["chromosome"]) or ""
        if not chrom:
            continue
        start = int(row["begin"])
        end = int(row["end"])
        if end >= start:
            by_chrom[chrom].append((start, end))

    merged = []
    for chrom, intervals in by_chrom.items():
        intervals.sort()
        current_start = None
        current_end = None
        for start, end in intervals:
            if current_start is None:
                current_start, current_end = start, end
                continue
            if start <= int(current_end) + 1:
                current_end = max(int(current_end), end)
            else:
                merged.append((chrom, int(current_start), int(current_end)))
                current_start, current_end = start, end
        if current_start is not None:
            merged.append((chrom, int(current_start), int(current_end)))
    return sorted(merged, key=lambda item: (chrom_sort_key(item[0]), item[1], item[2]))


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    if chrom.isdigit():
        return int(chrom), chrom
    order = {"X": 23, "Y": 24, "MT": 25}
    return order.get(chrom, 10**6), chrom


def write_regions_bed(path: Path, intervals: list[tuple[str, int, int]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for chrom, start1, end1 in intervals:
            writer.writerow([chrom, max(0, start1 - 1), end1])


def query_clinvar_variant_universe(
    clinvar_vcf: Path,
    regions_path: Path,
    context_index: dict[str, tuple[list[dict], list[int]]],
) -> tuple[list[dict[str, str]], Counter]:
    tabix = shutil.which("tabix")
    if tabix is None:
        raise FileNotFoundError("tabix executable not found; it is required for indexed ClinVar queries.")

    proc = subprocess.run(
        [tabix, "-R", str(regions_path), str(clinvar_vcf)],
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"tabix ClinVar query failed: {proc.stderr.strip()}")

    raw_by_key: dict[str, dict[str, object]] = {}
    counts = Counter()
    for line in proc.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        chrom, pos_text, rec_id, ref, alt_text, _qual, _filter, info_text = fields[:8]
        chrom = normalize_chrom(chrom) or ""
        pos = int(pos_text)
        ref = ref.upper()
        sig = info_value(info_text, "CLNSIG")
        label = clinvar_label(sig)
        for alt in alt_text.split(","):
            alt = alt.upper()
            vtype = variant_type(ref, alt)
            if vtype == "unsupported":
                counts["excluded_unsupported_allele_count"] += 1
                continue
            if vtype == "complex":
                counts["excluded_complex_allele_count"] += 1
                continue
            counts["raw_allele_count"] += 1
            counts[f"raw_{vtype}_allele_count"] += 1
            if label.startswith("excluded_"):
                counts[f"{label}_count"] += 1
                counts[f"{label}_{vtype}_count"] += 1
            if label not in {"benign", "pathogenic"}:
                continue

            normalized_items = normalize_clinvar_allele_for_targets(context_index, chrom, pos, ref, alt)
            if not normalized_items:
                counts[f"excluded_normalization_{vtype}_count"] += 1
                continue
            for key, key_type, context in normalized_items:
                key_text = variant_key_text(key)
                entry = raw_by_key.setdefault(
                    key_text,
                    {
                        "variant_key": key_text,
                        "variant_type": key_type,
                        "chrom": key[0],
                        "pos": key[1],
                        "ref": key[2],
                        "alt": key[3],
                        "labels": set(),
                        "clinvar_ids": set(),
                        "clinvar_sigs": set(),
                        "gene_ids": set(),
                    },
                )
                entry["labels"].add(label)
                if rec_id and rec_id != ".":
                    entry["clinvar_ids"].add(rec_id)
                if sig:
                    entry["clinvar_sigs"].add(sig)
                entry["gene_ids"].add(str(context["gene_id"]))

    rows = []
    for entry in raw_by_key.values():
        labels = set(entry["labels"])
        if labels == {"benign"}:
            label_class = "benign"
        elif labels == {"pathogenic"}:
            label_class = "pathogenic"
        else:
            counts["ambiguous_mixed_label_count"] += 1
            counts[f"ambiguous_mixed_label_{entry['variant_type']}_count"] += 1
            continue
        if len(entry["clinvar_ids"]) > 1 or len(entry["clinvar_sigs"]) > 1:
            counts["duplicate_usable_key_count"] += 1
        rows.append(
            {
                "variant_key": entry["variant_key"],
                "variant_type": entry["variant_type"],
                "chrom": entry["chrom"],
                "pos": entry["pos"],
                "ref": entry["ref"],
                "alt": entry["alt"],
                "label_class": label_class,
                "clinvar_ids": "|".join(sorted(entry["clinvar_ids"])),
                "clinvar_sigs": "|".join(sorted(entry["clinvar_sigs"])),
                "gene_ids": "|".join(sorted(entry["gene_ids"], key=gene_sort_key)),
            }
        )
    rows.sort(key=lambda row: (row["variant_type"], chrom_sort_key(str(row["chrom"])), int(row["pos"]), row["ref"], row["alt"]))
    return rows, counts


def normalize_clinvar_allele_for_targets(
    context_index: dict[str, tuple[list[dict], list[int]]],
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> list[tuple[tuple[str, int, str, str], str, dict]]:
    normalized = []
    seen = set()
    for context in contexts_for_variant(context_index, chrom, pos):
        key, status = normalize_vcf_key_for_context(context, chrom, pos, ref, alt)
        if status != "ok" or key is None:
            continue
        key_type = variant_type(key[2], key[3])
        if key_type not in VALIDATION_TYPES:
            continue
        dedupe_key = (key, context["gene_id"])
        if dedupe_key in seen:
            continue
        seen.add(dedupe_key)
        normalized.append((key, key_type, context))
    return normalized


def info_value(info_text: str, key: str) -> str:
    prefix = f"{key}="
    for item in info_text.split(";"):
        if item.startswith(prefix):
            return item[len(prefix) :]
    return ""


def clinvar_label(clnsig: str) -> str:
    text = str(clnsig or "").lower()
    if not text:
        return "excluded_missing"
    if "conflicting" in text:
        return "excluded_other"
    if "uncertain" in text or "vus" in text:
        return "excluded_vus"
    benign = "benign" in text
    pathogenic = "pathogenic" in text
    if benign and not pathogenic:
        return "benign"
    if pathogenic and not benign:
        return "pathogenic"
    return "excluded_other"


def gene_sort_key(value: str) -> tuple[int, str]:
    return (int(value), value) if str(value).isdigit() else (10**18, str(value))


def write_universe(path: Path, rows: list[dict[str, str]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=UNIVERSE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in UNIVERSE_FIELDS})


def collect_observed_keys_by_strategy_type(
    *,
    universe: pd.DataFrame,
    variant_annotations_tsv: Path,
    strategies: list[str],
    chunksize: int = 250_000,
) -> dict[tuple[str, str], set[str]]:
    if universe.empty:
        return {(strategy, variant_kind): set() for strategy in strategies for variant_kind in VALIDATION_TYPES}

    types_by_key = dict(zip(universe["variant_key"].astype(str), universe["variant_type"].astype(str)))
    universe_keys = set(types_by_key)
    observed_by_strategy_type: dict[tuple[str, str], set[str]] = {
        (strategy, variant_kind): set() for strategy in strategies for variant_kind in VALIDATION_TYPES
    }

    for chunk in pd.read_csv(
        variant_annotations_tsv,
        sep="\t",
        compression="gzip",
        usecols=["variant_key", "strategies"],
        keep_default_na=False,
        chunksize=chunksize,
    ):
        matched = chunk[chunk["variant_key"].astype(str).isin(universe_keys)]
        if matched.empty:
            continue
        for variant_key, strategy_text in zip(matched["variant_key"].astype(str), matched["strategies"].astype(str)):
            variant_kind = types_by_key.get(variant_key)
            if variant_kind not in VALIDATION_TYPES:
                continue
            for strategy in split_strategies(strategy_text):
                observed_by_strategy_type.setdefault((strategy, variant_kind), set()).add(variant_key)

    return observed_by_strategy_type


def compute_strategy_results(
    *,
    universe: pd.DataFrame,
    strategies: list[str],
    observed_by_strategy_type: dict[tuple[str, str], set[str]],
) -> pd.DataFrame:
    if universe.empty:
        return pd.DataFrame(columns=result_columns())

    labels_by_key = dict(zip(universe["variant_key"].astype(str), universe["label_class"].astype(str)))
    rows = []
    all_strategies = sorted(set(strategies) | {strategy for strategy, _variant_kind in observed_by_strategy_type})
    for variant_kind in VALIDATION_TYPES:
        subset = universe[universe["variant_type"].astype(str) == variant_kind]
        if subset.empty:
            continue
        benign_total = int((subset["label_class"] == "benign").sum())
        pathogenic_total = int((subset["label_class"] == "pathogenic").sum())
        for strategy in all_strategies:
            observed = observed_by_strategy_type.get((strategy, variant_kind), set())
            benign_observed = sum(1 for key in observed if labels_by_key.get(key) == "benign")
            pathogenic_observed = sum(1 for key in observed if labels_by_key.get(key) == "pathogenic")
            result = enrichment_result(
                strategy,
                benign_observed,
                pathogenic_observed,
                benign_total - benign_observed,
                pathogenic_total - pathogenic_observed,
            )
            rows.append(
                {
                    "variant_type": variant_kind,
                    "strategy": strategy,
                    "benign_observed": result.benign_observed,
                    "pathogenic_observed": result.pathogenic_observed,
                    "benign_not_observed": result.benign_not_observed,
                    "pathogenic_not_observed": result.pathogenic_not_observed,
                    "odds_ratio": result.odds_ratio,
                    "ci_low": result.ci_low,
                    "ci_high": result.ci_high,
                    "fisher_p": result.fisher_p,
                }
            )
    return pd.DataFrame(rows, columns=result_columns())


def result_columns() -> list[str]:
    return [
        "variant_type",
        "strategy",
        "benign_observed",
        "pathogenic_observed",
        "benign_not_observed",
        "pathogenic_not_observed",
        "odds_ratio",
        "ci_low",
        "ci_high",
        "fisher_p",
    ]


def split_strategies(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def path_metadata(path: Path) -> dict[str, object]:
    if not path.exists():
        raise FileNotFoundError(path)
    stat = path.stat()
    return {"path": str(path.resolve()), "size_bytes": stat.st_size, "mtime": int(stat.st_mtime)}


def directory_metadata(path: Path) -> dict[str, object]:
    if not path.exists():
        raise FileNotFoundError(path)
    files = sorted(item for item in path.glob("*.fa.gz") if item.is_file())
    sizes = 0
    mtimes = []
    for item in files:
        stat = item.stat()
        sizes += stat.st_size
        mtimes.append(int(stat.st_mtime))
    return {
        "path": str(path.resolve()),
        "file_count": len(files),
        "size_bytes": sizes,
        "max_mtime": max(mtimes) if mtimes else 0,
    }
