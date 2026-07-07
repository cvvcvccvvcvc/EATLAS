"""ClinVar ALT-observed validation for completed GAPH runs."""

from __future__ import annotations

import bisect
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


UNIVERSE_FIELDS = [
    "variant_key",
    "chrom",
    "pos",
    "ref",
    "alt",
    "label_class",
    "clinvar_ids",
    "clinvar_sigs",
    "gene_ids",
]
CACHE_VERSION = 1


@dataclass(frozen=True)
class ClinvarValidation:
    universe_path: Path
    manifest_path: Path
    universe: pd.DataFrame
    strategy_results: pd.DataFrame
    manifest: dict


def build_validation(
    *,
    run_dir: Path,
    variant_annotations_tsv: Path,
    genes_tsv: Path,
    clinvar_vcf: Path,
    strategies: list[str],
) -> ClinvarValidation:
    analytics_dir = run_dir / "analytics"
    analytics_dir.mkdir(parents=True, exist_ok=True)
    universe_path = analytics_dir / "clinvar_universe.snv.tsv.gz"
    manifest_path = analytics_dir / "clinvar_universe.snv.manifest.json"
    regions_path = analytics_dir / "clinvar_target_regions.bed"

    manifest = build_or_load_clinvar_universe(
        genes_tsv=genes_tsv,
        clinvar_vcf=clinvar_vcf,
        universe_path=universe_path,
        manifest_path=manifest_path,
        regions_path=regions_path,
    )
    universe = pd.read_csv(universe_path, sep="\t", compression="gzip", keep_default_na=False)
    strategy_results = compute_strategy_results(
        universe=universe,
        variant_annotations_tsv=variant_annotations_tsv,
        strategies=strategies,
    )
    return ClinvarValidation(universe_path, manifest_path, universe, strategy_results, manifest)


def build_or_load_clinvar_universe(
    *,
    genes_tsv: Path,
    clinvar_vcf: Path,
    universe_path: Path,
    manifest_path: Path,
    regions_path: Path,
) -> dict:
    expected_inputs = {
        "genes_tsv": path_metadata(genes_tsv),
        "clinvar_vcf": path_metadata(clinvar_vcf),
        "clinvar_tbi": path_metadata(Path(f"{clinvar_vcf}.tbi")),
        "mode": "snv_only",
        "cache_version": CACHE_VERSION,
    }
    if universe_path.exists() and manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("inputs") == expected_inputs:
            return manifest

    genes = read_genes(genes_tsv)
    intervals = merged_intervals(genes)
    write_regions_bed(regions_path, intervals)
    rows, counts = query_clinvar_snv_universe(clinvar_vcf, regions_path, genes)
    write_universe(universe_path, rows)
    manifest = {
        "inputs": expected_inputs,
        "target_gene_count": len(genes),
        "target_region_count": len(intervals),
        "raw_snv_allele_count": counts["raw_snv_allele_count"],
        "usable_snv_allele_count": len(rows),
        "benign_count": counts["benign"],
        "pathogenic_count": counts["pathogenic"],
        "excluded_vus_count": counts["excluded_vus"],
        "excluded_missing_count": counts["excluded_missing"],
        "excluded_other_count": counts["excluded_other"],
        "ambiguous_mixed_label_count": counts["ambiguous_mixed_label"],
        "duplicate_usable_key_count": counts["duplicate_usable_key"],
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


def normalize_chrom(value: str) -> str:
    chrom = str(value or "").strip()
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    return "MT" if chrom == "M" else chrom


def merged_intervals(genes: list[dict[str, str]]) -> list[tuple[str, int, int]]:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in genes:
        chrom = normalize_chrom(row["chromosome"])
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


def query_clinvar_snv_universe(
    clinvar_vcf: Path,
    regions_path: Path,
    genes: list[dict[str, str]],
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

    gene_index = build_gene_index(genes)
    raw_by_key: dict[str, dict[str, object]] = {}
    counts = Counter()
    for line in proc.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        chrom, pos_text, rec_id, ref, alt_text, _qual, _filter, info_text = fields[:8]
        chrom = normalize_chrom(chrom)
        pos = int(pos_text)
        ref = ref.upper()
        sig = info_value(info_text, "CLNSIG")
        label = clinvar_label(sig)
        for alt in alt_text.split(","):
            alt = alt.upper()
            if not is_snv(ref, alt):
                continue
            counts["raw_snv_allele_count"] += 1
            if label == "excluded_vus":
                counts["excluded_vus"] += 1
            elif label == "excluded_missing":
                counts["excluded_missing"] += 1
            elif label == "excluded_other":
                counts["excluded_other"] += 1
            if label not in {"benign", "pathogenic"}:
                continue
            key = f"{chrom}:{pos}:{ref}>{alt}"
            entry = raw_by_key.setdefault(
                key,
                {
                    "variant_key": key,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "labels": set(),
                    "clinvar_ids": set(),
                    "clinvar_sigs": set(),
                    "gene_ids": set(overlapping_gene_ids(gene_index, chrom, pos)),
                },
            )
            entry["labels"].add(label)
            if rec_id and rec_id != ".":
                entry["clinvar_ids"].add(rec_id)
            if sig:
                entry["clinvar_sigs"].add(sig)

    rows = []
    for entry in raw_by_key.values():
        labels = set(entry["labels"])
        if labels == {"benign"}:
            label_class = "benign"
            counts["benign"] += 1
        elif labels == {"pathogenic"}:
            label_class = "pathogenic"
            counts["pathogenic"] += 1
        else:
            counts["ambiguous_mixed_label"] += 1
            continue
        if len(entry["clinvar_ids"]) > 1 or len(entry["clinvar_sigs"]) > 1:
            counts["duplicate_usable_key"] += 1
        rows.append(
            {
                "variant_key": entry["variant_key"],
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
    rows.sort(key=lambda row: (chrom_sort_key(str(row["chrom"])), int(row["pos"]), row["ref"], row["alt"]))
    return rows, counts


def is_snv(ref: str, alt: str) -> bool:
    bases = {"A", "C", "G", "T"}
    return len(ref) == 1 and len(alt) == 1 and ref in bases and alt in bases


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


def build_gene_index(genes: list[dict[str, str]]) -> dict[str, tuple[list[dict[str, object]], list[int]]]:
    by_chrom: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in genes:
        chrom = normalize_chrom(row["chromosome"])
        if not chrom:
            continue
        by_chrom[chrom].append(
            {
                "gene_id": row["gene_id"],
                "begin": int(row["begin"]),
                "end": int(row["end"]),
            }
        )

    index = {}
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda row: (int(row["begin"]), int(row["end"]), str(row["gene_id"])))
        index[chrom] = (rows, [int(row["begin"]) for row in rows])
    return index


def overlapping_gene_ids(
    gene_index: dict[str, tuple[list[dict[str, object]], list[int]]],
    chrom: str,
    pos: int,
) -> list[str]:
    rows, starts = gene_index.get(chrom, ([], []))
    limit = bisect.bisect_right(starts, pos)
    return [str(row["gene_id"]) for row in rows[:limit] if int(row["end"]) >= pos]


def gene_sort_key(value: str) -> tuple[int, str]:
    return (int(value), value) if str(value).isdigit() else (10**18, str(value))


def write_universe(path: Path, rows: list[dict[str, str]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=UNIVERSE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in UNIVERSE_FIELDS})


def compute_strategy_results(
    *,
    universe: pd.DataFrame,
    variant_annotations_tsv: Path,
    strategies: list[str],
    chunksize: int = 250_000,
) -> pd.DataFrame:
    if universe.empty:
        return pd.DataFrame(columns=result_columns())

    labels_by_key = dict(zip(universe["variant_key"].astype(str), universe["label_class"].astype(str)))
    universe_keys = set(labels_by_key)
    observed_by_strategy: dict[str, set[str]] = {strategy: set() for strategy in strategies}

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
            for strategy in split_strategies(strategy_text):
                observed_by_strategy.setdefault(strategy, set()).add(variant_key)

    benign_total = int((universe["label_class"] == "benign").sum())
    pathogenic_total = int((universe["label_class"] == "pathogenic").sum())
    rows = []
    for strategy in sorted(observed_by_strategy):
        observed = observed_by_strategy[strategy]
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
