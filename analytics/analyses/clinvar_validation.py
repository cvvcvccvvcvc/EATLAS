"""ClinVar ALT-observed validation for completed GAPH runs."""

from __future__ import annotations

import csv
import gzip
import io
import json
import os
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import pandas as pd

from analytics.io.artifacts import (
    directory_metadata,
    path_metadata,
    write_json_atomic,
    write_text_atomic,
    write_tsv_atomic,
)
from analytics.io.performance import PerformanceProfile, profile_stage
from genomics.variants import (
    build_context_index,
    contexts_for_variant,
    load_target_contexts,
    normalize_chrom,
    normalize_vcf_key_for_context,
    parse_variant_key,
    variant_key_text,
    variant_type,
)
from analytics.vep.consequences import VEP_CONSEQUENCE_ORDER
from genomics.vep.annotator import annotate_vep_consequences
from genomics.vep.result_cache import DEFAULT_TILE_SIZE_BP
from analytics.analyses.observed_variant_store import ObservedVariantStore


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
    "clinvar_mc_so_ids",
    "clinvar_mc_terms",
    "gene_ids",
    "clinvar_disease_names",
    "clinvar_disease_ids",
]
CACHE_VERSION = 7
OBSERVED_MEMBERSHIP_CACHE_VERSION = 1
VALIDATION_TYPES = ["snv", "indel"]
OBSERVED_MEMBERSHIP_COLUMNS = ["strategy", "variant_type", "variant_key"]


@dataclass(frozen=True)
class ClinvarValidation:
    universe_path: Path
    manifest_path: Path
    universe: pd.DataFrame
    manifest: dict
    observed_by_strategy_type: dict[tuple[str, str], set[str]]
    consequence_column: str = "vep_consequence_terms"
    consequence_source: str = "Ensembl VEP"
    observed_memberships_path: Path | None = None
    observed_memberships_manifest_path: Path | None = None


def build_validation(
    *,
    analytics_dir: Path,
    genes_tsv: Path | Sequence[Path],
    target_sequences_dir: Path | Sequence[Path],
    clinvar_vcf: Path,
    strategies: list[str],
    observed_store: ObservedVariantStore,
    vep_backend: str = "rest",
    vep_release: str | None = None,
    vep_executable: str | Path = "vep",
    vep_cache_dir: Path | None = None,
    vep_forks: int = 1,
    vep_result_cache_dir: Path | None = None,
    vep_result_cache_tile_size_bp: int = DEFAULT_TILE_SIZE_BP,
    performance_profile: PerformanceProfile | None = None,
) -> ClinvarValidation:
    analytics_dir.mkdir(parents=True, exist_ok=True)
    universe_path = analytics_dir / "clinvar_universe.snv_indel.tsv.gz"
    manifest_path = analytics_dir / "clinvar_universe.snv_indel.manifest.json"
    regions_path = analytics_dir / "clinvar_target_regions.bed"

    with profile_stage(performance_profile, "ClinVar universe") as timing:
        manifest = build_or_load_clinvar_universe(
            genes_tsv=genes_tsv,
            target_sequences_dir=target_sequences_dir,
            clinvar_vcf=clinvar_vcf,
            universe_path=universe_path,
            manifest_path=manifest_path,
            regions_path=regions_path,
        )
        universe = pd.read_csv(
            universe_path,
            sep="\t",
            compression="gzip",
            keep_default_na=False,
        )
        timing["metrics"] = {"alleles": int(len(universe))}
    membership_universe_path = universe_path
    if not vep_release:
        raise ValueError("A pinned VEP release is required for ClinVar consequence annotation")
    with profile_stage(performance_profile, "ClinVar VEP consequences") as timing:
        universe_path, vep_manifest_path, universe, vep_manifest = build_or_load_vep_universe(
            universe=universe,
            universe_path=universe_path,
            analytics_dir=analytics_dir,
            backend=vep_backend,
            release=str(vep_release),
            vep_executable=vep_executable,
            vep_cache_dir=vep_cache_dir,
            vep_forks=vep_forks,
            vep_result_cache_dir=vep_result_cache_dir,
            vep_result_cache_tile_size_bp=vep_result_cache_tile_size_bp,
        )
        timing["details"] = "cache hit" if vep_manifest["cache_hit"] else "cache miss"
        timing["metrics"] = {
            "alleles": int(len(universe)),
            "requested": int(vep_manifest.get("request_count", 0)),
        }
    manifest = {**manifest, "consequence_source": "Ensembl VEP", "vep": vep_manifest}
    manifest_path = vep_manifest_path
    with profile_stage(performance_profile, "ClinVar observed memberships") as timing:
        (
            observed_by_strategy_type,
            membership_manifest,
            observed_memberships_path,
            observed_memberships_manifest_path,
        ) = build_or_load_observed_keys_by_strategy_type(
            universe=universe,
            universe_path=membership_universe_path,
            observed_store=observed_store,
            strategies=strategies,
            analytics_dir=analytics_dir,
        )
        timing["details"] = (
            "cache hit" if membership_manifest["cache_hit"] else "cache miss"
        )
        timing["metrics"] = {
            "memberships": int(membership_manifest["membership_count"]),
        }
    return ClinvarValidation(
        universe_path,
        manifest_path,
        universe,
        manifest,
        observed_by_strategy_type,
        "vep_consequence_terms",
        "Ensembl VEP",
        observed_memberships_path,
        observed_memberships_manifest_path,
    )


def build_or_load_vep_universe(
    *,
    universe: pd.DataFrame,
    universe_path: Path,
    analytics_dir: Path,
    backend: str,
    release: str,
    vep_executable: str | Path,
    vep_cache_dir: Path | None,
    vep_forks: int,
    vep_result_cache_dir: Path | None = None,
    vep_result_cache_tile_size_bp: int = DEFAULT_TILE_SIZE_BP,
) -> tuple[Path, Path, pd.DataFrame, dict[str, object]]:
    """Annotate the compact ClinVar validation universe once with RefSeq VEP."""

    output_path = analytics_dir / "clinvar_universe.snv_indel.vep.tsv.gz"
    manifest_path = analytics_dir / "clinvar_universe.snv_indel.vep.manifest.json"
    contract = {
        "schema_version": 2,
        "source": path_metadata(universe_path),
        "backend": backend,
        "release": release,
    }
    if output_path.exists() and manifest_path.exists():
        existing = json.loads(manifest_path.read_text())
        if (
            existing.get("contract") == contract
            and existing.get("output") == path_metadata(output_path)
        ):
            return (
                output_path,
                manifest_path,
                pd.read_csv(output_path, sep="\t", compression="gzip", keep_default_na=False),
                {**existing, "cache_hit": True},
            )

    requests = _clinvar_vep_requests(universe)
    with tempfile.TemporaryDirectory(prefix="gaph_clinvar_vep_") as temporary:
        annotations, summary = annotate_vep_consequences(
            requests,
            Path(temporary) / "clinvar_vep.sqlite",
            backend=backend,
            release=release,
            vep_executable=vep_executable,
            vep_cache_dir=vep_cache_dir,
            vep_forks=vep_forks,
            vep_result_cache_dir=vep_result_cache_dir,
            vep_result_cache_tile_size_bp=vep_result_cache_tile_size_bp,
        )
    aggregated = _aggregate_vep_by_variant(annotations)
    enriched = universe.merge(
        aggregated,
        on="variant_key",
        how="left",
        validate="one_to_one",
    )
    if enriched["vep_status"].isna().any():
        raise ValueError("VEP did not return a status for every ClinVar allele")
    disease_columns = ["clinvar_disease_names", "clinvar_disease_ids"]
    missing_disease_columns = set(disease_columns) - set(enriched.columns)
    if missing_disease_columns:
        raise ValueError(
            "ClinVar universe is missing condition columns: "
            + ", ".join(sorted(missing_disease_columns))
        )
    enriched = enriched[
        [column for column in enriched.columns if column not in disease_columns]
        + disease_columns
    ]
    write_tsv_atomic(output_path, enriched)
    manifest = {
        "status": "complete",
        "contract": contract,
        "allele_count": len(enriched),
        "request_count": len(requests),
        "status_counts": {
            str(status): int(count)
            for status, count in enriched["vep_status"].value_counts().sort_index().items()
        },
        "vep": summary,
        "output": path_metadata(output_path),
    }
    write_json_atomic(manifest_path, manifest)
    return output_path, manifest_path, enriched, {**manifest, "cache_hit": False}


def _clinvar_vep_requests(universe: pd.DataFrame) -> pd.DataFrame:
    requests = []
    for row in universe[["variant_key", "gene_ids"]].itertuples(index=False):
        key = str(row.variant_key)
        parsed = parse_variant_key(key)
        if parsed is None:
            raise ValueError(f"Invalid normalized ClinVar variant_key: {key}")
        chrom, pos, ref, alt = parsed
        for gene_id in str(row.gene_ids).split("|"):
            if gene_id:
                requests.append(
                    {
                        "variant_key": key,
                        "gene_id": gene_id,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                    }
                )
    return pd.DataFrame(
        requests,
        columns=["variant_key", "gene_id", "chrom", "pos", "ref", "alt"],
    ).drop_duplicates(["variant_key", "gene_id"])


def _aggregate_vep_by_variant(annotations: pd.DataFrame) -> pd.DataFrame:
    rank = {term: index for index, term in enumerate(VEP_CONSEQUENCE_ORDER)}
    rows = []
    for variant_key, group in annotations.groupby("variant_key", sort=False):
        successful = group[group["status"] == "ok"]
        terms = {
            term
            for value in successful["consequence_terms"]
            for term in str(value).split("&")
            if term
        }
        ordered_terms = sorted(terms, key=lambda term: (rank.get(term, len(rank)), term))
        statuses = sorted({str(value) for value in group["status"] if str(value)})
        status = (
            "ok"
            if ordered_terms
            else (statuses[0] if len(statuses) == 1 else "no_consequence")
        )
        rows.append(
            {
                "variant_key": str(variant_key),
                "vep_status": status,
                "vep_primary_consequence": ordered_terms[0] if ordered_terms else "",
                "vep_consequence_terms": "|".join(ordered_terms),
            }
        )
    return pd.DataFrame(
        rows,
        columns=[
            "variant_key",
            "vep_status",
            "vep_primary_consequence",
            "vep_consequence_terms",
        ],
    )


def build_or_load_clinvar_universe(
    *,
    genes_tsv: Path | Sequence[Path],
    target_sequences_dir: Path | Sequence[Path],
    clinvar_vcf: Path,
    universe_path: Path,
    manifest_path: Path,
    regions_path: Path,
) -> dict:
    expected_inputs = {
        "genes_tsv": [path_metadata(path) for path in _paths(genes_tsv)],
        "target_sequences_dir": [
            directory_metadata(path, "*.fa.gz")
            for path in _paths(target_sequences_dir)
        ],
        "clinvar_vcf": path_metadata(clinvar_vcf),
        "clinvar_tbi": path_metadata(Path(f"{clinvar_vcf}.tbi")),
        "mode": "snv_indel",
        "cache_version": CACHE_VERSION,
    }
    if universe_path.exists() and manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if (
            manifest.get("complete") is True
            and manifest.get("inputs") == expected_inputs
            and manifest.get("output") == path_metadata(universe_path)
        ):
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
        "complete": True,
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
        "missing_molecular_consequence_count": sum(
            1 for row in rows if not row["clinvar_mc_terms"]
        ),
        "multiple_molecular_consequence_count": sum(
            1 for row in rows if len(str(row["clinvar_mc_terms"]).split("|")) > 1
        ),
        "regions_bed": str(regions_path),
        "universe_tsv": str(universe_path),
        "output": path_metadata(universe_path),
    }
    write_json_atomic(manifest_path, manifest)
    return manifest


def read_genes(path: Path | Sequence[Path]) -> list[dict[str, str]]:
    rows = []
    seen: set[str] = set()
    for item in _paths(path):
        with gzip.open(item, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"gene_id", "chromosome", "begin", "end"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"Genes table missing required columns: {', '.join(sorted(missing))}"
                )
            for row in reader:
                gene_id = str(row["gene_id"])
                if gene_id in seen:
                    raise ValueError(
                        f"Duplicate target Gene ID across source runs: {gene_id}"
                    )
                seen.add(gene_id)
                rows.append(row)
    return rows


def _paths(value: Path | Sequence[Path]) -> tuple[Path, ...]:
    return (value,) if isinstance(value, Path) else tuple(value)


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
    lines = [
        f"{chrom}\t{max(0, start1 - 1)}\t{end1}\n"
        for chrom, start1, end1 in intervals
    ]
    write_text_atomic(path, "".join(lines))


def query_clinvar_variant_universe(
    clinvar_vcf: Path,
    regions_path: Path,
    context_index: dict[str, tuple[list[dict], list[int]]],
) -> tuple[list[dict[str, str]], Counter]:
    tabix = shutil.which("tabix")
    if tabix is None:
        raise FileNotFoundError("tabix executable not found; it is required for indexed ClinVar queries.")

    raw_by_key: dict[str, dict[str, object]] = {}
    counts = Counter()
    with tabix_output_lines(tabix, clinvar_vcf, regions_path) as output_lines:
        for line in output_lines:
            line = line.rstrip("\n")
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
            disease_names = info_value(info_text, "CLNDN")
            disease_ids = info_value(info_text, "CLNDISDB")
            molecular_consequences = parse_molecular_consequences(info_value(info_text, "MC"))
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

                normalized_items = normalize_clinvar_allele_for_targets(
                    context_index, chrom, pos, ref, alt
                )
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
                            "clinvar_mc_so_ids": set(),
                            "clinvar_mc_terms": set(),
                            "gene_ids": set(),
                            "clinvar_condition_records": set(),
                        },
                    )
                    entry["labels"].add(label)
                    if rec_id and rec_id != ".":
                        entry["clinvar_ids"].add(rec_id)
                    if sig:
                        entry["clinvar_sigs"].add(sig)
                    if (disease_names and disease_names != ".") or (
                        disease_ids and disease_ids != "."
                    ):
                        entry["clinvar_condition_records"].add(
                            (disease_names or ".", disease_ids or ".")
                        )
                    for so_id, term in molecular_consequences:
                        if so_id:
                            entry["clinvar_mc_so_ids"].add(so_id)
                        if term:
                            entry["clinvar_mc_terms"].add(term)
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
        condition_records = sorted(
            entry["clinvar_condition_records"], key=lambda item: (item[1], item[0])
        )
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
                "clinvar_mc_so_ids": "|".join(sorted(entry["clinvar_mc_so_ids"])),
                "clinvar_mc_terms": "|".join(sorted(entry["clinvar_mc_terms"])),
                "gene_ids": "|".join(sorted(entry["gene_ids"], key=gene_sort_key)),
                "clinvar_disease_names": ";".join(item[0] for item in condition_records),
                "clinvar_disease_ids": ";".join(item[1] for item in condition_records),
            }
        )
    rows.sort(key=lambda row: (row["variant_type"], chrom_sort_key(str(row["chrom"])), int(row["pos"]), row["ref"], row["alt"]))
    return rows, counts


@contextmanager
def tabix_output_lines(
    tabix: str,
    clinvar_vcf: Path,
    regions_path: Path,
) -> Iterator[TextIO]:
    with tempfile.TemporaryFile(mode="w+", encoding="utf-8") as stderr_handle:
        proc = subprocess.Popen(
            [tabix, "-R", str(regions_path), str(clinvar_vcf)],
            text=True,
            encoding="utf-8",
            stdout=subprocess.PIPE,
            stderr=stderr_handle,
        )
        if proc.stdout is None:
            proc.kill()
            proc.wait()
            raise RuntimeError("tabix ClinVar query did not provide an output stream.")
        try:
            yield proc.stdout
        except BaseException:
            proc.kill()
            proc.wait()
            raise
        finally:
            proc.stdout.close()
        returncode = proc.wait()
        stderr_handle.seek(0)
        stderr = stderr_handle.read().strip()
    if returncode not in (0, 1):
        raise RuntimeError(f"tabix ClinVar query failed: {stderr}")


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


def parse_molecular_consequences(value: str) -> list[tuple[str, str]]:
    """Parse ClinVar's comma-separated MC values into SO identifier/term pairs."""
    consequences = []
    for item in str(value or "").split(","):
        item = item.strip()
        if not item or item == ".":
            continue
        so_id, separator, term = item.partition("|")
        if not separator:
            so_id, term = "", so_id
        pair = (so_id.strip(), term.strip())
        if pair not in consequences:
            consequences.append(pair)
    return consequences


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
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
        with temporary_path.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                compresslevel=6,
                fileobj=raw_handle,
                mtime=0,
            ) as gzip_handle:
                with io.TextIOWrapper(
                    gzip_handle,
                    encoding="utf-8",
                    newline="",
                ) as text_handle:
                    writer = csv.DictWriter(
                        text_handle,
                        fieldnames=UNIVERSE_FIELDS,
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerows(rows)
            raw_handle.flush()
            os.fsync(raw_handle.fileno())
        temporary_path.chmod(0o644)
        temporary_path.replace(path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def build_or_load_observed_keys_by_strategy_type(
    *,
    universe: pd.DataFrame,
    universe_path: Path,
    observed_store: ObservedVariantStore,
    strategies: list[str],
    analytics_dir: Path,
) -> tuple[dict[tuple[str, str], set[str]], dict[str, object], Path, Path]:
    output_path = analytics_dir / "clinvar_observed_memberships.tsv.gz"
    manifest_path = analytics_dir / "clinvar_observed_memberships.manifest.json"
    expected_inputs = {
        "cache_version": OBSERVED_MEMBERSHIP_CACHE_VERSION,
        "observed_store": {
            "manifest": path_metadata(observed_store.manifest_path),
            "alleles": path_metadata(observed_store.allele_path),
        },
        "clinvar_universe": path_metadata(universe_path),
        "strategies": sorted(strategies),
    }
    if output_path.exists() and manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
            if (
                manifest.get("complete") is True
                and manifest.get("inputs") == expected_inputs
                and manifest.get("output") == path_metadata(output_path)
            ):
                frame = pd.read_csv(
                    output_path,
                    sep="\t",
                    compression="gzip",
                    keep_default_na=False,
                )
                _validate_observed_membership_columns(frame)
                return (
                    _observed_membership_sets(frame, strategies),
                    {**manifest, "cache_hit": True},
                    output_path,
                    manifest_path,
                )
        except (OSError, json.JSONDecodeError, ValueError):
            pass

    observed = collect_observed_keys_by_strategy_type(
        universe=universe,
        observed_store=observed_store,
        strategies=strategies,
    )
    rows = [
        {
            "strategy": strategy,
            "variant_type": variant_type,
            "variant_key": variant_key,
        }
        for (strategy, variant_type), keys in sorted(observed.items())
        for variant_key in sorted(keys)
    ]
    frame = pd.DataFrame(rows, columns=OBSERVED_MEMBERSHIP_COLUMNS)
    write_tsv_atomic(output_path, frame)
    manifest = {
        "complete": True,
        "inputs": expected_inputs,
        "membership_count": int(len(frame)),
        "output": path_metadata(output_path),
    }
    write_json_atomic(manifest_path, manifest)
    return observed, {**manifest, "cache_hit": False}, output_path, manifest_path


def _validate_observed_membership_columns(frame: pd.DataFrame) -> None:
    missing = set(OBSERVED_MEMBERSHIP_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(
            "ClinVar observed-membership cache missing columns: "
            f"{', '.join(sorted(missing))}"
        )


def _observed_membership_sets(
    frame: pd.DataFrame,
    strategies: list[str],
) -> dict[tuple[str, str], set[str]]:
    observed = {
        (strategy, variant_type): set()
        for strategy in strategies
        for variant_type in VALIDATION_TYPES
    }
    for row in frame.itertuples(index=False):
        key = (str(row.strategy), str(row.variant_type))
        observed.setdefault(key, set()).add(str(row.variant_key))
    return observed


def collect_observed_keys_by_strategy_type(
    *,
    universe: pd.DataFrame,
    observed_store: ObservedVariantStore,
    strategies: list[str],
) -> dict[tuple[str, str], set[str]]:
    if universe.empty:
        return {
            (strategy, variant_kind): set()
            for strategy in strategies
            for variant_kind in VALIDATION_TYPES
        }

    types_by_key = dict(
        zip(
            universe["variant_key"].astype(str),
            universe["variant_type"].astype(str),
        )
    )
    observed_by_strategy_type: dict[tuple[str, str], set[str]] = {
        (strategy, variant_kind): set()
        for strategy in strategies
        for variant_kind in VALIDATION_TYPES
    }

    for variant_key, strategy in observed_store.observed_strategy_keys(
        universe["variant_key"],
        strategies,
    ):
        variant_kind = types_by_key.get(variant_key)
        if variant_kind in VALIDATION_TYPES:
            observed_by_strategy_type[(strategy, variant_kind)].add(variant_key)

    return observed_by_strategy_type
