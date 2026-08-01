"""ClinVar and gnomAD evidence for consequence-matched target-space alleles."""

from __future__ import annotations

import bisect
import concurrent.futures
import json
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

from analytics.io.artifacts import path_metadata, write_json_atomic, write_tsv_atomic
from analytics.io.performance import PerformanceProfile, profile_stage
from genomics.clinvar import record_category, significance_class
from genomics.gnomad_cache import GnomadRegionCache
from genomics.gnomad import (
    GNOMAD_API_URL,
    fetch_region_variants_recursive,
    select_af_metrics,
)

from genomics.variants import normalize_chrom


CACHE_VERSION = 2
GNOMAD_DATASET = "gnomad_r4"
GNOMAD_CLUSTER_GAP_BP = 200_000
GNOMAD_WORKERS = 5
EVIDENCE_COLUMNS = [
    "variant_key",
    "clinvar_found",
    "clinvar_classified",
    "clinvar_class",
    "gnomad_status",
    "gnomad_found",
    "gnomad_af",
]


def build_external_evidence(
    *,
    matched: pd.DataFrame,
    matched_path: Path,
    clinvar_vcf: Path,
    output_path: Path,
    manifest_path: Path,
    gnomad_cache_dir: Path | None = None,
    performance_profile: PerformanceProfile | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Build or load exact-allele ClinVar and gnomAD annotations."""

    expected_inputs = {
        "cache_version": CACHE_VERSION,
        "matched_tsv": path_metadata(matched_path),
        "clinvar_vcf": path_metadata(clinvar_vcf),
        "clinvar_tbi": path_metadata(Path(f"{clinvar_vcf}.tbi")),
        "gnomad_api": GNOMAD_API_URL,
        "gnomad_dataset": GNOMAD_DATASET,
    }
    cached_evidence = None
    cached_manifest: dict[str, object] = {}
    if output_path.exists() and manifest_path.exists():
        try:
            cached_manifest = json.loads(manifest_path.read_text())
            if (
                cached_manifest.get("inputs") == expected_inputs
                and cached_manifest.get("output") == path_metadata(output_path)
            ):
                cached_evidence = _read_evidence(output_path)
        except (OSError, json.JSONDecodeError):
            cached_manifest = {}
        if cached_evidence is not None and cached_manifest.get("complete"):
            return cached_evidence, {**cached_manifest, "cache_hit": True}

    variants = (
        matched[["variant_key", "chrom", "pos", "ref", "alt"]]
        .drop_duplicates("variant_key")
        .sort_values(["chrom", "pos", "variant_key"], kind="mergesort")
        .reset_index(drop=True)
    )
    if cached_evidence is not None and set(cached_evidence["variant_key"].astype(str)) != set(
        variants["variant_key"].astype(str)
    ):
        cached_evidence = None

    if cached_evidence is None:
        with profile_stage(performance_profile, "External evidence ClinVar lookup") as timing:
            clinvar = _annotate_clinvar(variants, clinvar_vcf)
            timing["metrics"] = {
                "requested_alleles": int(len(variants)),
                "found_alleles": int(clinvar["clinvar_found"].sum()),
            }
        gnomad = _empty_gnomad_evidence(variants)
    else:
        clinvar = cached_evidence[["variant_key", "clinvar_found", "clinvar_classified", "clinvar_class"]]
        gnomad = cached_evidence[["variant_key", "gnomad_status", "gnomad_found", "gnomad_af"]].copy()

    unresolved_keys = set(
        gnomad.loc[~gnomad["gnomad_status"].eq("ok"), "variant_key"].astype(str)
    )
    unresolved = variants[variants["variant_key"].astype(str).isin(unresolved_keys)]
    with profile_stage(performance_profile, "External evidence gnomAD lookup") as timing:
        fetched, gnomad_summary = _annotate_gnomad(
            unresolved,
            gnomad_cache_dir=gnomad_cache_dir,
        )
        shared_cache = gnomad_summary.get("shared_cache", {})
        timing["metrics"] = {
            "requested_alleles": int(len(unresolved)),
            "regions": int(gnomad_summary.get("region_count", 0)),
            "raw_variants": int(gnomad_summary.get("raw_variant_count", 0)),
            "tile_hits": int(shared_cache.get("tile_hit_count", 0)),
            "tile_misses": int(shared_cache.get("tile_miss_count", 0)),
        }
    if not fetched.empty:
        fetched = fetched.copy()
        fetched["gnomad_af"] = pd.to_numeric(fetched["gnomad_af"], errors="coerce")
        gnomad = gnomad.set_index("variant_key")
        fetched = fetched.set_index("variant_key")
        gnomad.loc[fetched.index, ["gnomad_status", "gnomad_found", "gnomad_af"]] = fetched[
            ["gnomad_status", "gnomad_found", "gnomad_af"]
        ]
        gnomad = gnomad.reset_index()

    evidence = clinvar.merge(gnomad, on="variant_key", how="outer", validate="one_to_one")
    evidence = variants[["variant_key"]].merge(
        evidence,
        on="variant_key",
        how="left",
        validate="one_to_one",
    )
    write_tsv_atomic(output_path, evidence)

    failed_allele_count = int((~evidence["gnomad_status"].eq("ok")).sum())
    gnomad_summary.update(
        {
            "queried_allele_count": int(len(unresolved)),
            "cached_ok_allele_count": int(
                0
                if cached_evidence is None
                else cached_evidence["gnomad_status"].eq("ok").sum()
            ),
            "failed_allele_count": failed_allele_count,
        }
    )

    manifest = {
        "inputs": expected_inputs,
        "complete": failed_allele_count == 0,
        "unique_allele_count": int(len(evidence)),
        "clinvar_found_count": int(evidence["clinvar_found"].sum()),
        "clinvar_classified_count": int(evidence["clinvar_classified"].sum()),
        "gnomad_found_count": int(evidence["gnomad_found"].sum()),
        "gnomad": gnomad_summary,
        "evidence_tsv": str(output_path),
        "output": path_metadata(output_path),
    }
    write_json_atomic(manifest_path, manifest)
    return evidence, {**manifest, "cache_hit": False}


def _empty_gnomad_evidence(variants: pd.DataFrame) -> pd.DataFrame:
    frame = variants[["variant_key"]].copy()
    frame["gnomad_status"] = "error"
    frame["gnomad_found"] = False
    frame["gnomad_af"] = float("nan")
    return frame




def _read_evidence(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
    missing = set(EVIDENCE_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(f"External-evidence cache missing columns: {', '.join(sorted(missing))}")
    for column in ["clinvar_found", "clinvar_classified", "gnomad_found"]:
        frame[column] = frame[column].astype(str).str.lower().isin({"true", "1"})
    frame["gnomad_af"] = pd.to_numeric(frame["gnomad_af"], errors="coerce")
    return frame


def categorize_clinvar_sig(value: str) -> str:
    """Map unambiguous ClinVar significance text to a report class."""
    return significance_class(value) or ""


def _annotate_clinvar(variants: pd.DataFrame, clinvar_vcf: Path) -> pd.DataFrame:
    tabix = shutil.which("tabix")
    if tabix is None:
        raise FileNotFoundError("tabix executable not found; it is required for ClinVar lookup.")
    if not clinvar_vcf.exists() or not Path(f"{clinvar_vcf}.tbi").exists():
        raise FileNotFoundError(f"Indexed ClinVar VCF not found: {clinvar_vcf}")

    wanted = {
        (normalize_chrom(row.chrom) or "", int(row.pos), str(row.ref), str(row.alt)): str(row.variant_key)
        for row in variants.itertuples(index=False)
    }
    signatures: dict[str, set[str]] = defaultdict(set)
    found: set[str] = set()
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as handle:
        regions_path = Path(handle.name)
        for chrom, pos in sorted({(key[0], key[1]) for key in wanted}):
            handle.write(f"{chrom}\t{pos - 1}\t{pos}\n")
    try:
        proc = subprocess.run(
            [tabix, "-R", str(regions_path), str(clinvar_vcf)],
            text=True,
            capture_output=True,
            check=False,
        )
    finally:
        regions_path.unlink(missing_ok=True)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"tabix ClinVar query failed: {proc.stderr.strip()}")

    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        chrom = normalize_chrom(fields[0]) or ""
        pos = int(fields[1])
        ref = fields[3].upper()
        sig = _info_value(fields[7], "CLNSIG")
        for alt in fields[4].upper().split(","):
            variant_key = wanted.get((chrom, pos, ref, alt))
            if variant_key is None:
                continue
            found.add(variant_key)
            if sig:
                signatures[variant_key].add(sig)

    rows = []
    for variant_key in variants["variant_key"].astype(str):
        sig = "|".join(sorted(signatures.get(variant_key, set())))
        rows.append(
            {
                "variant_key": variant_key,
                "clinvar_found": variant_key in found,
                "clinvar_classified": bool(sig),
                "clinvar_class": record_category(sig, found=variant_key in found),
            }
        )
    return pd.DataFrame(rows)


def _info_value(info_text: str, key: str) -> str:
    prefix = f"{key}="
    for item in info_text.split(";"):
        if item.startswith(prefix):
            return item[len(prefix) :]
    return ""


def _cluster_positions(positions: list[int]) -> list[tuple[int, int]]:
    if not positions:
        return []
    ordered = sorted(set(positions))
    clusters = []
    start = previous = ordered[0]
    for position in ordered[1:]:
        if position - previous > GNOMAD_CLUSTER_GAP_BP:
            clusters.append((start, previous))
            start = position
        previous = position
    clusters.append((start, previous))
    return clusters


def _annotate_gnomad(
    variants: pd.DataFrame,
    gnomad_cache_dir: Path | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    region_cache = GnomadRegionCache(
        gnomad_cache_dir,
        fetcher=fetch_region_variants_recursive,
    )
    wanted = {
        (normalize_chrom(row.chrom) or "", int(row.pos), str(row.ref), str(row.alt)): str(row.variant_key)
        for row in variants.itertuples(index=False)
    }
    positions_by_chrom: dict[str, list[int]] = defaultdict(list)
    for chrom, pos, _ref, _alt in wanted:
        positions_by_chrom[chrom].append(pos)
    tasks = []
    keys_by_task: dict[tuple[str, int, int], set[str]] = {}
    for chrom, positions in positions_by_chrom.items():
        ordered = sorted(set(positions))
        keys_by_position: dict[int, set[str]] = defaultdict(set)
        for (key_chrom, pos, _ref, _alt), variant_key in wanted.items():
            if key_chrom == chrom:
                keys_by_position[pos].add(variant_key)
        for start, end in _cluster_positions(ordered):
            task = (chrom, start, end)
            tasks.append(task)
            left = bisect.bisect_left(ordered, start)
            right = bisect.bisect_right(ordered, end)
            keys_by_task[task] = {
                variant_key
                for position in ordered[left:right]
                for variant_key in keys_by_position[position]
            }

    status_by_key = {variant_key: "" for variant_key in wanted.values()}
    found_by_key: dict[str, float | None] = {}
    errors = []
    raw_variant_count = 0
    if tasks:
        task_iterator = iter(tasks)
        worker_count = min(GNOMAD_WORKERS, len(tasks))
        with concurrent.futures.ThreadPoolExecutor(max_workers=worker_count) as executor:
            pending = {}
            for _ in range(worker_count):
                chrom, start, end = next(task_iterator)
                future = executor.submit(
                    region_cache.fetch_region,
                    chrom,
                    max(1, start - 100),
                    end + 100,
                )
                pending[future] = (chrom, start, end)

            while pending:
                completed, _remaining = concurrent.futures.wait(
                    pending,
                    return_when=concurrent.futures.FIRST_COMPLETED,
                )
                for future in completed:
                    chrom, start, end = pending.pop(future)
                    region_keys = keys_by_task[(chrom, start, end)]
                    try:
                        records = future.result()
                    except Exception as exc:  # Failures are not interpreted as absence.
                        errors.append(
                            {
                                "chrom": chrom,
                                "start": start,
                                "end": end,
                                "error": f"{type(exc).__name__}: {exc}",
                            }
                        )
                        for variant_key in region_keys:
                            status_by_key[variant_key] = "error"
                    else:
                        raw_variant_count += len(records)
                        for variant_key in region_keys:
                            status_by_key[variant_key] = "ok"
                        for record in records:
                            key = (
                                normalize_chrom(record.get("chrom")) or "",
                                int(record.get("pos", 0)),
                                str(record.get("ref", "")).upper(),
                                str(record.get("alt", "")).upper(),
                            )
                            variant_key = wanted.get(key)
                            if variant_key is None:
                                continue
                            af, _source, *_rest = select_af_metrics(record)
                            found_by_key[variant_key] = af

                    try:
                        next_chrom, next_start, next_end = next(task_iterator)
                    except StopIteration:
                        continue
                    next_future = executor.submit(
                        region_cache.fetch_region,
                        next_chrom,
                        max(1, next_start - 100),
                        next_end + 100,
                    )
                    pending[next_future] = (next_chrom, next_start, next_end)

    rows = []
    for variant_key in variants["variant_key"].astype(str):
        found = variant_key in found_by_key
        rows.append(
            {
                "variant_key": variant_key,
                "gnomad_status": status_by_key[variant_key],
                "gnomad_found": found,
                "gnomad_af": found_by_key.get(variant_key),
            }
        )
    summary = {
        "dataset": GNOMAD_DATASET,
        "region_count": len(tasks),
        "successful_region_count": len(tasks) - len(errors),
        "failed_region_count": len(errors),
        "raw_variant_count": raw_variant_count,
        "errors": errors,
        "shared_cache": region_cache.snapshot(),
    }
    return pd.DataFrame(rows), summary
