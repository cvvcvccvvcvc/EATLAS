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

from bin.fetch_gnomad_variants import (
    GNOMAD_API_URL,
    fetch_region_variants_recursive,
    select_af_metrics,
)

from .clinvar_validation import path_metadata
from .variant_keys import normalize_chrom


CACHE_VERSION = 1
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
    if output_path.exists() and manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError):
            manifest = {}
        if manifest.get("complete") and manifest.get("inputs") == expected_inputs:
            return _read_evidence(output_path), manifest

    variants = (
        matched[["variant_key", "chrom", "pos", "ref", "alt"]]
        .drop_duplicates("variant_key")
        .sort_values(["chrom", "pos", "variant_key"], kind="mergesort")
        .reset_index(drop=True)
    )
    clinvar = _annotate_clinvar(variants, clinvar_vcf)
    gnomad, gnomad_summary = _annotate_gnomad(variants)
    evidence = clinvar.merge(gnomad, on="variant_key", how="outer", validate="one_to_one")
    evidence = variants[["variant_key"]].merge(
        evidence,
        on="variant_key",
        how="left",
        validate="one_to_one",
    )
    evidence.to_csv(output_path, sep="\t", index=False, compression="gzip", lineterminator="\n")

    manifest = {
        "inputs": expected_inputs,
        "complete": gnomad_summary["failed_region_count"] == 0,
        "unique_allele_count": int(len(evidence)),
        "clinvar_found_count": int(evidence["clinvar_found"].sum()),
        "clinvar_classified_count": int(evidence["clinvar_classified"].sum()),
        "gnomad_found_count": int(evidence["gnomad_found"].sum()),
        "gnomad": gnomad_summary,
        "evidence_tsv": str(output_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return evidence, manifest


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

    text = str(value or "").lower()
    if not text:
        return ""
    if "conflicting" in text:
        return "Other"
    if "uncertain" in text or "vus" in text:
        return "VUS"
    benign = "benign" in text
    pathogenic = "pathogenic" in text
    if benign and not pathogenic:
        return "B/LB"
    if pathogenic and not benign:
        return "P/LP"
    return "Other"


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
                "clinvar_class": categorize_clinvar_sig(sig),
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


def _annotate_gnomad(variants: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
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
    with concurrent.futures.ThreadPoolExecutor(max_workers=GNOMAD_WORKERS) as executor:
        futures = {
            executor.submit(
                fetch_region_variants_recursive,
                chrom,
                max(1, start - 100),
                end + 100,
            ): (chrom, start, end)
            for chrom, start, end in tasks
        }
        for future in concurrent.futures.as_completed(futures):
            chrom, start, end = futures[future]
            region_keys = keys_by_task[(chrom, start, end)]
            try:
                records = future.result()
            except Exception as exc:  # Network failures are recorded, not interpreted as absence.
                errors.append(
                    {"chrom": chrom, "start": start, "end": end, "error": f"{type(exc).__name__}: {exc}"}
                )
                for variant_key in region_keys:
                    status_by_key[variant_key] = "error"
                continue

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
    }
    return pd.DataFrame(rows), summary
