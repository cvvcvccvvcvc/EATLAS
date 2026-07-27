#!/usr/bin/env python3
"""Fetch gnomAD variants for a specific gene/region via GraphQL API."""

import argparse
import json
import logging
import random
import sys
import time
import urllib.request
from urllib.error import HTTPError, URLError
from pathlib import Path
from typing import Any

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_MAX_ATTEMPTS = 5
GNOMAD_RETRY_BASE_SECONDS = 5.0
GNOMAD_RETRY_MAX_SECONDS = 30.0
GNOMAD_TRANSIENT_HTTP_STATUSES = {408, 429, 500, 502, 503, 504}
GNOMAD_REGION_MIN_WINDOW_BP = 500

GNOMAD_REGION_QUERY = """
query VariantsInRegion($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      chrom
      pos
      ref
      alt
      consequence
      hgvsc
      hgvsp
      exome { af }
      genome { af }
      joint { an ac }
    }
  }
}
"""

def _to_float(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None

def _to_int(value) -> int | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None

def _joint_af_metrics(variant: dict) -> tuple[int | None, int | None, float | None]:
    joint = variant.get("joint")
    if not isinstance(joint, dict):
        return None, None, None

    an = _to_int(joint.get("an"))
    ac_raw = joint.get("ac")
    ac = _to_int(ac_raw[0]) if isinstance(ac_raw, list) and ac_raw else _to_int(ac_raw)

    if an is None or an <= 0 or ac is None or ac < 0:
        return an, ac, None
    return an, ac, ac / an

def select_af_metrics(variant: dict):
    exome = variant.get("exome")
    genome = variant.get("genome")
    af_exome = _to_float(exome.get("af")) if isinstance(exome, dict) else None
    af_genome = _to_float(genome.get("af")) if isinstance(genome, dict) else None
    an_joint, ac_joint, af_joint = _joint_af_metrics(variant)
    
    if af_joint is not None:
        return af_joint, "joint", af_exome, af_genome, af_joint, an_joint, ac_joint
    if af_exome is not None:
        return af_exome, "exome", af_exome, af_genome, af_joint, an_joint, ac_joint
    if af_genome is not None:
        return af_genome, "genome", af_exome, af_genome, af_joint, an_joint, ac_joint
    return None, None, af_exome, af_genome, af_joint, an_joint, ac_joint


_select_af_metrics = select_af_metrics

def execute_graphql(query: str, variables: dict) -> dict:
    req = urllib.request.Request(
        GNOMAD_API_URL,
        method="POST",
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36"
        },
        data=json.dumps({"query": query, "variables": variables}).encode()
    )
    with urllib.request.urlopen(req, timeout=60) as response:
        return json.loads(response.read())

def retry_sleep_seconds(attempt: int) -> float:
    delay = min(GNOMAD_RETRY_MAX_SECONDS, GNOMAD_RETRY_BASE_SECONDS * (2 ** (attempt - 1)))
    return delay + random.uniform(0.0, delay * 0.2)


def is_retryable_network_error(exc: Exception) -> bool:
    if isinstance(exc, HTTPError):
        return exc.code in GNOMAD_TRANSIENT_HTTP_STATUSES
    return isinstance(exc, (URLError, TimeoutError))


def fetch_region_variants_recursive(
    chrom: str,
    start: int,
    stop: int,
    max_attempts: int = GNOMAD_MAX_ATTEMPTS,
) -> list[dict]:
    if max_attempts < 1:
        raise ValueError("max_attempts must be >= 1")

    variables = {"chrom": chrom, "start": start, "stop": stop}
    for attempt in range(1, max_attempts + 1):
        retry_error: Exception | None = None
        try:
            data = execute_graphql(GNOMAD_REGION_QUERY, variables)
        except (HTTPError, URLError, TimeoutError) as exc:
            if not is_retryable_network_error(exc):
                raise
            retry_error = exc
        else:
            if "errors" not in data:
                region = data.get("data", {}).get("region")
                return region.get("variants", []) if region else []

            messages = [str(error.get("message", error)) for error in data["errors"]]
            joined = " | ".join(message.lower() for message in messages)
            if (
                ("select a smaller region" in joined or "too many variants" in joined)
                and (stop - start + 1) > GNOMAD_REGION_MIN_WINDOW_BP
            ):
                mid = (start + stop) // 2
                logger.info(f"Splitting gnomAD region {chrom}:{start}-{stop} due to API constraint.")
                left = fetch_region_variants_recursive(chrom, start, mid, max_attempts)
                right = fetch_region_variants_recursive(chrom, mid + 1, stop, max_attempts)
                return left + right
            if "rate limit" not in joined:
                raise RuntimeError(f"GraphQL errors: {joined}")
            retry_error = RuntimeError(f"GraphQL errors: {joined}")

        if attempt == max_attempts:
            assert retry_error is not None
            raise retry_error

        delay = retry_sleep_seconds(attempt)
        logger.warning(
            "gnomAD request failed for %s:%s-%s on attempt %s/%s (%s: %s); "
            "retrying in %.1fs",
            chrom,
            start,
            stop,
            attempt,
            max_attempts,
            type(retry_error).__name__,
            retry_error,
            delay,
        )
        time.sleep(delay)

def write_vcf(variants: list, out_path: Path) -> None:
    with open(out_path, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=gnomAD_GraphQL_API\n")
        vcf.write("##reference=GRCh38\n")
        vcf.write('##INFO=<ID=GNOMAD_VID,Number=1,Type=String,Description="gnomAD variant identifier">\n')
        vcf.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        vcf.write('##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor Allele Frequency">\n')
        vcf.write('##INFO=<ID=AF_SOURCE,Number=1,Type=String,Description="AF source priority">\n')
        vcf.write('##INFO=<ID=AF_EXOME,Number=1,Type=Float,Description="Exome AF">\n')
        vcf.write('##INFO=<ID=AF_GENOME,Number=1,Type=Float,Description="Genome AF">\n')
        vcf.write('##INFO=<ID=AF_JOINT,Number=1,Type=Float,Description="Joint AF">\n')
        vcf.write('##INFO=<ID=AN_JOINT,Number=1,Type=Integer,Description="Joint AN">\n')
        vcf.write('##INFO=<ID=AC_JOINT,Number=1,Type=Integer,Description="Joint AC">\n')
        vcf.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence">\n')
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for v in variants:
            chrom = v["chrom"]
            pos = v["pos"]
            ref = v["ref"]
            alt = v["alt"]
            vid = str(v.get("variant_id", "."))
            af, af_source, af_exome, af_genome, af_joint, an_joint, ac_joint = _select_af_metrics(v)
            maf = min(af, 1.0 - af) if af is not None else None

            info_fields = []
            if vid and vid != ".": info_fields.append(f"GNOMAD_VID={vid}")
            if af_source: info_fields.append(f"AF_SOURCE={af_source}")
            if af is not None: info_fields.append(f"AF={af:.6g}")
            if maf is not None: info_fields.append(f"MAF={maf:.6g}")
            if af_exome is not None: info_fields.append(f"AF_EXOME={af_exome:.6g}")
            if af_genome is not None: info_fields.append(f"AF_GENOME={af_genome:.6g}")
            if af_joint is not None: info_fields.append(f"AF_JOINT={af_joint:.6g}")
            if an_joint is not None: info_fields.append(f"AN_JOINT={an_joint}")
            if ac_joint is not None: info_fields.append(f"AC_JOINT={ac_joint}")
            if v.get("consequence"): info_fields.append(f"CSQ={v['consequence']}")

            info = ";".join(info_fields) if info_fields else "."
            vcf.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t{info}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--out-vcf", required=True, type=Path)
    args = parser.parse_args()

    chrom = args.chrom
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
        
    start, end = min(args.start, args.end), max(args.start, args.end)
    logger.info(f"Fetching gnomAD variants for region {chrom}:{start}-{end}")
    
    # Pad region slightly to catch edge variants
    start = max(1, start - 100)
    end = end + 100

    variants = fetch_region_variants_recursive(chrom, start, end)
    
    # Deduplicate variants
    seen = set()
    deduped = []
    for v in variants:
        key = (v.get("chrom"), v.get("pos"), v.get("ref"), v.get("alt"))
        if key not in seen:
            seen.add(key)
            deduped.append(v)
            
    deduped.sort(key=lambda v: (str(v.get("chrom")), int(v.get("pos"))))
    write_vcf(deduped, args.out_vcf)
    logger.info(f"Saved {len(deduped)} variants to {args.out_vcf}")

if __name__ == "__main__":
    main()
