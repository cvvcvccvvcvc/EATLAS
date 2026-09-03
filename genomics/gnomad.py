"""gnomAD GraphQL client."""

import json
import logging
import random
import time
import urllib.request
from urllib.error import HTTPError, URLError

logger = logging.getLogger(__name__)

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"
GNOMAD_MAX_ATTEMPTS = 10
GNOMAD_RETRY_BASE_SECONDS = 5.0
GNOMAD_RETRY_MAX_SECONDS = 60.0
GNOMAD_TRANSIENT_HTTP_STATUSES = {408, 429, 500, 502, 503, 504}
GNOMAD_REGION_MIN_WINDOW_BP = 500

GNOMAD_REGION_QUERY = """
query VariantsInRegion($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: %s) {
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
""" % GNOMAD_DATASET


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
    return min(GNOMAD_RETRY_MAX_SECONDS, delay + random.uniform(0.0, delay * 0.2))


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
