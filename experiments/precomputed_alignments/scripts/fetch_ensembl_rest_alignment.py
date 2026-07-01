#!/usr/bin/env python3
"""Fetch a small Ensembl Compara REST alignment region as JSON."""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", default="homo_sapiens")
    parser.add_argument("--region", required=True, help="Ensembl region, for example 4:122612500-122612650")
    parser.add_argument("--method", default="EPO")
    parser.add_argument("--species-set-group", default="mammals")
    parser.add_argument(
        "--display-species-set",
        action="append",
        default=[],
        help="Optional species filter. Can be supplied more than once.",
    )
    parser.add_argument("--server", default="https://rest.ensembl.org")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--retries", type=int, default=3)
    return parser.parse_args()


def build_url(args: argparse.Namespace) -> str:
    region = urllib.parse.quote(args.region, safe=":-")
    path = f"/alignment/region/{args.species}/{region}"
    query_items: list[tuple[str, str]] = [
        ("method", args.method),
        ("species_set_group", args.species_set_group),
        ("content-type", "application/json"),
    ]
    for species in args.display_species_set:
        query_items.append(("display_species_set", species))
    query = urllib.parse.urlencode(query_items)
    return args.server.rstrip("/") + path + "?" + query


def fetch_json(url: str, timeout: float, retries: int) -> object:
    headers = {
        "Accept": "application/json",
        "User-Agent": "gaph-precomputed-alignment-prototype/0.1",
    }
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code == 429 and attempt < retries:
                retry_after = exc.headers.get("Retry-After")
                delay = float(retry_after) if retry_after else 2.0 * attempt
                time.sleep(delay)
                continue
            raise RuntimeError(f"Ensembl REST request failed with HTTP {exc.code}: {exc.reason}") from exc
        except urllib.error.URLError as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(2.0 * attempt)
                continue
            raise RuntimeError(f"Ensembl REST request failed: {exc.reason}") from exc
    raise RuntimeError(f"Ensembl REST request failed: {last_error}")


def main() -> None:
    args = parse_args()
    url = build_url(args)
    payload = fetch_json(url, args.timeout, args.retries)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.output}")
    print(url)


if __name__ == "__main__":
    main()

