"""Compact, resumable Ensembl VEP consequence annotation for analytics."""

from __future__ import annotations

import concurrent.futures
import hashlib
import json
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd


VEP_BASE_URL = "https://rest.ensembl.org"
VEP_BATCH_SIZE = 200
VEP_OPTIONS = {
    "canonical": 1,
    "mane": 1,
    "pick_allele_gene": 1,
    "refseq": 1,
    "variant_class": 1,
}

# Ensembl's documented order, most to least severe. Unknown terms sort last.
VEP_CONSEQUENCE_ORDER = (
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
    "sequence_variant",
)
_CONSEQUENCE_RANK = {term: rank for rank, term in enumerate(VEP_CONSEQUENCE_ORDER)}


def annotate_vep_consequences(
    rows: pd.DataFrame,
    cache_path: Path,
    *,
    base_url: str = VEP_BASE_URL,
    release: str | None = None,
    max_workers: int = 2,
    retries: int = 4,
    timeout_seconds: float = 120.0,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Annotate unique SNV/gene pairs and persist each completed REST batch."""

    required = {"variant_key", "gene_id", "chrom", "pos", "ref", "alt"}
    missing = required - set(rows.columns)
    if missing:
        raise ValueError(f"VEP input missing columns: {', '.join(sorted(missing))}")
    if max_workers < 1:
        raise ValueError("VEP max_workers must be >= 1")

    unique = (
        rows[list(required)]
        .astype({"variant_key": str, "gene_id": str, "chrom": str, "ref": str, "alt": str})
        .drop_duplicates(["variant_key", "gene_id"])
        .sort_values(["chrom", "pos", "variant_key", "gene_id"], kind="mergesort")
        .reset_index(drop=True)
    )
    if unique.empty:
        return _empty_annotations(), {"status": "complete", "requested": 0, "cached": 0, "queried": 0}

    release = release or _fetch_release(base_url, retries, timeout_seconds)
    config = {
        "base_url": base_url.rstrip("/"),
        "options": VEP_OPTIONS,
        "release": str(release),
        "schema_version": 1,
    }
    config_hash = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(cache_path) as connection:
        _initialize_cache(connection)
        wanted = list(unique[["variant_key", "gene_id"]].itertuples(index=False, name=None))
        cached = _load_cached(connection, config_hash, wanted)
        missing_rows = [
            row
            for row in unique.to_dict(orient="records")
            if (str(row["variant_key"]), str(row["gene_id"])) not in cached
        ]

        batches = [
            missing_rows[index : index + VEP_BATCH_SIZE]
            for index in range(0, len(missing_rows), VEP_BATCH_SIZE)
        ]
        if batches:
            with concurrent.futures.ThreadPoolExecutor(max_workers=min(max_workers, len(batches))) as executor:
                futures = [
                    executor.submit(
                        _request_batch,
                        batch,
                        base_url,
                        retries,
                        timeout_seconds,
                    )
                    for batch in batches
                ]
                for future in concurrent.futures.as_completed(futures):
                    annotations = future.result()
                    _store_annotations(connection, config_hash, str(release), annotations)
                    cached.update(
                        {
                            (str(item["variant_key"]), str(item["gene_id"])): item
                            for item in annotations
                        }
                    )

    ordered = [cached[(variant_key, gene_id)] for variant_key, gene_id in wanted]
    result = pd.DataFrame(ordered, columns=_annotation_columns())
    status_counts = result["status"].value_counts().sort_index().to_dict()
    return result, {
        "status": "complete",
        "base_url": base_url.rstrip("/"),
        "release": str(release),
        "options": VEP_OPTIONS,
        "requested": len(wanted),
        "cached": len(wanted) - len(missing_rows),
        "queried": len(missing_rows),
        "batch_count": len(batches),
        "status_counts": {str(key): int(value) for key, value in status_counts.items()},
        "cache_path": str(cache_path),
    }


def _annotation_columns() -> list[str]:
    return [
        "variant_key",
        "gene_id",
        "status",
        "primary_consequence",
        "consequence_terms",
        "transcript_id",
        "mane_select",
        "canonical",
        "impact",
        "variant_class",
    ]


def _empty_annotations() -> pd.DataFrame:
    return pd.DataFrame(columns=_annotation_columns())


def _fetch_release(base_url: str, retries: int, timeout_seconds: float) -> str:
    payload = _request_json(
        f"{base_url.rstrip('/')}/info/software",
        data=None,
        retries=retries,
        timeout_seconds=timeout_seconds,
    )
    if isinstance(payload, dict) and payload.get("release") is not None:
        return str(payload["release"])
    raise RuntimeError("Ensembl REST /info/software did not return a release.")


def _request_batch(
    rows: list[dict[str, object]],
    base_url: str,
    retries: int,
    timeout_seconds: float,
) -> list[dict[str, object]]:
    request_ids: dict[str, dict[str, object]] = {}
    variants = []
    for index, row in enumerate(rows):
        request_id = f"gaph_{index:03d}"
        request_ids[request_id] = row
        variants.append(
            f"{str(row['chrom']).removeprefix('chr')} {int(row['pos'])} {request_id} "
            f"{str(row['ref']).upper()} {str(row['alt']).upper()} . . ."
        )
    query = urllib.parse.urlencode(VEP_OPTIONS)
    payload = _request_json(
        f"{base_url.rstrip('/')}/vep/homo_sapiens/region?{query}",
        data={"variants": variants},
        retries=retries,
        timeout_seconds=timeout_seconds,
    )
    if not isinstance(payload, list):
        raise RuntimeError("Ensembl VEP returned a non-list response.")

    responses: dict[str, dict[str, object]] = {}
    for record in payload:
        if not isinstance(record, dict):
            continue
        fields = str(record.get("input") or "").split()
        request_id = fields[2] if len(fields) >= 3 else str(record.get("id") or "")
        if request_id in request_ids:
            responses[request_id] = record

    annotations = []
    for request_id, row in request_ids.items():
        annotations.append(_parse_record(row, responses.get(request_id)))
    return annotations


def _parse_record(row: dict[str, object], record: dict[str, object] | None) -> dict[str, object]:
    base = {
        "variant_key": str(row["variant_key"]),
        "gene_id": str(row["gene_id"]),
        "status": "no_response",
        "primary_consequence": "",
        "consequence_terms": "",
        "transcript_id": "",
        "mane_select": "",
        "canonical": False,
        "impact": "",
        "variant_class": "",
    }
    if record is None:
        return base

    target = [
        item
        for item in record.get("transcript_consequences", [])
        if isinstance(item, dict) and str(item.get("gene_id") or "") == str(row["gene_id"])
    ]
    if not target:
        return {**base, "status": "no_target_gene"}
    item = target[0]
    terms = sorted(
        {str(term) for term in item.get("consequence_terms", []) if str(term)},
        key=lambda term: (_CONSEQUENCE_RANK.get(term, len(_CONSEQUENCE_RANK)), term),
    )
    if not terms:
        return {**base, "status": "no_consequence"}
    canonical = str(item.get("canonical") or "").lower() in {"1", "true", "yes"}
    return {
        **base,
        "status": "ok",
        "primary_consequence": terms[0],
        "consequence_terms": "&".join(terms),
        "transcript_id": str(item.get("transcript_id") or ""),
        "mane_select": str(item.get("mane_select") or ""),
        "canonical": canonical,
        "impact": str(item.get("impact") or ""),
        "variant_class": str(record.get("variant_class") or ""),
    }


def _request_json(
    url: str,
    *,
    data: dict[str, object] | None,
    retries: int,
    timeout_seconds: float,
) -> object:
    body = json.dumps(data).encode() if data is not None else None
    headers = {"Accept": "application/json", "User-Agent": "gaph-analytics/1.0"}
    if body is not None:
        headers["Content-Type"] = "application/json"
    for attempt in range(retries + 1):
        request = urllib.request.Request(url, data=body, headers=headers, method="POST" if body else "GET")
        try:
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                return json.loads(response.read().decode())
        except urllib.error.HTTPError as exc:
            retryable = exc.code == 429 or 500 <= exc.code < 600
            if not retryable or attempt >= retries:
                detail = exc.read().decode(errors="replace")[:500]
                raise RuntimeError(f"Ensembl REST request failed ({exc.code}): {detail}") from exc
            retry_after = exc.headers.get("Retry-After")
            try:
                delay = float(retry_after) if retry_after else min(30.0, 2.0**attempt)
            except ValueError:
                delay = min(30.0, 2.0**attempt)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            if attempt >= retries:
                raise RuntimeError(f"Ensembl REST request failed: {exc}") from exc
            delay = min(30.0, 2.0**attempt)
        time.sleep(delay)
    raise AssertionError("unreachable")


def _initialize_cache(connection: sqlite3.Connection) -> None:
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS consequences (
            config_hash TEXT NOT NULL,
            vep_release TEXT NOT NULL,
            variant_key TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            status TEXT NOT NULL,
            primary_consequence TEXT NOT NULL,
            consequence_terms TEXT NOT NULL,
            transcript_id TEXT NOT NULL,
            mane_select TEXT NOT NULL,
            canonical INTEGER NOT NULL,
            impact TEXT NOT NULL,
            variant_class TEXT NOT NULL,
            PRIMARY KEY (config_hash, variant_key, gene_id)
        )
        """
    )
    connection.commit()


def _load_cached(
    connection: sqlite3.Connection,
    config_hash: str,
    keys: list[tuple[str, str]],
) -> dict[tuple[str, str], dict[str, object]]:
    result: dict[tuple[str, str], dict[str, object]] = {}
    columns = _annotation_columns()
    selected = ", ".join(columns)
    for index in range(0, len(keys), 400):
        chunk = keys[index : index + 400]
        placeholders = ",".join(["(?, ?)"] * len(chunk))
        parameters = [config_hash, *(value for pair in chunk for value in pair)]
        query = (
            f"SELECT {selected} FROM consequences WHERE config_hash = ? "
            f"AND (variant_key, gene_id) IN ({placeholders})"
        )
        for values in connection.execute(query, parameters):
            item = dict(zip(columns, values))
            item["canonical"] = bool(item["canonical"])
            result[(str(item["variant_key"]), str(item["gene_id"]))] = item
    return result


def _store_annotations(
    connection: sqlite3.Connection,
    config_hash: str,
    release: str,
    rows: list[dict[str, object]],
) -> None:
    columns = _annotation_columns()
    connection.executemany(
        f"""
        INSERT OR REPLACE INTO consequences (
            config_hash, vep_release, {', '.join(columns)}
        ) VALUES ({', '.join(['?'] * (len(columns) + 2))})
        """,
        [
            (
                config_hash,
                release,
                *(int(row[column]) if column == "canonical" else row[column] for column in columns),
            )
            for row in rows
        ],
    )
    connection.commit()
