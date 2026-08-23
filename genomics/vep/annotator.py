"""Compact, resumable Ensembl VEP consequence annotation for analytics."""

from __future__ import annotations

import csv
import concurrent.futures
import hashlib
import json
import sqlite3
import subprocess
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

from .terms import VEP_CONSEQUENCE_ORDER
from .result_cache import DEFAULT_TILE_SIZE_BP, VepResultCache


VEP_BASE_URL = "https://rest.ensembl.org"
VEP_BATCH_SIZE = 200
VEP_ASSEMBLY = "GRCh38"
VEP_OPTIONS = {
    "canonical": 1,
    "mane": 1,
    "pick_allele_gene": 1,
    "refseq": 1,
    "variant_class": 1,
}

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
    backend: str = "rest",
    vep_executable: str | Path = "vep",
    vep_cache_dir: Path | None = None,
    vep_forks: int = 1,
    vep_result_cache_dir: Path | None = None,
    vep_result_cache_tile_size_bp: int = DEFAULT_TILE_SIZE_BP,
    publish_vep_result_cache: bool = True,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Annotate unique variant/gene pairs and persist completed work."""

    required = {"variant_key", "gene_id", "chrom", "pos", "ref", "alt"}
    missing = required - set(rows.columns)
    if missing:
        raise ValueError(f"VEP input missing columns: {', '.join(sorted(missing))}")
    if max_workers < 1:
        raise ValueError("VEP max_workers must be >= 1")
    if backend not in {"rest", "local"}:
        raise ValueError("VEP backend must be 'rest' or 'local'")
    if vep_forks < 1:
        raise ValueError("VEP forks must be >= 1")

    if backend == "rest":
        release = release or _fetch_release(base_url, retries, timeout_seconds)
        # Preserve the original REST hash so existing caches remain reusable.
        config = vep_result_cache_config(
            backend=backend,
            release=str(release),
            base_url=base_url,
        )
    else:
        if release is None:
            raise ValueError("Local VEP requires an explicit release")
        if vep_cache_dir is None:
            raise ValueError("Local VEP requires a cache directory")
        vep_cache_dir = Path(vep_cache_dir).expanduser()
        if not vep_cache_dir.is_dir():
            raise FileNotFoundError(f"Local VEP cache directory does not exist: {vep_cache_dir}")
        config = vep_result_cache_config(
            backend=backend,
            release=str(release),
            base_url=base_url,
        )

    unique = (
        rows[list(required)]
        .astype({"variant_key": str, "gene_id": str, "chrom": str, "ref": str, "alt": str})
        .drop_duplicates(["variant_key", "gene_id"])
        .sort_values(["chrom", "pos", "variant_key", "gene_id"], kind="mergesort")
        .reset_index(drop=True)
    )
    if unique.empty:
        return _empty_annotations(), {
            "status": "complete",
            "backend": backend,
            "release": str(release),
            "options": VEP_OPTIONS,
            "requested": 0,
            "cached": 0,
            "queried": 0,
            "status_counts": {},
        }
    config_hash = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()

    result_cache = (
        VepResultCache(
            vep_result_cache_dir,
            config=config,
            tile_size_bp=vep_result_cache_tile_size_bp,
        )
        if vep_result_cache_dir is not None
        else None
    )
    shared_cached: dict[tuple[str, str], dict[str, object]] = {}
    shared_lookup: dict[str, object] | None = None
    if result_cache is not None:
        shared_frame, shared_lookup = result_cache.lookup(unique)
        shared_cached = _annotation_map(shared_frame)

    wanted = list(unique[["variant_key", "gene_id"]].itertuples(index=False, name=None))
    shared_cached_count = len(shared_cached)
    cached = shared_cached
    local_cached_count = 0
    missing_rows: list[dict[str, object]] = []
    shared_missing_keys: list[tuple[str, str]] = []
    batch_count = 0

    # Shared cache entries are immutable and authoritative. Only consult the
    # per-run cache for shared misses; a complete shared hit therefore requires
    # no SQLite reads or writes.
    if shared_cached_count < len(wanted):
        shared_missing_rows = [
            row
            for row in unique.to_dict(orient="records")
            if (str(row["variant_key"]), str(row["gene_id"])) not in shared_cached
        ]
        shared_missing_keys = [
            (str(row["variant_key"]), str(row["gene_id"]))
            for row in shared_missing_rows
        ]
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with sqlite3.connect(cache_path) as connection:
            _initialize_cache(connection)
            local_cached = _load_cached(
                connection,
                config_hash,
                shared_missing_keys,
            )
            local_cached_count = len(local_cached)
            cached.update(local_cached)
            missing_rows = [
                row
                for row in shared_missing_rows
                if (str(row["variant_key"]), str(row["gene_id"])) not in local_cached
            ]

            if missing_rows and backend == "local":
                annotations = _run_local_vep(
                    missing_rows,
                    executable=vep_executable,
                    cache_dir=vep_cache_dir,
                    release=str(release),
                    forks=vep_forks,
                    temporary_root=cache_path.parent,
                )
                _store_annotations(connection, config_hash, str(release), annotations)
                cached.update(
                    {
                        (str(item["variant_key"]), str(item["gene_id"])): item
                        for item in annotations
                    }
                )
                batch_count = 1
            elif missing_rows:
                batches = [
                    missing_rows[index : index + VEP_BATCH_SIZE]
                    for index in range(0, len(missing_rows), VEP_BATCH_SIZE)
                ]
                batch_count = len(batches)
                worker_count = min(max_workers, len(batches))
                with concurrent.futures.ThreadPoolExecutor(max_workers=worker_count) as executor:
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
    shared_publish: dict[str, object] | None = None
    if result_cache is not None and publish_vep_result_cache and shared_missing_keys:
        publish_rows = pd.DataFrame(
            [cached[key] for key in shared_missing_keys],
            columns=_annotation_columns(),
        )
        shared_publish = result_cache.publish(publish_rows)
    status_counts = result["status"].value_counts().sort_index().to_dict()
    summary = {
        "status": "complete",
        "backend": backend,
        "release": str(release),
        "options": VEP_OPTIONS,
        "requested": len(wanted),
        "cached": len(wanted) - len(missing_rows),
        "shared_cached": shared_cached_count,
        "local_cached": local_cached_count,
        "queried": len(missing_rows),
        "batch_count": batch_count,
        "status_counts": {str(key): int(value) for key, value in status_counts.items()},
        "cache_path": str(cache_path),
    }
    if result_cache is not None:
        summary["shared_cache"] = {
            "lookup": shared_lookup,
            "publish": shared_publish,
        }
    if backend == "rest":
        summary["base_url"] = base_url.rstrip("/")
    else:
        summary.update(
            {
                "assembly": VEP_ASSEMBLY,
                "vep_cache_dir": str(vep_cache_dir),
                "vep_executable": str(vep_executable),
                "vep_forks": vep_forks,
            }
        )
    return result, summary


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


def vep_result_cache_config(
    *,
    backend: str,
    release: str,
    base_url: str = VEP_BASE_URL,
) -> dict[str, object]:
    """Return the semantic VEP configuration used to namespace shared results."""

    if backend == "rest":
        return {
            "base_url": base_url.rstrip("/"),
            "options": VEP_OPTIONS,
            "release": str(release),
            "schema_version": 1,
        }
    if backend == "local":
        return {
            "assembly": VEP_ASSEMBLY,
            "backend": "local",
            "options": VEP_OPTIONS,
            "refseq": True,
            "release": str(release),
            "schema_version": 1,
            "use_given_ref": True,
        }
    raise ValueError("VEP backend must be 'rest' or 'local'")


def _empty_annotations() -> pd.DataFrame:
    return pd.DataFrame(columns=_annotation_columns())


def _annotation_map(
    frame: pd.DataFrame,
) -> dict[tuple[str, str], dict[str, object]]:
    return {
        (str(item["variant_key"]), str(item["gene_id"])): item
        for item in frame[_annotation_columns()].to_dict(orient="records")
    }


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


def _run_local_vep(
    rows: list[dict[str, object]],
    *,
    executable: str | Path,
    cache_dir: Path,
    release: str,
    forks: int,
    temporary_root: Path,
) -> list[dict[str, object]]:
    request_ids = {f"gaph_{index:08d}": row for index, row in enumerate(rows)}
    temporary_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".vep_", dir=temporary_root) as directory:
        directory_path = Path(directory)
        input_path = directory_path / "input.vcf"
        output_path = directory_path / "output.tsv"
        with input_path.open("w", newline="") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
            for request_id, row in request_ids.items():
                ref = _validate_vcf_allele(row["ref"], "REF")
                alt = _validate_vcf_allele(row["alt"], "ALT")
                writer.writerow(
                    [
                        str(row["chrom"]).removeprefix("chr"),
                        int(row["pos"]),
                        request_id,
                        ref,
                        alt,
                        ".",
                        ".",
                        ".",
                    ]
                )

        command = [
            str(executable),
            "--offline",
            "--cache",
            "--refseq",
            "--use_given_ref",
            "--species",
            "homo_sapiens",
            "--assembly",
            VEP_ASSEMBLY,
            "--cache_version",
            release,
            "--dir_cache",
            str(cache_dir),
            "--input_file",
            str(input_path),
            "--output_file",
            str(output_path),
            "--format",
            "vcf",
            "--tab",
            "--fields",
            "Uploaded_variation,Gene,Feature,Consequence,IMPACT,CANONICAL,MANE_SELECT,VARIANT_CLASS",
            "--pick_allele_gene",
            "--canonical",
            "--mane",
            "--variant_class",
            "--fork",
            str(forks),
            "--no_stats",
            "--force_overwrite",
        ]
        try:
            process = subprocess.run(command, capture_output=True, text=True, check=False)
        except OSError as exc:
            raise RuntimeError(f"Could not start local VEP executable {executable}: {exc}") from exc
        if process.returncode != 0:
            detail = (process.stderr or process.stdout).strip()[-2_000:]
            raise RuntimeError(f"Local VEP failed with exit code {process.returncode}: {detail}")
        if not output_path.exists():
            raise RuntimeError("Local VEP completed without producing its output table")
        records = _read_local_output(output_path, request_ids)
    return [
        _parse_local_record(row, records.get(request_id))
        for request_id, row in request_ids.items()
    ]


def _validate_vcf_allele(value: object, label: str) -> str:
    allele = str(value).upper()
    if not allele or any(character.isspace() for character in allele):
        raise ValueError(f"Local VEP {label} allele must be non-empty and contain no whitespace")
    return allele


def _read_local_output(
    path: Path,
    request_ids: dict[str, dict[str, object]],
) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        reader = csv.DictReader(
            (line for line in handle if not line.startswith("##")),
            delimiter="\t",
        )
        if reader.fieldnames is None or "#Uploaded_variation" not in reader.fieldnames:
            raise RuntimeError("Local VEP output is missing #Uploaded_variation")
        for record in reader:
            request_id = str(record.get("#Uploaded_variation") or "")
            source = request_ids.get(request_id)
            if source is None or str(record.get("Gene") or "") != str(source["gene_id"]):
                continue
            records.setdefault(request_id, record)
    return records


def _parse_local_record(
    row: dict[str, object],
    record: dict[str, str] | None,
) -> dict[str, object]:
    base = _empty_annotation(row)
    if record is None:
        return {**base, "status": "no_target_gene"}
    terms = sorted(
        {
            term
            for term in str(record.get("Consequence") or "").split(",")
            if term and term != "-"
        },
        key=lambda term: (_CONSEQUENCE_RANK.get(term, len(_CONSEQUENCE_RANK)), term),
    )
    if not terms:
        return {**base, "status": "no_consequence"}
    return {
        **base,
        "status": "ok",
        "primary_consequence": terms[0],
        "consequence_terms": "&".join(terms),
        "transcript_id": _local_value(record, "Feature"),
        "mane_select": _local_value(record, "MANE_SELECT"),
        "canonical": _local_value(record, "CANONICAL").lower() in {"1", "true", "yes"},
        "impact": _local_value(record, "IMPACT"),
        "variant_class": _local_value(record, "VARIANT_CLASS"),
    }


def _local_value(record: dict[str, str], field: str) -> str:
    value = str(record.get(field) or "")
    return "" if value == "-" else value


def _parse_record(row: dict[str, object], record: dict[str, object] | None) -> dict[str, object]:
    base = _empty_annotation(row)
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


def _empty_annotation(row: dict[str, object]) -> dict[str, object]:
    return {
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
