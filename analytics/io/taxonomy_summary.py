"""Resolve and cache the run-level taxonomy summary for analytics."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import tempfile
from itertools import chain
from pathlib import Path
from typing import Sequence

from analytics.io.artifacts import content_identity, write_json_atomic
from analytics.derivations.taxonomy import (
    TAXONOMY_SUMMARY_FIELDS,
    build_taxonomy_summary_rows,
    load_taxonomy_profiles,
    write_taxonomy_summary,
)


CACHE_SCHEMA_VERSION = 3
CACHE_DIRNAME = "taxonomy_summary"
CACHE_FILENAME = "taxonomy_summary.tsv.gz"


def resolve_taxonomy_summary_path(
    run_dir: Path,
    *,
    analytics_dir: Path,
) -> Path:
    """Require the canonical Stage 1 taxonomy contract and derive its summary."""

    fetch_dir = run_dir / "fetch"
    taxonomy_tsv = fetch_dir / "taxonomy.tsv.gz"
    orthologs_tsv = fetch_dir / "orthologs.selected.tsv.gz"
    missing = [path for path in (taxonomy_tsv, orthologs_tsv) if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Incomplete Stage 1 taxonomy contract required for analytics; missing: "
            + ", ".join(str(path) for path in missing)
        )
    return build_or_load_taxonomy_summary(
        taxonomy_tsv=taxonomy_tsv,
        orthologs_tsv=orthologs_tsv,
        analytics_dir=analytics_dir,
    )


def build_or_load_taxonomy_summary(
    *,
    taxonomy_tsv: Path,
    orthologs_tsv: Path,
    analytics_dir: Path,
) -> Path:
    """Materialize one run's taxonomy summary under analytics."""

    return build_or_load_taxonomy_summary_many(
        taxonomy_tsvs=(taxonomy_tsv,),
        orthologs_tsvs=(orthologs_tsv,),
        analytics_dir=analytics_dir,
    )


def build_or_load_taxonomy_summary_many(
    *,
    taxonomy_tsvs: Sequence[Path],
    orthologs_tsvs: Sequence[Path],
    analytics_dir: Path,
) -> Path:
    """Derive one exact summary from compatible, disjoint source runs."""

    taxonomy_tsvs = tuple(taxonomy_tsvs)
    orthologs_tsvs = tuple(orthologs_tsvs)
    if not taxonomy_tsvs or len(taxonomy_tsvs) != len(orthologs_tsvs):
        raise ValueError(
            "Taxonomy summary requires equal non-empty taxonomy and ortholog inputs"
        )

    cache_dir = analytics_dir / CACHE_DIRNAME
    output_path = cache_dir / CACHE_FILENAME
    manifest_path = cache_dir / "manifest.json"
    inputs = {
        "taxonomy": [content_identity(path) for path in taxonomy_tsvs],
        "orthologs_selected": [content_identity(path) for path in orthologs_tsvs],
    }
    fingerprint = _fingerprint(inputs)
    if _cache_is_valid(manifest_path, output_path, inputs, fingerprint):
        return output_path

    profiles = {}
    for path in taxonomy_tsvs:
        for tax_id, profile in load_taxonomy_profiles(path).items():
            previous = profiles.setdefault(tax_id, profile)
            if previous != profile:
                raise ValueError(
                    f"Source taxonomy tables disagree for tax_id {tax_id}: {path}"
                )

    handles = [gzip.open(path, "rt", newline="") for path in orthologs_tsvs]
    try:
        rows = build_taxonomy_summary_rows(
            chain.from_iterable(
                csv.DictReader(handle, delimiter="\t") for handle in handles
            ),
            profiles,
        )
    finally:
        for handle in handles:
            handle.close()

    cache_dir.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=cache_dir,
            prefix=f".{CACHE_FILENAME}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
        write_taxonomy_summary(temporary_path, rows)
        temporary_path.chmod(0o644)
        temporary_path.replace(output_path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    write_json_atomic(
        manifest_path,
        {
            "schema_version": CACHE_SCHEMA_VERSION,
            "status": "complete",
            "fingerprint": fingerprint,
            "inputs": inputs,
            "columns": TAXONOMY_SUMMARY_FIELDS,
            "row_count": len(rows),
            "output": content_identity(output_path),
        },
    )
    return output_path


def _fingerprint(inputs: dict[str, object]) -> str:
    payload = {
        "schema_version": CACHE_SCHEMA_VERSION,
        "inputs": inputs,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _cache_is_valid(
    manifest_path: Path,
    output_path: Path,
    inputs: dict[str, object],
    fingerprint: str,
) -> bool:
    if not manifest_path.is_file() or not output_path.is_file():
        return False
    try:
        manifest = json.loads(manifest_path.read_text())
        return (
            manifest.get("schema_version") == CACHE_SCHEMA_VERSION
            and manifest.get("status") == "complete"
            and manifest.get("fingerprint") == fingerprint
            and manifest.get("inputs") == inputs
            and manifest.get("columns") == TAXONOMY_SUMMARY_FIELDS
            and manifest.get("output") == content_identity(output_path)
        )
    except (OSError, ValueError, json.JSONDecodeError):
        return False
