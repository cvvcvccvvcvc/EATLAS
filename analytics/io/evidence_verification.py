"""Reuse full evidence verification only while the source tree is unchanged."""

from __future__ import annotations

import json
from pathlib import Path

from analytics.io.artifacts import write_json_atomic
from provenance.evidence_inventory import evidence_stat_fingerprint, verify_run_evidence


def verify_source_evidence(
    *,
    run_dir: Path,
    source_id: str,
    inventory: dict[str, object],
    inventory_descriptor: dict[str, object],
    cache_dir: Path,
) -> bool:
    """Return True for a previously verified, unchanged source. Caller owns root lock."""

    marker = cache_dir / "evidence_inventory.verified.json"
    expected = {
        "schema_version": 1,
        "source_id": source_id,
        "run_dir": str(run_dir.resolve()),
        "inventory": inventory_descriptor,
        "stat_fingerprint": evidence_stat_fingerprint(run_dir, inventory),
    }
    if marker.is_file():
        try:
            if json.loads(marker.read_text()) == expected:
                return True
        except (OSError, json.JSONDecodeError):
            pass  # An unreadable marker cannot attest to a previous verification.
    expected["stat_fingerprint"] = verify_run_evidence(run_dir, inventory)
    write_json_atomic(marker, expected)
    return False
