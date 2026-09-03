from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from analytics.io.evidence_verification import verify_source_evidence
from provenance import evidence_inventory as inventory_module

from provenance.evidence_inventory import (
    EvidenceInventoryError,
    assert_inventory_matches_records,
    build_evidence_inventory,
    build_run_evidence_inventory,
    combine_evidence_inventories,
    evidence_stat_fingerprint,
    inventory_file_descriptor,
    load_bound_evidence_inventory,
    load_evidence_inventory,
    validate_evidence_inventory,
    verify_run_evidence,
    write_evidence_inventory,
)


def _run(tmp_path: Path) -> Path:
    run = tmp_path / "run"
    for scope in ("fetch", "alignment", "annotation"):
        directory = run / scope
        directory.mkdir(parents=True)
        (directory / "manifest.json").write_text(f"{scope}\n")
    (run / "fetch" / "empty.tsv").write_bytes(b"")
    return run


def test_inventory_is_deterministic_and_verifies_every_evidence_file(
    tmp_path: Path,
) -> None:
    run = _run(tmp_path)
    path = run / "evidence_inventory.json"

    first = write_evidence_inventory(
        path,
        {
            "annotation": run / "annotation",
            "fetch": run / "fetch",
            "alignment": run / "alignment",
        },
    )
    second = build_run_evidence_inventory(run)

    assert first == second == load_evidence_inventory(path)
    assert [row["path"] for row in first["files"]] == [
        "alignment/manifest.json",
        "annotation/manifest.json",
        "fetch/empty.tsv",
        "fetch/manifest.json",
    ]
    assert first["file_count"] == 4
    assert verify_run_evidence(run, first) == evidence_stat_fingerprint(run, first)
    assert_inventory_matches_records(
        first,
        {
            str(row["path"]): (int(row["size_bytes"]), str(row["sha256"]))
            for row in first["files"]
        },
    )


@pytest.mark.parametrize("mutation", ["changed", "added", "removed"])
def test_inventory_rejects_evidence_tree_changes(
    tmp_path: Path,
    mutation: str,
) -> None:
    run = _run(tmp_path)
    inventory = build_run_evidence_inventory(run)
    target = run / "fetch" / "manifest.json"
    if mutation == "changed":
        target.write_text("changed\n")
    elif mutation == "added":
        (run / "fetch" / "extra.txt").write_text("extra\n")
    else:
        target.unlink()

    with pytest.raises(EvidenceInventoryError, match="differs from its inventory"):
        verify_run_evidence(run, inventory)


def test_inventory_rejects_unsafe_paths_and_symlinks(tmp_path: Path) -> None:
    run = _run(tmp_path)
    inventory = build_run_evidence_inventory(run)
    inventory["files"][0]["path"] = "../outside"
    with pytest.raises(EvidenceInventoryError, match="Unsafe evidence path"):
        validate_evidence_inventory(inventory)

    (run / "fetch" / "linked").symlink_to(run / "annotation")
    with pytest.raises(EvidenceInventoryError, match="symlink"):
        build_run_evidence_inventory(run)


def test_combines_exactly_one_strict_fragment_per_stage(tmp_path: Path) -> None:
    fragments = {}
    for scope in ("fetch", "alignment", "annotation"):
        source = tmp_path / f"{scope}.txt"
        source.write_text(f"{scope}\n")
        fragments[scope] = build_evidence_inventory(
            {f"{scope}/manifest.json": source}, scopes=(scope,)
        )

    combined = combine_evidence_inventories(fragments)

    assert combined["scope"] == ["fetch", "alignment", "annotation"]
    assert [row["path"] for row in combined["files"]] == [
        "alignment/manifest.json",
        "annotation/manifest.json",
        "fetch/manifest.json",
    ]
    with pytest.raises(EvidenceInventoryError, match="must cover"):
        combine_evidence_inventories({"fetch": fragments["fetch"]})


def test_load_rejects_summary_tampering(tmp_path: Path) -> None:
    path = tmp_path / "inventory.json"
    payload = build_run_evidence_inventory(_run(tmp_path))
    payload["total_bytes"] = int(payload["total_bytes"]) + 1
    path.write_text(json.dumps(payload))

    with pytest.raises(EvidenceInventoryError, match="total_bytes"):
        load_evidence_inventory(path)


def test_bound_load_parses_only_inventory_bytes_authenticated_by_manifest(
    tmp_path: Path,
) -> None:
    run = _run(tmp_path)
    path = tmp_path / "evidence_inventory.json"
    expected = write_evidence_inventory(
        path,
        {scope: run / scope for scope in inventory_module.EVIDENCE_SCOPES},
    )
    descriptor = inventory_file_descriptor(path)

    observed, observed_descriptor = load_bound_evidence_inventory(path, descriptor)

    assert observed == expected
    assert observed_descriptor == descriptor
    path.write_text("{}\n")
    with pytest.raises(EvidenceInventoryError, match="does not match"):
        load_bound_evidence_inventory(path, descriptor)


def test_inventory_walk_fails_when_a_directory_cannot_be_read(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    def failed_walk(root, *, topdown, onerror, followlinks):
        assert topdown is True
        assert followlinks is False
        onerror(PermissionError(13, "permission denied", str(root)))

    monkeypatch.setattr(inventory_module.os, "walk", failed_walk)

    with pytest.raises(EvidenceInventoryError, match="Cannot read evidence directory"):
        build_run_evidence_inventory(_run(tmp_path))


def test_verified_marker_skips_hashes_and_detects_same_size_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    run = _run(tmp_path)
    path = run / "evidence_inventory.json"
    inventory = write_evidence_inventory(
        path, {scope: run / scope for scope in inventory_module.EVIDENCE_SCOPES}
    )
    kwargs = dict(
        run_dir=run, source_id="test-source", inventory=inventory,
        inventory_descriptor=inventory_file_descriptor(path), cache_dir=tmp_path / "cache",
    )
    original = inventory_module._sha256_file
    calls = []

    def tracked_hash(path):
        calls.append(path)
        return original(path)

    monkeypatch.setattr(inventory_module, "_sha256_file", tracked_hash)
    assert verify_source_evidence(**kwargs) is False
    assert len(calls) == inventory["file_count"]
    calls.clear()
    assert verify_source_evidence(**kwargs) is True
    assert not calls

    target = run / "fetch" / "manifest.json"
    previous = target.stat()
    target.write_bytes(b"FETCH\n")
    os.utime(target, ns=(previous.st_atime_ns, previous.st_mtime_ns))
    with pytest.raises(EvidenceInventoryError, match="differs from its inventory"):
        verify_source_evidence(**kwargs)
    assert calls
