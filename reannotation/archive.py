"""Verified small archive inputs and append-only supplement transfers."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
from pathlib import Path

from genomics.gnomad_cache import tiles_for_region
from provenance.evidence_inventory import load_bound_evidence_inventory
from run_archiving.archive import read_archive_manifest
from run_archiving.rclone import remote_join
from .store import atomic_json, digest, load_snapshot, tile_key

INPUT_FILES = ("run_manifest.json", "evidence_inventory.json", "annotation/manifest.json", "annotation/failures.tsv.gz")


def source_from_archive(client, remote_root: str, run_id: str, work: Path, budget: int) -> dict:
    archive = read_archive_manifest(client, remote_root=remote_root, run_id=run_id)
    records = {x["path"]: x for x in archive["files"]}
    needed = sum(records[name]["size"] for name in INPUT_FILES)
    if needed > budget or needed + 1024**3 > shutil.disk_usage(work.parent).free:
        raise RuntimeError(f"Insufficient space for {run_id}: needs {needed} metadata bytes plus 1 GiB reserve")
    work.mkdir(parents=True, exist_ok=True)
    try:
        for name in INPUT_FILES:
            destination = work / name
            client.download_file(remote_join(remote_root, "runs", run_id, "data", name), destination)
            if destination.stat().st_size != records[name]["size"] or digest(destination) != records[name]["sha256"]:
                raise ValueError(f"Archive input checksum mismatch: {run_id}/{name}")
        root = json.loads((work / "run_manifest.json").read_text())
        if root.get("status") != "complete" or root.get("success") is not True:
            raise ValueError(f"Source run is not complete: {run_id}")
        inventory, _ = load_bound_evidence_inventory(work / "evidence_inventory.json", root.get("evidence_inventory"))
        evidence = {item["path"]: item for item in inventory["files"]}
        for name in INPUT_FILES[2:]:
            if evidence[name]["sha256"] != digest(work / name):
                raise ValueError(f"Source evidence checksum mismatch: {run_id}/{name}")
        annotation = json.loads((work / "annotation/manifest.json").read_text())
        if annotation.get("gnomad_dataset") != "gnomad_r4":
            raise ValueError(f"Unexpected gnomAD source contract: {run_id}")
        regions = []
        with gzip.open(work / "annotation/failures.tsv.gz", "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["source"] == "gnomad" and row["scope"] == "region":
                    regions.append([row["chrom"], int(row["start"]), int(row["end"])])
        tiles = sorted({tile_key(t) for chrom, start, end in regions
                        for t in tiles_for_region(chrom, max(1, start - 100), end + 100)})
        return {"run_id": run_id, "archive_tree_sha256": archive["tree_sha256"],
                "root_manifest_sha256": digest(work / "run_manifest.json"),
                "evidence_tree_sha256": inventory["tree_sha256"], "regions": regions, "tiles": tiles}
    finally:
        shutil.rmtree(work)


def upload_snapshot(client, path: Path, remote_root: str) -> None:
    value = load_snapshot(path)
    root = path.parent.parent
    checksums = root / "work" / "upload.md5"
    checksums.parent.mkdir(parents=True, exist_ok=True)
    receipt = root / "work" / (hashlib.sha256(remote_root.encode()).hexdigest() + ".uploaded.json")
    uploaded = set(json.loads(receipt.read_text())) if receipt.exists() else set()
    fresh = [v for sha, v in value["objects"].items() if sha not in uploaded]
    for item in fresh:
        client.copy_file(root / item["path"], remote_join(remote_root, item["path"]))
    checksums.write_text("".join(f"{v['md5']}  {v['path']}\n" for v in fresh))
    if fresh:
        client.verify_checksum("MD5", checksums, remote_root)
        uploaded.update(v["sha256"] for v in fresh)
        atomic_json(receipt, sorted(uploaded))
    client.copy_file(path, remote_join(remote_root, "snapshots", path.name))
    checksums.write_text(f"{digest(path, 'md5')}  snapshots/{path.name}\n")
    client.verify_checksum("MD5", checksums, remote_root)
    checksums.unlink()


def download_snapshot(client, remote_root: str, snapshot_id: str, destination: Path) -> Path:
    if len(snapshot_id) != 64 or any(c not in "0123456789abcdef" for c in snapshot_id):
        raise ValueError("Snapshot ID must be a SHA-256 digest")
    path = destination / "snapshots" / f"{snapshot_id}.json"
    client.download_file(remote_join(remote_root, "snapshots", path.name), path)
    value = load_snapshot(path, verify_objects=False)
    required = sum(item["size_bytes"] for item in value["objects"].values()
                   if not (destination / item["path"]).is_file())
    if required + 1024**3 > shutil.disk_usage(destination).free:
        raise RuntimeError(f"Insufficient free space for supplement: needs {required} bytes plus 1 GiB reserve")
    for item in value["objects"].values():
        member = destination / item["path"]
        if not member.is_file() or digest(member) != item["sha256"]:
            client.download_file(remote_join(remote_root, item["path"]), member)
    load_snapshot(path)
    return path
