"""Content-addressed fragments and flattened, immutable snapshot manifests."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from functools import lru_cache
from pathlib import Path

import duckdb
import pandas as pd

from genomics.gnomad import GNOMAD_API_URL, GNOMAD_DATASET, select_af_metrics, validate_observation_window
from genomics.gnomad_cache import Tile

SCHEMA = "gaph_gnomad_supplement_v1"
COLUMNS = ["tile", "chrom", "pos", "ref", "alt", "gnomad_af", "gnomad_af_source",
           "gnomad_csq", "af_exome", "af_genome", "af_joint", "an_joint", "ac_joint"]


def canonical(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(path: Path, algorithm="sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


@lru_cache(maxsize=4096)
def _verified_digest(path: str, signature: tuple) -> str:
    return digest(Path(path))


def atomic_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(canonical(payload))
    os.replace(temporary, path)


def tile_key(tile: Tile) -> str:
    return f"{tile.chrom}:{tile.start}:{tile.end}"


def parse_tile(value: str) -> Tile:
    chrom, start, end = value.split(":")
    tile = Tile(chrom, int(start), int(end))
    if tile.start < 1 or tile.end != tile.start + 24999 or (tile.start - 1) % 25000:
        raise ValueError(f"Invalid supplement tile: {value}")
    return tile


def empty_state() -> dict:
    return {"schema": SCHEMA, "dataset": GNOMAD_DATASET, "api_url": GNOMAD_API_URL, "reference_genome": "GRCh38",
            "sources": {}, "tiles": {}, "objects": {}, "pending": {}}


def write_fragment(root: Path, observations: dict, state: dict) -> str:
    rows = []
    for key, (records, observed) in observations.items():
        parse_tile(key)
        validate_observation_window(observed)
        for record in records:
            af, source, exome, genome, joint, an, ac = select_af_metrics(record)
            rows.append([key, record["chrom"], int(record["pos"]), record["ref"], record["alt"],
                         af, source, record.get("consequence") or "", exome, genome, joint, an, ac])
    frame = pd.DataFrame(rows, columns=COLUMNS)
    for name in ["tile", "chrom", "ref", "alt", "gnomad_af_source", "gnomad_csq"]:
        frame[name] = frame[name].astype("string")
    for name in ["pos", "an_joint", "ac_joint"]:
        frame[name] = frame[name].astype("Int64")
    for name in ["gnomad_af", "af_exome", "af_genome", "af_joint"]:
        frame[name] = frame[name].astype("Float64")
    directory = root / "objects"
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=directory, suffix=".parquet", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with duckdb.connect() as connection:
            connection.register("fragment", frame)
            connection.sql("SELECT * FROM fragment ORDER BY chrom, pos, ref, alt").write_parquet(
                str(temporary), compression="zstd")
        sha = digest(temporary)
        destination = directory / f"{sha}.parquet"
        if destination.exists():
            if digest(destination) != sha:
                raise ValueError(f"Changed supplement object: {destination}")
        else:
            os.replace(temporary, destination)
        state["objects"][sha] = {"path": f"objects/{sha}.parquet", "sha256": sha,
                                  "size_bytes": destination.stat().st_size,
                                  "md5": digest(destination, "md5")}
        for key, (records, observed) in observations.items():
            state["tiles"][key] = {"object": sha, "observation_window": observed,
                                     "row_count": len(records)}
            state["pending"].pop(key, None)
        return sha
    finally:
        temporary.unlink(missing_ok=True)


def publish(root: Path, state: dict, provenance: dict) -> Path:
    needed = {key for source in state["sources"].values() for key in source["tiles"]}
    active = {key: value for key, value in state["tiles"].items() if key in needed}
    objects = {value["object"] for value in active.values()}
    payload = {**state, "tiles": active,
               "objects": {key: state["objects"][key] for key in sorted(objects)},
               "pending": {key: state["pending"].get(key, "not fetched") for key in sorted(needed - active.keys())},
               "provenance": provenance, "complete": needed <= active.keys()}
    sha = hashlib.sha256(canonical(payload)).hexdigest()
    path = root / "snapshots" / f"{sha}.json"
    if path.exists() and path.read_bytes() != canonical(payload):
        raise ValueError(f"Changed immutable snapshot: {path}")
    if not path.exists():
        atomic_json(path, payload)
    return path


def load_snapshot(path: Path, *, verify_objects=True) -> dict:
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != path.stem:
        raise ValueError(f"Snapshot digest mismatch: {path}")
    value = json.loads(raw)
    if value.get("schema") != SCHEMA or value.get("dataset") != GNOMAD_DATASET or value.get("reference_genome") != "GRCh38":
        raise ValueError(f"Unsupported supplement contract: {path}")
    root = path.parent.parent
    for key, item in value["tiles"].items():
        parse_tile(key)
        validate_observation_window(item["observation_window"])
        if item["object"] not in value["objects"]:
            raise ValueError(f"Missing object for tile {key}")
    for sha, item in value["objects"].items():
        if item["path"] != f"objects/{sha}.parquet" or item["sha256"] != sha:
            raise ValueError("Invalid supplement object identity")
        if verify_objects:
            member = root / item["path"]
            stat = member.stat()
            signature = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
            if stat.st_size != item["size_bytes"] or _verified_digest(str(member.resolve()), signature) != sha:
                raise ValueError(f"Supplement object changed: {member}")
    needed = {key for source in value["sources"].values() for key in source["tiles"]}
    if set(value["pending"]) != needed - value["tiles"].keys() or not value["tiles"].keys() <= needed:
        raise ValueError("Supplement pending set does not match source coverage")
    if bool(value["complete"]) != (needed <= value["tiles"].keys()):
        raise ValueError("Supplement completeness claim does not match tile coverage")
    return value
