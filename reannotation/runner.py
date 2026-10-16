"""Sequential archive planning, resumable tile acquisition, and publication."""

from __future__ import annotations

import copy
import hashlib
import json
import shutil
import time
from http.client import IncompleteRead, RemoteDisconnected
from urllib.error import HTTPError, URLError
from pathlib import Path

from genomics.gnomad_cache import GnomadRegionCache, tiles_for_region
from .archive import source_from_archive, upload_snapshot
from .store import atomic_json, canonical, empty_state, load_snapshot, parse_tile, publish, tile_key, write_fragment


def require_capacity(root: Path, budget: int) -> None:
    used = sum(p.stat().st_size for p in root.rglob("*") if p.is_file())
    if used >= budget or shutil.disk_usage(root).free < 1024**3:
        raise RuntimeError(f"Reannotation space limit reached: {used} bytes used; budget={budget}. Progress retained.")


def run_reannotation(*, root: Path, archive_root: str, output_remote: str | None,
                     run_ids: list[str], client, cache: GnomadRegionCache,
                     provenance: dict, parent: Path | None = None,
                     refresh: list[str] = (), budget: int = 20 * 1024**3,
                     max_fetch_tiles: int | None = None, passes: int = 0) -> Path:
    if not run_ids or len(set(run_ids)) != len(run_ids):
        raise ValueError("Run IDs must be nonempty and unique")
    run_ids = sorted(run_ids)
    root.mkdir(parents=True, exist_ok=True)
    (root / "work").mkdir(exist_ok=True)
    parent_id = parent.stem if parent else None
    plan = {"archive_root": archive_root, "run_ids": run_ids, "parent": parent_id, "refresh": list(refresh)}
    plan_id = hashlib.sha256(canonical(plan)).hexdigest()
    progress = root / "work" / f"{plan_id}.json"
    if progress.exists():
        state = json.loads(progress.read_text())
        if state["plan"] != plan:
            raise ValueError("Reannotation progress plan changed")
    else:
        state = copy.deepcopy(load_snapshot(parent)) if parent else empty_state()
        if not state["sources"].keys() <= set(run_ids):
            raise ValueError("Run IDs must include all source runs bound by the parent snapshot")
        if parent and parent.parent.parent.resolve() != root.resolve():
            raise ValueError("Parent snapshot must belong to the same object store")
        state["plan"] = plan
        state.pop("provenance", None)
        for spec in refresh:
            chrom, interval = spec.split(":")
            start, end = map(int, interval.split("-"))
            for tile in tiles_for_region(chrom, start, end):
                state["tiles"].pop(tile_key(tile), None)
        atomic_json(progress, state)
    state["producer_commits"] = sorted(set(state.get("producer_commits", [])) | {provenance["git_commit"]})
    # Only the four small, checksum-bound metadata files are downloaded per run.
    for run_id in run_ids:
        require_capacity(root, budget)
        source = source_from_archive(client, archive_root, run_id, root / "work" / run_id, budget)
        previous = state["sources"].get(run_id)
        if previous and previous != source:
            raise ValueError(f"Archive identity changed for {run_id}")
        state["sources"][run_id] = source
        atomic_json(progress, state)
        print(f"Planned {run_id}: failed_regions={len(source['regions'])}; tiles={len(source['tiles'])}", flush=True)
    force = {tile_key(tile) for spec in refresh
             for tile in tiles_for_region(spec.split(":")[0], *map(int, spec.split(":")[1].split("-")))}
    # Refreshed observations use a separate resume cache; old observations survive.
    refresh_cache = GnomadRegionCache(root / "work" / f"refresh-{plan_id}") if force else cache
    attempts = 0
    wave = 0
    while True:
        observations = {}
        for run_id in run_ids:
            for key in state["sources"][run_id]["tiles"]:
                if key in state["tiles"] or key in observations:
                    continue
                tile = parse_tile(key)
                selected_cache = refresh_cache if key in force else cache
                cached = selected_cache.read_cached_tile_with_observation(tile)
                if cached is None:
                    if max_fetch_tiles is not None and attempts >= max_fetch_tiles:
                        continue
                    attempts += 1
                    try:
                        selected_cache.fetch_tiles([tile])
                        cached = selected_cache.read_cached_tile_with_observation(tile)
                    except (HTTPError, URLError, TimeoutError, IncompleteRead, RemoteDisconnected, RuntimeError) as exc:
                        state["pending"][key] = f"{type(exc).__name__}: {exc}"[:4500]
                        atomic_json(progress, state)
                        print(f"Pending {run_id} {key}: {type(exc).__name__}: {exc}", flush=True)
                        continue
                if cached is None:
                    raise ValueError(f"Fetcher did not publish a valid tile: {key}")
                observations[key] = cached
                if len(observations) >= 32:
                    require_capacity(root, budget)
                    write_fragment(root, observations, state)
                    observations.clear()
                    atomic_json(progress, state)
            if observations:
                require_capacity(root, budget)
                write_fragment(root, observations, state)
                observations.clear()
                atomic_json(progress, state)
            snapshot = publish(root, state, provenance)
            if output_remote:
                upload_snapshot(client, snapshot, output_remote)
            remaining = len(set(state["sources"][run_id]["tiles"]) - state["tiles"].keys())
            print(f"{run_id}: remaining_tiles={remaining}; snapshot={snapshot.stem}", flush=True)
        value = load_snapshot(snapshot)
        wave += 1
        if value["complete"] or (passes and wave >= passes) or (max_fetch_tiles is not None and attempts >= max_fetch_tiles):
            return snapshot
        print(f"Reannotation waiting: {len(value['pending'])} tiles remain; successful tiles retained", flush=True)
        time.sleep(min(900, 60 * 2 ** min(wave - 1, 4)))
