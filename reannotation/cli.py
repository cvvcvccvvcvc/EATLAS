from __future__ import annotations

import argparse
import fcntl
import json
import os
import subprocess
from pathlib import Path

from genomics.gnomad_cache import GnomadRegionCache
from run_archiving.rclone import RcloneClient
from .archive import download_snapshot
from .runner import run_reannotation


def main():
    parser = argparse.ArgumentParser(description="Incremental, immutable gnomAD supplements; source runs are read-only.")
    parser.add_argument("--rclone", default="rclone")
    parser.add_argument("--rclone-config", type=Path, default=Path(os.environ.get("RCLONE_CONFIG", Path.home() / ".config/rclone/rclone.conf")))
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run")
    run.add_argument("--root", type=Path, required=True)
    run.add_argument("--archive-root", required=True)
    run.add_argument("--output-remote", required=True)
    run.add_argument("--run-id", action="append", required=True)
    run.add_argument("--gnomad-cache-dir", type=Path, required=True)
    run.add_argument("--parent", type=Path)
    run.add_argument("--refresh-region", action="append", default=[])
    run.add_argument("--disk-budget-gb", type=int, default=20)
    run.add_argument("--max-fetch-tiles", type=int)
    run.add_argument("--passes", type=int, default=0, help="0 continues until complete; positive values bound retry waves.")
    fetch = commands.add_parser("download")
    fetch.add_argument("--remote-root", required=True)
    fetch.add_argument("--snapshot", required=True)
    fetch.add_argument("--destination", required=True, type=Path)
    args = parser.parse_args()
    client = RcloneClient(executable=args.rclone, config_path=args.rclone_config, transfers=2, checkers=2)
    if args.command == "download":
        client.require_google_drive(args.remote_root)
        print(download_snapshot(client, args.remote_root, args.snapshot, args.destination))
        return
    if args.disk_budget_gb < 1 or args.passes < 0 or (args.max_fetch_tiles is not None and args.max_fetch_tiles < 0):
        parser.error("Budgets must be positive; passes and max-fetch-tiles must be nonnegative")
    if len(set(args.run_id)) != len(args.run_id):
        parser.error("Duplicate run IDs")
    code_root = Path(__file__).resolve().parents[1]
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=code_root, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=code_root, text=True).strip()
    if dirty:
        raise ValueError("Production reannotation requires a clean Git checkout")
    for remote in [args.archive_root, args.output_remote]:
        client.require_google_drive(remote)
    args.root = args.root.expanduser().resolve()
    args.root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("GAPH_GNOMAD_CACHE_DIR", str(args.gnomad_cache_dir.resolve()))
    with (args.root / ".writer.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = run_reannotation(root=args.root, archive_root=args.archive_root, output_remote=args.output_remote,
                                  run_ids=args.run_id, client=client, cache=GnomadRegionCache(args.gnomad_cache_dir),
                                  provenance={"git_commit": commit, "git_dirty": False}, parent=args.parent,
                                  refresh=args.refresh_region, budget=args.disk_budget_gb * 1024**3,
                                  max_fetch_tiles=args.max_fetch_tiles, passes=args.passes)
    print(json.dumps({"snapshot": str(result)}))
