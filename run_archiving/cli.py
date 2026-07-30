from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Sequence

from run_archiving.archive import (
    ArchiveError,
    archive_run,
    list_archives,
    remove_local_run,
    restore_run,
    verify_remote,
)
from run_archiving.rclone import RcloneClient, RcloneError, validate_remote_root


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return parsed


def _remote_default() -> str | None:
    return os.environ.get("GAPH_ARCHIVE_REMOTE")


def _config_default() -> Path | None:
    configured = os.environ.get("RCLONE_CONFIG")
    if configured:
        return Path(configured)
    standard = Path.home() / ".config" / "rclone" / "rclone.conf"
    return standard if standard.exists() else None


def _format_bytes(value: int) -> str:
    size = float(value)
    units = ("B", "KiB", "MiB", "GiB", "TiB", "PiB")
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{int(size)} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    raise AssertionError("unreachable")


def _add_remote_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--remote",
        default=_remote_default(),
        help=(
            "rclone archive root, for example gdrive:GAPH "
            "(default: GAPH_ARCHIVE_REMOTE)"
        ),
    )
    parser.add_argument(
        "--rclone-config",
        type=Path,
        default=_config_default(),
        help="rclone config path (default: RCLONE_CONFIG or standard user config)",
    )
    parser.add_argument("--rclone-bin", default="rclone")
    parser.add_argument("--transfers", type=_positive_int, default=4)
    parser.add_argument("--checkers", type=_positive_int, default=8)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m run_archiving",
        description="Archive complete pipeline run directories with rclone.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    archive_parser = subparsers.add_parser(
        "archive", help="upload, verify, and mark one run archive complete"
    )
    archive_parser.add_argument("--run-dir", type=Path, required=True)
    archive_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="hash the run and execute rclone copy with --dry-run",
    )
    archive_parser.add_argument(
        "--allow-legacy-run",
        action="store_true",
        help=(
            "allow a historical run without root run_manifest.json; "
            "failed or running manifests are never accepted"
        ),
    )
    _add_remote_arguments(archive_parser)

    list_parser = subparsers.add_parser(
        "list", help="list archives marked complete without rechecking their data"
    )
    _add_remote_arguments(list_parser)

    verify_parser = subparsers.add_parser(
        "verify", help="verify an existing remote archive"
    )
    verify_parser.add_argument("--run-id", required=True)
    _add_remote_arguments(verify_parser)

    restore_parser = subparsers.add_parser(
        "restore", help="restore and verify one archived run"
    )
    restore_parser.add_argument("--run-id", required=True)
    restore_parser.add_argument("--destination", type=Path, required=True)
    _add_remote_arguments(restore_parser)

    remove_parser = subparsers.add_parser(
        "remove-local",
        help="remove a local run only after a fresh remote and local verification",
    )
    remove_parser.add_argument("--run-dir", type=Path, required=True)
    remove_parser.add_argument(
        "--results-root",
        type=Path,
        default=(
            Path(os.environ["GAPH_ROOT"]) / "results"
            if os.environ.get("GAPH_ROOT")
            else None
        ),
        help="required parent of run-dir (default: $GAPH_ROOT/results)",
    )
    remove_parser.add_argument(
        "--confirm-run-id",
        required=True,
        help="must exactly equal the run directory name",
    )
    _add_remote_arguments(remove_parser)
    return parser


def _client(args: argparse.Namespace) -> RcloneClient:
    if not args.remote:
        raise ArchiveError(
            "Archive remote is required; pass --remote or set GAPH_ARCHIVE_REMOTE."
        )
    args.remote = validate_remote_root(args.remote)
    return RcloneClient(
        executable=args.rclone_bin,
        config_path=args.rclone_config,
        transfers=args.transfers,
        checkers=args.checkers,
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        client = _client(args)
        if args.command == "archive":
            result = archive_run(
                client,
                run_dir=args.run_dir,
                remote_root=args.remote,
                dry_run=args.dry_run,
                allow_legacy_run=args.allow_legacy_run,
            )
        elif args.command == "list":
            archives = list_archives(client, remote_root=args.remote)
            total_bytes = sum(int(item["total_bytes"]) for item in archives)
            result = {
                "status": "listed",
                "remote": args.remote,
                "archive_count": len(archives),
                "total_bytes": total_bytes,
                "total_size": _format_bytes(total_bytes),
                "archives": [
                    {
                        **item,
                        "size": _format_bytes(int(item["total_bytes"])),
                    }
                    for item in archives
                ],
            }
        elif args.command == "verify":
            manifest = verify_remote(
                client, remote_root=args.remote, run_id=args.run_id
            )
            result = {
                "status": "verified",
                "run_id": args.run_id,
                "file_count": manifest["file_count"],
                "total_bytes": manifest["total_bytes"],
                "tree_sha256": manifest["tree_sha256"],
            }
        elif args.command == "restore":
            result = restore_run(
                client,
                remote_root=args.remote,
                run_id=args.run_id,
                destination=args.destination,
            )
        elif args.command == "remove-local":
            if args.results_root is None:
                raise ArchiveError(
                    "--results-root is required when GAPH_ROOT is not set."
                )
            result = remove_local_run(
                client,
                run_dir=args.run_dir,
                results_root=args.results_root,
                remote_root=args.remote,
                confirmation=args.confirm_run_id,
            )
        else:
            parser.error(f"Unsupported command: {args.command}")
            return 2
    except (ArchiveError, RcloneError, FileNotFoundError, PermissionError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0
