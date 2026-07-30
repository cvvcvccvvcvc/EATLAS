from __future__ import annotations

import json
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path
from typing import Any


class RcloneError(RuntimeError):
    """Raised when an rclone operation fails."""


def validate_remote_root(remote_root: str) -> str:
    value = remote_root.strip().rstrip("/")
    remote_name, separator, _ = value.partition(":")
    if not separator or not remote_name or any(char.isspace() for char in remote_name):
        raise ValueError(
            "Remote root must use rclone syntax, for example 'gdrive:GAPH'."
        )
    return value


def remote_join(remote_root: str, *parts: str) -> str:
    root = validate_remote_root(remote_root)
    clean_parts = [part.strip("/") for part in parts if part.strip("/")]
    if not clean_parts:
        return root
    separator = "" if root.endswith(":") else "/"
    return f"{root}{separator}{'/'.join(clean_parts)}"


def check_config_permissions(config_path: Path | None) -> None:
    if config_path is None:
        return
    path = config_path.expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"rclone config does not exist: {path}")
    mode = stat.S_IMODE(path.stat().st_mode)
    if mode & 0o077:
        raise PermissionError(
            f"rclone config must not be accessible by group or others: "
            f"{path} has mode {mode:o}; expected 600"
        )


class RcloneClient:
    def __init__(
        self,
        *,
        executable: str = "rclone",
        config_path: Path | None = None,
        transfers: int = 4,
        checkers: int = 8,
    ) -> None:
        resolved = shutil.which(executable)
        if resolved is None and executable == "rclone":
            environment_binary = Path(sys.executable).with_name("rclone")
            if environment_binary.is_file() and os.access(environment_binary, os.X_OK):
                resolved = str(environment_binary)
        if resolved is None:
            raise FileNotFoundError(
                f"rclone executable was not found on PATH: {executable}"
            )
        check_config_permissions(config_path)
        self.executable = resolved
        self.config_path = config_path.expanduser() if config_path else None
        self.transfers = transfers
        self.checkers = checkers

    def _command(self, *arguments: str) -> list[str]:
        command = [self.executable]
        if self.config_path is not None:
            command.extend(["--config", str(self.config_path)])
        command.extend(arguments)
        return command

    def _run(
        self,
        *arguments: str,
        capture: bool = False,
        check: bool = True,
    ) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            self._command(*arguments),
            check=False,
            text=True,
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.PIPE if capture else None,
        )
        if check and result.returncode != 0:
            detail = (result.stderr or result.stdout or "").strip()
            raise RcloneError(
                f"rclone command failed with exit code {result.returncode}: "
                f"{detail or 'no diagnostic output'}"
            )
        return result

    def preflight(self, remote_root: str) -> None:
        validated = validate_remote_root(remote_root)
        remote_name = f"{validated.split(':', 1)[0]}:"
        self._run("version", capture=True)
        self._run("about", remote_name, "--json", capture=True)

    def read_text_optional(self, remote_path: str) -> str | None:
        result = self._run("cat", remote_path, capture=True, check=False)
        if result.returncode == 0:
            return result.stdout
        detail = (result.stderr or "").lower()
        missing_markers = (
            "not found",
            "directory not found",
            "object not found",
            "couldn't find",
        )
        if any(marker in detail for marker in missing_markers):
            return None
        raise RcloneError(
            f"Could not read remote object {remote_path}: "
            f"{(result.stderr or result.stdout or '').strip()}"
        )

    def list_files(self, remote_path: str, *, include: str) -> tuple[str, ...]:
        result = self._run(
            "lsjson",
            remote_path,
            "--recursive",
            "--files-only",
            "--include",
            include,
            "--no-modtime",
            "--no-mimetype",
            capture=True,
        )
        try:
            payload = json.loads(result.stdout)
        except json.JSONDecodeError as exc:
            raise RcloneError(
                f"Invalid JSON returned by 'rclone lsjson' for {remote_path}"
            ) from exc
        if not isinstance(payload, list):
            raise RcloneError(
                f"Invalid result returned by 'rclone lsjson' for {remote_path}"
            )
        paths: list[str] = []
        for item in payload:
            if not isinstance(item, dict) or not isinstance(item.get("Path"), str):
                raise RcloneError(
                    f"Invalid file row returned by 'rclone lsjson' for {remote_path}"
                )
            paths.append(item["Path"])
        return tuple(sorted(paths))

    def copy_tree(
        self,
        source: str | Path,
        destination: str | Path,
        *,
        dry_run: bool = False,
    ) -> None:
        arguments = [
            "copy",
            str(source),
            str(destination),
            "--checksum",
            "--immutable",
            "--transfers",
            str(self.transfers),
            "--checkers",
            str(self.checkers),
            "--stats",
            "30s",
            "--stats-one-line-date",
        ]
        if dry_run:
            arguments.append("--dry-run")
        self._run(*arguments)

    def copy_file(self, source: str | Path, destination: str) -> None:
        self._run(
            "copyto",
            str(source),
            destination,
            "--checksum",
            "--immutable",
        )

    def download_file(self, source: str, destination: Path) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        self._run("copyto", source, str(destination), "--checksum")

    def verify_checksum(
        self,
        algorithm: str,
        checksum_file: Path,
        destination: str,
    ) -> None:
        self._run(
            "checksum",
            algorithm,
            str(checksum_file),
            destination,
            "--checkers",
            str(self.checkers),
        )

    def size(self, remote_path: str) -> tuple[int, int]:
        result = self._run("size", remote_path, "--json", capture=True)
        try:
            payload: dict[str, Any] = json.loads(result.stdout)
            return int(payload["count"]), int(payload["bytes"])
        except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
            raise RcloneError(
                f"Invalid JSON returned by 'rclone size' for {remote_path}"
            ) from exc
