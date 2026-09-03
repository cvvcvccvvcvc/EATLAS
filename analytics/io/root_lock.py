"""Process-level ownership of an analytics output root."""

from __future__ import annotations

import fcntl
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path


LOCK_FILENAME = ".strategy_report.lock"


@contextmanager
def analytics_root_lock(analytics_root: Path) -> Iterator[Path]:
    """Hold exclusive nonblocking ownership while a report may write the root."""

    root = analytics_root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    lock_path = root / LOCK_FILENAME
    with lock_path.open("a+b") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                "Another strategy report is already writing analytics root: "
                f"{root}"
            ) from exc
        try:
            yield lock_path
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
