"""Shared request pacing for all gnomAD clients on a shared filesystem."""

from __future__ import annotations

import fcntl
import json
import os
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path


REQUEST_INTERVAL_SECONDS = 7.5  # Eight requests/minute, below the provider's ten.


def state_path() -> Path:
    configured = os.environ.get("GAPH_GNOMAD_RATE_STATE")
    if configured:
        return Path(configured).expanduser()
    cache = os.environ.get("GAPH_GNOMAD_CACHE_DIR")
    if not cache:
        root = os.environ.get("GAPH_ROOT")
        cache = str(Path(root) / "cache/gnomad") if root else str(Path.home() / ".cache/gaph/gnomad")
    return Path(cache) / "request_pacing.json"


@contextmanager
def _locked_state(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with os.fdopen(os.open(path, os.O_CREAT | os.O_RDWR, 0o600), "r+") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        handle.seek(0)
        raw = handle.read()
        state = json.loads(raw) if raw else {"next_request": 0.0}
        yield state
        handle.seek(0)
        json.dump(state, handle)
        handle.truncate()
        handle.flush()
        os.fsync(handle.fileno())


def wait_for_request_slot() -> None:
    with _locked_state(state_path()) as state:
        delay = max(0.0, float(state["next_request"]) - time.time())
        if delay:
            time.sleep(delay)
        state["next_request"] = time.time() + REQUEST_INTERVAL_SECONDS


def defer_requests(retry_after: str | None, *, minimum: float = 60.0) -> None:
    delay = minimum
    if retry_after:
        try:
            delay = max(delay, float(retry_after))
        except ValueError:
            try:
                stamp = parsedate_to_datetime(retry_after)
                delay = max(delay, (stamp - datetime.now(timezone.utc)).total_seconds())
            except (TypeError, ValueError):
                pass
    with _locked_state(state_path()) as state:
        state["next_request"] = max(float(state["next_request"]), time.time() + delay)
