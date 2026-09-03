from __future__ import annotations

from pathlib import Path

import pytest

from analytics.io.root_lock import LOCK_FILENAME, analytics_root_lock


def test_analytics_root_lock_is_nonblocking_and_released(tmp_path: Path) -> None:
    root = tmp_path / "analytics"

    with analytics_root_lock(root) as lock_path:
        assert lock_path == root / LOCK_FILENAME
        assert lock_path.is_file()
        with pytest.raises(RuntimeError, match="already writing analytics root"):
            with analytics_root_lock(root):
                pass

    with analytics_root_lock(root):
        pass


def test_analytics_root_lock_is_released_after_failure(tmp_path: Path) -> None:
    root = tmp_path / "analytics"

    with pytest.raises(ValueError, match="report failed"):
        with analytics_root_lock(root):
            raise ValueError("report failed")

    with analytics_root_lock(root):
        pass
