from __future__ import annotations

import sys
from pathlib import Path


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

import check_runtime  # noqa: E402


def test_map_ont_strategy_checks_minimap2_dependency(
    monkeypatch,
    tmp_path: Path,
) -> None:
    requested: list[str] = []
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "check_runtime.py",
            "--stage",
            "align",
            "--alignment-strategies",
            "minimap2_map_ont_pseudoreads_30000_15000",
            "--out-json",
            str(tmp_path / "runtime.json"),
        ],
    )
    monkeypatch.setattr(
        check_runtime,
        "require_executable",
        lambda name, _raw, _errors: requested.append(name),
    )

    check_runtime.main()

    assert requested == ["bedtools", "minimap2"]
