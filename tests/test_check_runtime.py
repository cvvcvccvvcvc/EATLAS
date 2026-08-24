from __future__ import annotations

import sys
from pathlib import Path


from bin import check_runtime


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
    monkeypatch.setattr(
        check_runtime,
        "require_python_module",
        lambda _name, _errors: None,
    )

    check_runtime.main()

    assert requested == ["minimap2"]
