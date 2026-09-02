from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


from bin import check_runtime


def test_explicit_relative_executable_resolves_to_absolute_path(
    monkeypatch,
    tmp_path: Path,
) -> None:
    executable = tmp_path / "vep"
    executable.touch(mode=0o755)
    monkeypatch.chdir(tmp_path)

    assert check_runtime.resolve_executable("./vep") == str(executable)


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
            "--vep-backend",
            "rest",
            "--vep-release",
            "",
            "--vep-executable",
            "vep",
            "--vep-cache-dir",
            "",
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


def test_local_vep_probe_uses_annotation_cache_contract(
    monkeypatch,
    tmp_path: Path,
) -> None:
    executable = tmp_path / "vep"
    executable.touch(mode=0o755)
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    commands: list[list[str]] = []

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(command)
        assert kwargs["timeout"] == check_runtime.LOCAL_VEP_PROBE_TIMEOUT_SECONDS
        return subprocess.CompletedProcess(command, 0, stdout="cache ok", stderr="")

    monkeypatch.setattr(check_runtime.subprocess, "run", fake_run)
    errors: list[str] = []

    result = check_runtime.probe_local_vep(
        release="116",
        executable=str(executable),
        cache_dir=str(cache_dir),
        errors=errors,
    )

    assert errors == []
    assert result["probe"] == "passed"
    assert commands == [
        [
            str(executable),
            "--offline",
            "--cache",
            "--refseq",
            "--use_given_ref",
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh38",
            "--cache_version",
            "116",
            "--dir_cache",
            str(cache_dir),
            "--show_cache_info",
        ]
    ]


def test_local_vep_probe_failure_is_reported_in_runtime_json(
    monkeypatch,
    tmp_path: Path,
) -> None:
    executable = tmp_path / "vep"
    executable.touch(mode=0o755)
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    output = tmp_path / "runtime.json"
    monkeypatch.setattr(
        check_runtime.subprocess,
        "run",
        lambda command, **_kwargs: subprocess.CompletedProcess(
            command,
            2,
            stdout="",
            stderr="cache release unavailable",
        ),
    )
    monkeypatch.setattr(check_runtime, "require_python_module", lambda _name, _errors: None)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "check_runtime.py",
            "--alignment-strategies",
            "",
            "--vep-backend",
            "local",
            "--vep-release",
            "116",
            "--vep-executable",
            str(executable),
            "--vep-cache-dir",
            str(cache_dir),
            "--out-json",
            str(output),
        ],
    )

    try:
        check_runtime.main()
    except SystemExit as exc:
        assert "cache release unavailable" in str(exc)
    else:
        raise AssertionError("Expected a failed runtime dependency check")

    payload = json.loads(output.read_text())
    assert payload["ok"] is False
    assert payload["vep"]["backend"] == "local"
    assert payload["vep"]["probe"] == "failed"
    assert "cache release unavailable" in payload["errors"][0]
