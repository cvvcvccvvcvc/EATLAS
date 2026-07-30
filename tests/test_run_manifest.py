from __future__ import annotations

import json
import os
import shutil
import stat
import subprocess
import time
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "run_manifest_workflow.nf"


def _command(tmp_path: Path, outdir: Path, *, fail: bool, hold_seconds: int) -> list[str]:
    nextflow = shutil.which("nextflow")
    if nextflow is None:
        pytest.skip("nextflow is not installed")
    empty_config = tmp_path / "empty.config"
    empty_config.write_text("")
    schema_path = tmp_path / "parameters.schema.json"
    schema_path.write_text(
        json.dumps(
            {
                "definitions": {
                    "test": {
                        "properties": {
                            name: {}
                            for name in [
                                "api_token",
                                "endpoint",
                                "fail",
                                "hold_seconds",
                                "outdir",
                                "source_dir",
                                "stage",
                            ]
                        }
                    }
                }
            }
        )
    )
    return [
        nextflow,
        "-C",
        str(empty_config),
        "run",
        str(FIXTURE),
        "-lib",
        str(ROOT / "lib"),
        "-ansi-log",
        "false",
        "-w",
        str(tmp_path / "work"),
        "--outdir",
        str(outdir),
        "--source_dir",
        str(ROOT),
        "--schema_path",
        str(schema_path),
        "--api_token",
        "super-secret",
        "--endpoint",
        "https://user:password@example.test/data?token=query-secret",
        "--fail",
        str(fail).lower(),
        "--hold_seconds",
        str(hold_seconds),
    ]


def _environment(tmp_path: Path) -> dict[str, str]:
    return {
        **os.environ,
        "NXF_OFFLINE": "true",
        "NXF_HOME": str(tmp_path / "nextflow"),
    }


def _wait_for_running_manifest(path: Path) -> dict:
    deadline = time.monotonic() + 15
    while time.monotonic() < deadline:
        if path.exists():
            payload = json.loads(path.read_text())
            if payload.get("status") == "running":
                return payload
        time.sleep(0.05)
    raise AssertionError("Run manifest did not reach running state")


def _assert_provenance(payload: dict) -> None:
    expected_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()
    serialized = json.dumps(payload)

    assert payload["schema_version"] == 1
    assert payload["pipeline"] == "gaph_v2"
    assert payload["git_commit"] == expected_commit
    assert isinstance(payload["git_dirty"], bool)
    assert payload["profiles"] == ["standard"]
    assert payload["stage"] == "test"
    assert payload["nextflow_version"]
    assert payload["parameters"]["api_token"] == "<redacted>"
    assert payload["parameters"]["endpoint"] == (
        "https://<redacted>@example.test/data?token=<redacted>"
    )
    assert set(payload["parameters"]) == {
        "api_token",
        "endpoint",
        "fail",
        "hold_seconds",
        "outdir",
        "source_dir",
        "stage",
    }
    assert "super-secret" not in serialized
    assert "query-secret" not in serialized
    assert "user:password" not in serialized


def test_run_manifest_records_running_and_complete_states(tmp_path: Path) -> None:
    outdir = tmp_path / "run"
    manifest_path = outdir / "run_manifest.json"
    process = subprocess.Popen(
        _command(tmp_path, outdir, fail=False, hold_seconds=2),
        cwd=tmp_path,
        env=_environment(tmp_path),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    running = _wait_for_running_manifest(manifest_path)
    stdout, stderr = process.communicate(timeout=30)

    assert process.returncode == 0, stdout + stderr
    assert running["status"] == "running"
    assert running["completed_at"] is None
    assert running["success"] is None
    completed = json.loads(manifest_path.read_text())
    _assert_provenance(completed)
    assert completed["status"] == "complete"
    assert completed["success"] is True
    assert completed["exit_status"] == 0
    assert completed["completed_at"]
    assert stat.S_IMODE(manifest_path.stat().st_mode) == 0o644


def test_run_manifest_records_failed_completion(tmp_path: Path) -> None:
    outdir = tmp_path / "run"
    completed = subprocess.run(
        _command(tmp_path, outdir, fail=True, hold_seconds=0),
        cwd=tmp_path,
        env=_environment(tmp_path),
        text=True,
        capture_output=True,
        timeout=30,
    )

    assert completed.returncode != 0
    payload = json.loads((outdir / "run_manifest.json").read_text())
    _assert_provenance(payload)
    assert payload["status"] == "failed"
    assert payload["success"] is False
    assert payload["exit_status"] != 0
    assert payload["completed_at"]
