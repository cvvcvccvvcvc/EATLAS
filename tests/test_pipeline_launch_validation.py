from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_unknown_parameter_stops_pipeline_before_execution(tmp_path: Path) -> None:
    nextflow = shutil.which("nextflow")
    if nextflow is None:
        pytest.skip("nextflow is not installed")

    completed = subprocess.run(
        [
            nextflow,
            "run",
            str(ROOT),
            "-ansi-log",
            "false",
            "-w",
            str(tmp_path / "work"),
            "--outdir",
            str(tmp_path / "out"),
            "--definitely_unknown",
            "value",
        ],
        cwd=tmp_path,
        env={**os.environ, "NXF_ANSI_LOG": "false"},
        text=True,
        capture_output=True,
        timeout=30,
    )

    output = completed.stdout + completed.stderr
    assert completed.returncode != 0
    assert "definitely_unknown" in output
    assert "invalid input values" in output
    assert "Target annotation GFF3 not found" not in output
