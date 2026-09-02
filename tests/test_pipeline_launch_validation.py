from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _run_nextflow(
    tmp_path: Path,
    *arguments: str,
    extra_env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    nextflow = shutil.which("nextflow")
    if nextflow is None:
        pytest.skip("nextflow is not installed")
    return subprocess.run(
        [
            nextflow,
            "run",
            str(ROOT),
            "-ansi-log",
            "false",
            "-w",
            str(tmp_path / "work"),
            *arguments,
        ],
        cwd=tmp_path,
        env={**os.environ, "NXF_ANSI_LOG": "false", **(extra_env or {})},
        text=True,
        capture_output=True,
        timeout=30,
    )


def test_unknown_parameter_stops_pipeline_before_execution(tmp_path: Path) -> None:
    completed = _run_nextflow(
        tmp_path,
        "--outdir",
        str(tmp_path / "out"),
        "--definitely_unknown",
        "value",
    )

    output = completed.stdout + completed.stderr
    assert completed.returncode != 0
    assert "definitely_unknown" in output
    assert "invalid input values" in output
    assert "Target annotation GFF3 not found" not in output


def test_help_uses_nf_schema_and_stops_before_pipeline_validation(tmp_path: Path) -> None:
    completed = _run_nextflow(tmp_path, "--help")

    output = completed.stdout + completed.stderr
    assert completed.returncode == 0
    assert "Typical pipeline command" in output
    assert "--help" in output
    assert "--ids_file" in output
    assert "Target annotation GFF3 not found" not in output


def test_full_help_stops_before_pipeline_validation(tmp_path: Path) -> None:
    completed = _run_nextflow(tmp_path, "--helpFull")

    output = completed.stdout + completed.stderr
    assert completed.returncode == 0
    assert "--vep_result_cache_tile_size_bp" in output
    assert "Target annotation GFF3 not found" not in output


def test_named_parameter_help_is_detailed(tmp_path: Path) -> None:
    completed = _run_nextflow(tmp_path, "--help", "ids_file")

    output = completed.stdout + completed.stderr
    assert completed.returncode == 0
    assert "--ids_file" in output
    assert "Path to the input file containing gene IDs" in output
    assert "--vep_backend" not in output
    assert "Target annotation GFF3 not found" not in output


def test_numeric_vep_release_passes_schema_validation(tmp_path: Path) -> None:
    completed = _run_nextflow(
        tmp_path,
        "--ids_file",
        str(ROOT / "assets/inputs/gene_ids/smoke_5_genes.txt"),
        "--vep_backend",
        "local",
        "--vep_release",
        "116",
        "--vep_cache_dir",
        str(tmp_path),
        "--target_annotation_gff3",
        str(tmp_path / "missing.gff.gz"),
    )

    output = completed.stdout + completed.stderr
    assert completed.returncode != 0
    assert "expected type" not in output
    assert "invalid input values" not in output
    assert "Target annotation GFF3 not found" in output


def test_numeric_vep_release_from_environment_passes_schema_validation(
    tmp_path: Path,
) -> None:
    completed = _run_nextflow(
        tmp_path,
        "--ids_file",
        str(ROOT / "assets/inputs/gene_ids/smoke_5_genes.txt"),
        "--vep_backend",
        "local",
        "--vep_cache_dir",
        str(tmp_path),
        "--target_annotation_gff3",
        str(tmp_path / "missing.gff.gz"),
        extra_env={"GAPH_VEP_RELEASE": "116"},
    )

    output = completed.stdout + completed.stderr
    assert completed.returncode != 0
    assert "expected type" not in output
    assert "invalid input values" not in output
    assert "Target annotation GFF3 not found" in output
