from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = PROJECT_ROOT / "analytics" / "slurm" / "submit_strategy_report.sh"
BATCH = PROJECT_ROOT / "analytics" / "slurm" / "strategy_report.sbatch"


def _read_argv(path: Path) -> list[str]:
    payload = path.read_bytes()
    assert payload.endswith(b"\0")
    return [part.decode() for part in payload[:-1].split(b"\0")]


def _fixture(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    project = tmp_path / "project"
    launcher = project / "analytics" / "slurm" / LAUNCHER.name
    batch = launcher.parent / BATCH.name
    launcher.parent.mkdir(parents=True)
    shutil.copy2(LAUNCHER, launcher)
    shutil.copy2(BATCH, batch)

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    (fake_bin / "git").write_text(
        "#!/usr/bin/env bash\n"
        "if [[ \"$*\" == *'rev-parse HEAD'* ]]; then printf 'abc123\\n'; fi\n"
        "exit 0\n"
    )
    (fake_bin / "sbatch").write_text(
        "#!/usr/bin/env bash\n"
        "printf '%s\\0' \"$@\" > \"$SBATCH_CAPTURE\"\n"
        "printf '99123\\n'\n"
    )
    (fake_bin / "git").chmod(0o755)
    (fake_bin / "sbatch").chmod(0o755)
    environment = os.environ.copy()
    environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
    environment["SBATCH_CAPTURE"] = str(tmp_path / "sbatch.args")
    return launcher, environment


def _run(path: Path) -> Path:
    (path / "annotation" / "variant_annotations").mkdir(parents=True)
    (path / "run_manifest.json").write_text("{}\n")
    (path / "annotation" / "manifest.json").write_text("{}\n")
    (path / "annotation" / "variant_annotations" / "manifest.json").write_text(
        "{}\n"
    )
    return path


def test_report_launcher_forwards_one_root_and_repeated_source_runs(
    tmp_path: Path,
) -> None:
    launcher, environment = _fixture(tmp_path)
    analytics_root = tmp_path / "analytics root"
    first = _run(tmp_path / "run one")
    second = _run(tmp_path / "run two")

    completed = subprocess.run(
        [
            "bash",
            str(launcher),
            "--analytics-root",
            str(analytics_root),
            "--run-dir",
            str(first),
            "--run-dir",
            str(second),
            "--report-name",
            "combined",
            "--",
            "--target-space-null",
            "--target-space-null-seed",
            "7",
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 0, completed.stderr
    arguments = _read_argv(Path(environment["SBATCH_CAPTURE"]))
    batch_index = arguments.index(str(launcher.parent / BATCH.name))
    assert arguments[batch_index + 1 :] == [
        str(analytics_root),
        "combined",
        "abc123",
        str(launcher.parents[2]),
        "2",
        str(first),
        str(second),
        "--target-space-null",
        "--target-space-null-seed",
        "7",
    ]
    assert (analytics_root / "slurm").is_dir()


def test_report_launcher_requires_external_analytics_root(tmp_path: Path) -> None:
    launcher, environment = _fixture(tmp_path)
    run = _run(tmp_path / "run")

    completed = subprocess.run(
        ["bash", str(launcher), "--run-dir", str(run), "--report-name", "report"],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 2
    assert "--analytics-root is required" in completed.stderr
    assert not Path(environment["SBATCH_CAPTURE"]).exists()


def test_report_launcher_does_not_write_inside_source_run(tmp_path: Path) -> None:
    launcher, environment = _fixture(tmp_path)
    run = _run(tmp_path / "run")

    completed = subprocess.run(
        [
            "bash",
            str(launcher),
            "--analytics-root",
            str(run / "analytics"),
            "--run-dir",
            str(run),
            "--report-name",
            "report",
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 2
    assert "must be outside source run" in completed.stderr
    assert not (run / "analytics").exists()
    assert not Path(environment["SBATCH_CAPTURE"]).exists()
