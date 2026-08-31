from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = PROJECT_ROOT / "analytics" / "slurm" / "submit_strategy_report.sh"
BATCH = PROJECT_ROOT / "analytics" / "slurm" / "strategy_report.sbatch"
COMMIT = "a" * 40


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
        f"if [[ \"$*\" == *'rev-parse origin/main'* ]]; then printf '{COMMIT}\\n';\n"
        f"elif [[ \"$*\" == *'rev-parse HEAD'* ]]; then printf '{COMMIT}\\n'; fi\n"
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


def _worker_fixture(tmp_path: Path) -> tuple[Path, Path, dict[str, str]]:
    project = tmp_path / "project with spaces"
    worker = project / "analytics" / "slurm" / BATCH.name
    worker.parent.mkdir(parents=True)
    shutil.copy2(BATCH, worker)

    fake_bin = tmp_path / "worker bin"
    fake_bin.mkdir()
    git = fake_bin / "git"
    git.write_text(
        "#!/usr/bin/env bash\n"
        "if [[ \"$*\" == *'rev-parse HEAD'* ]]; then\n"
        f"  printf '{COMMIT}\\n'\n"
        "  exit 0\n"
        "fi\n"
        "if [[ \"$*\" == *'diff --cached --quiet'* ]]; then\n"
        "  exit \"${GIT_CACHED_DIFF_STATUS:-0}\"\n"
        "fi\n"
        "if [[ \"$*\" == *'diff --quiet'* ]]; then\n"
        "  exit \"${GIT_DIFF_STATUS:-0}\"\n"
        "fi\n"
        "exit 0\n"
    )
    git.chmod(0o755)

    gaph_root = tmp_path / "gaph root"
    analytics_python = gaph_root / "envs" / "analytics" / "bin" / "python"
    analytics_python.parent.mkdir(parents=True)
    analytics_python.write_text(
        "#!/usr/bin/env bash\n"
        "printf '%s\\0' \"$@\" > \"$PYTHON_CAPTURE\"\n"
    )
    analytics_python.chmod(0o755)

    home = tmp_path / "home"
    home.mkdir()
    (home / ".gaph_v2_cluster_env.sh").write_text(
        f"export GAPH_ROOT={shlex.quote(str(gaph_root))}\n"
    )
    environment = os.environ.copy()
    environment.update(
        {
            "HOME": str(home),
            "PATH": f"{fake_bin}{os.pathsep}{environment['PATH']}",
            "PYTHON_CAPTURE": str(tmp_path / "python.args"),
        }
    )
    return worker, project, environment


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
            "--expected-commit",
            COMMIT,
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
        COMMIT,
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
        [
            "bash",
            str(launcher),
            "--run-dir",
            str(run),
            "--report-name",
            "report",
            "--expected-commit",
            COMMIT,
        ],
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
            "--expected-commit",
            COMMIT,
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 2
    assert "must be outside source run" in completed.stderr
    assert not (run / "analytics").exists()
    assert not Path(environment["SBATCH_CAPTURE"]).exists()


def test_report_worker_forwards_sources_and_report_arguments(tmp_path: Path) -> None:
    worker, project, environment = _worker_fixture(tmp_path)
    analytics_root = tmp_path / "analytics root"
    first = tmp_path / "run one"
    second = tmp_path / "run two"

    completed = subprocess.run(
        [
            "bash",
            str(worker),
            str(analytics_root),
            "combined",
            COMMIT,
            str(project),
            "2",
            str(first),
            str(second),
            "--target-space-null",
            "--target-space-null-seed",
            "7",
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 0, completed.stderr
    assert _read_argv(Path(environment["PYTHON_CAPTURE"])) == [
        "-m",
        "analytics.strategy_report",
        "--analytics-root",
        str(analytics_root),
        "--report-name",
        "combined",
        "--run-dir",
        str(first),
        "--run-dir",
        str(second),
        "--target-space-null",
        "--target-space-null-seed",
        "7",
    ]
    assert (analytics_root / "slurm" / "combined.git_commit").read_text() == (
        f"{COMMIT}\n"
    )


@pytest.mark.parametrize(
    ("status_variable", "message"),
    [
        ("GIT_DIFF_STATUS", "tracked changes"),
        ("GIT_CACHED_DIFF_STATUS", "staged changes"),
    ],
)
def test_report_worker_rejects_changes_after_submission(
    tmp_path: Path,
    status_variable: str,
    message: str,
) -> None:
    worker, project, environment = _worker_fixture(tmp_path)
    environment[status_variable] = "1"

    completed = subprocess.run(
        [
            "bash",
            str(worker),
            str(tmp_path / "analytics"),
            "report",
            COMMIT,
            str(project),
            "1",
            str(tmp_path / "run"),
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 1
    assert message in completed.stderr
    assert not Path(environment["PYTHON_CAPTURE"]).exists()
