from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = PROJECT_ROOT / "scripts" / "slurm" / "run_and_report.sh"


def _read_argv(path: Path) -> list[str]:
    payload = path.read_bytes()
    assert payload.endswith(b"\0")
    return [value.decode() for value in payload[:-1].split(b"\0")]


def _launcher_fixture(tmp_path: Path, pipeline_status: int = 0) -> tuple[Path, dict[str, str]]:
    project = tmp_path / "project with spaces"
    launcher = project / "scripts" / "slurm" / "run_and_report.sh"
    launcher.parent.mkdir(parents=True)
    shutil.copy2(LAUNCHER, launcher)

    report_launcher = project / "analytics" / "slurm" / "submit_strategy_report.sh"
    report_launcher.parent.mkdir(parents=True)
    report_launcher.write_text(
        "#!/usr/bin/env bash\n"
        "printf '%s\\0' \"$@\" > \"$REPORT_ARGS_CAPTURE\"\n"
        "printf 'report\\n' >> \"$CALL_ORDER_CAPTURE\"\n"
    )

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    micromamba = fake_bin / "micromamba"
    micromamba.write_text(
        "#!/usr/bin/env bash\n"
        "printf '%s' \"$PWD\" > \"$PIPELINE_CWD_CAPTURE\"\n"
        "printf '%s\\0' \"$@\" > \"$PIPELINE_ARGS_CAPTURE\"\n"
        "printf 'pipeline\\n' >> \"$CALL_ORDER_CAPTURE\"\n"
        "exit \"$PIPELINE_STATUS\"\n"
    )
    micromamba.chmod(0o755)

    gaph_root = tmp_path / "gaph root"
    home = tmp_path / "home"
    home.mkdir()
    (home / ".gaph_v2_cluster_env.sh").write_text(
        f"export GAPH_ROOT={shlex.quote(str(gaph_root))}\n"
    )

    env = os.environ.copy()
    env.update(
        {
            "HOME": str(home),
            "PATH": f"{fake_bin}{os.pathsep}{env['PATH']}",
            "PIPELINE_STATUS": str(pipeline_status),
            "PIPELINE_ARGS_CAPTURE": str(tmp_path / "pipeline.args"),
            "PIPELINE_CWD_CAPTURE": str(tmp_path / "pipeline.cwd"),
            "REPORT_ARGS_CAPTURE": str(tmp_path / "report.args"),
            "CALL_ORDER_CAPTURE": str(tmp_path / "calls.txt"),
        }
    )
    return launcher, env


def _required_args(tmp_path: Path) -> tuple[list[str], Path, Path, Path]:
    ids = tmp_path / "ids with spaces.txt"
    ids.write_text("1\n")
    run = tmp_path / "results" / "run one"
    work = tmp_path / "work" / "run one"
    analytics = tmp_path / "analytics root"
    return (
        [
            "--ids-file",
            str(ids),
            "--run-dir",
            str(run),
            "--work-dir",
            str(work),
            "--analytics-root",
            str(analytics),
            "--report-name",
            "strategy.compare",
        ],
        ids,
        run,
        work,
    )


def test_launcher_runs_pipeline_then_forwards_report_arguments_exactly(tmp_path: Path) -> None:
    launcher, env = _launcher_fixture(tmp_path)
    required, ids, run, work = _required_args(tmp_path)
    report_args = [
        "--target-space-null",
        "--target-space-null-sample-size",
        "5000",
        "--label",
        "value with spaces",
        "β*?[x]",
        "",
    ]
    command = [
        "bash",
        str(launcher),
        *required,
        "--alignment-strategies",
        "minimap2_asm20,nucmer",
        "--slurm-cpus",
        "12",
        "--slurm-memory",
        "96G",
        "--slurm-time",
        "08:00:00",
        "--slurm-partition",
        "main-long",
        "--",
        *report_args,
    ]

    completed = subprocess.run(command, env=env, text=True, capture_output=True)

    assert completed.returncode == 0, completed.stderr
    project = launcher.parents[2]
    assert _read_argv(Path(env["PIPELINE_ARGS_CAPTURE"])) == [
        "run",
        "-p",
        str(tmp_path / "gaph root" / "envs" / "controller"),
        "nextflow",
        "run",
        str(project),
        "-profile",
        "slurm",
        "--ids_file",
        str(ids),
        "--outdir",
        str(run),
        "-work-dir",
        str(work),
        "-resume",
        "--alignment_strategies",
        "minimap2_asm20,nucmer",
    ]
    assert Path(env["PIPELINE_CWD_CAPTURE"]).read_text() == str(project)
    assert _read_argv(Path(env["REPORT_ARGS_CAPTURE"])) == [
        "--analytics-root",
        str(tmp_path / "analytics root"),
        "--run-dir",
        str(run),
        "--report-name",
        "strategy.compare",
        "--slurm-cpus",
        "12",
        "--slurm-memory",
        "96G",
        "--slurm-time",
        "08:00:00",
        "--slurm-partition",
        "main-long",
        "--",
        *report_args,
    ]
    assert Path(env["CALL_ORDER_CAPTURE"]).read_text().splitlines() == [
        "pipeline",
        "report",
    ]


def test_launcher_uses_existing_defaults_without_redeclaring_them(tmp_path: Path) -> None:
    launcher, env = _launcher_fixture(tmp_path)
    required, _ids, _run, _work = _required_args(tmp_path)

    completed = subprocess.run(
        ["bash", str(launcher), *required],
        env=env,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 0, completed.stderr
    pipeline_args = _read_argv(Path(env["PIPELINE_ARGS_CAPTURE"]))
    report_args = _read_argv(Path(env["REPORT_ARGS_CAPTURE"]))
    assert "-resume" in pipeline_args
    assert "--alignment_strategies" not in pipeline_args
    assert all(not value.startswith("--slurm-") for value in report_args)


def test_pipeline_failure_prevents_report_submission(tmp_path: Path) -> None:
    launcher, env = _launcher_fixture(tmp_path, pipeline_status=17)
    required, _ids, _run, _work = _required_args(tmp_path)

    completed = subprocess.run(
        ["bash", str(launcher), *required, "--", "--target-space-null"],
        env=env,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 17
    assert not Path(env["REPORT_ARGS_CAPTURE"]).exists()
    assert Path(env["CALL_ORDER_CAPTURE"]).read_text().splitlines() == ["pipeline"]


def test_launcher_rejects_invalid_inputs_before_starting_pipeline(tmp_path: Path) -> None:
    launcher, env = _launcher_fixture(tmp_path)
    required, _ids, _run, _work = _required_args(tmp_path)
    report_name_index = required.index("--report-name") + 1
    required[report_name_index] = "bad/name"

    completed = subprocess.run(
        ["bash", str(launcher), *required],
        env=env,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 2
    assert "--report-name may contain" in completed.stderr
    assert not Path(env["PIPELINE_ARGS_CAPTURE"]).exists()
