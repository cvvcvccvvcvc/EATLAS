from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = PROJECT_ROOT / "scripts" / "slurm" / "run_pipelines.sh"
COMMIT = "a" * 40


def _fixture(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    project = tmp_path / "project with spaces"
    launcher = project / "scripts" / "slurm" / LAUNCHER.name
    launcher.parent.mkdir(parents=True)
    shutil.copy2(LAUNCHER, launcher)

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    git = fake_bin / "git"
    git.write_text(
        "#!/usr/bin/env bash\n"
        "repo=''\n"
        "if [[ \"${1:-}\" == '-C' ]]; then repo=$2; shift 2; fi\n"
        "if [[ \"$*\" == 'status --porcelain=v1 --untracked-files=normal' ]]; then\n"
        "  if [[ -n \"${PIPELINE_ROOT:-}\" && \"$repo\" == \"$PIPELINE_ROOT\" ]]; then\n"
        "    printf '%s' \"${PIPELINE_GIT_STATUS_OUTPUT:-}\"\n"
        "  else\n"
        "    printf '%s' \"${GIT_STATUS_OUTPUT:-}\"\n"
        "  fi\n"
        "elif [[ \"$*\" == 'rev-parse origin/main' ]]; then\n"
        "  printf '%s\\n' \"${ORIGIN_COMMIT}\"\n"
        "elif [[ \"$*\" == 'rev-parse HEAD' ]]; then\n"
        "  if [[ -n \"${PIPELINE_ROOT:-}\" && \"$repo\" == \"$PIPELINE_ROOT\" ]]; then\n"
        "    printf '%s\\n' \"${PIPELINE_HEAD_COMMIT}\"\n"
        "  else\n"
        "    printf '%s\\n' \"${LAUNCHER_HEAD_COMMIT}\"\n"
        "  fi\n"
        "elif [[ \"$*\" == merge-base\\ --is-ancestor* ]]; then\n"
        "  exit \"${ANCESTOR_STATUS:-0}\"\n"
        "fi\n"
    )
    git.chmod(0o755)

    micromamba = fake_bin / "micromamba"
    micromamba.write_text(
        "#!/usr/bin/env python3\n"
        "import hashlib, json, os, shutil, sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "nextflow_index = args.index('nextflow')\n"
        "nextflow_args = args[nextflow_index + 1:]\n"
        "cleanup_state = Path(os.environ['CLEANUP_STATE'])\n"
        "if nextflow_args[0] == 'clean':\n"
        "    session = nextflow_args[-1]\n"
        "    entries = [json.loads(line) for line in cleanup_state.read_text().splitlines()] if cleanup_state.exists() else []\n"
        "    paths = [Path(entry['path']) for entry in entries if entry['session'] == session and Path(entry['path']).exists()]\n"
        "    outside = os.environ.get('CLEANUP_OUTSIDE_PATH')\n"
        "    if outside and '-n' in nextflow_args:\n"
        "        paths.append(Path(outside))\n"
        "    if '-n' in nextflow_args:\n"
        "        for path in paths:\n"
        "            print(f'Would remove {path}')\n"
        "    else:\n"
        "        for path in paths:\n"
        "            if path.exists():\n"
        "                shutil.rmtree(path)\n"
        "        with Path(os.environ['CLEANUP_CAPTURE']).open('a') as handle:\n"
        "            handle.write(session + '\\n')\n"
        "    raise SystemExit(0)\n"
        "def value(flag, default=None):\n"
        "    return args[args.index(flag) + 1] if flag in args else default\n"
        "ids_file = value('--ids_file')\n"
        "outdir = Path(value('--outdir'))\n"
        "run_name = outdir.name\n"
        "manifest_path = outdir / 'run_manifest.json'\n"
        "old = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}\n"
        "session = value('-resume', old.get('session_id', f'session-{run_name}'))\n"
        "failed = run_name == os.environ.get('FAIL_RUN')\n"
        "outdir.mkdir(parents=True, exist_ok=True)\n"
        "work_dir = Path(value('-work-dir'))\n"
        "attempt = sum(1 for entry in (cleanup_state.read_text().splitlines() if cleanup_state.exists() else []) if json.loads(entry)['session'] == session)\n"
        "task_dir = work_dir / 'aa' / f'{run_name}-{attempt}'\n"
        "task_dir.mkdir(parents=True, exist_ok=True)\n"
        "(task_dir / 'task-output').write_text(run_name + '\\n')\n"
        "with cleanup_state.open('a') as handle:\n"
        "    handle.write(json.dumps({'session': session, 'path': str(task_dir)}) + '\\n')\n"
        "inventory_descriptor = None\n"
        "if not failed:\n"
        "    inventory_bytes = (json.dumps({'schema_version': 1}) + '\\n').encode()\n"
        "    inventory_path = outdir / 'evidence_inventory.json'\n"
        "    inventory_path.write_bytes(inventory_bytes)\n"
        "    inventory_descriptor = {\n"
        "        'path': inventory_path.name, 'schema_version': 1,\n"
        "        'size_bytes': len(inventory_bytes),\n"
        "        'sha256': hashlib.sha256(inventory_bytes).hexdigest(),\n"
        "    }\n"
        "manifest = {\n"
        "    'schema_version': 3, 'pipeline': 'gaph_v2',\n"
        "    'status': 'failed' if failed else 'complete',\n"
        "    'success': not failed, 'exit_status': 17 if failed else 0,\n"
        "    'session_id': session, 'git_commit': os.environ['PIPELINE_HEAD_COMMIT'],\n"
        "    'git_dirty': False,\n"
        "    'evidence_inventory': inventory_descriptor,\n"
        "    'parameters': {\n"
        "        'ids_file': ids_file, 'outdir': str(outdir),\n"
        "        'alignment_strategies': value('--alignment_strategies', 'default'),\n"
        "        'fetch_max_forks': int(value('--fetch_max_forks', 2)),\n"
        "        'alignment_max_forks': int(value('--alignment_max_forks', 4)),\n"
        "        'annotation_max_forks': int(value('--annotation_max_forks', 4)),\n"
        "    },\n"
        "}\n"
        "manifest_path.write_text(json.dumps(manifest) + '\\n')\n"
        "with Path(os.environ['PIPELINE_CAPTURE']).open('a') as handle:\n"
        "    handle.write(json.dumps(args) + '\\n')\n"
        "raise SystemExit(17 if failed else 0)\n"
    )
    micromamba.chmod(0o755)

    gaph_root = tmp_path / "gaph root"
    work_root = tmp_path / "work root"
    home = tmp_path / "home"
    home.mkdir()
    (home / ".gaph_v2_cluster_env.sh").write_text(
        f"export GAPH_ROOT={shlex.quote(str(gaph_root))}\n"
        f"export GAPH_WORK_DIR={shlex.quote(str(work_root))}\n"
    )
    environment = os.environ.copy()
    environment.update(
        {
            "HOME": str(home),
            "PATH": f"{fake_bin}{os.pathsep}{environment['PATH']}",
            "LAUNCHER_HEAD_COMMIT": COMMIT,
            "PIPELINE_HEAD_COMMIT": COMMIT,
            "ORIGIN_COMMIT": COMMIT,
            "PIPELINE_CAPTURE": str(tmp_path / "pipeline.jsonl"),
            "CLEANUP_CAPTURE": str(tmp_path / "cleanup.txt"),
            "CLEANUP_STATE": str(tmp_path / "cleanup.jsonl"),
        }
    )
    return launcher, environment


def _ids(tmp_path: Path, *names: str) -> list[Path]:
    paths = []
    for index, name in enumerate(names, start=1):
        path = tmp_path / "inputs" / name
        path.parent.mkdir(exist_ok=True)
        path.write_text(f"{index}\n")
        paths.append(path)
    return paths


def _calls(environment: dict[str, str]) -> list[list[str]]:
    return [
        json.loads(line)
        for line in Path(environment["PIPELINE_CAPTURE"]).read_text().splitlines()
    ]


def test_runs_inputs_sequentially_with_derived_result_and_work_paths(tmp_path: Path) -> None:
    launcher, environment = _fixture(tmp_path)
    first, second = _ids(tmp_path, "batch_001.txt", "batch_002.ids")
    results = tmp_path / "results" / "all genes"

    completed = subprocess.run(
        [
            "bash",
            str(launcher),
            "--results-root",
            str(results),
            "--expected-commit",
            COMMIT,
            "--alignment-strategies",
            "minimap2_asm20,nucmer",
            "--alignment-max-forks",
            "7",
            str(first),
            str(second),
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 2
    assert "--results-root basename" in completed.stderr
    assert not Path(environment["PIPELINE_CAPTURE"]).exists()

    results = tmp_path / "results" / "all_genes"
    completed = subprocess.run(
        [
            "bash",
            str(launcher),
            "--results-root",
            str(results),
            "--expected-commit",
            COMMIT,
            "--alignment-strategies",
            "minimap2_asm20,nucmer",
            "--alignment-max-forks",
            "7",
            str(first),
            str(second),
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 0, completed.stderr
    calls = _calls(environment)
    assert [Path(call[call.index("--outdir") + 1]).name for call in calls] == [
        "batch_001",
        "batch_002",
    ]
    assert all("-resume" not in call for call in calls)
    assert calls[0][calls[0].index("-work-dir") + 1] == str(
        tmp_path / "work root" / "all_genes"
    )
    assert calls[1][calls[1].index("--alignment_strategies") + 1] == (
        "minimap2_asm20,nucmer"
    )
    assert calls[1][calls[1].index("--alignment_max_forks") + 1] == "7"
    assert Path(environment["CLEANUP_CAPTURE"]).read_text().splitlines() == [
        "session-batch_001",
        "session-batch_002",
    ]
    assert not any((tmp_path / "work root" / "all_genes").rglob("task-output"))


def test_stops_on_failure_then_skips_complete_and_resumes_exact_session(tmp_path: Path) -> None:
    launcher, environment = _fixture(tmp_path)
    inputs = _ids(tmp_path, "batch_001.txt", "batch_002.txt", "batch_003.txt")
    results = tmp_path / "results" / "group"
    command = [
        "bash",
        str(launcher),
        "--results-root",
        str(results),
        "--expected-commit",
        COMMIT,
        *map(str, inputs),
    ]
    environment["FAIL_RUN"] = "batch_002"

    failed = subprocess.run(command, env=environment, text=True, capture_output=True)

    assert failed.returncode == 17
    assert [Path(call[call.index("--outdir") + 1]).name for call in _calls(environment)] == [
        "batch_001",
        "batch_002",
    ]

    environment.pop("FAIL_RUN")
    resumed = subprocess.run(command, env=environment, text=True, capture_output=True)

    assert resumed.returncode == 0, resumed.stderr
    calls = _calls(environment)
    assert [Path(call[call.index("--outdir") + 1]).name for call in calls] == [
        "batch_001",
        "batch_002",
        "batch_002",
        "batch_003",
    ]
    resumed_call = calls[2]
    assert resumed_call[resumed_call.index("-resume") + 1] == "session-batch_002"
    assert resumed_call[resumed_call.index("--alignment_strategies") + 1] == "default"
    assert "Skipping completed run batch_001" in resumed.stdout
    assert Path(environment["CLEANUP_CAPTURE"]).read_text().splitlines() == [
        "session-batch_001",
        "session-batch_002",
        "session-batch_003",
    ]
    assert not any((tmp_path / "work root" / "group").rglob("task-output"))


def test_resume_may_reduce_but_not_increase_alignment_concurrency(tmp_path: Path) -> None:
    launcher, environment = _fixture(tmp_path)
    [ids_file] = _ids(tmp_path, "batch.txt")
    results = tmp_path / "results" / "group"
    base_command = [
        "bash",
        str(launcher),
        "--results-root",
        str(results),
        "--expected-commit",
        COMMIT,
    ]
    environment["FAIL_RUN"] = "batch"

    failed = subprocess.run(
        [*base_command, "--alignment-max-forks", "4", str(ids_file)],
        env=environment,
        text=True,
        capture_output=True,
    )
    assert failed.returncode == 17

    environment.pop("FAIL_RUN")
    refused = subprocess.run(
        [*base_command, "--alignment-max-forks", "5", str(ids_file)],
        env=environment,
        text=True,
        capture_output=True,
    )
    assert refused.returncode == 2
    assert "may only reduce concurrency" in refused.stderr
    assert len(_calls(environment)) == 1

    resumed = subprocess.run(
        [*base_command, "--alignment-max-forks", "2", str(ids_file)],
        env=environment,
        text=True,
        capture_output=True,
    )
    assert resumed.returncode == 0, resumed.stderr
    resumed_call = _calls(environment)[1]
    assert resumed_call[resumed_call.index("-resume") + 1] == "session-batch"
    assert resumed_call[resumed_call.index("--alignment_max_forks") + 1] == "2"
    assert "Reducing alignment concurrency for run batch from 4 to 2" in resumed.stdout


def test_refuses_completed_session_cleanup_outside_group_work(tmp_path: Path) -> None:
    launcher, environment = _fixture(tmp_path)
    [ids_file] = _ids(tmp_path, "batch.txt")
    outside = tmp_path / "outside-work"
    outside.mkdir()
    environment["CLEANUP_OUTSIDE_PATH"] = str(outside)

    completed = subprocess.run(
        [
            "bash",
            str(launcher),
            "--results-root",
            str(tmp_path / "results" / "group"),
            "--expected-commit",
            COMMIT,
            str(ids_file),
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 2
    assert "outside the group work directory" in completed.stderr
    assert outside.is_dir()
    assert not Path(environment["CLEANUP_CAPTURE"]).exists()


def test_rejects_completed_run_with_modified_evidence_inventory(tmp_path: Path) -> None:
    launcher, environment = _fixture(tmp_path)
    [ids_file] = _ids(tmp_path, "batch.txt")
    results = tmp_path / "results" / "group"
    command = [
        "bash",
        str(launcher),
        "--results-root",
        str(results),
        "--expected-commit",
        COMMIT,
        str(ids_file),
    ]

    first = subprocess.run(command, env=environment, text=True, capture_output=True)
    assert first.returncode == 0, first.stderr
    (results / "batch" / "evidence_inventory.json").write_text("{}\n")

    repeated = subprocess.run(command, env=environment, text=True, capture_output=True)

    assert repeated.returncode == 2
    assert "cannot read run manifest" in repeated.stderr
    assert len(_calls(environment)) == 1


def test_rejects_duplicate_run_names_before_launch(tmp_path: Path) -> None:
    launcher, environment = _fixture(tmp_path)
    first = tmp_path / "one" / "batch.txt"
    second = tmp_path / "two" / "batch.txt"
    first.parent.mkdir()
    second.parent.mkdir()
    first.write_text("1\n")
    second.write_text("2\n")

    completed = subprocess.run(
        [
            "bash",
            str(launcher),
            "--results-root",
            str(tmp_path / "results" / "group"),
            "--expected-commit",
            COMMIT,
            str(first),
            str(second),
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 2
    assert "duplicate run name: batch" in completed.stderr
    assert not Path(environment["PIPELINE_CAPTURE"]).exists()


def test_revision_gate_runs_before_pipeline(tmp_path: Path) -> None:
    launcher, environment = _fixture(tmp_path)
    [ids_file] = _ids(tmp_path, "batch.txt")
    environment["ANCESTOR_STATUS"] = "1"

    completed = subprocess.run(
        [
            "bash",
            str(launcher),
            "--results-root",
            str(tmp_path / "results" / "group"),
            "--expected-commit",
            COMMIT,
            str(ids_file),
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 2
    assert "is not reachable from fetched origin/main" in completed.stderr
    assert not Path(environment["PIPELINE_CAPTURE"]).exists()


def test_revision_gate_requires_current_clean_launcher(tmp_path: Path) -> None:
    launcher, environment = _fixture(tmp_path)
    [ids_file] = _ids(tmp_path, "batch.txt")
    environment["ORIGIN_COMMIT"] = "b" * 40

    completed = subprocess.run(
        [
            "bash",
            str(launcher),
            "--results-root",
            str(tmp_path / "results" / "group"),
            "--expected-commit",
            COMMIT,
            str(ids_file),
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 2
    assert "launcher HEAD" in completed.stderr
    assert not Path(environment["PIPELINE_CAPTURE"]).exists()


def test_runs_historical_commit_from_separate_clean_checkout(tmp_path: Path) -> None:
    launcher, environment = _fixture(tmp_path)
    [ids_file] = _ids(tmp_path, "batch.txt")
    pipeline_root = tmp_path / "historical checkout"
    pipeline_root.mkdir()
    historical_commit = "b" * 40
    environment.update(
        {
            "PIPELINE_ROOT": str(pipeline_root),
            "PIPELINE_HEAD_COMMIT": historical_commit,
        }
    )

    completed = subprocess.run(
        [
            "bash",
            str(launcher),
            "--results-root",
            str(tmp_path / "results" / "group"),
            "--pipeline-root",
            str(pipeline_root),
            "--expected-commit",
            historical_commit,
            str(ids_file),
        ],
        env=environment,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 0, completed.stderr
    call = _calls(environment)[0]
    nextflow_index = call.index("nextflow")
    pipeline_run_index = call.index("run", nextflow_index)
    assert call[pipeline_run_index + 1] == str(pipeline_root)
    assert f"Pipeline commit: {historical_commit}" in completed.stdout
