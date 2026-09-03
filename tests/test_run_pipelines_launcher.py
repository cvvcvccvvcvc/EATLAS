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
        "if [[ \"$*\" == *'status --porcelain=v1 --untracked-files=normal'* ]]; then\n"
        "  printf '%s' \"${GIT_STATUS_OUTPUT:-}\"\n"
        "elif [[ \"$*\" == *'rev-parse origin/main'* ]]; then\n"
        "  printf '%s\\n' \"${ORIGIN_COMMIT}\"\n"
        "elif [[ \"$*\" == *'rev-parse HEAD'* ]]; then\n"
        "  printf '%s\\n' \"${HEAD_COMMIT}\"\n"
        "fi\n"
    )
    git.chmod(0o755)

    micromamba = fake_bin / "micromamba"
    micromamba.write_text(
        "#!/usr/bin/env python3\n"
        "import hashlib, json, os, sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
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
        "    'session_id': session, 'git_commit': os.environ['HEAD_COMMIT'],\n"
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
            "HEAD_COMMIT": COMMIT,
            "ORIGIN_COMMIT": COMMIT,
            "PIPELINE_CAPTURE": str(tmp_path / "pipeline.jsonl"),
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
    assert "fetched origin/main" in completed.stderr
    assert not Path(environment["PIPELINE_CAPTURE"]).exists()
