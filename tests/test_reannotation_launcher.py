from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from pathlib import Path

import pytest

COMMIT = "a" * 40


@pytest.mark.parametrize("dirty,head,expected", [(False, COMMIT, 0), (True, COMMIT, 2), (False, "b" * 40, 2)])
def test_reannotation_launcher_revision_gate_and_exact_arguments(tmp_path, dirty, head, expected):
    project = tmp_path / 'project with spaces'
    scripts = project / 'reannotation/slurm'
    shutil.copytree(Path(__file__).resolve().parents[1] / 'reannotation/slurm', scripts)
    scratch = tmp_path / 'scratch'
    python = scratch / 'envs/reannotation/bin/python'
    python.parent.mkdir(parents=True)
    python.write_text('#!/bin/sh\nexit 0\n')
    python.chmod(0o755)
    home = tmp_path / 'home'
    home.mkdir()
    (home / '.gaph_v2_cluster_env.sh').write_text(f'export GAPH_ROOT={shlex.quote(str(scratch))}\n')
    binaries = tmp_path / 'bin'
    binaries.mkdir()
    (binaries / 'git').write_text('#!/bin/bash\n'
        f'if [[ "$*" == *"rev-parse HEAD"* ]]; then echo {head}; '
        f'elif [[ "$*" == *"rev-parse origin/main"* ]]; then echo {COMMIT}; '
        f'elif [[ "$*" == *"status --porcelain"* ]]; then printf %s "{"dirty" if dirty else ""}"; fi\n')
    (binaries / 'sbatch').write_text('#!/bin/bash\nprintf "%s\\0" "$@" > "$CAPTURE"\necho 12345\n')
    for file in binaries.iterdir():
        file.chmod(0o755)
    capture = tmp_path / 'args'
    env = {**os.environ, 'HOME':str(home), 'PATH':str(binaries)+os.pathsep+os.environ['PATH'], 'CAPTURE':str(capture)}
    args = ['--root', str(scratch / 'supplement with spaces'), '--run-id', 'batch_001', '--run-id', 'batch_002']
    result = subprocess.run(['bash',str(scripts / 'submit_reannotation.sh'), '--expected-commit', COMMIT, '--', *args],
                            env=env, capture_output=True, text=True)
    assert result.returncode == expected, result.stderr
    if expected == 0:
        forwarded = capture.read_bytes().decode().rstrip('\0').split('\0')
        assert forwarded[-len(args):] == args
    else:
        assert not capture.exists()
