from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
WRAPPER = PROJECT_ROOT / "bin" / "gaph-vep116"


@pytest.mark.parametrize("runtime_name", ["singularity", "apptainer"])
def test_vep_wrapper_discovers_runtime_and_forwards_arguments(
    tmp_path: Path,
    runtime_name: str,
) -> None:
    gaph_root = tmp_path / "gaph"
    image = gaph_root / "containers" / "ensembl-vep_release_116.0.sif"
    image.parent.mkdir(parents=True)
    image.write_bytes(b"test image")

    runtime_dir = tmp_path / "runtime"
    runtime_dir.mkdir()
    capture_path = tmp_path / "arguments.txt"
    runtime = runtime_dir / runtime_name
    runtime.write_text(
        "#!/usr/bin/env bash\nprintf '%s\\n' \"$@\" > \"$CAPTURE_PATH\"\n",
        encoding="utf-8",
    )
    runtime.chmod(0o755)

    env = os.environ.copy()
    env.pop("GAPH_CONTAINER_RUNTIME", None)
    env.update(
        {
            "CAPTURE_PATH": str(capture_path),
            "GAPH_ROOT": str(gaph_root),
            "PATH": f"{runtime_dir}:{env['PATH']}",
        }
    )
    subprocess.run(
        [str(WRAPPER), "--offline", "--cache"],
        env=env,
        check=True,
    )

    assert capture_path.read_text(encoding="utf-8").splitlines() == [
        "exec",
        "--bind",
        f"{gaph_root}:{gaph_root}",
        str(image),
        "vep",
        "--offline",
        "--cache",
    ]
    assert (gaph_root / "singularity" / "cache").is_dir()
    assert (gaph_root / "singularity" / "tmp").is_dir()
