"""Actual software versions used by alignment producer tasks."""

from __future__ import annotations

import platform
import re
import subprocess

import pysam


VERSION_PROBE_TIMEOUT_SECONDS = 30


def minimap2_software(executable: str) -> dict[str, object]:
    return _runtime({"minimap2": _first_line([executable, "--version"])})


def nucmer_software(executable: str) -> dict[str, object]:
    return _runtime({"nucmer": _first_line([executable, "--version"])})


def bwa_software(bwa_executable: str, samtools_executable: str) -> dict[str, object]:
    bwa_output = _probe([bwa_executable], allow_nonzero=True)
    match = re.search(r"(?m)^\s*Version:\s*(\S+)", bwa_output)
    if match is None:
        raise RuntimeError(f"Could not parse BWA version from {bwa_executable}")
    return _runtime(
        {
            "bwa": match.group(1),
            "samtools": _first_line([samtools_executable, "--version"]),
        }
    )


def _runtime(tools: dict[str, str]) -> dict[str, object]:
    return {
        "python": platform.python_version(),
        "pysam": pysam.__version__,
        "tools": dict(sorted(tools.items())),
    }


def _first_line(command: list[str]) -> str:
    output = _probe(command)
    line = next((value.strip() for value in output.splitlines() if value.strip()), "")
    if not line:
        raise RuntimeError(f"Version probe returned no output: {command[0]}")
    return line


def _probe(command: list[str], *, allow_nonzero: bool = False) -> str:
    try:
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
            timeout=VERSION_PROBE_TIMEOUT_SECONDS,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(
            f"Could not probe software version: {command[0]}: {exc}"
        ) from exc
    output = "\n".join(value for value in (completed.stdout, completed.stderr) if value)
    if completed.returncode != 0 and not allow_nonzero:
        detail = output.strip()[-500:] or "no output"
        raise RuntimeError(
            f"Software version probe failed with exit {completed.returncode}: "
            f"{command[0]}: {detail}"
        )
    return output
