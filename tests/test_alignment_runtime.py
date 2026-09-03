from types import SimpleNamespace

import pytest

from bin import alignment_runtime


def test_alignment_software_versions_parse_each_selected_tool(monkeypatch) -> None:
    responses = {
        "minimap2": SimpleNamespace(returncode=0, stdout="2.30-r1287\n", stderr=""),
        "nucmer": SimpleNamespace(
            returncode=0,
            stdout="nucmer (NUCmer) 4.0.1\n",
            stderr="",
        ),
        "bwa": SimpleNamespace(
            returncode=1,
            stdout="",
            stderr="Program: bwa\nVersion: 0.7.19-r1273\n",
        ),
        "samtools": SimpleNamespace(
            returncode=0,
            stdout="samtools 1.22.1\nUsing htslib 1.22.1\n",
            stderr="",
        ),
    }
    monkeypatch.setattr(
        alignment_runtime.subprocess,
        "run",
        lambda command, **_kwargs: responses[command[0]],
    )

    minimap2 = alignment_runtime.minimap2_software("minimap2")
    nucmer = alignment_runtime.nucmer_software("nucmer")
    bwa = alignment_runtime.bwa_software("bwa", "samtools")

    assert minimap2["tools"] == {"minimap2": "2.30-r1287"}
    assert nucmer["tools"] == {"nucmer": "nucmer (NUCmer) 4.0.1"}
    assert bwa["tools"] == {
        "bwa": "0.7.19-r1273",
        "samtools": "samtools 1.22.1",
    }
    assert minimap2["python"]
    assert minimap2["pysam"]


def test_version_probe_rejects_failed_standard_version_command(monkeypatch) -> None:
    monkeypatch.setattr(
        alignment_runtime.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(
            returncode=2,
            stdout="",
            stderr="bad invocation",
        ),
    )

    with pytest.raises(RuntimeError, match="exit 2"):
        alignment_runtime.minimap2_software("minimap2")
