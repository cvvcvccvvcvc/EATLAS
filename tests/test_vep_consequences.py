from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import pytest

from analytics.annotation import vep


def test_parse_record_selects_target_gene_and_most_severe_term() -> None:
    row = {"variant_key": "1:10:A>G", "gene_id": "25"}
    record = {
        "variant_class": "SNV",
        "transcript_consequences": [
            {
                "gene_id": "999",
                "transcript_id": "NM_other",
                "consequence_terms": ["stop_gained"],
            },
            {
                "gene_id": "25",
                "transcript_id": "NM_target",
                "consequence_terms": ["splice_region_variant", "missense_variant"],
                "impact": "MODERATE",
                "canonical": 1,
                "mane_select": "NM_target.1",
            },
        ],
    }

    result = vep._parse_record(row, record)

    assert result["status"] == "ok"
    assert result["primary_consequence"] == "missense_variant"
    assert result["consequence_terms"] == "missense_variant&splice_region_variant"
    assert result["transcript_id"] == "NM_target"


def test_sqlite_cache_reuses_completed_annotations(tmp_path: Path, monkeypatch) -> None:
    rows = pd.DataFrame(
        [
            {"variant_key": "1:10:A>G", "gene_id": "25", "chrom": "1", "pos": 10, "ref": "A", "alt": "G"},
            {"variant_key": "1:10:A>G", "gene_id": "25", "chrom": "1", "pos": 10, "ref": "A", "alt": "G"},
        ]
    )
    calls = []

    def fake_request(batch, *_args):
        calls.append(batch)
        return [
            {
                "variant_key": item["variant_key"],
                "gene_id": item["gene_id"],
                "status": "ok",
                "primary_consequence": "missense_variant",
                "consequence_terms": "missense_variant",
                "transcript_id": "NM_1",
                "mane_select": "NM_1",
                "canonical": True,
                "impact": "MODERATE",
                "variant_class": "SNV",
            }
            for item in batch
        ]

    monkeypatch.setattr(vep, "_request_batch", fake_request)
    cache = tmp_path / "vep.sqlite"
    first, first_summary = vep.annotate_vep_consequences(rows, cache, release="116")
    second, second_summary = vep.annotate_vep_consequences(rows, cache, release="116")

    assert len(calls) == 1
    assert len(calls[0]) == 1
    assert first_summary["queried"] == 1
    assert second_summary["cached"] == 1
    pd.testing.assert_frame_equal(first, second)


def test_local_vep_uses_offline_refseq_cache_and_reuses_sqlite(
    tmp_path: Path,
    monkeypatch,
) -> None:
    rows = pd.DataFrame(
        [
            {
                "variant_key": "17:43044295:T>A",
                "gene_id": "672",
                "chrom": "17",
                "pos": 43044295,
                "ref": "T",
                "alt": "A",
            }
        ]
    )
    commands = []

    def fake_run(command, **kwargs):
        commands.append((command, kwargs))
        output_path = Path(command[command.index("--output_file") + 1])
        output_path.write_text(
            "## ENSEMBL VARIANT EFFECT PREDICTOR v116.0\n"
            "#Uploaded_variation\tGene\tFeature\tConsequence\tIMPACT\tCANONICAL\tMANE_SELECT\tVARIANT_CLASS\n"
            "gaph_00000000\t672\tNM_007294.4\t3_prime_UTR_variant\tMODIFIER\tYES\t-\tSNV\n"
        )
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(vep.subprocess, "run", fake_run)
    cache_dir = tmp_path / "vep-cache"
    cache_dir.mkdir()
    sqlite_path = tmp_path / "vep.sqlite"

    first, first_summary = vep.annotate_vep_consequences(
        rows,
        sqlite_path,
        backend="local",
        release="116",
        vep_executable="gaph-vep",
        vep_cache_dir=cache_dir,
        vep_forks=4,
    )
    second, second_summary = vep.annotate_vep_consequences(
        rows,
        sqlite_path,
        backend="local",
        release="116",
        vep_executable="gaph-vep",
        vep_cache_dir=cache_dir,
        vep_forks=4,
    )

    assert len(commands) == 1
    command, kwargs = commands[0]
    assert command[0] == "gaph-vep"
    assert "--offline" in command
    assert "--refseq" in command
    assert "--use_given_ref" in command
    assert command[command.index("--fork") + 1] == "4"
    assert kwargs == {"capture_output": True, "text": True, "check": False}
    assert first.loc[0, "status"] == "ok"
    assert first.loc[0, "primary_consequence"] == "3_prime_UTR_variant"
    assert first.loc[0, "transcript_id"] == "NM_007294.4"
    assert first.loc[0, "mane_select"] == ""
    assert bool(first.loc[0, "canonical"])
    assert first_summary["backend"] == "local"
    assert first_summary["queried"] == 1
    assert second_summary["cached"] == 1
    pd.testing.assert_frame_equal(first, second)


def test_local_vep_requires_release_and_cache(tmp_path: Path) -> None:
    rows = pd.DataFrame(
        [
            {
                "variant_key": "1:10:A>G",
                "gene_id": "25",
                "chrom": "1",
                "pos": 10,
                "ref": "A",
                "alt": "G",
            }
        ]
    )

    with pytest.raises(ValueError, match="explicit release"):
        vep.annotate_vep_consequences(rows, tmp_path / "one.sqlite", backend="local")
    with pytest.raises(ValueError, match="cache directory"):
        vep.annotate_vep_consequences(
            rows,
            tmp_path / "two.sqlite",
            backend="local",
            release="116",
        )
