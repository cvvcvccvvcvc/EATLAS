from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from analytics import vep_annotation as bulk


def test_partitioned_vep_annotation_resumes_and_finalizes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    run_dir = tmp_path / "run"
    annotation_dir = run_dir / "annotation"
    annotation_dir.mkdir(parents=True)
    source = annotation_dir / "variant_annotations.tsv.gz"
    pd.DataFrame(
        [
            {"variant_key": "1:10:A>G", "gene_id": "1", "strategies": "s1"},
            {"variant_key": "1:11:C>CT", "gene_id": "1", "strategies": "s1|s2"},
            {"variant_key": "1:12:A>", "gene_id": "1", "strategies": "s2"},
        ]
    ).to_csv(source, sep="\t", index=False, compression="gzip")
    (annotation_dir / "manifest.json").write_text(
        json.dumps({"variant_context_count": 3}) + "\n"
    )
    outdir = run_dir / "analytics" / "vep_consequences"

    plan = bulk.prepare_partitions(
        annotation_tsv=source,
        outdir=outdir,
        partition_size=2,
    )
    cached_plan = bulk.prepare_partitions(
        annotation_tsv=source,
        outdir=outdir,
        partition_size=2,
    )

    assert plan["partition_count"] == 2
    assert plan["row_count"] == 3
    assert not plan["cache_hit"]
    assert cached_plan["cache_hit"]

    first_input = outdir / plan["partitions"][0]["path"]
    second_input = outdir / plan["partitions"][1]["path"]
    first_identity = bulk._file_identity(first_input)
    second_input.unlink()
    resumed_plan = bulk.prepare_partitions(
        annotation_tsv=source,
        outdir=outdir,
        partition_size=2,
    )
    assert not resumed_plan["cache_hit"]
    assert bulk._file_identity(first_input) == first_identity
    assert second_input.exists()

    calls = []

    def fake_annotate(frame, _cache_path, **kwargs):
        calls.append((frame.copy(), kwargs))
        result = frame[["variant_key", "gene_id"]].drop_duplicates().copy()
        result["status"] = "ok"
        result["primary_consequence"] = "intron_variant"
        result["consequence_terms"] = "intron_variant"
        result["transcript_id"] = "NM_1"
        result["mane_select"] = ""
        result["canonical"] = False
        result["impact"] = "MODIFIER"
        result["variant_class"] = "SNV"
        return result, {
            "status": "complete",
            "backend": kwargs["backend"],
            "release": kwargs["release"],
            "requested": len(result),
        }

    monkeypatch.setattr(bulk, "annotate_vep_consequences", fake_annotate)
    first = bulk.annotate_partition(
        outdir=outdir,
        partition_index=1,
        backend="local",
        release="116",
        vep_executable="vep",
        vep_cache_dir=tmp_path / "cache",
        vep_forks=4,
    )
    second = bulk.annotate_partition(
        outdir=outdir,
        partition_index=2,
        backend="local",
        release="116",
        vep_executable="vep",
        vep_cache_dir=tmp_path / "cache",
        vep_forks=4,
    )
    cached_first = bulk.annotate_partition(
        outdir=outdir,
        partition_index=1,
        backend="local",
        release="116",
        vep_executable="vep",
        vep_cache_dir=tmp_path / "cache",
        vep_forks=4,
    )

    assert first["status_counts"] == {"ok": 2}
    assert second["status_counts"] == {"invalid_variant_key": 1}
    assert cached_first["cache_hit"]
    assert len(calls) == 1
    assert calls[0][1]["vep_forks"] == 4

    manifest = bulk.finalize_annotations(outdir=outdir)
    cached_manifest = bulk.finalize_annotations(outdir=outdir)
    output = pd.read_csv(
        outdir / "variant_annotations.vep.tsv.gz",
        sep="\t",
        compression="gzip",
        dtype=str,
        keep_default_na=False,
    )

    assert manifest["row_count"] == 3
    assert manifest["partition_count"] == 2
    assert manifest["status_counts"] == {"invalid_variant_key": 1, "ok": 2}
    assert cached_manifest["cache_hit"]
    assert output["variant_key"].tolist() == ["1:10:A>G", "1:11:C>CT", "1:12:A>"]
    assert output["vep_status"].tolist() == ["ok", "ok", "invalid_variant_key"]
    assert output.loc[0, "vep_primary_consequence"] == "intron_variant"


def test_prepare_refuses_to_reuse_a_different_partition_contract(tmp_path: Path) -> None:
    source = tmp_path / "variant_annotations.tsv.gz"
    pd.DataFrame([{"variant_key": "1:10:A>G", "gene_id": "1"}]).to_csv(
        source,
        sep="\t",
        index=False,
        compression="gzip",
    )
    outdir = tmp_path / "vep"
    bulk.prepare_partitions(annotation_tsv=source, outdir=outdir, partition_size=1)

    with pytest.raises(ValueError, match="different input contract"):
        bulk.prepare_partitions(annotation_tsv=source, outdir=outdir, partition_size=2)
