from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

from bin import annotate_vep_partition as vep_partition
from bin.annotate_events import (
    VARIANT_ANNOTATION_FIELDS,
    write_variant_annotation_shards,
)


def test_source_annotation_shards_are_bounded_headered_and_ordered(
    tmp_path: Path,
) -> None:
    rows = [
        {"variant_key": f"1:{position}:A>G", "gene_id": "1"}
        for position in range(1, 6)
    ]

    dataset = write_variant_annotation_shards(
        tmp_path / "variant_annotation_shards",
        rows,
        shard_size=2,
    )

    assert dataset["row_count"] == 5
    assert dataset["shard_count"] == 3
    assert [item["row_count"] for item in dataset["shards"]] == [2, 2, 1]
    observed = []
    for item in dataset["shards"]:
        frame = pd.read_csv(
            tmp_path / "variant_annotation_shards" / item["path"],
            sep="\t",
            compression="gzip",
            dtype=str,
            keep_default_na=False,
        )
        assert list(frame.columns) == VARIANT_ANNOTATION_FIELDS
        observed.extend(frame["variant_key"].tolist())
    assert observed == [row["variant_key"] for row in rows]


def test_empty_source_annotation_dataset_has_one_runnable_shard(tmp_path: Path) -> None:
    dataset = write_variant_annotation_shards(
        tmp_path / "variant_annotation_shards",
        [],
        shard_size=2,
    )

    assert dataset["row_count"] == 0
    assert dataset["shard_count"] == 1
    frame = pd.read_csv(
        tmp_path / "variant_annotation_shards" / dataset["shards"][0]["path"],
        sep="\t",
        compression="gzip",
        dtype=str,
        keep_default_na=False,
    )
    assert frame.empty
    assert list(frame.columns) == VARIANT_ANNOTATION_FIELDS


def test_vep_partition_preserves_rows_and_marks_invalid_keys(
    tmp_path: Path,
    monkeypatch,
) -> None:
    source = tmp_path / "source.tsv.gz"
    pd.DataFrame(
        [
            {"variant_key": "1:10:A>G", "gene_id": "1", "strategies": "s1"},
            {"variant_key": "", "gene_id": "1", "strategies": "s2"},
            {"variant_key": "1:11:C>T", "gene_id": "2", "strategies": "s1,s2"},
        ]
    ).to_csv(source, sep="\t", index=False, compression="gzip")

    calls = []

    def fake_annotate(requests, _cache_path, **kwargs):
        calls.append((requests.copy(), kwargs))
        result = requests[["variant_key", "gene_id"]].copy()
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
            "options": {"pick_allele_gene": 1},
            "requested": len(result),
            "queried": len(result),
            "cached": 0,
            "status_counts": {"ok": len(result)},
        }

    monkeypatch.setattr(vep_partition, "annotate_vep_consequences", fake_annotate)
    outdir = tmp_path / "out"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "annotate_vep_partition.py",
            "--input-tsv",
            str(source),
            "--outdir",
            str(outdir),
            "--partition-id",
            "partition_000001",
            "--shard-id",
            "shard_000001",
            "--expected-row-count",
            "3",
            "--vep-backend",
            "rest",
            "--vep-release",
            "116",
        ],
    )

    vep_partition.main()

    output = pd.read_csv(
        outdir / "variant_annotations.tsv.gz",
        sep="\t",
        compression="gzip",
        dtype=str,
        keep_default_na=False,
    )
    manifest = json.loads((outdir / "manifest.json").read_text())
    assert output["variant_key"].tolist() == ["1:10:A>G", "", "1:11:C>T"]
    assert output["vep_status"].tolist() == ["ok", "invalid_variant_key", "ok"]
    assert manifest["schema"] == vep_partition.SCHEMA
    assert manifest["config"]["backend"] == "rest"
    assert manifest["config"]["release"] == "116"
    assert manifest["status_counts"] == {"invalid_variant_key": 1, "ok": 2}
    assert calls[0][0]["variant_key"].tolist() == ["1:10:A>G", "1:11:C>T"]
