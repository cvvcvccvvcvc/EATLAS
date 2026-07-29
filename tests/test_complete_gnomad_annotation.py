from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

import complete_gnomad_annotation as completion  # noqa: E402


FAILURE_FIELDS = ["source", "scope", "chrom", "start", "end", "failure_type", "message"]
VARIANT_ANNOTATION_FIELDS = [
    "variant_key",
    "gene_id",
    "event_type",
    "ref",
    "alt",
    "lookup_status",
    "support_row_count",
    "support_ortholog_count",
    "strategies",
    "clinvar_sig",
    "gnomad_af",
    "gnomad_af_source",
    "gnomad_csq",
]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_completion_retries_failed_regions_and_preserves_other_columns(
    tmp_path: Path,
    monkeypatch,
) -> None:
    run_dir = tmp_path / "run"
    source = run_dir / "annotation"
    targets = run_dir / "fetch" / "sequences" / "targets"
    targets.mkdir(parents=True)
    with gzip.open(targets / "1.fa.gz", "wt") as handle:
        handle.write(">1\nACCG\n")
    write_tsv(
        run_dir / "fetch" / "genes.tsv.gz",
        ["gene_id", "genomic_accession", "chromosome", "begin", "end"],
        [
            {
                "gene_id": "1",
                "genomic_accession": "NC_000001.11",
                "chromosome": "1",
                "begin": 100,
                "end": 103,
            }
        ],
    )

    rows = [
        {
            "variant_key": "1:101:C>T",
            "gene_id": "1",
            "event_type": "snv",
            "ref": "C",
            "alt": "T",
            "lookup_status": "ok",
            "support_row_count": "3",
            "support_ortholog_count": "2",
            "strategies": "s1",
        },
        {
            "variant_key": "1:102:C>CT",
            "gene_id": "1",
            "event_type": "ins",
            "ref": "",
            "alt": "T",
            "lookup_status": "ok",
            "support_row_count": "1",
            "support_ortholog_count": "1",
            "strategies": "s1",
        },
    ]
    write_tsv(source / "variant_annotations.tsv.gz", VARIANT_ANNOTATION_FIELDS, rows)
    write_tsv(
        source / "variant_strategy_support.tsv.gz",
        ["variant_key", "gene_id", "strategy"],
        [{"variant_key": "1:101:C>T", "gene_id": "1", "strategy": "s1"}],
    )
    write_tsv(
        source / "failures.tsv.gz",
        FAILURE_FIELDS,
        [
            {
                "source": "gnomad",
                "scope": "region",
                "chrom": "1",
                "start": "101",
                "end": "102",
                "failure_type": "TimeoutError",
                "message": "timed out",
            },
            {
                "source": "clinvar",
                "scope": "region",
                "chrom": "1",
                "start": "101",
                "end": "102",
                "failure_type": "test",
                "message": "keep me",
            },
        ],
    )
    (source / "manifest.json").write_text(
        json.dumps(
            {
                "annotated_variant_context_count": 2,
                "failure_count": 2,
                "gnomad_region_count": 1,
                "gnomad_region_success_count": 0,
                "gnomad_region_failure_count": 1,
                "gnomad_raw_variant_count": 0,
                "gnomad_cached_variant_count": 0,
                "gnomad_key_status_counts": {},
            }
        )
    )

    records = [
        {
            "chrom": "1",
            "pos": 101,
            "ref": "C",
            "alt": "T",
            "consequence": "missense_variant",
            "joint": {"an": 100, "ac": [2]},
        },
        {
            "chrom": "1",
            "pos": 102,
            "ref": "C",
            "alt": "CT",
            "consequence": "inframe_insertion",
            "joint": {"an": 100, "ac": [1]},
        },
    ]
    monkeypatch.setattr(
        completion,
        "fetch_region_variants_recursive",
        lambda *_args, **_kwargs: records,
    )

    shared_cache_dir = tmp_path / "gnomad_cache"
    manifest = completion.complete_gnomad_annotation(
        run_dir,
        workers=1,
        gnomad_cache_dir=shared_cache_dir,
    )
    completed = read_tsv(run_dir / "annotation_gnomad_complete" / "variant_annotations.tsv.gz")

    assert [row["gnomad_af"] for row in completed] == ["0.02", "0.01"]
    assert [row["gnomad_csq"] for row in completed] == [
        "missense_variant",
        "inframe_insertion",
    ]
    for before, after in zip(rows, completed):
        for column, value in before.items():
            assert after[column] == value
    assert read_tsv(run_dir / "annotation_gnomad_complete" / "failures.tsv.gz") == [
        {
            "source": "clinvar",
            "scope": "region",
            "chrom": "1",
            "start": "101",
            "end": "102",
            "failure_type": "test",
            "message": "keep me",
        }
    ]
    assert manifest["gnomad_region_success_count"] == 1
    assert manifest["gnomad_region_failure_count"] == 0
    assert manifest["gnomad_completion"]["updated_variant_context_count"] == 2
    assert manifest["gnomad_completion"]["shared_cache"]["enabled"] is True
    assert manifest["gnomad_completion"]["shared_cache"]["tile_write_count"] == 1
    assert not (
        run_dir / "annotation_gnomad_complete" / "variant_strategy_support.tsv.gz"
    ).exists()

    def unexpected_fetch(*_args):
        raise AssertionError("A completed output must not refetch gnomAD regions")

    monkeypatch.setattr(completion, "fetch_region_variants_recursive", unexpected_fetch)
    resumed = completion.complete_gnomad_annotation(
        run_dir,
        workers=1,
        gnomad_cache_dir=shared_cache_dir,
    )

    assert resumed["gnomad_completion"]["attempted_region_count"] == 0
