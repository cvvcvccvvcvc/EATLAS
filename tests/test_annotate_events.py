from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path

import duckdb
import pytest


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
PROJECT_DIR = BIN_DIR.parent
FINALIZE_SCRIPT = BIN_DIR / "finalize_annotation_partitions.py"

from bin.annotate_events import (
    EVENT_VARIANT_MAP_FIELDS,
    VARIANT_ANNOTATION_FIELDS,
    event_variant_map_row,
    load_alignment_manifest,
    path_metadata,
    write_tsv_gz,
    write_variant_annotation_shards,
)
from bin.annotate_vep_partition import SCHEMA as VEP_SHARD_SCHEMA, VEP_FIELDS
from analytics.derivations.support import (  # noqa: E402
    VARIANT_STRATEGY_SUPPORT_FIELDS,
    EventOrthologSupportStream,
    ExactSupportSpool,
    StrategySupport,
    aggregate_exact_support,
    build_variant_strategy_support,
    load_snv_alt_family_support,
)
from analytics.io.artifacts import content_identity  # noqa: E402
from genomics.variants import event_vcf_key, variant_aggregate_key  # noqa: E402
from bin.finalize_annotation_partitions import (
    COUNTER_FIELDS,
    COUNT_FIELDS,
    VARIANT_DATASET_SCHEMA,
    copy_event_variant_map_dataset,
    copy_variant_annotation_dataset,
    event_variant_map_manifest,
    merge_gnomad_shared_cache,
    merge_partition_timings,
    validate_partition_manifests,
)


def canonical_partition_manifest(partition_id: str) -> dict:
    return {
        "partition_id": partition_id,
        "stage": "annotation",
        "schema": "normalized_annotation_evidence_partition_v2",
        **{field: 0 for field in COUNT_FIELDS},
        **{field: {} for field in COUNTER_FIELDS},
        "failure_count": 0,
        "clinvar_vcf": {"path": "clinvar.vcf.gz", "size_bytes": 1, "mtime": 1},
        "clinvar_tbi": {"path": "clinvar.vcf.gz.tbi", "size_bytes": 1, "mtime": 1},
        "gnomad_api_url": "https://gnomad.example/api",
        "gnomad_dataset": "gnomad_r4",
    }


def canonical_alignment_manifest(partition_id: str = "") -> dict:
    return {
        "stage": "alignment",
        "partition_id": partition_id,
        "schema": "normalized_alignment_evidence_partition_v2",
        "alignment_event_mode": "compact_support",
        "event_ortholog_support_format": "event_group_id_v2",
        "alignment_event_count": 0,
        "event_ortholog_support_count": 0,
    }


def test_alignment_manifest_is_required_and_partition_bound(tmp_path: Path) -> None:
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(canonical_alignment_manifest("partition_000001")) + "\n")

    assert load_alignment_manifest(path, "partition_000001")["alignment_event_count"] == 0
    with pytest.raises(ValueError, match="partition mismatch"):
        load_alignment_manifest(path, "partition_000002")


def test_variant_annotation_schema_is_analysis_ready() -> None:
    assert VARIANT_ANNOTATION_FIELDS == [
        "variant_key",
        "gene_id",
        "event_type",
        "ref",
        "alt",
        "lookup_status",
        "strategies",
        "clinvar_sig",
        "clinvar_revstat",
        "clinvar_review_stars",
        "clinvar_review_stars_status",
        "clinvar_id",
        "clinvar_allele_id",
        "clinvar_scv_count",
        "clinvar_hgvs",
        "clinvar_disease",
        "clinvar_variant_type",
        "gnomad_af",
        "gnomad_af_source",
        "gnomad_csq",
    ]


@pytest.mark.parametrize(
    "entrypoint",
    [
        "annotate_events.py",
        "annotate_vep_partition.py",
        "finalize_annotation_partitions.py",
    ],
)
def test_annotation_entrypoints_accept_large_tsv_fields(
    tmp_path: Path,
    entrypoint: str,
) -> None:
    source = tmp_path / "large_field.tsv.gz"
    large_field = "A" * 165_969
    with gzip.open(source, "wt", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerow([large_field])

    probe = (
        "import csv,gzip,runpy,sys;"
        "runpy.run_path(sys.argv[1],run_name='csv_limit_probe');"
        "handle=gzip.open(sys.argv[2],'rt',newline='');"
        "print(len(next(csv.reader(handle,delimiter=chr(9)))[0]))"
    )
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            probe,
            str(BIN_DIR / entrypoint),
            str(source),
        ],
        cwd=BIN_DIR.parent,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(len(large_field))


def test_partitioned_manifest_keeps_durable_lineage_counts_only() -> None:
    assert "excluded_non_concrete_event_count" not in COUNT_FIELDS
    assert "event_variant_map_count" in COUNT_FIELDS


def write_event_variant_map(
    path: Path,
    rows: list[dict[str, object]],
    fields: list[str] | None = None,
) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields or EVENT_VARIANT_MAP_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_vep_shard(
    root: Path,
    *,
    partition_id: str,
    shard: dict[str, object],
    source_fields: list[str],
    rows: list[dict[str, object]],
) -> Path:
    shard_id = str(shard["shard_id"])
    directory = root / f"vep_{partition_id}_{shard_id}"
    directory.mkdir(parents=True)
    output = directory / "variant_annotations.tsv.gz"
    output_fields = [*source_fields, *VEP_FIELDS]
    write_tsv_gz(output, output_fields, rows)
    status_counts = {}
    for row in rows:
        status = str(row.get("vep_status") or "")
        status_counts[status] = status_counts.get(status, 0) + 1
    manifest = {
        "stage": "annotation",
        "schema": VEP_SHARD_SCHEMA,
        "partition_id": partition_id,
        "shard_id": shard_id,
        "row_count": len(rows),
        "status_counts": status_counts,
        "input": {
            "name": str(shard["path"]),
            "size_bytes": int(shard["size_bytes"]),
            "fields": source_fields,
        },
        "output": {
            "name": output.name,
            "size_bytes": output.stat().st_size,
            "fields": output_fields,
        },
        "config": {"backend": "local", "release": "116"},
    }
    (directory / "manifest.json").write_text(json.dumps(manifest) + "\n")
    return directory


def test_event_variant_map_preserves_collapsed_and_non_concrete_lineage() -> None:
    contexts = {
        "1": {
            "chrom": "1",
            "begin": 100,
            "end": 104,
            "seq": "AAAAA",
        }
    }
    raw_deletions = [
        {
            "gene_id": "1",
            "event_type": "del",
            "target_start0": target_start0,
            "genomic_accession": "NC_000001.11",
            "genomic_start1": str(100 + target_start0),
            "ref": "A",
            "alt": "",
        }
        for target_start0 in [2, 3]
    ]
    normalized = [event_vcf_key(row, contexts) for row in raw_deletions]
    assert normalized == [
        (("1", 100, "AA", "A"), "ok"),
        (("1", 100, "AA", "A"), "ok"),
    ]

    rows = [
        event_variant_map_row(1, *normalized[0]),
        event_variant_map_row(2, *normalized[1]),
        event_variant_map_row(3, None, "non_concrete_allele"),
    ]
    assert rows == [
        {
            "event_group_id": 1,
            "variant_key": "1:100:AA>A",
            "normalization_status": "ok",
        },
        {
            "event_group_id": 2,
            "variant_key": "1:100:AA>A",
            "normalization_status": "ok",
        },
        {
            "event_group_id": 3,
            "variant_key": "",
            "normalization_status": "non_concrete_allele",
        },
    ]


def test_event_variant_map_manifest_declares_exact_partitioned_schema() -> None:
    assert event_variant_map_manifest(3, 2) == {
        "layout": "partitioned",
        "format": "tsv_gzip_v1",
        "path": "event_variant_map/partitions",
        "partition_count": 2,
        "row_count": 3,
        "fields": [
            "event_group_id",
            "variant_key",
            "normalization_status",
        ],
        "event_group_id_scope": "partition",
    }


def test_event_variant_map_dataset_is_validated_and_copied_byte_for_byte(
    tmp_path: Path,
) -> None:
    partitions = []
    source_bytes = {}
    partition_rows = [
        [
            event_variant_map_row(1, ("1", 100, "AA", "A"), "ok"),
            event_variant_map_row(2, ("1", 100, "AA", "A"), "ok"),
        ],
        [event_variant_map_row(1, None, "non_concrete_allele")],
    ]
    for index, rows in enumerate(partition_rows, start=1):
        partition_id = f"partition_{index:06d}"
        partition = tmp_path / f"annotation_{partition_id}"
        partition.mkdir()
        source = partition / "event_variant_map.tsv.gz"
        write_event_variant_map(source, rows)
        source_bytes[partition_id] = source.read_bytes()
        manifest = canonical_partition_manifest(partition_id)
        manifest["event_row_count"] = len(rows)
        manifest["event_variant_map_count"] = len(rows)
        partitions.append((partition, manifest))

    output = tmp_path / "event_variant_map"
    row_count = copy_event_variant_map_dataset(partitions, output)

    assert row_count == 3
    for partition_id, expected_bytes in source_bytes.items():
        copied = output / "partitions" / partition_id / "event_variant_map.tsv.gz"
        assert copied.read_bytes() == expected_bytes


@pytest.mark.parametrize(
    ("rows", "map_count", "event_count", "error"),
    [
        pytest.param(
            [event_variant_map_row(2, ("1", 100, "A", "G"), "ok")],
            1,
            1,
            "event_group_id values must be consecutive",
            id="foreign-key-gap",
        ),
        pytest.param(
            [event_variant_map_row(1, ("1", 100, "A", "G"), "ok")],
            2,
            1,
            "row count does not match partition manifest",
            id="manifest-count",
        ),
        pytest.param(
            [event_variant_map_row(1, ("1", 100, "", "G"), "ok")],
            1,
            1,
            "invalid canonical variant_key",
            id="invalid-canonical-key",
        ),
        pytest.param(
            [event_variant_map_row(1, ("1", 100, "A", "G"), "unsupported_allele")],
            1,
            1,
            "Unresolved event has a canonical variant_key",
            id="unresolved-key",
        ),
    ],
)
def test_event_variant_map_rejects_invalid_count_or_event_fk(
    tmp_path: Path,
    rows: list[dict[str, object]],
    map_count: int,
    event_count: int,
    error: str,
) -> None:
    partition_id = "partition_000001"
    partition = tmp_path / f"annotation_{partition_id}"
    partition.mkdir()
    write_event_variant_map(partition / "event_variant_map.tsv.gz", rows)
    manifest = canonical_partition_manifest(partition_id)
    manifest["event_variant_map_count"] = map_count
    manifest["event_row_count"] = event_count

    with pytest.raises(ValueError, match=error):
        copy_event_variant_map_dataset(
            [(partition, manifest)],
            tmp_path / "event_variant_map",
        )


def test_event_variant_map_rejects_noncanonical_schema(tmp_path: Path) -> None:
    partition_id = "partition_000001"
    partition = tmp_path / f"annotation_{partition_id}"
    partition.mkdir()
    write_event_variant_map(
        partition / "event_variant_map.tsv.gz",
        [event_variant_map_row(1, ("1", 100, "A", "G"), "ok")],
        fields=list(reversed(EVENT_VARIANT_MAP_FIELDS)),
    )
    manifest = canonical_partition_manifest(partition_id)
    manifest["event_variant_map_count"] = 1
    manifest["event_row_count"] = 1

    with pytest.raises(ValueError, match="Unexpected event-variant map fields"):
        copy_event_variant_map_dataset(
            [(partition, manifest)],
            tmp_path / "event_variant_map",
        )


def test_partition_manifest_validation_requires_current_contract(tmp_path: Path) -> None:
    partition = tmp_path / "annotation_partition_000001"
    manifest = canonical_partition_manifest("partition_000001")

    validate_partition_manifests([(partition, manifest)])

    del manifest["event_row_count"]
    with pytest.raises(ValueError, match="missing event_row_count"):
        validate_partition_manifests([(partition, manifest)])


def test_finalizer_publishes_only_source_annotation_evidence(tmp_path: Path) -> None:
    partition_root = tmp_path / "partitions"
    partition = partition_root / "annotation_partition_000001"
    partition.mkdir(parents=True)
    clinvar_vcf = tmp_path / "clinvar.vcf.gz"
    clinvar_tbi = tmp_path / "clinvar.vcf.gz.tbi"
    clinvar_vcf.write_bytes(b"clinvar")
    clinvar_tbi.write_bytes(b"index")
    manifest = canonical_partition_manifest("partition_000001")
    manifest.update(
        {
            "event_row_count": 1,
            "event_variant_map_count": 1,
            "variant_context_count": 1,
            "annotated_variant_context_count": 1,
            "clinvar_vcf": path_metadata(clinvar_vcf),
            "clinvar_tbi": path_metadata(clinvar_tbi),
        }
    )
    source_rows = [
        {"variant_key": "1:100:A>G", "gene_id": "1", "lookup_status": "ok"}
    ]
    source_dataset = write_variant_annotation_shards(
        partition / "variant_annotation_shards",
        source_rows,
    )
    manifest["variant_annotations"] = source_dataset
    (partition / "manifest.json").write_text(json.dumps(manifest) + "\n")
    write_event_variant_map(
        partition / "event_variant_map.tsv.gz",
        [event_variant_map_row(1, ("1", 100, "A", "G"), "ok")],
    )
    with gzip.open(partition / "failures.tsv.gz", "wt", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerow(
            ["source", "scope", "chrom", "start", "end", "failure_type", "message"]
        )
    vep_root = tmp_path / "vep"
    vep_directory = write_vep_shard(
        vep_root,
        partition_id="partition_000001",
        shard=source_dataset["shards"][0],
        source_fields=VARIANT_ANNOTATION_FIELDS,
        rows=[
            {
                **source_rows[0],
                "vep_status": "ok",
                "vep_primary_consequence": "intron_variant",
            }
        ],
    )
    outdir = tmp_path / "annotation"

    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "bin.finalize_annotation_partitions",
            "--partition-root",
            str(partition_root),
            "--vep-root",
            str(vep_root),
            "--clinvar-vcf",
            str(clinvar_vcf),
            "--clinvar-tbi",
            str(clinvar_tbi),
            "--outdir",
            str(outdir),
        ],
        cwd=PROJECT_DIR,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert {path.name for path in outdir.iterdir()} == {
        "variant_annotations",
        "event_variant_map",
        "failures.tsv.gz",
        "manifest.json",
    }
    final_manifest = json.loads((outdir / "manifest.json").read_text())
    assert final_manifest["stage"] == "annotation"
    assert final_manifest["schema"] == "normalized_annotation_evidence_v4"
    assert final_manifest["clinvar_vcf"] == content_identity(clinvar_vcf)
    assert final_manifest["clinvar_tbi"] == content_identity(clinvar_tbi)
    dataset_manifest = json.loads(
        (outdir / "variant_annotations" / "manifest.json").read_text()
    )
    assert dataset_manifest["schema"] == VARIANT_DATASET_SCHEMA
    assert dataset_manifest["row_count"] == 1
    assert dataset_manifest["vep_status_counts"] == {"ok": 1}
    durable_shard = (
        outdir
        / "variant_annotations"
        / dataset_manifest["partitions"][0]["shards"][0]["path"]
    )
    assert durable_shard.read_bytes() == (
        vep_directory / "variant_annotations.tsv.gz"
    ).read_bytes()
    forbidden_keys = {
        "variant_strategy_support_count",
        "variant_ortholog_support_count",
        "variant_ortholog_support_format",
        "ortholog_evidence_summary_count",
    }
    assert forbidden_keys.isdisjoint(final_manifest)


def test_partition_timings_are_preserved_and_summed(tmp_path: Path) -> None:
    partitions = [
        (
            tmp_path / "partition_000001",
            {
                "partition_id": "partition_000001",
                "timings_seconds": {"collapse_events": 1.25, "gnomad_lookup": 0.5},
            },
        ),
        (
            tmp_path / "partition_000002",
            {
                "partition_id": "partition_000002",
                "timings_seconds": {"collapse_events": 2.75},
            },
        ),
    ]

    by_partition = merge_partition_timings(partitions)

    assert by_partition == {
        "partition_000001": {"collapse_events": 1.25, "gnomad_lookup": 0.5},
        "partition_000002": {"collapse_events": 2.75},
    }


def test_analytics_support_schema_is_not_part_of_stage3_manifest() -> None:
    assert VARIANT_STRATEGY_SUPPORT_FIELDS[-1] == "site_aligned_ortholog_count"
    assert "alt_support_family_count" in VARIANT_STRATEGY_SUPPORT_FIELDS
    assert "variant_strategy_site_depth_count" not in COUNT_FIELDS


def test_exact_support_schema_is_not_part_of_stage3_manifest() -> None:
    assert "variant_ortholog_support_count" not in COUNT_FIELDS


def test_partitioned_manifest_aggregates_shared_gnomad_cache_metrics(tmp_path: Path) -> None:
    identity = {
        "enabled": True,
        "directory": "/cache/gnomad",
        "schema_version": 1,
        "dataset": "gnomad_r4",
        "reference_genome": "GRCh38",
        "tile_size_bp": 25_000,
    }
    partitions = [
        (
            tmp_path / "one",
            {
                "gnomad_shared_cache": {
                    **identity,
                    "tile_hit_count": 2,
                    "tile_miss_count": 3,
                    "tile_write_count": 3,
                    "corrupt_tile_count": 0,
                    "fetch_batch_count": 1,
                    "split_count": 0,
                }
            },
        ),
        (
            tmp_path / "two",
            {
                "gnomad_shared_cache": {
                    **identity,
                    "tile_hit_count": 5,
                    "tile_miss_count": 0,
                    "tile_write_count": 0,
                    "corrupt_tile_count": 0,
                    "fetch_batch_count": 0,
                    "split_count": 0,
                }
            },
        ),
    ]

    merged = merge_gnomad_shared_cache(partitions)

    assert merged is not None
    assert merged["tile_hit_count"] == 7
    assert merged["tile_miss_count"] == 3
    assert merged["tile_write_count"] == 3
    assert merged["fetch_batch_count"] == 1


def test_exact_support_spool_counts_distinct_orthologs(tmp_path: Path) -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "_support_by_strategy": {},
    }
    aggregates_by_id = [None, aggregate]
    spool = ExactSupportSpool(tmp_path / "exact_support.tsv")
    for strategy, ortholog in [
        ("s1", "101"),
        ("s1", "101"),
        ("s1", "102"),
        ("s2", "101"),
    ]:
        spool.add_group(
            variant_context_id=1,
            gene_id="1",
            strategy=strategy,
            support_rows=[
                {"ortholog_gene_id": ortholog, "support_row_count": "1"}
            ],
        )
    with duckdb.connect() as connection:
        exact_edge_count = aggregate_exact_support(
            connection,
            spool,
            aggregates_by_id,
        )

    rows, missing_key_count = build_variant_strategy_support([aggregate])

    assert exact_edge_count == 3
    assert not spool.path.exists()
    assert missing_key_count == 0
    assert rows == [
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "strategy": "s1",
            "alt_support_row_count": 3,
            "alt_support_ortholog_count": 2,
            "alt_support_family_count": "",
            "site_aligned_ortholog_count": "",
        },
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "strategy": "s2",
            "alt_support_row_count": 1,
            "alt_support_ortholog_count": 1,
            "alt_support_family_count": "",
            "site_aligned_ortholog_count": "",
        },
    ]


def test_variant_strategy_support_loads_exact_alt_family_count(tmp_path: Path) -> None:
    path = tmp_path / "snv_alt_taxonomic_support.tsv.gz"
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "gene_id",
                "strategy",
                "target_start0",
                "ref",
                "alt",
                "known_family_count",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "gene_id": "1",
                "strategy": "s1",
                "target_start0": 4,
                "ref": "a",
                "alt": "g",
                "known_family_count": 2,
            }
        )
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "event_type": "snv",
        "target_start0": 4,
        "ref": "A",
        "alt": "G",
        "_support_by_strategy": {
            "s1": StrategySupport(row_count=3, ortholog_count=3)
        },
    }

    family_supports = load_snv_alt_family_support(path)
    rows, _missing = build_variant_strategy_support(
        [aggregate],
        {("1", "s1", 4): 3},
        family_supports,
    )

    assert rows[0]["alt_support_family_count"] == 2


def test_event_ortholog_support_stream_reads_consecutive_compact_groups(
    tmp_path: Path,
) -> None:
    support_tsv = tmp_path / "event_ortholog_support.tsv.gz"
    write_tsv_gz(
        support_tsv,
        [
            "event_group_id",
            "ortholog_gene_id",
            "tax_id",
            "mapq",
            "native_alignment_type",
            "support_row_count",
        ],
        [
            {
                "event_group_id": "1",
                "ortholog_gene_id": "101",
                "tax_id": "10090",
                "mapq": "60",
                "native_alignment_type": "P",
                "support_row_count": "2",
            },
            {
                "event_group_id": "2",
                "ortholog_gene_id": "201",
                "tax_id": "10116",
                "mapq": "10",
                "native_alignment_type": "supplementary",
                "support_row_count": "1",
            },
        ],
    )

    with EventOrthologSupportStream(support_tsv) as stream:
        assert [row["ortholog_gene_id"] for row in stream.take(1)] == ["101"]
        assert [row["ortholog_gene_id"] for row in stream.take(2)] == ["201"]
        stream.finish()


def test_event_ortholog_support_stream_rejects_unmatched_group(tmp_path: Path) -> None:
    support_tsv = tmp_path / "event_ortholog_support.tsv.gz"
    write_tsv_gz(
        support_tsv,
        [
            "event_group_id",
            "ortholog_gene_id",
            "tax_id",
            "mapq",
            "native_alignment_type",
            "support_row_count",
        ],
        [
            {
                "event_group_id": "2",
                "ortholog_gene_id": "201",
                "tax_id": "10116",
                "mapq": "",
                "native_alignment_type": "",
                "support_row_count": "1",
            }
        ],
    )

    with EventOrthologSupportStream(support_tsv) as stream:
        assert stream.take(1) == []
        with pytest.raises(ValueError, match="no matching compact event"):
            stream.finish()


def test_snv_support_uses_site_aligned_depth() -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "event_type": "snv",
        "target_start0": "9",
        "_support_by_strategy": {
            "s1": StrategySupport(row_count=1, ortholog_count=1)
        },
    }
    rows, _missing_key_count = build_variant_strategy_support(
        [aggregate],
        {("1", "s1", 9): 4},
    )
    assert rows[0]["site_aligned_ortholog_count"] == 4


def test_snv_support_rejects_alt_count_above_site_depth() -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "event_type": "snv",
        "_support_by_strategy": {
            "s1": StrategySupport(row_count=2, ortholog_count=2)
        },
    }

    with pytest.raises(ValueError, match="exceeds site-aligned"):
        build_variant_strategy_support(
            [aggregate],
            {("1", "s1", 0): 1},
        )


def test_variant_strategy_support_uses_exact_edge_multiplicity() -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "_support_by_strategy": {
            "s1": StrategySupport(row_count=10, ortholog_count=3)
        },
    }

    rows, _missing_key_count = build_variant_strategy_support([aggregate])

    assert rows[0]["alt_support_row_count"] == 10
    assert rows[0]["alt_support_ortholog_count"] == 3


def test_exact_support_spool_requires_one_strategy(tmp_path: Path) -> None:
    with ExactSupportSpool(tmp_path / "exact_support.tsv") as spool:
        with pytest.raises(ValueError, match="requires one alignment strategy"):
            spool.add_group(
                variant_context_id=1,
                gene_id="1",
                strategy="",
                support_rows=[
                    {"ortholog_gene_id": "101", "support_row_count": "1"}
                ],
            )


def test_canonical_variant_key_collapses_raw_indel_representations() -> None:
    left = {
        "gene_id": "1",
        "event_type": "del",
        "target_start0": "10",
        "target_end0": "11",
        "ref": "A",
        "alt": "",
    }
    right = {
        "gene_id": "1",
        "event_type": "del",
        "target_start0": "11",
        "target_end0": "12",
        "ref": "A",
        "alt": "",
    }

    assert variant_aggregate_key(left, "1:100:AA>A") == variant_aggregate_key(
        right,
        "1:100:AA>A",
    )
    assert variant_aggregate_key(left, "") != variant_aggregate_key(right, "")


@pytest.mark.parametrize(
    ("ref", "alt"),
    [
        ("G", "Y"),
        ("T", "R"),
        ("A", "N"),
        ("", "."),
    ],
)
def test_event_vcf_key_rejects_non_concrete_alleles(ref: str, alt: str) -> None:
    key, status = event_vcf_key(
        {
            "gene_id": "1",
            "event_type": "snv",
            "genomic_accession": "NC_000001.11",
            "genomic_start1": "100",
            "target_start0": "0",
            "ref": ref,
            "alt": alt,
        },
        {},
    )

    assert key is None
    assert status == "non_concrete_allele"


def test_event_vcf_key_requires_target_context() -> None:
    key, status = event_vcf_key(
        {
            "gene_id": "1",
            "event_type": "snv",
            "genomic_accession": "NC_000001.11",
            "genomic_start1": "100",
            "target_start0": "0",
            "ref": "A",
            "alt": "G",
        },
        {},
    )

    assert key is None
    assert status == "raw_no_context"


@pytest.mark.parametrize(
    ("event_type", "ref", "alt"),
    [
        ("ins", "", "A"),
        ("del", "A", ""),
    ],
)
def test_event_vcf_key_leaves_ambiguous_anchor_unresolved(
    event_type: str,
    ref: str,
    alt: str,
) -> None:
    key, status = event_vcf_key(
        {
            "gene_id": "1",
            "event_type": event_type,
            "genomic_accession": "NC_000001.11",
            "genomic_start1": "102",
            "target_start0": "2",
            "ref": ref,
            "alt": alt,
        },
        {
            "1": {
                "chrom": "1",
                "begin": 100,
                "end": 104,
                "seq": "ANAAA",
            }
        },
    )

    assert key is None
    assert status == "unsupported_allele"
