from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
PROJECT_DIR = BIN_DIR.parent
FINALIZE_SCRIPT = BIN_DIR / "finalize_annotation_partitions.py"
sys.path.insert(0, str(BIN_DIR))

from annotate_events import (  # noqa: E402
    EVENT_VARIANT_MAP_FIELDS,
    EventOrthologSupportStream,
    PARTITION_TSV_SHARD_FORMAT,
    VARIANT_ANNOTATION_FIELDS,
    VARIANT_STRATEGY_SUPPORT_FIELDS,
    add_strategy_support,
    build_variant_strategy_support,
    event_vcf_key,
    event_variant_map_row,
    iter_variant_strategy_snv_sites,
    load_snv_alt_genus_support,
    load_alignment_manifest,
    resolve_target_feature_paths,
    variant_aggregate_key,
    write_tsv_gz,
)
from finalize_annotation_partitions import (  # noqa: E402
    COUNTER_FIELDS,
    COUNT_FIELDS,
    concatenate_tsv_gz_members,
    copy_event_variant_map_dataset,
    event_variant_map_manifest,
    merge_gnomad_shared_cache,
    merge_partition_timings,
    validate_partition_manifests,
)


def canonical_partition_manifest(partition_id: str) -> dict:
    return {
        "partition_id": partition_id,
        "stage": "annotation",
        "schema": "normalized_annotation_evidence_partition_v1",
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
        "schema": "normalized_alignment_evidence_partition_v1",
        "alignment_event_mode": "compact_support",
        "event_ortholog_support_format": "event_group_id_v1",
        "alignment_event_count": 0,
        "event_ortholog_support_count": 0,
    }


def test_alignment_manifest_is_required_and_partition_bound(tmp_path: Path) -> None:
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(canonical_alignment_manifest("partition_000001")) + "\n")

    assert load_alignment_manifest(path, "partition_000001")["alignment_event_count"] == 0
    with pytest.raises(ValueError, match="partition mismatch"):
        load_alignment_manifest(path, "partition_000002")


def test_target_features_accept_one_canonical_table_or_partition_directory(
    tmp_path: Path,
) -> None:
    table = tmp_path / "target_features.tsv.gz"
    table.write_bytes(b"")
    assert resolve_target_feature_paths(table) == [table]

    directory = tmp_path / "target_features"
    directory.mkdir()
    partition = directory / "1.tsv.gz"
    partition.write_bytes(b"")
    assert resolve_target_feature_paths(directory) == [partition]


def test_variant_annotation_schema_is_analysis_ready() -> None:
    assert VARIANT_ANNOTATION_FIELDS == [
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
    ["annotate_events.py", "finalize_annotation_partitions.py"],
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
        "sys.path.insert(0,sys.argv[3]);"
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
            str(BIN_DIR),
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
    manifest = canonical_partition_manifest("partition_000001")
    manifest.update(
        {
            "event_row_count": 1,
            "event_variant_map_count": 1,
            "variant_context_count": 1,
            "annotated_variant_context_count": 1,
            "partition_tsv_shard_format": PARTITION_TSV_SHARD_FORMAT,
            "partition_tsv_shard_fields": {
                "variant_annotations.tsv.gz": VARIANT_ANNOTATION_FIELDS,
            },
        }
    )
    (partition / "manifest.json").write_text(json.dumps(manifest) + "\n")
    write_tsv_gz(
        partition / "variant_annotations.tsv.gz",
        VARIANT_ANNOTATION_FIELDS,
        [{"variant_key": "1:100:A>G", "gene_id": "1", "lookup_status": "ok"}],
        include_header=False,
    )
    write_event_variant_map(
        partition / "event_variant_map.tsv.gz",
        [event_variant_map_row(1, ("1", 100, "A", "G"), "ok")],
    )
    with gzip.open(partition / "failures.tsv.gz", "wt", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerow(
            ["source", "scope", "chrom", "start", "end", "failure_type", "message"]
        )
    outdir = tmp_path / "annotation"

    completed = subprocess.run(
        [
            sys.executable,
            str(FINALIZE_SCRIPT),
            "--partition-root",
            str(partition_root),
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
        "variant_annotations.tsv.gz",
        "event_variant_map",
        "failures.tsv.gz",
        "manifest.json",
    }
    final_manifest = json.loads((outdir / "manifest.json").read_text())
    assert final_manifest["stage"] == "annotation"
    assert final_manifest["schema"] == "normalized_annotation_evidence_v1"
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


def test_partition_tsv_members_are_concatenated_without_recompression(
    tmp_path: Path,
) -> None:
    filename = "large_table.tsv.gz"
    fields = ["variant_key", "support_count"]
    partition_rows = [
        [{"variant_key": "1:1:A>G", "support_count": 2}],
        [],
        [{"variant_key": "2:2:C>T", "support_count": 3}],
    ]
    partitions = []
    source_bytes = []
    for index, rows in enumerate(partition_rows, start=1):
        partition = tmp_path / f"partition_{index:06d}"
        partition.mkdir()
        source = partition / filename
        write_tsv_gz(source, fields, rows, include_header=False)
        source_bytes.append(source.read_bytes())
        partitions.append(
            (
                partition,
                {
                    "partition_tsv_shard_format": PARTITION_TSV_SHARD_FORMAT,
                    "partition_tsv_shard_fields": {filename: fields},
                    "row_count": len(rows),
                },
            )
        )
    output = tmp_path / "merged.tsv.gz"

    row_count = concatenate_tsv_gz_members(
        partitions,
        filename,
        "row_count",
        output,
    )

    with gzip.open(output, "rt", newline="") as handle:
        assert list(csv.reader(handle, delimiter="\t")) == [
            fields,
            ["1:1:A>G", "2"],
            ["2:2:C>T", "3"],
        ]
    output_bytes = output.read_bytes()
    member_positions = [output_bytes.find(member) for member in source_bytes]
    assert row_count == 2
    assert all(position >= 0 for position in member_positions)
    assert member_positions == sorted(member_positions)
    assert output_bytes.endswith(source_bytes[-1])


def test_partition_tsv_member_rejects_embedded_header(tmp_path: Path) -> None:
    filename = "large_table.tsv.gz"
    fields = ["variant_key", "support_count"]
    partition = tmp_path / "partition_000001"
    partition.mkdir()
    write_tsv_gz(
        partition / filename,
        fields,
        [{"variant_key": "1:1:A>G", "support_count": 2}],
    )

    with pytest.raises(ValueError, match="unexpectedly contains a header"):
        concatenate_tsv_gz_members(
            [
                (
                    partition,
                    {
                        "partition_tsv_shard_format": PARTITION_TSV_SHARD_FORMAT,
                        "partition_tsv_shard_fields": {filename: fields},
                        "row_count": 1,
                    },
                )
            ],
            filename,
            "row_count",
            tmp_path / "merged.tsv.gz",
        )


def test_analytics_support_schema_is_not_part_of_stage3_manifest() -> None:
    assert VARIANT_STRATEGY_SUPPORT_FIELDS[-1] == "site_aligned_ortholog_count"
    assert "alt_support_genus_count" in VARIANT_STRATEGY_SUPPORT_FIELDS
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


def test_variant_strategy_support_counts_distinct_orthologs() -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "_support_by_strategy": {},
    }
    for strategy, ortholog in [
        ("s1", "101"),
        ("s1", "101"),
        ("s1", "102"),
        ("s2", "101"),
    ]:
        add_strategy_support(
            aggregate,
            {
                "strategy": strategy,
                "ortholog_gene_id": ortholog,
            },
        )

    rows, missing_key_count = build_variant_strategy_support([aggregate])

    assert missing_key_count == 0
    assert rows == [
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "strategy": "s1",
            "alt_support_row_count": 3,
            "alt_support_ortholog_count": 2,
            "alt_support_genus_count": "",
            "site_aligned_ortholog_count": "",
        },
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "strategy": "s2",
            "alt_support_row_count": 1,
            "alt_support_ortholog_count": 1,
            "alt_support_genus_count": "",
            "site_aligned_ortholog_count": "",
        },
    ]


def test_variant_strategy_support_loads_exact_alt_genus_count(tmp_path: Path) -> None:
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
                "all__genus",
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
                "all__genus": 2,
            }
        )
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "event_type": "snv",
        "target_start0": 4,
        "ref": "A",
        "alt": "G",
        "_support_by_strategy": {},
    }
    for ortholog in ("101", "102", "103"):
        add_strategy_support(
            aggregate,
            {"strategy": "s1", "ortholog_gene_id": ortholog},
        )

    genus_supports = load_snv_alt_genus_support(path)
    rows, _missing = build_variant_strategy_support(
        [aggregate],
        {("1", "s1", 4): 3},
        genus_supports,
    )

    assert rows[0]["alt_support_genus_count"] == 2


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
            "taxname",
            "mapq",
            "native_alignment_type",
            "support_row_count",
        ],
        [
            {
                "event_group_id": "1",
                "ortholog_gene_id": "101",
                "tax_id": "10090",
                "taxname": "Mus musculus",
                "mapq": "60",
                "native_alignment_type": "P",
                "support_row_count": "2",
            },
            {
                "event_group_id": "2",
                "ortholog_gene_id": "201",
                "tax_id": "10116",
                "taxname": "Rattus norvegicus",
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
            "taxname",
            "mapq",
            "native_alignment_type",
            "support_row_count",
        ],
        [
            {
                "event_group_id": "2",
                "ortholog_gene_id": "201",
                "tax_id": "10116",
                "taxname": "Rattus norvegicus",
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
        "_support_by_strategy": {},
    }
    add_strategy_support(
        aggregate,
        {"strategy": "s1", "ortholog_gene_id": "101"},
    )

    assert list(iter_variant_strategy_snv_sites([aggregate])) == [
        {
            "gene_id": "1",
            "strategy": "s1",
            "target_start0": "9",
        }
    ]
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
        "_support_by_strategy": {},
    }
    add_strategy_support(
        aggregate,
        {"strategy": "s1", "support_ortholog_count": "2"},
    )

    with pytest.raises(ValueError, match="exceeds site-aligned"):
        build_variant_strategy_support(
            [aggregate],
            {("1", "s1", 0): 1},
        )


def test_variant_strategy_support_accepts_single_strategy_compact_counts() -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "_support_by_strategy": {},
    }
    add_strategy_support(
        aggregate,
        {
            "strategy": "s1",
            "support_row_count": "10",
            "support_ortholog_count": "3",
        },
    )

    rows, _missing_key_count = build_variant_strategy_support([aggregate])

    assert rows[0]["alt_support_row_count"] == 10
    assert rows[0]["alt_support_ortholog_count"] == 3


def test_variant_strategy_support_requires_one_strategy() -> None:
    aggregate = {
        "variant_key": "1:100:A>G",
        "gene_id": "1",
        "_support_by_strategy": {},
    }

    with pytest.raises(ValueError, match="requires one alignment strategy"):
        add_strategy_support(
            aggregate,
            {
                "strategies": "s1,s2",
                "support_row_count": "10",
                "support_ortholog_count": "3",
            },
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


def test_event_vcf_key_keeps_concrete_alleles() -> None:
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

    assert key == ("1", 100, "A", "G")
    assert status == "raw_no_context"
