from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "bin"))

from fetch_taxonomy import fetch_taxonomy_records, taxonomy_row
from finalize_annotation_partitions import merge_ortholog_evidence
from merge_alignment_results import write_compact_events
from ortholog_evidence_summary import write_ortholog_evidence_summary
from taxonomic_evidence import (
    COUNT_KEYS,
    build_taxonomy_summary_rows,
    count_member_groups,
    load_taxonomy_profiles,
)


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def taxonomy_fixture(tmp_path: Path) -> Path:
    path = tmp_path / "taxonomy.tsv.gz"
    fields = ["tax_id", "species_id", "genus_id", "family_id", "order_id", "parent_tax_ids"]
    write_tsv_gz(
        path,
        fields,
        [
            {
                "tax_id": "9598",
                "species_id": "9598",
                "genus_id": "9596",
                "family_id": "9604",
                "order_id": "9443",
                "parent_tax_ids": "2759,33208,7742,32523,32524,40674,9443,9598",
            },
            {
                "tax_id": "10090",
                "species_id": "10090",
                "genus_id": "10088",
                "family_id": "10066",
                "order_id": "9989",
                "parent_tax_ids": "2759,33208,7742,32523,32524,40674,10090",
            },
        ],
    )
    return path


def test_taxonomy_row_reads_ncbi_taxonomy_lineage() -> None:
    row = taxonomy_row(
        "9598",
        {
            "tax_id": 9598,
            "rank": "SPECIES",
            "current_scientific_name": {"name": "Pan troglodytes"},
            "group_name": "primates",
            "classification": {
                "class": {"id": 40674, "name": "Mammalia"},
                "order": {"id": 9443, "name": "Primates"},
                "family": {"id": 9604, "name": "Hominidae"},
                "genus": {"id": 9596, "name": "Pan"},
                "species": {"id": 9598, "name": "Pan troglodytes"},
            },
            "parents": [9443, 40674, 7742, 33208, 2759],
        },
    )

    assert row["scientific_name"] == "Pan troglodytes"
    assert row["species_id"] == "9598"
    assert row["genus_id"] == "9596"
    assert row["is_primate"] == "true"
    assert row["is_mammal"] == "true"
    assert row["is_vertebrate"] == "true"


def test_taxonomy_batch_request_does_not_use_single_taxon_parents_flag(
    tmp_path: Path,
    monkeypatch,
) -> None:
    def fake_run(command, *, text, stdout, stderr):
        assert command[:4] == ["datasets", "summary", "taxonomy", "taxon"]
        assert "--inputfile" in command
        assert "--parents" not in command
        assert text is True
        assert stderr is not None
        stdout.write(json.dumps({"taxonomy": {"tax_id": 9598}}) + "\n")
        return type("Result", (), {"returncode": 0, "stderr": ""})()

    monkeypatch.setattr("fetch_taxonomy.subprocess.run", fake_run)

    records = fetch_taxonomy_records(["9598"], "datasets", tmp_path)

    assert records == {"9598": {"tax_id": 9598}}


def test_scope_and_unit_counts_use_any_member_semantics(tmp_path: Path) -> None:
    profiles = load_taxonomy_profiles(taxonomy_fixture(tmp_path))
    members = [("chimp_gene", "9598"), ("mouse_gene", "10090")]
    counts = count_member_groups(members, profiles)

    assert counts["all__ortholog"] == 2
    assert counts["all__species"] == 2
    assert counts["mammalia__family"] == 2
    assert counts["primates__ortholog"] == 1
    assert counts["primates__order"] == 1

    summary = build_taxonomy_summary_rows(
        [
            {"query_gene_id": "1", "ortholog_gene_id": gene_id, "tax_id": tax_id}
            for gene_id, tax_id in members
        ],
        profiles,
    )
    by_key = {
        (row["taxonomic_scope"], row["evidence_unit"]): row
        for row in summary
    }
    assert by_key[("all", "ortholog")]["orthologs_per_gene_median"] == "2.0"
    assert by_key[("primates", "species")]["units_per_gene_median"] == "1.0"


def test_compact_events_write_taxonomic_alt_counts_in_sqlite(tmp_path: Path) -> None:
    events = tmp_path / "alignment_events.tsv.gz"
    compact = tmp_path / "compact.tsv.gz"
    ortholog_support = tmp_path / "event_ortholog_support.tsv.gz"
    support = tmp_path / "support.tsv.gz"
    fields = [
        "gene_id",
        "event_type",
        "target_start0",
        "target_end0",
        "genomic_accession",
        "genomic_start1",
        "genomic_end1",
        "ref",
        "alt",
        "ortholog_gene_id",
        "strategy",
        "tool",
        "preset",
        "tax_id",
        "taxname",
        "qc_flags",
    ]
    write_tsv_gz(
        events,
        fields,
        [
            {
                "gene_id": "1",
                "event_type": "snv",
                "target_start0": 4,
                "target_end0": 5,
                "genomic_accession": "NC_000001.11",
                "genomic_start1": 5,
                "genomic_end1": 5,
                "ref": "A",
                "alt": "G",
                "ortholog_gene_id": ortholog_gene_id,
                "strategy": "s1",
                "tool": "tool",
                "tax_id": tax_id,
            }
            for ortholog_gene_id, tax_id in [("chimp_gene", "9598"), ("mouse_gene", "10090")]
        ],
    )

    compact_count, raw_count, support_count, ortholog_support_count = write_compact_events(
        [events],
        compact,
        ortholog_support,
        taxonomy_fixture(tmp_path),
        support,
    )

    assert (compact_count, raw_count, support_count, ortholog_support_count) == (1, 2, 1, 2)
    with gzip.open(support, "rt", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["all__ortholog"] == "2"
    assert row["mammalia__species"] == "2"
    assert row["primates__species"] == "1"
    with gzip.open(compact, "rt", newline="") as handle:
        compact_row = next(csv.DictReader(handle, delimiter="\t"))
    with gzip.open(ortholog_support, "rt", newline="") as handle:
        ortholog_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert compact_row["event_group_id"] == "1"
    assert {row["event_group_id"] for row in ortholog_rows} == {"1"}
    assert {row["ortholog_gene_id"] for row in ortholog_rows} == {
        "chimp_gene",
        "mouse_gene",
    }


def test_compact_taxonomic_alt_counts_use_numeric_site_order(tmp_path: Path) -> None:
    events = tmp_path / "alignment_events.tsv.gz"
    compact = tmp_path / "compact.tsv.gz"
    ortholog_support = tmp_path / "event_ortholog_support.tsv.gz"
    support = tmp_path / "support.tsv.gz"
    fields = [
        "gene_id",
        "event_type",
        "target_start0",
        "target_end0",
        "genomic_accession",
        "genomic_start1",
        "genomic_end1",
        "ref",
        "alt",
        "ortholog_gene_id",
        "strategy",
        "tool",
        "preset",
        "tax_id",
        "taxname",
        "qc_flags",
    ]
    event_sites = [
        ("s1", 10_080, "chimp_10080"),
        ("s2", 10, "chimp_s2"),
        ("s1", 1_016, "chimp_1016"),
    ]
    write_tsv_gz(
        events,
        fields,
        [
            {
                "gene_id": "1",
                "event_type": "snv",
                "target_start0": position,
                "target_end0": position + 1,
                "genomic_accession": "NC_000001.11",
                "genomic_start1": position + 1,
                "genomic_end1": position + 1,
                "ref": "A",
                "alt": "G",
                "ortholog_gene_id": ortholog_gene_id,
                "strategy": strategy,
                "tool": "tool",
                "tax_id": "9598",
            }
            for strategy, position, ortholog_gene_id in event_sites
        ],
    )

    compact_count, raw_count, support_count, ortholog_support_count = write_compact_events(
        [events],
        compact,
        ortholog_support,
        taxonomy_fixture(tmp_path),
        support,
    )

    assert (compact_count, raw_count, support_count, ortholog_support_count) == (3, 3, 3, 3)
    with gzip.open(support, "rt", newline="") as handle:
        support_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [
        (row["gene_id"], row["strategy"], int(row["target_start0"]))
        for row in support_rows
    ] == [
        ("1", "s1", 1_016),
        ("1", "s1", 10_080),
        ("1", "s2", 10),
    ]


def test_compact_evidence_summary_preserves_scope_unit_and_gnomad_status(tmp_path: Path) -> None:
    profiles = load_taxonomy_profiles(taxonomy_fixture(tmp_path))
    depth_counts = count_member_groups(
        [("chimp_gene", "9598"), ("mouse_gene", "10090")],
        profiles,
    )
    chimp_alt = count_member_groups([("chimp_gene", "9598")], profiles)
    mouse_alt = count_member_groups([("mouse_gene", "10090")], profiles)
    depth = tmp_path / "snv_taxonomic_depth.tsv.gz"
    alt = tmp_path / "snv_alt_taxonomic_support.tsv.gz"
    features = tmp_path / "1.tsv.gz"
    output = tmp_path / "ortholog_evidence_summary.tsv.gz"
    write_tsv_gz(
        depth,
        ["gene_id", "strategy", "target_start0", *COUNT_KEYS],
        [{"gene_id": "1", "strategy": "s1", "target_start0": 4, **depth_counts}],
    )
    write_tsv_gz(
        alt,
        ["gene_id", "strategy", "target_start0", "ref", "alt", *COUNT_KEYS],
        [
            {"gene_id": "1", "strategy": "s1", "target_start0": 4, "ref": "A", "alt": "G", **chimp_alt},
            {"gene_id": "1", "strategy": "s1", "target_start0": 4, "ref": "A", "alt": "T", **mouse_alt},
        ],
    )
    write_tsv_gz(
        features,
        ["gene_id", "feature_type", "target_start0", "target_end0"],
        [
            {"gene_id": "1", "feature_type": "gene", "target_start0": 0, "target_end0": 10},
            {"gene_id": "1", "feature_type": "cds", "target_start0": 0, "target_end0": 10},
        ],
    )

    row_count = write_ortholog_evidence_summary(
        depth,
        alt,
        [features],
        {
            ("1", 4, "A", "G"): "found",
            ("1", 4, "A", "T"): "lookup_failed",
        },
        output,
    )

    assert row_count > 0
    with gzip.open(output, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    all_ortholog = next(
        row
        for row in rows
        if row["taxonomic_scope"] == "all"
        and row["evidence_unit"] == "ortholog"
        and row["alt_support_count"] == "1"
    )
    assert all_ortholog["site_aligned_count"] == "2"
    assert all_ortholog["gnomad_found_count"] == "1"
    assert all_ortholog["gnomad_lookup_failed_count"] == "1"
    primate = next(
        row
        for row in rows
        if row["taxonomic_scope"] == "primates"
        and row["evidence_unit"] == "species"
        and row["alt_support_count"] == "1"
    )
    assert primate["site_aligned_count"] == "1"
    assert primate["alt_support_count"] == "1"
    assert primate["gnomad_found_count"] == "1"


def test_finalizer_sums_matching_partition_histograms(tmp_path: Path) -> None:
    fields = [
        "strategy",
        "target_context",
        "taxonomic_scope",
        "evidence_unit",
        "site_aligned_count",
        "alt_support_count",
        "gnomad_found_count",
        "gnomad_not_found_count",
        "gnomad_lookup_failed_count",
    ]
    partitions = []
    for index, found in enumerate((2, 3), start=1):
        partition = tmp_path / f"partition_{index:06d}"
        partition.mkdir()
        write_tsv_gz(
            partition / "ortholog_evidence_summary.tsv.gz",
            fields,
            [
                {
                    "strategy": "s1",
                    "target_context": "cds",
                    "taxonomic_scope": "mammalia",
                    "evidence_unit": "species",
                    "site_aligned_count": 10,
                    "alt_support_count": 2,
                    "gnomad_found_count": found,
                    "gnomad_not_found_count": 1,
                    "gnomad_lookup_failed_count": 0,
                }
            ],
        )
        partitions.append((partition, {"ortholog_evidence_summary_count": 1}))

    output = tmp_path / "merged.tsv.gz"
    assert merge_ortholog_evidence(partitions, output) == 1
    with gzip.open(output, "rt", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["gnomad_found_count"] == "5"
    assert row["gnomad_not_found_count"] == "2"
