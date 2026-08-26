from __future__ import annotations

import csv
import gzip
import json
import random
from pathlib import Path

import pytest

from bin.fetch_taxonomy import fetch_taxonomy_records, taxonomy_row
from bin.merge_alignment_results import write_compact_events
from analytics.derivations.ortholog_evidence import write_ortholog_evidence_summary
from analytics.derivations.support import merge_ortholog_evidence
from analytics.derivations.taxonomy import (
    COUNT_KEYS,
    SCOPE_ANCESTORS,
    UNIT_ORDER,
    TaxonomyProfile,
    build_taxonomy_summary_rows,
    count_member_groups,
    load_taxonomy_profiles,
)
from genomics.taxonomy import TAXONOMY_FIELDS


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def taxonomy_fixture(tmp_path: Path) -> Path:
    path = tmp_path / "taxonomy.tsv.gz"
    write_tsv_gz(
        path,
        TAXONOMY_FIELDS,
        [
            {
                "tax_id": "9598",
                "taxonomy_status": "resolved",
                "species_id": "9598",
                "genus_id": "9596",
                "family_id": "9604",
                "order_id": "9443",
                "lineage_tax_ids": "2759,33208,7742,32523,32524,40674,9443,9598",
            },
            {
                "tax_id": "10090",
                "taxonomy_status": "resolved",
                "species_id": "10090",
                "genus_id": "10088",
                "family_id": "10066",
                "order_id": "9989",
                "lineage_tax_ids": "2759,33208,7742,32523,32524,40674,10090",
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
                "domain": {"id": 2759, "name": "Eukaryota"},
                "kingdom": {"id": 33208, "name": "Metazoa"},
                "phylum": {"id": 7711, "name": "Chordata"},
                "class": {"id": 40674, "name": "Mammalia"},
                "order": {"id": 9443, "name": "Primates"},
                "family": {"id": 9604, "name": "Hominidae"},
                "genus": {"id": 9596, "name": "Pan"},
                "species": {"id": 9598, "name": "Pan troglodytes"},
            },
            "parents": [1, 2759, 33208, 7711, 7742, 40674, 9443, 9443],
        },
    )

    assert list(row) == TAXONOMY_FIELDS
    assert row["taxonomy_status"] == "resolved"
    assert row["scientific_name"] == "Pan troglodytes"
    assert row["domain_id"] == "2759"
    assert row["kingdom_name"] == "Metazoa"
    assert row["phylum_id"] == "7711"
    assert row["species_id"] == "9598"
    assert row["genus_id"] == "9596"
    assert row["lineage_tax_ids"] == "1,2759,33208,7711,7742,40674,9443,9598"
    assert "is_primate" not in row
    assert "parent_tax_ids" not in row


def test_taxonomy_row_marks_missing_response_without_inventing_lineage() -> None:
    row = taxonomy_row("12345", None)

    assert list(row) == TAXONOMY_FIELDS
    assert row["tax_id"] == "12345"
    assert row["taxonomy_status"] == "not_returned"
    assert row["scientific_name"] == ""
    assert row["domain_id"] == ""
    assert row["species_id"] == ""
    assert row["lineage_tax_ids"] == ""


def test_taxonomy_loader_rejects_noncanonical_lineage_column(tmp_path: Path) -> None:
    path = tmp_path / "taxonomy.tsv.gz"
    write_tsv_gz(
        path,
        [
            "tax_id",
            "species_id",
            "genus_id",
            "family_id",
            "order_id",
            "parent_tax_ids",
        ],
        [
            {
                "tax_id": "9598",
                "species_id": "9598",
                "genus_id": "9596",
                "family_id": "9604",
                "order_id": "9443",
                "parent_tax_ids": "2759, 33208; 7742,9443",
            }
        ],
    )

    with pytest.raises(ValueError, match="exact canonical fields"):
        load_taxonomy_profiles(path)


def test_taxonomy_loader_rejects_extra_alias_and_invalid_status(tmp_path: Path) -> None:
    path = tmp_path / "taxonomy.tsv.gz"
    row = taxonomy_row("9598", None)
    write_tsv_gz(
        path,
        [*TAXONOMY_FIELDS, "parent_tax_ids"],
        [{**row, "parent_tax_ids": "2759,9443"}],
    )
    with pytest.raises(ValueError, match="exact canonical fields"):
        load_taxonomy_profiles(path)

    row["taxonomy_status"] = "unknown"
    write_tsv_gz(path, TAXONOMY_FIELDS, [row])
    with pytest.raises(ValueError, match="invalid taxonomy_status"):
        load_taxonomy_profiles(path)


def test_taxonomic_counts_reject_tax_id_absent_from_canonical_taxonomy(
    tmp_path: Path,
) -> None:
    profiles = load_taxonomy_profiles(taxonomy_fixture(tmp_path))

    with pytest.raises(ValueError, match="tax_id absent from canonical taxonomy"):
        count_member_groups([("unknown_gene", "999999")], profiles)


def test_taxonomy_summary_requires_exact_selected_tax_id_coverage(tmp_path: Path) -> None:
    profiles = load_taxonomy_profiles(taxonomy_fixture(tmp_path))
    rows = [
        {"query_gene_id": "1", "ortholog_gene_id": "chimp", "tax_id": "9598"},
        {"query_gene_id": "1", "ortholog_gene_id": "mouse", "tax_id": "10090"},
        {"query_gene_id": "1", "ortholog_gene_id": "unknown", "tax_id": "999999"},
    ]

    with pytest.raises(ValueError, match=r"missing taxonomy tax_id\(s\): 999999"):
        build_taxonomy_summary_rows(rows, profiles)

    with pytest.raises(ValueError, match=r"unexpected taxonomy tax_id\(s\): 10090"):
        build_taxonomy_summary_rows(rows[:1], profiles)


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

    monkeypatch.setattr("bin.fetch_taxonomy.subprocess.run", fake_run)

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


def test_missing_rank_fallback_does_not_collide_with_real_taxon_id() -> None:
    profiles = {
        "missing": TaxonomyProfile(
            tax_id="123",
            ancestor_ids=frozenset({"123"}),
            species_id="",
            genus_id="",
            family_id="",
            order_id="",
        ),
        "resolved": TaxonomyProfile(
            tax_id="456",
            ancestor_ids=frozenset({"456"}),
            species_id="123",
            genus_id="123",
            family_id="123",
            order_id="123",
        ),
    }

    counts = count_member_groups(
        [("missing_gene", "missing"), ("resolved_gene", "resolved")],
        profiles,
    )

    assert counts["all__species"] == 2
    assert counts["all__genus"] == 2


def test_taxonomic_counter_matches_set_reference_for_arbitrary_profiles() -> None:
    randomizer = random.Random(20_260_826)
    ancestor_ids = tuple(SCOPE_ANCESTORS.values())[1:]

    def reference_count(members, profiles):
        groups = [set() for _ in COUNT_KEYS]
        for ortholog_gene_id, tax_id in members:
            if not ortholog_gene_id:
                continue
            profile = profiles[tax_id]
            fallback_id = ("taxon", profile.tax_id)
            unit_ids = (
                ortholog_gene_id,
                profile.species_id or fallback_id,
                profile.genus_id or fallback_id,
                profile.family_id or fallback_id,
                profile.order_id or fallback_id,
            )
            for scope_index, ancestor_id in enumerate(SCOPE_ANCESTORS.values()):
                if ancestor_id and ancestor_id not in profile.ancestor_ids:
                    continue
                offset = scope_index * len(UNIT_ORDER)
                for unit_index, unit_id in enumerate(unit_ids):
                    groups[offset + unit_index].add(unit_id)
        return dict(
            zip(
                COUNT_KEYS,
                (len(group) for group in groups),
                strict=True,
            )
        )

    for _case in range(250):
        profiles = {}
        for tax_index in range(randomizer.randint(1, 16)):
            tax_id = str(10_000 + tax_index)
            profiles[tax_id] = TaxonomyProfile(
                tax_id=tax_id,
                ancestor_ids=frozenset(
                    ancestor
                    for ancestor in ancestor_ids
                    if randomizer.random() < 0.65
                ),
                species_id=randomizer.choice(("", tax_id, f"s{tax_index % 5}")),
                genus_id=randomizer.choice(("", f"g{tax_index % 4}")),
                family_id=randomizer.choice(("", f"f{tax_index % 3}")),
                order_id=randomizer.choice(("", f"o{tax_index % 2}")),
            )
        tax_ids = tuple(profiles)
        members = [
            (
                randomizer.choice(("", f"ortholog_{randomizer.randrange(12)}")),
                randomizer.choice(tax_ids),
            )
            for _member in range(randomizer.randint(0, 50))
        ]

        assert count_member_groups(members, profiles) == reference_count(
            members,
            profiles,
        )


def test_compact_events_preserve_exact_taxonomic_identities(tmp_path: Path) -> None:
    events = tmp_path / "alignment_events.tsv.gz"
    compact = tmp_path / "compact.tsv.gz"
    ortholog_support = tmp_path / "event_ortholog_support.tsv.gz"
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
        "mapq",
        "native_alignment_type",
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

    compact_count, raw_count, ortholog_support_count = write_compact_events(
        [events],
        compact,
        ortholog_support,
    )

    assert (compact_count, raw_count, ortholog_support_count) == (1, 2, 2)
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
    assert {row["tax_id"] for row in ortholog_rows} == {"9598", "10090"}


def test_compact_events_use_numeric_site_order(tmp_path: Path) -> None:
    events = tmp_path / "alignment_events.tsv.gz"
    compact = tmp_path / "compact.tsv.gz"
    ortholog_support = tmp_path / "event_ortholog_support.tsv.gz"
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
        "mapq",
        "native_alignment_type",
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

    compact_count, raw_count, ortholog_support_count = write_compact_events(
        [events],
        compact,
        ortholog_support,
    )

    assert (compact_count, raw_count, ortholog_support_count) == (3, 3, 3)
    with gzip.open(compact, "rt", newline="") as handle:
        compact_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [
        (row["gene_id"], row["strategy"], int(row["target_start0"]))
        for row in compact_rows
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
