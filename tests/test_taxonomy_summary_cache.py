from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import pytest

from analytics.io.taxonomy_summary import (
    build_or_load_taxonomy_summary,
    build_or_load_taxonomy_summary_many,
    resolve_taxonomy_summary_path,
)
from bin.fetch_taxonomy import TAXONOMY_FIELDS
from analytics.derivations.taxonomy import (
    build_taxonomy_summary_rows,
    load_taxonomy_profiles,
    write_taxonomy_summary,
)


ORTHOLOG_FIELDS = ["query_gene_id", "ortholog_gene_id", "tax_id"]


def _write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
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


def _write_taxonomy_inputs(fetch_dir: Path) -> tuple[Path, Path]:
    taxonomy = fetch_dir / "taxonomy.tsv.gz"
    orthologs = fetch_dir / "orthologs.selected.tsv.gz"
    _write_tsv_gz(
        taxonomy,
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
    _write_tsv_gz(
        orthologs,
        ORTHOLOG_FIELDS,
        [
            {"query_gene_id": "1", "ortholog_gene_id": "chimp_1", "tax_id": "9598"},
            {"query_gene_id": "1", "ortholog_gene_id": "mouse_1", "tax_id": "10090"},
            {"query_gene_id": "2", "ortholog_gene_id": "chimp_2", "tax_id": "9598"},
        ],
    )
    return taxonomy, orthologs


def test_analytics_cache_exactly_matches_current_taxonomy_summary_builder(
    tmp_path: Path,
) -> None:
    taxonomy, orthologs = _write_taxonomy_inputs(tmp_path / "fetch")
    actual = build_or_load_taxonomy_summary(
        taxonomy_tsv=taxonomy,
        orthologs_tsv=orthologs,
        analytics_dir=tmp_path / "analytics",
    )

    expected = tmp_path / "alignment" / "taxonomy_summary.tsv.gz"
    profiles = load_taxonomy_profiles(taxonomy)
    with gzip.open(orthologs, "rt", newline="") as handle:
        rows = build_taxonomy_summary_rows(
            csv.DictReader(handle, delimiter="\t"),
            profiles,
        )
    write_taxonomy_summary(expected, rows)

    with gzip.open(actual, "rt", newline="") as handle:
        actual_text = handle.read()
    with gzip.open(expected, "rt", newline="") as handle:
        expected_text = handle.read()
    assert actual_text == expected_text


def test_taxonomy_summary_cache_is_fingerprinted_by_both_inputs(tmp_path: Path) -> None:
    taxonomy, orthologs = _write_taxonomy_inputs(tmp_path / "fetch")
    analytics_dir = tmp_path / "analytics"
    output = build_or_load_taxonomy_summary(
        taxonomy_tsv=taxonomy,
        orthologs_tsv=orthologs,
        analytics_dir=analytics_dir,
    )
    manifest_path = analytics_dir / "taxonomy_summary" / "manifest.json"
    first_manifest = json.loads(manifest_path.read_text())
    first_mtime = output.stat().st_mtime_ns

    assert build_or_load_taxonomy_summary(
        taxonomy_tsv=taxonomy,
        orthologs_tsv=orthologs,
        analytics_dir=analytics_dir,
    ) == output
    assert output.stat().st_mtime_ns == first_mtime

    _write_tsv_gz(
        orthologs,
        ORTHOLOG_FIELDS,
        [
            {"query_gene_id": "1", "ortholog_gene_id": "chimp_1", "tax_id": "9598"},
            {"query_gene_id": "1", "ortholog_gene_id": "mouse_2", "tax_id": "10090"},
        ],
    )
    build_or_load_taxonomy_summary(
        taxonomy_tsv=taxonomy,
        orthologs_tsv=orthologs,
        analytics_dir=analytics_dir,
    )
    second_manifest = json.loads(manifest_path.read_text())

    assert second_manifest["fingerprint"] != first_manifest["fingerprint"]
    assert second_manifest["inputs"] != first_manifest["inputs"]


def test_taxonomy_summary_combines_disjoint_source_rows_exactly(tmp_path: Path) -> None:
    first_taxonomy, first_orthologs = _write_taxonomy_inputs(tmp_path / "first")
    second_taxonomy, second_orthologs = _write_taxonomy_inputs(tmp_path / "second")
    _write_tsv_gz(
        first_orthologs,
        ORTHOLOG_FIELDS,
        [
            {"query_gene_id": "1", "ortholog_gene_id": "chimp_1", "tax_id": "9598"},
            {"query_gene_id": "1", "ortholog_gene_id": "mouse_1", "tax_id": "10090"},
        ],
    )
    _write_tsv_gz(
        second_orthologs,
        ORTHOLOG_FIELDS,
        [
            {"query_gene_id": "2", "ortholog_gene_id": "chimp_2", "tax_id": "9598"}
        ],
    )

    output = build_or_load_taxonomy_summary_many(
        taxonomy_tsvs=(first_taxonomy, second_taxonomy),
        orthologs_tsvs=(first_orthologs, second_orthologs),
        analytics_dir=tmp_path / "analytics",
    )

    with gzip.open(output, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    all_orthologs = next(
        row
        for row in rows
        if row["taxonomic_scope"] == "all"
        and row["evidence_unit"] == "ortholog"
    )
    assert all_orthologs["gene_count"] == "2"
    assert all_orthologs["ortholog_count"] == "3"
    assert all_orthologs["taxon_count"] == "2"
    assert all_orthologs["unit_count"] == "3"
    assert all_orthologs["orthologs_per_gene_median"] == "1.5"


def test_taxonomy_summary_resolution_requires_canonical_fetch_contract(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "run"
    _write_tsv_gz(
        run_dir / "fetch" / "orthologs.selected.tsv.gz",
        ORTHOLOG_FIELDS,
        [],
    )
    with pytest.raises(FileNotFoundError, match="Incomplete Stage 1 taxonomy contract"):
        resolve_taxonomy_summary_path(run_dir, analytics_dir=tmp_path / "analytics")

    (run_dir / "fetch" / "orthologs.selected.tsv.gz").unlink()
    _write_tsv_gz(run_dir / "fetch" / "taxonomy.tsv.gz", TAXONOMY_FIELDS, [])
    with pytest.raises(FileNotFoundError, match="Incomplete Stage 1 taxonomy contract"):
        resolve_taxonomy_summary_path(run_dir, analytics_dir=tmp_path / "analytics")
