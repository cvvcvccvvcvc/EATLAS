from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from analytics.io import run_inputs as run_inputs_module
from analytics.io.taxonomy_summary import (
    build_or_load_taxonomy_summary,
    resolve_taxonomy_summary_path,
)
from bin.taxonomic_evidence import (
    build_taxonomy_summary_rows,
    load_taxonomy_profiles,
    write_taxonomy_summary,
)


TAXONOMY_FIELDS = [
    "tax_id",
    "species_id",
    "genus_id",
    "family_id",
    "order_id",
    "parent_tax_ids",
]
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


def test_taxonomy_summary_resolution_uses_legacy_when_fetch_taxonomy_is_absent(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "run"
    legacy = run_dir / "alignment" / "taxonomy_summary.tsv.gz"
    _write_tsv_gz(
        run_dir / "fetch" / "orthologs.selected.tsv.gz",
        ORTHOLOG_FIELDS,
        [],
    )
    assert resolve_taxonomy_summary_path(run_dir) == legacy

    (run_dir / "fetch" / "orthologs.selected.tsv.gz").unlink()
    _write_tsv_gz(run_dir / "fetch" / "taxonomy.tsv.gz", TAXONOMY_FIELDS, [])
    with pytest.raises(FileNotFoundError, match="Incomplete Stage 1 taxonomy inputs"):
        resolve_taxonomy_summary_path(run_dir)


def test_run_inputs_exposes_new_taxonomy_summary_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    run_dir = tmp_path / "run"
    fetch_dir = run_dir / "fetch"
    annotation_dir = run_dir / "annotation"
    (fetch_dir / "sequences" / "targets").mkdir(parents=True)
    annotation_dir.mkdir()
    _write_taxonomy_inputs(fetch_dir)
    pd.DataFrame(columns=["gene_id"]).to_csv(
        fetch_dir / "genes.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame(columns=["gene_id"]).to_csv(
        fetch_dir / "target_features.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    source_annotations = annotation_dir / "variant_annotations.tsv.gz"
    pd.DataFrame(columns=["variant_key"]).to_csv(
        source_annotations, sep="\t", index=False, compression="gzip"
    )
    monkeypatch.setattr(
        run_inputs_module,
        "resolve_vep_variant_annotations",
        lambda _run_dir, source: source,
    )

    inputs = run_inputs_module.resolve_run_inputs(run_dir)

    assert inputs.taxonomy_summary_tsv == (
        run_dir / "analytics" / "taxonomy_summary" / "taxonomy_summary.tsv.gz"
    )
    assert inputs.taxonomy_summary_tsv.is_file()
