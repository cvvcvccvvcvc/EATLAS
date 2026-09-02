from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

from bin import build_fetch_dataset as fetch_dataset


def write_tsv_gz(path: Path, header: list[str] | None, rows: list[list[str]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if header is not None:
            writer.writerow(header)
        writer.writerows(rows)


def test_merge_tsv_gz_requires_every_chunk_table(tmp_path) -> None:
    existing = tmp_path / "existing.tsv.gz"
    write_tsv_gz(existing, ["gene_id"], [["1"]])

    with pytest.raises(FileNotFoundError, match="Missing chunk table"):
        fetch_dataset.merge_tsv_gz(
            [existing, tmp_path / "missing.tsv.gz"],
            tmp_path / "merged.tsv.gz",
        )


def test_merge_tsv_gz_rejects_missing_header(tmp_path) -> None:
    source = tmp_path / "empty.tsv.gz"
    write_tsv_gz(source, None, [])

    with pytest.raises(ValueError, match="has no header"):
        fetch_dataset.merge_tsv_gz([source], tmp_path / "merged.tsv.gz")


def test_merge_tsv_gz_rejects_header_mismatch(tmp_path) -> None:
    first = tmp_path / "first.tsv.gz"
    second = tmp_path / "second.tsv.gz"
    write_tsv_gz(first, ["gene_id", "status"], [["1", "ok"]])
    write_tsv_gz(second, ["gene_id", "result"], [["2", "ok"]])

    with pytest.raises(ValueError, match="header mismatch"):
        fetch_dataset.merge_tsv_gz([first, second], tmp_path / "merged.tsv.gz")


def test_merge_tsv_gz_rejects_row_width_mismatch(tmp_path) -> None:
    source = tmp_path / "source.tsv.gz"
    write_tsv_gz(source, ["gene_id", "status"], [["1"]])

    with pytest.raises(ValueError, match="row width mismatch"):
        fetch_dataset.merge_tsv_gz([source], tmp_path / "merged.tsv.gz")


def test_merge_tsv_gz_accepts_header_only_tables(tmp_path) -> None:
    first = tmp_path / "first.tsv.gz"
    second = tmp_path / "second.tsv.gz"
    output = tmp_path / "merged.tsv.gz"
    write_tsv_gz(first, ["gene_id", "status"], [])
    write_tsv_gz(second, ["gene_id", "status"], [])

    assert fetch_dataset.merge_tsv_gz([first, second], output) == 0
    with gzip.open(output, "rt") as handle:
        assert handle.read() == "gene_id\tstatus\n"


@pytest.mark.parametrize(
    ("target_gene_count", "selected_ortholog_count", "message"),
    [
        (0, 1, "no target genes"),
        (1, 0, "no selected orthologs"),
    ],
)
def test_validate_fetch_counts_rejects_unusable_dataset(
    target_gene_count: int,
    selected_ortholog_count: int,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        fetch_dataset.validate_fetch_counts(target_gene_count, selected_ortholog_count)


def test_validate_chunk_manifests_requires_every_expected_chunk() -> None:
    with pytest.raises(ValueError, match="missing=.*chunk_000002"):
        fetch_dataset.validate_chunk_manifests(
            {"chunk_000001", "chunk_000002"},
            [{"chunk_id": "chunk_000001"}],
        )


def test_consistent_manifest_value_returns_shared_value() -> None:
    assert fetch_dataset.consistent_manifest_value(
        [
            {"target_tax_id": "9606"},
            {"target_tax_id": "9606"},
        ],
        "target_tax_id",
    ) == "9606"


def test_consistent_manifest_value_rejects_mismatch() -> None:
    with pytest.raises(ValueError, match="must agree"):
        fetch_dataset.consistent_manifest_value(
            [
                {"target_assembly_name": "GRCh38.p14"},
                {"target_assembly_name": "other"},
            ],
            "target_assembly_name",
        )


def test_validate_gene_outcomes_accepts_success_or_failure() -> None:
    fetch_dataset.validate_gene_outcomes(
        {"1", "2"},
        ["1"],
        ["2"],
    )


def test_validate_gene_outcomes_rejects_unaccounted_gene() -> None:
    with pytest.raises(ValueError, match="missing=.*2"):
        fetch_dataset.validate_gene_outcomes(
            {"1", "2"},
            ["1"],
            [],
        )


def test_validate_sequence_gene_ids_accepts_exact_identity() -> None:
    fetch_dataset.validate_sequence_gene_ids(
        "Ortholog",
        {"1", "2"},
        {"2", "1"},
    )


@pytest.mark.parametrize(
    ("expected", "observed", "message"),
    [
        ({"1", "2"}, {"1"}, "missing_count=1, missing_sample=\\['2'\\]"),
        ({"1"}, {"1", "2"}, "unexpected_count=1, unexpected_sample=\\['2'\\]"),
    ],
)
def test_copied_sequence_gene_ids_reject_metadata_mismatch(
    tmp_path: Path,
    expected: set[str],
    observed: set[str],
    message: str,
) -> None:
    chunk_dir = tmp_path / "chunk"
    ortholog_dir = chunk_dir / "sequences" / "orthologs"
    ortholog_dir.mkdir(parents=True)
    for gene_id in observed:
        (ortholog_dir / f"{gene_id}.fa.gz").write_bytes(b"sequence")

    _, copied_gene_ids = fetch_dataset.copy_sequences(
        [chunk_dir],
        tmp_path / "fetch",
    )

    with pytest.raises(ValueError, match=message):
        fetch_dataset.validate_sequence_gene_ids(
            "Ortholog",
            expected,
            copied_gene_ids,
        )
