from __future__ import annotations

import pytest

from bin import build_fetch_dataset as fetch_dataset


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
