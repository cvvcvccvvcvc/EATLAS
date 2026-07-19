from __future__ import annotations

import sys
from pathlib import Path

import pytest


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

import build_fetch_dataset as fetch_dataset  # noqa: E402


def test_validate_chunk_manifests_requires_every_expected_chunk() -> None:
    with pytest.raises(ValueError, match="missing=.*chunk_000002"):
        fetch_dataset.validate_chunk_manifests(
            {"chunk_000001", "chunk_000002"},
            [{"chunk_id": "chunk_000001"}],
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
