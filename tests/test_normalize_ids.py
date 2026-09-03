from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "normalize_ids.py"


def run_normalizer(
    input_path: Path,
    output_dir: Path,
    *,
    chunk_size: int,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--ids-file",
            str(input_path),
            "--chunk-size",
            str(chunk_size),
            "--outdir",
            str(output_dir),
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def test_normalizer_preserves_first_occurrence_and_writes_deterministic_chunks(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "ids.txt"
    input_path.write_text("# panel\n672, 7157\n672\n5728 1956,4609\n")
    output_dir = tmp_path / "normalized"

    completed = run_normalizer(input_path, output_dir, chunk_size=2)

    assert completed.returncode == 0, completed.stderr
    rows = pd.read_csv(output_dir / "input.ids.tsv", sep="\t", dtype=str).fillna("")
    assert rows.to_dict(orient="records") == [
        {
            "input_position": "1",
            "line_number": "2",
            "raw_value": "672",
            "gene_id": "672",
            "accepted": "true",
            "accepted_index": "1",
            "duplicate_of_index": "",
        },
        {
            "input_position": "2",
            "line_number": "2",
            "raw_value": "7157",
            "gene_id": "7157",
            "accepted": "true",
            "accepted_index": "2",
            "duplicate_of_index": "",
        },
        {
            "input_position": "3",
            "line_number": "3",
            "raw_value": "672",
            "gene_id": "672",
            "accepted": "false",
            "accepted_index": "",
            "duplicate_of_index": "1",
        },
        {
            "input_position": "4",
            "line_number": "4",
            "raw_value": "5728",
            "gene_id": "5728",
            "accepted": "true",
            "accepted_index": "3",
            "duplicate_of_index": "",
        },
        {
            "input_position": "5",
            "line_number": "4",
            "raw_value": "1956",
            "gene_id": "1956",
            "accepted": "true",
            "accepted_index": "4",
            "duplicate_of_index": "",
        },
        {
            "input_position": "6",
            "line_number": "4",
            "raw_value": "4609",
            "gene_id": "4609",
            "accepted": "true",
            "accepted_index": "5",
            "duplicate_of_index": "",
        },
    ]
    chunks = pd.read_csv(output_dir / "chunks.tsv", sep="\t", dtype=str)
    assert chunks[["chunk_id", "chunk_file", "gene_count"]].to_dict(
        orient="records"
    ) == [
        {
            "chunk_id": "chunk_000001",
            "chunk_file": "chunks/chunk_000001.ids.txt",
            "gene_count": "2",
        },
        {
            "chunk_id": "chunk_000002",
            "chunk_file": "chunks/chunk_000002.ids.txt",
            "gene_count": "2",
        },
        {
            "chunk_id": "chunk_000003",
            "chunk_file": "chunks/chunk_000003.ids.txt",
            "gene_count": "1",
        },
    ]
    assert [
        path.read_text()
        for path in sorted((output_dir / "chunks").glob("*.ids.txt"))
    ] == ["672\n7157\n", "5728\n1956\n", "4609\n"]


@pytest.mark.parametrize(
    ("contents", "error"),
    [
        ("672\nnot-an-id\n", "expected a positive integer"),
        ("# comments only\n\n", "No Entrez Gene IDs found"),
    ],
)
def test_normalizer_rejects_invalid_or_empty_input(
    tmp_path: Path,
    contents: str,
    error: str,
) -> None:
    input_path = tmp_path / "ids.txt"
    input_path.write_text(contents)

    completed = run_normalizer(input_path, tmp_path / "normalized", chunk_size=2)

    assert completed.returncode != 0
    assert error in completed.stderr


def test_normalizer_rejects_nonpositive_chunk_size(tmp_path: Path) -> None:
    input_path = tmp_path / "ids.txt"
    input_path.write_text("672\n")

    completed = run_normalizer(input_path, tmp_path / "normalized", chunk_size=0)

    assert completed.returncode != 0
    assert "--chunk-size must be a positive integer" in completed.stderr
