from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
MERGE_SCRIPT = PROJECT_DIR / "bin" / "merge_alignment_results.py"

TABLE_HEADERS = {
    "ortholog_alignment_summary.tsv.gz": [
        "gene_id",
        "strategy",
        "status",
        "event_count",
        "aligned_target_bp",
    ],
    "alignment_segments.tsv.gz": ["gene_id"],
    "feature_coverage.tsv.gz": ["gene_id"],
    "alignment_events.tsv.gz": [
        "gene_id",
        "event_type",
        "target_start0",
        "target_end0",
        "genomic_accession",
        "genomic_start1",
        "genomic_end1",
        "ref",
        "alt",
    ],
    "failures.tsv.gz": ["gene_id"],
}


def write_tsv_gz(path: Path, header: list[str], rows: list[list[str]] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows or [])


def write_result_dir(
    root: Path,
    name: str,
    manifest: dict,
    *,
    missing_table: str | None = None,
) -> Path:
    result_dir = root / name
    result_dir.mkdir(parents=True)
    (result_dir / "manifest.json").write_text(json.dumps(manifest) + "\n")

    gene_ids = manifest.get("gene_ids") or [manifest.get("gene_id", "")]
    strategies = manifest.get("strategies") or [manifest.get("strategy", "")]
    summary_rows = [
        [gene_id, strategy, "aligned", "0", "1"]
        for gene_id in gene_ids
        for strategy in strategies
    ]
    for filename, header in TABLE_HEADERS.items():
        if filename == missing_table:
            continue
        rows = summary_rows if filename == "ortholog_alignment_summary.tsv.gz" else []
        write_tsv_gz(result_dir / filename, header, rows)
    return result_dir


def run_merge(arguments: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(MERGE_SCRIPT), *arguments],
        cwd=PROJECT_DIR,
        text=True,
        capture_output=True,
        check=False,
    )


def partition_arguments(
    result_dirs: list[Path],
    outdir: Path,
    *,
    gene_ids: str,
    strategies: str,
) -> list[str]:
    arguments = [
        "--partition-id",
        "partition_000001",
        "--expected-gene-ids",
        gene_ids,
        "--expected-strategies",
        strategies,
        "--outdir",
        str(outdir),
    ]
    for result_dir in result_dirs:
        arguments.extend(["--result-dir", str(result_dir)])
    return arguments


def write_final_inputs(
    root: Path,
    task_rows: list[list[str]],
) -> tuple[Path, Path, Path, Path]:
    alignment_tasks = root / "alignment_tasks.tsv.gz"
    taxonomy_presets = root / "taxonomy_presets.tsv.gz"
    taxonomy_failures = root / "taxonomy_failures.tsv.gz"
    target_features = root / "target_features.tsv.gz"
    write_tsv_gz(alignment_tasks, ["gene_id", "status"], task_rows)
    write_tsv_gz(taxonomy_presets, ["tax_id"])
    write_tsv_gz(taxonomy_failures, ["tax_id"])
    write_tsv_gz(target_features, ["gene_id"])
    return alignment_tasks, taxonomy_presets, taxonomy_failures, target_features


def final_arguments(
    result_dirs: list[Path],
    outdir: Path,
    inputs: tuple[Path, Path, Path, Path],
    *,
    strategies: str,
) -> list[str]:
    alignment_tasks, taxonomy_presets, taxonomy_failures, target_features = inputs
    arguments = [
        "--alignment-tasks",
        str(alignment_tasks),
        "--taxonomy-presets",
        str(taxonomy_presets),
        "--taxonomy-failures",
        str(taxonomy_failures),
        "--target-features",
        str(target_features),
        "--expected-strategies",
        strategies,
        "--outdir",
        str(outdir),
    ]
    for result_dir in result_dirs:
        arguments.extend(["--result-dir", str(result_dir)])
    return arguments


def test_partition_merge_writes_exact_gene_manifest(tmp_path: Path) -> None:
    result_dirs = []
    for gene_id in ["1", "2"]:
        result_dirs.append(
            write_result_dir(
                tmp_path,
                f"gene_{gene_id}_s1",
                {"gene_id": gene_id, "strategy": "s1"},
            )
        )
        result_dirs.append(
            write_result_dir(
                tmp_path,
                f"gene_{gene_id}_combined",
                {"gene_id": gene_id, "strategies": ["s2", "s3"]},
            )
        )
    outdir = tmp_path / "merged"

    completed = run_merge(
        partition_arguments(
            result_dirs,
            outdir,
            gene_ids="1,2",
            strategies="s1,s2,s3",
        )
    )

    assert completed.returncode == 0, completed.stderr
    manifest = json.loads((outdir / "manifest.json").read_text())
    assert manifest["gene_count"] == 2
    assert manifest["gene_ids"] == ["1", "2"]
    assert manifest["strategies"] == ["s1", "s2", "s3"]


def test_partition_merge_rejects_missing_strategy(tmp_path: Path) -> None:
    result_dir = write_result_dir(
        tmp_path,
        "gene_1_s1",
        {"gene_id": "1", "strategy": "s1"},
    )

    completed = run_merge(
        partition_arguments(
            [result_dir],
            tmp_path / "merged",
            gene_ids="1",
            strategies="s1,s2",
        )
    )

    assert completed.returncode != 0
    assert "missing=1:s2" in completed.stderr


def test_partition_merge_supports_compact_events(tmp_path: Path) -> None:
    result_dirs = [
        write_result_dir(
            tmp_path,
            f"gene_1_{strategy}",
            {"gene_id": "1", "strategy": strategy},
        )
        for strategy in ["s1", "s2"]
    ]
    arguments = partition_arguments(
        result_dirs,
        tmp_path / "merged",
        gene_ids="1",
        strategies="s1,s2",
    )
    arguments.append("--compact-events")

    completed = run_merge(arguments)

    assert completed.returncode == 0, completed.stderr
    manifest = json.loads((tmp_path / "merged" / "manifest.json").read_text())
    assert manifest["alignment_event_mode"] == "compact_support"
    assert manifest["alignment_event_count"] == 0


def test_partition_merge_rejects_duplicate_gene_strategy(tmp_path: Path) -> None:
    result_dirs = [
        write_result_dir(
            tmp_path,
            name,
            {"gene_id": "1", "strategy": "s1"},
        )
        for name in ["first", "second"]
    ]

    completed = run_merge(
        partition_arguments(
            result_dirs,
            tmp_path / "merged",
            gene_ids="1",
            strategies="s1",
        )
    )

    assert completed.returncode != 0
    assert "Duplicate gene-strategy alignment result" in completed.stderr


def test_partition_merge_rejects_missing_required_table(tmp_path: Path) -> None:
    result_dir = write_result_dir(
        tmp_path,
        "gene_1_s1",
        {"gene_id": "1", "strategy": "s1"},
        missing_table="alignment_events.tsv.gz",
    )

    completed = run_merge(
        partition_arguments(
            [result_dir],
            tmp_path / "merged",
            gene_ids="1",
            strategies="s1",
        )
    )

    assert completed.returncode != 0
    assert "Missing required alignment table" in completed.stderr
    assert "alignment_events.tsv.gz" in completed.stderr


def test_final_merge_writes_exact_gene_manifest(tmp_path: Path) -> None:
    partition_dirs = [
        write_result_dir(
            tmp_path,
            f"partition_{index:06d}",
            {
                "partition_id": f"partition_{index:06d}",
                "gene_count": 1,
                "gene_ids": [gene_id],
                "strategies": ["s1"],
            },
        )
        for index, gene_id in enumerate(["1", "2"], start=1)
    ]
    inputs = write_final_inputs(tmp_path, [["1", "ready"], ["2", "ready"]])
    outdir = tmp_path / "merged"

    completed = run_merge(
        final_arguments(
            partition_dirs,
            outdir,
            inputs,
            strategies="s1",
        )
    )

    assert completed.returncode == 0, completed.stderr
    manifest = json.loads((outdir / "manifest.json").read_text())
    assert manifest["gene_count"] == 2
    assert manifest["gene_ids"] == ["1", "2"]


def test_final_merge_rejects_missing_ready_gene(tmp_path: Path) -> None:
    partition_dir = write_result_dir(
        tmp_path,
        "partition_000001",
        {
            "partition_id": "partition_000001",
            "gene_count": 1,
            "gene_ids": ["1"],
            "strategies": ["s1"],
        },
    )
    inputs = write_final_inputs(tmp_path, [["1", "ready"], ["2", "ready"]])

    completed = run_merge(
        final_arguments(
            [partition_dir],
            tmp_path / "merged",
            inputs,
            strategies="s1",
        )
    )

    assert completed.returncode != 0
    assert "Final alignment gene coverage mismatch" in completed.stderr
    assert "missing=['2']" in completed.stderr
