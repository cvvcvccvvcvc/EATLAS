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
    "strategy_summary.tsv.gz": [
        "strategy",
        "summary_row_count",
        "gene_count",
        "aligned_summary_row_count",
        "event_count",
        "aligned_target_bp",
    ],
    "alignment_segments.tsv.gz": [
        "gene_id",
        "strategy",
        "ortholog_gene_id",
        "target_start0",
        "target_end0",
    ],
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
        "ortholog_gene_id",
        "strategy",
        "tool",
        "preset",
        "tax_id",
        "taxname",
        "qc_flags",
    ],
    "snv_site_depth.tsv.gz": [
        "gene_id",
        "strategy",
        "target_start0",
        "site_aligned_ortholog_count",
    ],
    "failures.tsv.gz": ["gene_id"],
}
EVENT_ORTHOLOG_SUPPORT_HEADER = [
    "gene_id",
    "event_type",
    "target_start0",
    "target_end0",
    "genomic_accession",
    "genomic_start1",
    "genomic_end1",
    "ref",
    "alt",
    "strategy",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "support_row_count",
]


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
    strategy_rows = [
        [strategy, str(len(gene_ids)), str(len(gene_ids)), str(len(gene_ids)), "0", str(len(gene_ids))]
        for strategy in strategies
    ]
    for filename, header in TABLE_HEADERS.items():
        if filename == missing_table:
            continue
        if filename == "ortholog_alignment_summary.tsv.gz":
            rows = summary_rows
        elif filename == "strategy_summary.tsv.gz":
            rows = strategy_rows
        else:
            rows = []
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
    alignment_tasks: Path | None = None,
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
    if alignment_tasks is not None:
        arguments.extend(["--alignment-tasks", str(alignment_tasks)])
    for result_dir in result_dirs:
        arguments.extend(["--result-dir", str(result_dir)])
    return arguments


def write_final_inputs(
    root: Path,
    task_rows: list[list[str]],
    task_header: list[str] | None = None,
) -> tuple[Path, Path, Path, Path]:
    alignment_tasks = root / "alignment_tasks.tsv.gz"
    taxonomy_presets = root / "taxonomy_presets.tsv.gz"
    taxonomy_failures = root / "taxonomy_failures.tsv.gz"
    target_features = root / "target_features.tsv.gz"
    write_tsv_gz(alignment_tasks, task_header or ["gene_id", "status"], task_rows)
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


def test_partition_merge_uses_strategy_specific_gene_eligibility(tmp_path: Path) -> None:
    ensembl = "precomputed_ensembl_92_mammals_epo_extended"
    alignment_tasks = tmp_path / "alignment_tasks.tsv.gz"
    write_tsv_gz(
        alignment_tasks,
        ["gene_id", "status", "target_ready", "ortholog_ready"],
        [
            ["1", "ready", "true", "true"],
            ["2", "missing_ortholog_fasta", "true", "false"],
        ],
    )
    result_dirs = [
        write_result_dir(tmp_path, "gene_1_local", {"gene_id": "1", "strategy": "s1"}),
        write_result_dir(tmp_path, "gene_1_ensembl", {"gene_id": "1", "strategy": ensembl}),
        write_result_dir(tmp_path, "gene_2_ensembl", {"gene_id": "2", "strategy": ensembl}),
    ]
    outdir = tmp_path / "merged"

    completed = run_merge(
        partition_arguments(
            result_dirs,
            outdir,
            gene_ids="1,2",
            strategies=f"s1,{ensembl}",
            alignment_tasks=alignment_tasks,
        )
    )

    assert completed.returncode == 0, completed.stderr
    manifest = json.loads((outdir / "manifest.json").read_text())
    assert manifest["gene_ids"] == ["1", "2"]
    assert manifest["strategy_eligible_gene_counts"] == {"s1": 1, ensembl: 2}


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


def test_partition_annotation_input_keeps_compact_annotation_tables(tmp_path: Path) -> None:
    result_dir = write_result_dir(
        tmp_path,
        "gene_1_s1",
        {
            "gene_id": "1",
            "strategy": "s1",
            "alignment_segment_count": 2,
        },
    )
    write_tsv_gz(
        result_dir / "alignment_segments.tsv.gz",
        TABLE_HEADERS["alignment_segments.tsv.gz"],
        [
            ["1", "s1", "101", "0", "10"],
            ["1", "s1", "102", "0", "10"],
        ],
    )
    write_tsv_gz(
        result_dir / "alignment_events.tsv.gz",
        TABLE_HEADERS["alignment_events.tsv.gz"],
        [
            ["1", "snv", "4", "5", "NC_1", "5", "5", "A", "G", "101", "s1", "tool", "", "1", "species", ""],
            ["1", "snv", "4", "5", "NC_1", "5", "5", "A", "T", "102", "s1", "tool", "", "2", "species", ""],
        ],
    )
    arguments = partition_arguments(
        [result_dir],
        tmp_path / "merged",
        gene_ids="1",
        strategies="s1",
    )
    arguments.extend(["--output-profile", "annotation-input"])

    completed = run_merge(arguments)

    assert completed.returncode == 0, completed.stderr
    outdir = tmp_path / "merged"
    assert {
        path.name
        for path in outdir.iterdir()
    } == {
        "alignment_events.tsv.gz",
        "failures.tsv.gz",
        "feature_coverage.tsv.gz",
        "manifest.json",
        "snv_site_depth.tsv.gz",
        "strategy_summary.tsv.gz",
    }
    manifest = json.loads((outdir / "manifest.json").read_text())
    assert manifest["output_profile"] == "annotation-input"
    assert manifest["alignment_segment_count"] == 2
    assert manifest["snv_site_depth_count"] == 1
    with gzip.open(outdir / "snv_site_depth.tsv.gz", "rt", newline="") as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == [
            {
                "gene_id": "1",
                "strategy": "s1",
                "target_start0": "4",
                "site_aligned_ortholog_count": "2",
            }
        ]


def test_compact_events_preserve_strategy_specific_support(tmp_path: Path) -> None:
    result_dirs = [
        write_result_dir(
            tmp_path,
            f"gene_1_{strategy}",
            {"gene_id": "1", "strategy": strategy},
        )
        for strategy in ["s1", "s2"]
    ]
    event_header = TABLE_HEADERS["alignment_events.tsv.gz"]
    for result_dir, strategy in zip(result_dirs, ["s1", "s2"]):
        write_tsv_gz(
            result_dir / "alignment_segments.tsv.gz",
            TABLE_HEADERS["alignment_segments.tsv.gz"],
            [["1", strategy, "101", "0", "1"]],
        )
        write_tsv_gz(
            result_dir / "alignment_events.tsv.gz",
            event_header,
            [["1", "snv", "0", "1", "NC_1", "1", "1", "A", "G", "101", strategy, "tool", "", "1", "species", ""]],
        )

    arguments = partition_arguments(
        result_dirs,
        tmp_path / "merged",
        gene_ids="1",
        strategies="s1,s2",
    )
    arguments.append("--compact-events")

    completed = run_merge(arguments)

    assert completed.returncode == 0, completed.stderr
    with gzip.open(tmp_path / "merged" / "alignment_events.tsv.gz", "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["strategy"] for row in rows] == ["s1", "s2"]
    assert [row["support_ortholog_count"] for row in rows] == ["1", "1"]
    with gzip.open(
        tmp_path / "merged" / "event_ortholog_support.tsv.gz",
        "rt",
        newline="",
    ) as handle:
        ortholog_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["strategy"] for row in ortholog_rows] == ["s1", "s2"]
    assert [row["ortholog_gene_id"] for row in ortholog_rows] == ["101", "101"]
    assert [row["support_row_count"] for row in ortholog_rows] == ["1", "1"]


def test_compact_events_accept_large_allele_fields(tmp_path: Path) -> None:
    strategy = "precomputed_ensembl_92_mammals_epo_extended"
    result_dir = write_result_dir(
        tmp_path,
        "partition_000041",
        {
            "partition_id": "partition_000041",
            "gene_count": 1,
            "gene_ids": ["3492"],
            "strategies": [strategy],
        },
    )
    large_alt = "A" * 165_969
    write_tsv_gz(
        result_dir / "alignment_events.tsv.gz",
        TABLE_HEADERS["alignment_events.tsv.gz"],
        [
            [
                "3492",
                "ins",
                "431282",
                "431282",
                "NC_000014.9",
                "106017719",
                "106017719",
                "",
                large_alt,
                "pongo_abelii",
                strategy,
                "ensembl_compara",
                "",
                "9601",
                "Pongo abelii",
                "",
            ]
        ],
    )

    completed = run_merge(
        [
            *final_arguments(
                [result_dir],
                tmp_path / "merged",
                write_final_inputs(tmp_path, [["3492", "ready"]]),
                strategies=strategy,
            ),
            "--compact-events",
        ]
    )

    assert completed.returncode == 0, completed.stderr
    with gzip.open(tmp_path / "merged" / "alignment_events.tsv.gz", "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        row = handle.readline().rstrip("\n").split("\t")
    assert len(row[header.index("alt")]) == len(large_alt)


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


def test_final_merge_reports_strategy_specific_eligible_gene_counts(tmp_path: Path) -> None:
    ensembl = "precomputed_ensembl_92_mammals_epo_extended"
    partition_dirs = [
        write_result_dir(
            tmp_path,
            "partition_000001",
            {
                "partition_id": "partition_000001",
                "gene_count": 1,
                "gene_ids": ["1"],
                "strategies": ["s1", ensembl],
            },
        ),
        write_result_dir(
            tmp_path,
            "partition_000002",
            {
                "partition_id": "partition_000002",
                "gene_count": 1,
                "gene_ids": ["2"],
                "strategies": [ensembl],
            },
        ),
    ]
    inputs = write_final_inputs(
        tmp_path,
        [
            ["1", "ready", "true", "true"],
            ["2", "missing_ortholog_fasta", "true", "false"],
        ],
        ["gene_id", "status", "target_ready", "ortholog_ready"],
    )
    outdir = tmp_path / "merged"

    completed = run_merge(
        final_arguments(
            partition_dirs,
            outdir,
            inputs,
            strategies=f"s1,{ensembl}",
        )
    )

    assert completed.returncode == 0, completed.stderr
    manifest = json.loads((outdir / "manifest.json").read_text())
    assert manifest["gene_ids"] == ["1", "2"]
    assert manifest["strategy_eligible_gene_counts"] == {"s1": 1, ensembl: 2}


def test_final_merge_preserves_precompacted_ortholog_support(tmp_path: Path) -> None:
    partition_dirs = []
    for index, (gene_id, ortholog_gene_id) in enumerate([("1", "101"), ("2", "201")], 1):
        partition = write_result_dir(
            tmp_path,
            f"partition_{index:06d}",
            {
                "partition_id": f"partition_{index:06d}",
                "gene_count": 1,
                "gene_ids": [gene_id],
                "strategies": ["s1"],
                "alignment_event_mode": "compact_support",
                "alignment_event_count": 1,
                "raw_alignment_event_count": 1,
                "event_ortholog_support_count": 1,
            },
        )
        write_tsv_gz(
            partition / "event_ortholog_support.tsv.gz",
            EVENT_ORTHOLOG_SUPPORT_HEADER,
            [
                [
                    gene_id,
                    "snv",
                    "0",
                    "1",
                    "NC_1",
                    "1",
                    "1",
                    "A",
                    "G",
                    "s1",
                    ortholog_gene_id,
                    "10090",
                    "Mus musculus",
                    "1",
                ]
            ],
        )
        partition_dirs.append(partition)
    inputs = write_final_inputs(tmp_path, [["1", "ready"], ["2", "ready"]])
    outdir = tmp_path / "merged"
    arguments = final_arguments(partition_dirs, outdir, inputs, strategies="s1")
    arguments.extend(["--compact-events", "--events-already-compacted"])

    completed = run_merge(arguments)

    assert completed.returncode == 0, completed.stderr
    with gzip.open(outdir / "event_ortholog_support.tsv.gz", "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["ortholog_gene_id"] for row in rows] == ["101", "201"]
    manifest = json.loads((outdir / "manifest.json").read_text())
    assert manifest["event_ortholog_support_count"] == 2


def test_final_report_input_omits_handoff_tables(tmp_path: Path) -> None:
    partition_dirs = [
        write_result_dir(
            tmp_path,
            f"partition_{index:06d}",
            {
                "partition_id": f"partition_{index:06d}",
                "gene_count": 1,
                "gene_ids": [gene_id],
                "strategies": ["s1"],
                "ortholog_alignment_summary_count": 4,
                "alignment_segment_count": 5,
                "alignment_event_count": 6,
                "raw_alignment_event_count": 6,
                "alignment_event_mode": "raw",
            },
        )
        for index, gene_id in enumerate(["1", "2"], start=1)
    ]
    inputs = write_final_inputs(tmp_path, [["1", "ready"], ["2", "ready"]])
    outdir = tmp_path / "merged"
    arguments = final_arguments(
        partition_dirs,
        outdir,
        inputs,
        strategies="s1",
    )
    arguments.extend(["--output-profile", "report-input"])

    completed = run_merge(arguments)

    assert completed.returncode == 0, completed.stderr
    assert {
        path.name
        for path in outdir.iterdir()
    } == {
        "failures.tsv.gz",
        "feature_coverage.tsv.gz",
        "manifest.json",
        "strategy_summary.tsv.gz",
    }
    manifest = json.loads((outdir / "manifest.json").read_text())
    assert manifest["output_profile"] == "report-input"
    assert manifest["ortholog_alignment_summary_count"] == 8
    assert manifest["alignment_segment_count"] == 10
    assert manifest["alignment_event_count"] == 12


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


def test_bwa_parameters_survive_partition_and_final_merge(tmp_path: Path) -> None:
    bwa_parameters = {
        "pseudoread_len": 150,
        "pseudoread_step": 75,
        "pseudoread_phred": 30,
        "task_cpus": 3,
        "bwa_threads": 2,
    }
    result_dir = write_result_dir(
        tmp_path,
        "gene_1_bwa",
        {
            "gene_id": "1",
            "strategy": "bwa_pseudoreads",
            **bwa_parameters,
        },
    )
    partition_dir = tmp_path / "partition"

    partition = run_merge(
        partition_arguments(
            [result_dir],
            partition_dir,
            gene_ids="1",
            strategies="bwa_pseudoreads",
        )
    )

    assert partition.returncode == 0, partition.stderr
    partition_manifest = json.loads((partition_dir / "manifest.json").read_text())
    assert partition_manifest["strategy_parameters"] == {
        "bwa_pseudoreads": bwa_parameters
    }

    inputs = write_final_inputs(tmp_path, [["1", "ready"]])
    final_dir = tmp_path / "final"
    final = run_merge(
        final_arguments(
            [partition_dir],
            final_dir,
            inputs,
            strategies="bwa_pseudoreads",
        )
    )

    assert final.returncode == 0, final.stderr
    final_manifest = json.loads((final_dir / "manifest.json").read_text())
    assert final_manifest["strategy_parameters"] == {
        "bwa_pseudoreads": bwa_parameters
    }


def test_partition_merge_rejects_inconsistent_bwa_parameters(tmp_path: Path) -> None:
    result_dirs = [
        write_result_dir(
            tmp_path,
            f"gene_{gene_id}_bwa",
            {
                "gene_id": gene_id,
                "strategy": "bwa_pseudoreads",
                "pseudoread_len": 150,
                "pseudoread_step": step,
                "pseudoread_phred": 30,
            },
        )
        for gene_id, step in [("1", 75), ("2", 100)]
    ]

    completed = run_merge(
        partition_arguments(
            result_dirs,
            tmp_path / "merged",
            gene_ids="1,2",
            strategies="bwa_pseudoreads",
        )
    )

    assert completed.returncode != 0
    assert "Inconsistent strategy parameters for bwa_pseudoreads" in completed.stderr
