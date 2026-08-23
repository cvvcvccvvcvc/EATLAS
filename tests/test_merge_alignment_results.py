from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest


PROJECT_DIR = Path(__file__).resolve().parents[1]

from bin.alignment_table_schema import (
    EVENT_FIELDS,
    FAILURE_FIELDS,
    SEGMENT_FIELDS,
    SUMMARY_FIELDS,
)
TABLE_HEADERS = {
    "ortholog_alignment_summary.tsv.gz": SUMMARY_FIELDS,
    "strategy_summary.tsv.gz": [
        "strategy",
        "summary_row_count",
        "gene_count",
        "aligned_summary_row_count",
        "event_count",
        "aligned_target_bp",
    ],
    "alignment_segments.tsv.gz": SEGMENT_FIELDS,
    "alignment_events.tsv.gz": EVENT_FIELDS,
    "snv_site_depth.tsv.gz": [
        "gene_id",
        "strategy",
        "target_start0",
        "site_aligned_ortholog_count",
    ],
    "failures.tsv.gz": FAILURE_FIELDS,
}
EVENT_ORTHOLOG_SUPPORT_HEADER = [
    "event_group_id",
    "ortholog_gene_id",
    "tax_id",
    "mapq",
    "native_alignment_type",
    "support_row_count",
]
COMPACT_EVENT_HEADER = [
    "event_group_id",
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
    "qc_flags",
]
def write_tsv_gz(path: Path, header: list[str], rows: list[list[str]] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows or [])


def schema_row(fields: list[str], **values: object) -> list[str]:
    return [str(values.get(field, "")) for field in fields]


def write_result_dir(
    root: Path,
    name: str,
    manifest: dict,
    *,
    missing_table: str | None = None,
) -> Path:
    result_dir = root / name
    result_dir.mkdir(parents=True)
    manifest = dict(manifest)
    gene_ids = manifest.pop("gene_ids", None)
    if gene_ids is None:
        gene_ids = [manifest.pop("gene_id")]
    else:
        manifest.pop("gene_id", None)
    strategies = manifest.pop("strategies", None)
    if strategies is None:
        strategies = [manifest.pop("strategy")]
    else:
        manifest.pop("strategy", None)
    manifest["gene_ids"] = gene_ids
    manifest["strategies"] = strategies
    manifest.setdefault(
        "strategy_parameters",
        {strategy: {} for strategy in strategies},
    )
    manifest.setdefault("ortholog_alignment_summary_count", len(gene_ids) * len(strategies))
    manifest.setdefault("alignment_segment_count", 0)
    manifest.setdefault("alignment_event_mode", "raw")
    manifest.setdefault("raw_alignment_event_count", manifest.get("alignment_event_count", 0))
    manifest.setdefault("alignment_event_count", 0)
    manifest.setdefault("event_ortholog_support_count", 0)
    manifest.setdefault("snv_site_depth_count", 0)
    manifest.setdefault("snv_taxonomic_depth_count", 0)
    manifest.setdefault("snv_alt_taxonomic_support_count", 0)
    manifest.setdefault("failure_count", 0)
    (result_dir / "manifest.json").write_text(json.dumps(manifest) + "\n")

    summary_rows = [
        schema_row(
            SUMMARY_FIELDS,
            gene_id=gene_id,
            strategy=strategy,
            status="aligned",
            event_count=0,
            aligned_target_bp=1,
        )
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


def write_compact_result_dir(root: Path, name: str, manifest: dict) -> Path:
    compact_manifest = {
        **manifest,
        "alignment_event_mode": "compact_support",
    }
    result_dir = write_result_dir(root, name, compact_manifest)
    write_tsv_gz(
        result_dir / "alignment_events.tsv.gz",
        COMPACT_EVENT_HEADER,
    )
    write_tsv_gz(
        result_dir / "event_ortholog_support.tsv.gz",
        EVENT_ORTHOLOG_SUPPORT_HEADER,
    )
    partition_manifest = json.loads((result_dir / "manifest.json").read_text())
    partition_manifest["schema"] = "normalized_alignment_evidence_partition_v2"
    for field in (
        "strategy_summary_count",
        "snv_site_depth_count",
        "snv_taxonomic_depth_count",
        "snv_alt_taxonomic_support_count",
    ):
        partition_manifest.pop(field, None)
    (result_dir / "manifest.json").write_text(json.dumps(partition_manifest) + "\n")
    for filename in (
        "strategy_summary.tsv.gz",
        "snv_site_depth.tsv.gz",
    ):
        (result_dir / filename).unlink(missing_ok=True)
    return result_dir


def run_merge(arguments: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "bin.merge_alignment_results", *arguments],
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
    if alignment_tasks is None:
        alignment_tasks = outdir.parent / f"{outdir.name}.alignment_tasks.tsv.gz"
        write_tsv_gz(
            alignment_tasks,
            ["gene_id", "status", "target_ready", "ortholog_ready"],
            [
                [gene_id, "ready", "true", "true"]
                for gene_id in gene_ids.split(",")
            ],
        )
    arguments = [
        "--partition-id",
        "partition_000001",
        "--expected-gene-ids",
        gene_ids,
        "--expected-strategies",
        strategies,
        "--outdir",
        str(outdir),
        "--alignment-tasks",
        str(alignment_tasks),
    ]
    for result_dir in result_dirs:
        arguments.extend(["--result-dir", str(result_dir)])
    return arguments


def write_final_inputs(
    root: Path,
    task_rows: list[list[str]],
    task_header: list[str] | None = None,
) -> tuple[Path, Path, Path]:
    alignment_tasks = root / "alignment_tasks.tsv.gz"
    source_genes = root / "genes.tsv.gz"
    source_target_features = root / "target_features.tsv.gz"
    if task_header is None:
        task_header = ["gene_id", "status", "target_ready", "ortholog_ready"]
        task_rows = [
            [*row, "true", "true"] if len(row) == 2 else row
            for row in task_rows
        ]
    write_tsv_gz(alignment_tasks, task_header, task_rows)
    write_tsv_gz(source_genes, ["gene_id"])
    write_tsv_gz(source_target_features, ["gene_id"])
    return (
        alignment_tasks,
        source_genes,
        source_target_features,
    )


def final_arguments(
    result_dirs: list[Path],
    outdir: Path,
    inputs: tuple[Path, Path, Path],
    *,
    strategies: str,
) -> list[str]:
    (
        alignment_tasks,
        source_genes,
        source_target_features,
    ) = inputs
    arguments = [
        "--alignment-tasks",
        str(alignment_tasks),
        "--source-genes",
        str(source_genes),
        "--source-target-features",
        str(source_target_features),
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


def test_partition_merge_rejects_legacy_singular_manifest(tmp_path: Path) -> None:
    result_dir = write_result_dir(
        tmp_path,
        "gene_1_s1",
        {"gene_id": "1", "strategy": "s1"},
    )
    (result_dir / "manifest.json").write_text(
        json.dumps({"gene_id": "1", "strategy": "s1"}) + "\n"
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
    assert "invalid gene_ids" in completed.stderr


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
    assert manifest["strategies"] == [ensembl, "s1"]


def test_partition_merge_always_publishes_compact_events(tmp_path: Path) -> None:
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
    completed = run_merge(arguments)

    assert completed.returncode == 0, completed.stderr
    manifest = json.loads((tmp_path / "merged" / "manifest.json").read_text())
    assert manifest["alignment_event_mode"] == "compact_support"
    assert manifest["alignment_event_count"] == 0


def test_partition_merge_publishes_only_normalized_evidence(tmp_path: Path) -> None:
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
            schema_row(
                SEGMENT_FIELDS,
                gene_id="1",
                strategy="s1",
                ortholog_gene_id="101",
                target_start0=0,
                target_end0=10,
            ),
            schema_row(
                SEGMENT_FIELDS,
                gene_id="1",
                strategy="s1",
                ortholog_gene_id="102",
                target_start0=0,
                target_end0=10,
            ),
        ],
    )
    write_tsv_gz(
        result_dir / "alignment_events.tsv.gz",
        TABLE_HEADERS["alignment_events.tsv.gz"],
        [
            schema_row(
                EVENT_FIELDS,
                gene_id="1",
                ortholog_gene_id="101",
                tax_id="1",
                taxname="species",
                strategy="s1",
                tool="tool",
                event_type="snv",
                target_start0=4,
                target_end0=5,
                genomic_accession="NC_1",
                genomic_start1=5,
                genomic_end1=5,
                ref="A",
                alt="G",
            ),
            schema_row(
                EVENT_FIELDS,
                gene_id="1",
                ortholog_gene_id="102",
                tax_id="2",
                taxname="species",
                strategy="s1",
                tool="tool",
                event_type="snv",
                target_start0=4,
                target_end0=5,
                genomic_accession="NC_1",
                genomic_start1=5,
                genomic_end1=5,
                ref="A",
                alt="T",
            ),
        ],
    )
    arguments = partition_arguments(
        [result_dir],
        tmp_path / "merged",
        gene_ids="1",
        strategies="s1",
    )
    completed = run_merge(arguments)

    assert completed.returncode == 0, completed.stderr
    outdir = tmp_path / "merged"
    assert {
        path.name
        for path in outdir.iterdir()
    } == {
        "alignment_segments.tsv.gz",
        "alignment_events.tsv.gz",
        "event_ortholog_support.tsv.gz",
        "failures.tsv.gz",
        "manifest.json",
        "ortholog_alignment_summary.tsv.gz",
    }
    manifest = json.loads((outdir / "manifest.json").read_text())
    assert manifest["schema"] == "normalized_alignment_evidence_partition_v2"
    assert manifest["alignment_segment_count"] == 2
    assert {
        "strategy_summary_count",
        "feature_coverage_count",
        "snv_site_depth_count",
        "snv_taxonomic_depth_count",
        "snv_alt_taxonomic_support_count",
    }.isdisjoint(manifest)


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
            [
                schema_row(
                    SEGMENT_FIELDS,
                    gene_id="1",
                    strategy=strategy,
                    ortholog_gene_id="101",
                    target_start0=0,
                    target_end0=1,
                )
            ],
        )
        write_tsv_gz(
            result_dir / "alignment_events.tsv.gz",
            event_header,
            [
                schema_row(
                    EVENT_FIELDS,
                    gene_id="1",
                    ortholog_gene_id="101",
                    tax_id="1",
                    taxname="species",
                    strategy=strategy,
                    tool="tool",
                    event_type="snv",
                    target_start0=0,
                    target_end0=1,
                    genomic_accession="NC_1",
                    genomic_start1=1,
                    genomic_end1=1,
                    ref="A",
                    alt="G",
                    mapq="20" if strategy == "s1" else "30",
                    native_alignment_type="P" if strategy == "s1" else "primary",
                )
            ],
        )

    arguments = partition_arguments(
        result_dirs,
        tmp_path / "merged",
        gene_ids="1",
        strategies="s1,s2",
    )
    completed = run_merge(arguments)

    assert completed.returncode == 0, completed.stderr
    with gzip.open(tmp_path / "merged" / "alignment_events.tsv.gz", "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == COMPACT_EVENT_HEADER
        rows = list(reader)
    assert [row["strategy"] for row in rows] == ["s1", "s2"]
    with gzip.open(
        tmp_path / "merged" / "event_ortholog_support.tsv.gz",
        "rt",
        newline="",
    ) as handle:
        ortholog_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["event_group_id"] for row in rows] == ["1", "2"]
    assert [row["event_group_id"] for row in ortholog_rows] == ["1", "2"]
    assert [row["ortholog_gene_id"] for row in ortholog_rows] == ["101", "101"]
    assert [row["mapq"] for row in ortholog_rows] == ["20", "30"]
    assert [row["native_alignment_type"] for row in ortholog_rows] == ["P", "primary"]
    assert [row["support_row_count"] for row in ortholog_rows] == ["1", "1"]
    manifest = json.loads((tmp_path / "merged" / "manifest.json").read_text())
    assert set(manifest["timings_seconds"]) >= {
        "load_events_sqlite",
        "build_event_index",
        "stream_event_groups",
    }


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
            schema_row(
                EVENT_FIELDS,
                gene_id="3492",
                ortholog_gene_id="pongo_abelii",
                tax_id="9601",
                taxname="Pongo abelii",
                strategy=strategy,
                tool="ensembl_compara",
                event_type="ins",
                target_start0=431282,
                target_end0=431282,
                genomic_accession="NC_000014.9",
                genomic_start1=106017719,
                genomic_end1=106017719,
                alt=large_alt,
            )
        ],
    )

    completed = run_merge(
        partition_arguments(
            [result_dir],
            tmp_path / "merged",
            gene_ids="3492",
            strategies=strategy,
        )
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


@pytest.mark.parametrize(
    ("filename", "invalid_header", "expected_fragment"),
    [
        pytest.param(
            "alignment_events.tsv.gz",
            [field for field in EVENT_FIELDS if field != "ortholog_gene_id"],
            "ortholog_gene_id",
            id="events-missing-ortholog-gene-id",
        ),
        pytest.param(
            "alignment_events.tsv.gz",
            [field for field in EVENT_FIELDS if field != "tax_id"],
            "tax_id",
            id="events-missing-tax-id",
        ),
        pytest.param(
            "alignment_events.tsv.gz",
            [field for field in EVENT_FIELDS if field != "strategy"],
            "strategy",
            id="events-missing-strategy",
        ),
        pytest.param(
            "alignment_segments.tsv.gz",
            [SEGMENT_FIELDS[1], SEGMENT_FIELDS[0], *SEGMENT_FIELDS[2:]],
            "observed",
            id="segments-reordered",
        ),
        pytest.param(
            "ortholog_alignment_summary.tsv.gz",
            [*SUMMARY_FIELDS, "unexpected"],
            "unexpected",
            id="summary-extra-field",
        ),
        pytest.param(
            "failures.tsv.gz",
            [field for field in FAILURE_FIELDS if field != "message"],
            "message",
            id="failures-missing-message",
        ),
    ],
)
def test_partition_merge_rejects_noncanonical_aligner_header(
    tmp_path: Path,
    filename: str,
    invalid_header: list[str],
    expected_fragment: str,
) -> None:
    result_dir = write_result_dir(
        tmp_path,
        "gene_1_s1",
        {"gene_id": "1", "strategy": "s1"},
    )
    write_tsv_gz(result_dir / filename, invalid_header)

    completed = run_merge(
        partition_arguments(
            [result_dir],
            tmp_path / "merged",
            gene_ids="1",
            strategies="s1",
        )
    )

    assert completed.returncode != 0
    assert (
        f"Alignment table {result_dir / filename} has invalid header"
        in completed.stderr
    )
    assert "expected" in completed.stderr
    assert "observed" in completed.stderr
    assert expected_fragment in completed.stderr


@pytest.mark.parametrize(
    ("filename", "invalid_header", "expected_fragment"),
    [
        pytest.param(
            "alignment_events.tsv.gz",
            [
                COMPACT_EVENT_HEADER[1],
                COMPACT_EVENT_HEADER[0],
                *COMPACT_EVENT_HEADER[2:],
            ],
            "event_group_id",
            id="compact-events-reordered",
        ),
        pytest.param(
            "event_ortholog_support.tsv.gz",
            [
                field
                for field in EVENT_ORTHOLOG_SUPPORT_HEADER
                if field != "ortholog_gene_id"
            ],
            "ortholog_gene_id",
            id="compact-support-missing-ortholog-gene-id",
        ),
    ],
)
def test_final_merge_rejects_noncanonical_compact_event_header(
    tmp_path: Path,
    filename: str,
    invalid_header: list[str],
    expected_fragment: str,
) -> None:
    partition = write_compact_result_dir(
        tmp_path,
        "partition_000001",
        {
            "partition_id": "partition_000001",
            "gene_count": 1,
            "gene_ids": ["1"],
            "strategies": ["s1"],
        },
    )
    write_tsv_gz(partition / filename, invalid_header)
    inputs = write_final_inputs(tmp_path, [["1", "ready"]])

    completed = run_merge(
        final_arguments(
            [partition],
            tmp_path / "merged",
            inputs,
            strategies="s1",
        )
    )

    assert completed.returncode != 0
    assert (
        f"Alignment table {partition / filename} has invalid header"
        in completed.stderr
    )
    assert "expected" in completed.stderr
    assert "observed" in completed.stderr
    assert expected_fragment in completed.stderr


def test_final_merge_writes_exact_gene_manifest(tmp_path: Path) -> None:
    partition_dirs = [
        write_compact_result_dir(
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
        write_compact_result_dir(
            tmp_path,
            "partition_000001",
            {
                "partition_id": "partition_000001",
                "gene_count": 1,
                "gene_ids": ["1"],
                "strategies": ["s1", ensembl],
            },
        ),
        write_compact_result_dir(
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
    assert manifest["strategies"] == [ensembl, "s1"]


def test_final_merge_preserves_precompacted_ortholog_support(tmp_path: Path) -> None:
    partition_dirs = []
    for index, (gene_id, ortholog_gene_id) in enumerate([("1", "101"), ("2", "201")], 1):
        partition = write_compact_result_dir(
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
            partition / "alignment_events.tsv.gz",
            COMPACT_EVENT_HEADER,
            [
                schema_row(
                    COMPACT_EVENT_HEADER,
                    event_group_id="1",
                    gene_id=gene_id,
                    event_type="snv",
                    target_start0="0",
                    target_end0="1",
                    genomic_accession="NC_1",
                    genomic_start1="1",
                    genomic_end1="1",
                    ref="A",
                    alt="G",
                    strategy="s1",
                )
            ],
        )
        write_tsv_gz(
            partition / "event_ortholog_support.tsv.gz",
            EVENT_ORTHOLOG_SUPPORT_HEADER,
            [
                [
                    "1",
                    ortholog_gene_id,
                    "10090",
                    "",
                    "",
                    "1",
                ]
            ],
        )
        partition_dirs.append(partition)
    inputs = write_final_inputs(tmp_path, [["1", "ready"], ["2", "ready"]])
    outdir = tmp_path / "merged"
    completed = run_merge(final_arguments(partition_dirs, outdir, inputs, strategies="s1"))

    assert completed.returncode == 0, completed.stderr
    copied_root = outdir / "evidence" / "partitions"
    copied_rows = []
    for partition in partition_dirs:
        with gzip.open(
            copied_root / partition.name / "event_ortholog_support.tsv.gz",
            "rt",
            newline="",
        ) as handle:
            copied_rows.append(list(csv.DictReader(handle, delimiter="\t")))
    assert [rows[0]["event_group_id"] for rows in copied_rows] == ["1", "1"]
    assert [rows[0]["ortholog_gene_id"] for rows in copied_rows] == ["101", "201"]
    manifest = json.loads((outdir / "manifest.json").read_text())
    assert manifest["event_ortholog_support_count"] == 2


def test_final_merge_publishes_partitioned_normalized_evidence(
    tmp_path: Path,
) -> None:
    partition_dirs = [
        write_compact_result_dir(
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
            },
        )
        for index, gene_id in enumerate(["1", "2"], start=1)
    ]
    for partition_dir, gene_id in zip(partition_dirs, ["1", "2"]):
        write_tsv_gz(
            partition_dir / "alignment_events.tsv.gz",
            COMPACT_EVENT_HEADER,
            [
                [
                    "1",
                    gene_id,
                    "snv",
                    "4",
                    "5",
                    f"NC_{gene_id}",
                    "5",
                    "5",
                    "A",
                    "G",
                    "s1",
                    "",
                ]
            ],
        )
        write_tsv_gz(
            partition_dir / "event_ortholog_support.tsv.gz",
            EVENT_ORTHOLOG_SUPPORT_HEADER,
            [
                [
                    "1",
                    f"ortholog_{gene_id}",
                    gene_id,
                    "60",
                    "P",
                    "1",
                ]
            ],
        )
    evidence_files = {
        "manifest.json",
        "ortholog_alignment_summary.tsv.gz",
        "alignment_segments.tsv.gz",
        "alignment_events.tsv.gz",
        "event_ortholog_support.tsv.gz",
    }
    source_bytes = {
        (partition_dir.name, path.name): path.read_bytes()
        for partition_dir in partition_dirs
        for path in partition_dir.iterdir()
        if path.name in evidence_files
    }
    inputs = write_final_inputs(tmp_path, [["1", "ready"], ["2", "ready"]])
    outdir = tmp_path / "merged"
    arguments = final_arguments(
        partition_dirs,
        outdir,
        inputs,
        strategies="s1",
    )
    completed = run_merge(arguments)

    assert completed.returncode == 0, completed.stderr
    assert {
        path.name
        for path in outdir.iterdir()
    } == {
        "evidence",
        "failures.tsv.gz",
        "manifest.json",
    }
    manifest = json.loads((outdir / "manifest.json").read_text())
    assert manifest["schema"] == "normalized_alignment_evidence_v2"
    assert manifest["ortholog_alignment_summary_count"] == 8
    assert manifest["alignment_segment_count"] == 10
    assert manifest["alignment_event_count"] == 12
    assert {
        "strategy_summary_count",
        "feature_coverage_count",
        "snv_site_depth_count",
        "snv_taxonomic_depth_count",
        "snv_alt_taxonomic_support_count",
    }.isdisjoint(manifest)
    assert manifest["normalized_evidence"] == {
        "layout": "partitioned",
        "format": "tsv_gzip_v1",
        "path": "evidence/partitions",
        "partition_count": 2,
        "partition_files": [
            "manifest.json",
            "ortholog_alignment_summary.tsv.gz",
            "alignment_segments.tsv.gz",
            "alignment_events.tsv.gz",
            "event_ortholog_support.tsv.gz",
        ],
        "event_group_id_scope": "partition",
    }
    evidence_partitions = outdir / "evidence" / "partitions"
    assert {path.name for path in evidence_partitions.iterdir()} == {
        "partition_000001",
        "partition_000002",
    }
    for partition_dir in partition_dirs:
        copied_dir = evidence_partitions / partition_dir.name
        assert {path.name for path in copied_dir.iterdir()} == evidence_files
        for copied_path in copied_dir.iterdir():
            assert copied_path.read_bytes() == source_bytes[
                (partition_dir.name, copied_path.name)
            ]
        with gzip.open(copied_dir / "alignment_events.tsv.gz", "rt", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert [row["event_group_id"] for row in rows] == ["1"]


def test_final_merge_rejects_raw_partition_results(tmp_path: Path) -> None:
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
    inputs = write_final_inputs(tmp_path, [["1", "ready"]])

    completed = run_merge(
        final_arguments(
            [partition_dir],
            tmp_path / "merged",
            inputs,
            strategies="s1",
        )
    )

    assert completed.returncode != 0
    assert "Final alignment merge requires compact_support inputs" in completed.stderr


def test_final_merge_rejects_missing_ready_gene(tmp_path: Path) -> None:
    partition_dir = write_compact_result_dir(
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
    }
    result_dir = write_result_dir(
        tmp_path,
        "gene_1_bwa",
        {
            "gene_id": "1",
            "strategy": "bwa_pseudoreads_150_75",
            "strategy_parameters": {"bwa_pseudoreads_150_75": bwa_parameters},
        },
    )
    partition_dir = tmp_path / "partition"

    partition = run_merge(
        partition_arguments(
            [result_dir],
            partition_dir,
            gene_ids="1",
            strategies="bwa_pseudoreads_150_75",
        )
    )

    assert partition.returncode == 0, partition.stderr
    partition_manifest = json.loads((partition_dir / "manifest.json").read_text())
    assert partition_manifest["strategy_parameters"] == {
        "bwa_pseudoreads_150_75": bwa_parameters
    }

    inputs = write_final_inputs(tmp_path, [["1", "ready"]])
    final_dir = tmp_path / "final"
    final = run_merge(
        final_arguments(
            [partition_dir],
            final_dir,
            inputs,
            strategies="bwa_pseudoreads_150_75",
        )
    )

    assert final.returncode == 0, final.stderr
    final_manifest = json.loads((final_dir / "manifest.json").read_text())
    assert final_manifest["strategy_parameters"] == {
        "bwa_pseudoreads_150_75": bwa_parameters
    }


def test_partition_merge_rejects_inconsistent_bwa_parameters(tmp_path: Path) -> None:
    result_dirs = [
        write_result_dir(
            tmp_path,
            f"gene_{gene_id}_bwa",
            {
                "gene_id": gene_id,
                "strategy": "bwa_pseudoreads_150_75",
                "strategy_parameters": {
                    "bwa_pseudoreads_150_75": {
                        "pseudoread_len": 150,
                        "pseudoread_step": step,
                        "pseudoread_phred": 30,
                    }
                },
            },
        )
        for gene_id, step in [("1", 75), ("2", 100)]
    ]

    completed = run_merge(
        partition_arguments(
            result_dirs,
            tmp_path / "merged",
            gene_ids="1,2",
            strategies="bwa_pseudoreads_150_75",
        )
    )

    assert completed.returncode != 0
    assert "Inconsistent strategy parameters for bwa_pseudoreads_150_75" in completed.stderr
