from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from analytics.analyses.variant_summary import (
    _summary_from_grouped_aggregation,
    build_variant_summary,
)
from analytics.analyses.variant_summary_aggregation import (
    DUCKDB_MEMORY_LIMIT_ENV,
    _configure_duckdb_memory,
    _parse_memory_setting,
    _slurm_memory_bytes,
    aggregate_strategy_masks,
    aggregate_variant_groups,
    available_cpu_count,
    resolve_variant_aggregation_source,
)
from analytics.io.artifacts import file_identity


COLUMNS = ["variant_key", "gene_id", "event_type", "ref", "alt", "strategies"]


def test_available_cpu_count_prefers_slurm_allocation(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "12")
    assert available_cpu_count() == 12


def test_slurm_memory_budget_uses_node_or_per_cpu_allocation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "16000")
    assert _slurm_memory_bytes(4) == (16000 * 1024**2, "SLURM_MEM_PER_NODE")

    monkeypatch.delenv("SLURM_MEM_PER_NODE")
    monkeypatch.setenv("SLURM_MEM_PER_CPU", "2000")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "6")
    assert _slurm_memory_bytes(4) == (12000 * 1024**2, "SLURM_MEM_PER_CPU")


def test_duckdb_memory_uses_half_the_slurm_allocation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    duckdb = pytest.importorskip("duckdb")
    monkeypatch.delenv(DUCKDB_MEMORY_LIMIT_ENV, raising=False)
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "2048")
    connection = duckdb.connect()
    try:
        diagnostics = _configure_duckdb_memory(connection, thread_count=4)
    finally:
        connection.close()

    assert diagnostics["memory_limit_source"] == "SLURM_MEM_PER_NODE"
    observed = _parse_memory_setting(str(diagnostics["memory_limit"]))
    assert 1000 * 1024**2 <= observed <= 1050 * 1024**2


def test_strategy_masks_preserve_gene_context_and_union_global_alleles(tmp_path: Path) -> None:
    path = tmp_path / "variant_annotations.tsv.gz"
    pd.DataFrame(
        [
            ["1:100:A>G", "gene_a", "snv", "A", "G", "s1,s2"],
            ["1:100:A>G", "gene_a", "snv", "A", "G", "s1, s2"],
            ["1:100:A>G", "gene_b", "snv", "A", "G", "s1,s3"],
            ["1:200:C>T", "gene_a", "snv", "C", "T", "s2"],
        ],
        columns=COLUMNS,
    ).to_csv(path, sep="\t", index=False, compression="gzip")

    source = resolve_variant_aggregation_source(path)
    result = aggregate_strategy_masks(source, threads=2)

    assert not source.partitioned
    assert result.input_row_count == 4
    assert result.strategies == ("s1", "s2", "s3")
    assert sum(result.allele_gene_mask_counts.values()) == 3
    assert result.unique_variant_count == 2
    assert result.strategy_counts() == {"s1": 1, "s2": 2, "s3": 1}
    assert result.strategy_record_count == 4
    assert result.unique_strategy_counts() == {"s1": 0, "s2": 1, "s3": 0}
    assert result.all_strategy_variant_count == 1
    assert result.intersections() == [[1, 1, 1], [1, 2, 1], [1, 1, 1]]

    with pytest.raises(ValueError, match="thread count"):
        aggregate_strategy_masks(source, threads=0)


def test_partitioned_source_validates_manifests_and_matches_legacy(tmp_path: Path) -> None:
    artifact = tmp_path / "vep_consequences"
    partitions = artifact / "partitions"
    partitions.mkdir(parents=True)
    rows = [
        ["1:100:A>G", "gene_a", "snv", "A", "G", "s1,s2"],
        ["1:200:C>T", "gene_a", "snv", "C", "T", "s2"],
    ]
    entries = []
    manifest_files = []
    for index, row in enumerate(rows, start=1):
        partition_id = f"partition_{index:06d}"
        path = partitions / f"{partition_id}.tsv.gz"
        pd.DataFrame([row], columns=COLUMNS).to_csv(
            path,
            sep="\t",
            index=False,
            header=False,
            compression="gzip",
        )
        entry = {
            "partition_id": partition_id,
            "path": f"inputs/{partition_id}.tsv.gz",
            "row_count": 1,
            "file": {"size_bytes": 1, "mtime_ns": 1},
        }
        entries.append(entry)
        partition_manifest = {
            "status": "complete",
            "input": entry,
            "row_count": 1,
            "output_columns": COLUMNS,
            "output": file_identity(path),
        }
        manifest_path = partitions / f"{partition_id}.json"
        manifest_path.write_text(json.dumps(partition_manifest))
        manifest_files.append({"partition_id": partition_id, "file": file_identity(manifest_path)})

    plan = {
        "status": "complete",
        "row_count": 2,
        "output_columns": COLUMNS,
        "partitions": entries,
    }
    (artifact / "plan.json").write_text(json.dumps(plan))
    merged = artifact / "variant_annotations.vep.tsv.gz"
    with gzip.open(merged, "wt") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
    (artifact / "manifest.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "row_count": 2,
                "columns": COLUMNS,
                "partition_manifests": manifest_files,
                "output": file_identity(merged),
            }
        )
    )

    source = resolve_variant_aggregation_source(merged)
    result = aggregate_strategy_masks(source, threads=2)

    assert source.partitioned
    assert source.row_count == 2
    assert result.strategy_counts() == {"s1": 1, "s2": 2}

    first_manifest = partitions / "partition_000001.json"
    payload = json.loads(first_manifest.read_text())
    payload["row_count"] = 2
    first_manifest.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="row count changed"):
        resolve_variant_aggregation_source(merged)


def test_variant_groups_keep_gene_specific_context_and_consequence(tmp_path: Path) -> None:
    path = tmp_path / "variant_annotations.tsv.gz"
    genes = tmp_path / "genes.tsv.gz"
    features = tmp_path / "target_features.tsv.gz"
    columns = [
        *COLUMNS,
        "lookup_status",
        "clinvar_id",
        "clinvar_sig",
        "clinvar_review_stars",
        "clinvar_scv_count",
        "gnomad_af",
        "gnomad_csq",
        "vep_status",
        "vep_primary_consequence",
    ]
    pd.DataFrame(
        [
            ["1:100:A>G", "gene_a", "snv", "A", "G", "s1,s2", "ok", "VCV1", "Benign", "2", "3", "0.01", "", "ok", "missense_variant"],
            ["1:100:A>G", "gene_b", "snv", "A", "G", "s1,s3", "failed", "VCV1", "Benign", "2", "3", "", "", "ok", "intron_variant"],
            ["1:200:C>T", "gene_a", "snv", "C", "T", "s2", "ok", "", "", "", "", "0", "", "ok", "synonymous_variant"],
        ],
        columns=columns,
    ).to_csv(path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(
        [["gene_a", 1], ["gene_b", 1]], columns=["gene_id", "begin"]
    ).to_csv(genes, sep="\t", index=False, compression="gzip")
    pd.DataFrame(
        [
            ["gene_a", "gene", 0, 300],
            ["gene_a", "cds", 0, 150],
            ["gene_a", "intron", 150, 300],
            ["gene_b", "gene", 0, 300],
            ["gene_b", "intron", 0, 300],
        ],
        columns=["gene_id", "feature_type", "target_start0", "target_end0"],
    ).to_csv(features, sep="\t", index=False, compression="gzip")

    result = aggregate_variant_groups(
        resolve_variant_aggregation_source(path),
        genes_path=genes,
        target_features_path=features,
        threads=2,
    )

    assert result.masks.unique_variant_count == 2
    assert int(result.global_groups["variant_count"].sum()) == 2
    assert int(result.allele_gene_groups["variant_count"].sum()) == 3
    assert result.diagnostics["allele_gene_row_count"] == 3
    assert result.diagnostics["global_allele_row_count"] == 2
    assert _parse_memory_setting(str(result.diagnostics["memory_limit"])) > 0
    assert result.timings["compacted_relations"] >= 0
    observed = set(
        result.allele_gene_groups[
            ["gene_id", "target_context", "consequence"]
        ].itertuples(index=False, name=None)
    )
    assert ("gene_a", "cds", "missense_variant") in observed
    assert ("gene_b", "intron", "intron_variant") in observed
    assert int(
        result.global_groups.loc[
            result.global_groups["gnomad_status"].eq("lookup_failed"),
            "variant_count",
        ].sum()
    ) == 0
    assert int(
        result.allele_gene_groups.loc[
            result.allele_gene_groups["gnomad_status"].eq("lookup_failed"),
            "variant_count",
        ].sum()
    ) == 1
    assert result.gnomad_af_summary.iloc[0]["Median gnomAD AF"] == 0.01
    s2_af = result.gnomad_af_summary.set_index("strategy").loc["s2"]
    assert s2_af["Count"] == 1
    assert s2_af["Median gnomAD AF"] == 0.005

    duckdb_summary = _summary_from_grouped_aggregation(result, str, None, None)
    built_summary = build_variant_summary(
        path,
        tmp_path / "legacy_work",
        str,
        target_features_path=features,
        genes_path=genes,
    )
    assert duckdb_summary.unique_variant_count == built_summary.unique_variant_count
    assert duckdb_summary.strategy_record_count == built_summary.strategy_record_count
    assert duckdb_summary.strategy_stats.set_index("Strategy").loc["s1", "Genes"] == 2
    pd.testing.assert_frame_equal(
        duckdb_summary.strategy_stats.drop(columns="Genes"),
        built_summary.strategy_stats.drop(columns="Genes"),
        check_dtype=False,
    )
    for name in (
        "event_counts",
        "clinvar_counts",
        "gnomad_event_counts",
        "gnomad_af_summary",
    ):
        left = getattr(duckdb_summary, name).sort_values(
            list(getattr(duckdb_summary, name).columns)
        ).reset_index(drop=True)
        right = getattr(built_summary, name).sort_values(
            list(getattr(built_summary, name).columns)
        ).reset_index(drop=True)
        pd.testing.assert_frame_equal(left, right, check_dtype=False)


def test_duckdb_memory_override_takes_precedence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "variant_annotations.tsv.gz"
    pd.DataFrame(
        [["1:100:A>G", "gene_a", "snv", "A", "G", "s1"]],
        columns=COLUMNS,
    ).to_csv(path, sep="\t", index=False, compression="gzip")
    monkeypatch.setenv(DUCKDB_MEMORY_LIMIT_ENV, "512MB")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "16000")

    result = aggregate_variant_groups(
        resolve_variant_aggregation_source(path),
        threads=1,
        temp_dir=tmp_path / "spill",
    )

    assert result.diagnostics["memory_limit_source"] == DUCKDB_MEMORY_LIMIT_ENV
    observed = _parse_memory_setting(str(result.diagnostics["memory_limit"]))
    assert 450 * 1024**2 <= observed <= 550 * 1024**2
