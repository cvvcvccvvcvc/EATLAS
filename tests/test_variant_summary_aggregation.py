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
from genomics.clinvar import record_category
from genomics.variants import ALLELE_ANNOTATION_FIELDS


CORE_COLUMNS = ["variant_key", "gene_id", "event_type", "ref", "alt", "strategies"]
COLUMNS = [*CORE_COLUMNS, "vep_status", "vep_primary_consequence"]
EVIDENCE_COLUMNS = [
    *CORE_COLUMNS,
    "lookup_status",
    *ALLELE_ANNOTATION_FIELDS,
    "vep_status",
    "vep_primary_consequence",
]


def _evidence_row(gene_id: str, **overrides: str) -> dict[str, str]:
    row = dict.fromkeys(EVIDENCE_COLUMNS, "")
    row.update(
        {
            "variant_key": "1:100:A>G",
            "gene_id": gene_id,
            "event_type": "snv",
            "ref": "A",
            "alt": "G",
            "strategies": "s1",
            "lookup_status": "ok",
            "vep_status": "ok",
            "vep_primary_consequence": "missense_variant",
        }
    )
    row.update(overrides)
    return row


def _write_evidence_rows(path: Path, rows: list[dict[str, str]]) -> None:
    pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).to_csv(
        path,
        sep="\t",
        index=False,
        compression="gzip",
    )


def test_available_cpu_count_prefers_slurm_allocation(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "12")
    assert available_cpu_count() == 12


def test_variant_source_requires_vep_columns(tmp_path: Path) -> None:
    path = tmp_path / "variant_annotations.tsv.gz"
    pd.DataFrame(
        [["1:100:A>G", "gene_a", "snv", "A", "G", "s1"]],
        columns=CORE_COLUMNS,
    ).to_csv(path, sep="\t", index=False, compression="gzip")

    with pytest.raises(ValueError, match="vep_primary_consequence, vep_status"):
        resolve_variant_aggregation_source(path)


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
            ["1:100:A>G", "gene_a", "snv", "A", "G", "s1,s2", "ok", "missense_variant"],
            ["1:100:A>G", "gene_a", "snv", "A", "G", "s1, s2", "ok", "missense_variant"],
            ["1:100:A>G", "gene_b", "snv", "A", "G", "s1,s3", "ok", "intron_variant"],
            ["1:200:C>T", "gene_a", "snv", "C", "T", "s2", "ok", "synonymous_variant"],
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


def test_partitioned_source_validates_pipeline_dataset(tmp_path: Path) -> None:
    artifact = tmp_path / "variant_annotations"
    partitions = artifact / "partitions"
    rows = [
        ["1:100:A>G", "gene_a", "snv", "A", "G", "s1,s2", "ok", "missense_variant"],
        ["1:200:C>T", "gene_a", "snv", "C", "T", "s2", "ok", "synonymous_variant"],
    ]
    partition_entries = []
    for index, row in enumerate(rows, start=1):
        partition_id = f"partition_{index:06d}"
        path = partitions / partition_id / "shard_000001.tsv.gz"
        path.parent.mkdir(parents=True)
        pd.DataFrame([row], columns=COLUMNS).to_csv(
            path,
            sep="\t",
            index=False,
            compression="gzip",
        )
        partition_entries.append({
            "partition_id": partition_id,
            "shard_count": 1,
            "row_count": 1,
            "shards": [{
                "shard_id": "shard_000001",
                "path": f"partitions/{partition_id}/shard_000001.tsv.gz",
                "row_count": 1,
                "size_bytes": path.stat().st_size,
            }],
        })
    manifest_path = artifact / "manifest.json"
    manifest_path.write_text(json.dumps({
        "schema": "gaph_variant_annotation_dataset_v1",
        "status": "complete",
        "layout": "partitioned",
        "format": "tsv_gzip_v1",
        "partition_count": 2,
        "shard_count": 2,
        "row_count": 2,
        "fields": COLUMNS,
        "partitions": partition_entries,
    }))

    source = resolve_variant_aggregation_source(manifest_path)
    result = aggregate_strategy_masks(source, threads=2)

    assert source.partitioned
    assert source.row_count == 2
    assert result.strategy_counts() == {"s1": 1, "s2": 2}

    payload = json.loads(manifest_path.read_text())
    payload["partitions"][0]["row_count"] = 2
    manifest_path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="row count changed"):
        resolve_variant_aggregation_source(manifest_path)


def test_variant_groups_keep_gene_specific_context_and_consequence(tmp_path: Path) -> None:
    path = tmp_path / "variant_annotations.tsv.gz"
    genes = tmp_path / "genes.tsv.gz"
    features = tmp_path / "target_features.tsv.gz"
    columns = [
        *CORE_COLUMNS,
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

    ortholog_evidence = tmp_path / "ortholog_evidence_summary.tsv.gz"
    pd.DataFrame(
        columns=[
            "strategy",
            "target_context",
            "taxonomic_scope",
            "evidence_unit",
            "site_aligned_count",
            "alt_support_count",
            "gnomad_found_count",
            "gnomad_not_found_count",
            "gnomad_lookup_failed_count",
        ]
    ).to_csv(ortholog_evidence, sep="\t", index=False, compression="gzip")

    duckdb_summary = _summary_from_grouped_aggregation(
        result,
        str,
        None,
        ortholog_evidence,
    )
    built_summary = build_variant_summary(
        path,
        tmp_path / "summary_work",
        str,
        ortholog_evidence_summary_path=ortholog_evidence,
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
        "gnomad_consequence_counts",
        "gnomad_af_summary",
    ):
        left = getattr(duckdb_summary, name).sort_values(
            list(getattr(duckdb_summary, name).columns)
        ).reset_index(drop=True)
        right = getattr(built_summary, name).sort_values(
            list(getattr(built_summary, name).columns)
        ).reset_index(drop=True)
        pd.testing.assert_frame_equal(left, right, check_dtype=False)


def test_allele_evidence_is_reconciled_independently_of_row_order(
    tmp_path: Path,
) -> None:
    missing = _evidence_row(
        "gene_a",
        lookup_status="failed",
        vep_primary_consequence="intron_variant",
    )
    found = _evidence_row(
        "gene_b",
        clinvar_id="VCV1",
        clinvar_sig="Pathogenic",
        clinvar_review_stars="3",
        clinvar_scv_count="2",
        gnomad_af="0.01",
        gnomad_af_source="exomes",
        gnomad_csq="missense_variant",
    )
    results = []
    for index, rows in enumerate(([missing, found], [found, missing])):
        path = tmp_path / f"variants_{index}.tsv.gz"
        _write_evidence_rows(path, rows)
        result = aggregate_variant_groups(
            resolve_variant_aggregation_source(path),
            threads=1,
        )
        results.append(result)

        global_row = result.global_groups.iloc[0]
        assert bool(global_row["clinvar_found"])
        assert global_row["clinvar_category"] == "P/LP"
        assert global_row["review_stars"] == "3"
        assert global_row["gnomad_status"] == "found"

        contexts = result.allele_gene_groups.set_index("gene_id")
        assert set(contexts["clinvar_category"]) == {"P/LP"}
        assert contexts.loc["gene_a", "gnomad_status"] == "lookup_failed"
        assert contexts.loc["gene_b", "gnomad_status"] == "found"

        pathogenic = result.pathogenic_rows.iloc[0]
        assert pathogenic["gene_id"] == "gene_a, gene_b"
        assert pathogenic["lookup_status"] == "failed|ok"
        assert pathogenic["clinvar_sig"] == "Pathogenic"
        assert pathogenic["vep_primary_consequence"] == (
            "intron_variant|missense_variant"
        )

    pd.testing.assert_frame_equal(results[0].global_groups, results[1].global_groups)
    pd.testing.assert_frame_equal(
        results[0].allele_gene_groups,
        results[1].allele_gene_groups,
    )
    pd.testing.assert_frame_equal(results[0].pathogenic_rows, results[1].pathogenic_rows)


def test_duckdb_clinvar_categories_match_shared_semantics(tmp_path: Path) -> None:
    signatures = [
        ("", False),
        ("", True),
        ("Likely_benign", True),
        ("Likely_pathogenic", True),
        ("Uncertain_significance", True),
        ("Conflicting_classifications_of_pathogenicity", True),
    ]
    rows = [
        _evidence_row(
            f"gene_{index}",
            variant_key=f"1:{100 + index}:A>G",
            clinvar_id=f"VCV{index}" if found else "",
            clinvar_sig=significance,
        )
        for index, (significance, found) in enumerate(signatures)
    ]
    path = tmp_path / "variant_annotations.tsv.gz"
    _write_evidence_rows(path, rows)

    result = aggregate_variant_groups(
        resolve_variant_aggregation_source(path),
        threads=1,
    )
    observed = result.global_groups.set_index("clinvar_category")[
        "variant_count"
    ].to_dict()
    expected = pd.Series(
        [
            record_category(significance, found=found)
            for significance, found in signatures
        ]
    ).value_counts().to_dict()

    assert observed == expected


def test_same_gene_gnomad_status_reduction_is_order_independent(tmp_path: Path) -> None:
    found = _evidence_row("gene_a", gnomad_af="0.01")
    failed = _evidence_row("gene_a", lookup_status="failed")
    statuses = []
    for index, rows in enumerate(([found, failed], [failed, found])):
        path = tmp_path / f"same_gene_{index}.tsv.gz"
        _write_evidence_rows(path, rows)
        result = aggregate_variant_groups(
            resolve_variant_aggregation_source(path),
            threads=1,
        )
        statuses.append(result.allele_gene_groups.iloc[0]["gnomad_status"])

    assert statuses == ["found", "found"]


@pytest.mark.parametrize(
    ("field", "left", "right"),
    [
        ("clinvar_sig", "Pathogenic", "Benign"),
        ("gnomad_af", "0.01", "0.02"),
    ],
)
def test_allele_evidence_rejects_conflicting_nonempty_values(
    tmp_path: Path,
    field: str,
    left: str,
    right: str,
) -> None:
    path = tmp_path / f"conflicting_{field}.tsv.gz"
    _write_evidence_rows(
        path,
        [
            _evidence_row("gene_a", **{field: left}),
            _evidence_row("gene_b", **{field: right}),
        ],
    )

    with pytest.raises(ValueError, match=rf"1:100:A>G \({field}\)"):
        aggregate_variant_groups(
            resolve_variant_aggregation_source(path),
            threads=1,
        )


def test_allele_evidence_rejects_non_numeric_gnomad_af(tmp_path: Path) -> None:
    path = tmp_path / "invalid_af.tsv.gz"
    _write_evidence_rows(path, [_evidence_row("gene_a", gnomad_af="invalid")])

    with pytest.raises(ValueError, match=r"1:100:A>G \(gnomad_af invalid\)"):
        aggregate_variant_groups(
            resolve_variant_aggregation_source(path),
            threads=1,
        )


def test_allele_evidence_accepts_equivalent_numeric_af_values(tmp_path: Path) -> None:
    path = tmp_path / "equivalent_af.tsv.gz"
    _write_evidence_rows(
        path,
        [
            _evidence_row("gene_a", gnomad_af="0.1"),
            _evidence_row("gene_b", gnomad_af="0.10"),
        ],
    )

    result = aggregate_variant_groups(
        resolve_variant_aggregation_source(path),
        threads=1,
    )

    assert result.gnomad_af_summary.iloc[0]["Median gnomAD AF"] == 0.1


def test_duckdb_memory_override_takes_precedence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "variant_annotations.tsv.gz"
    pd.DataFrame(
        [["1:100:A>G", "gene_a", "snv", "A", "G", "s1", "ok", "missense_variant"]],
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
