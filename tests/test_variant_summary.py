from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pandas as pd

from analytics.annotation.consequences import UNANNOTATED_CONSEQUENCE
from analytics.analyses.variant_summary import (
    VARIANT_USECOLS,
    _categorize_clinvar,
    build_variant_summary,
    read_taxonomic_ortholog_evidence,
)
from analytics.io.performance import PerformanceProfile


def test_clinvar_categories_distinguish_unclassified_records_from_absence() -> None:
    categories = _categorize_clinvar(
        pd.Series(["", "", "Benign"]),
        pd.Series(["", "VCV1", "VCV2"]),
    )

    assert list(categories.astype(str)) == ["Not in ClinVar", "Unclassified", "B/LB"]


def test_taxonomic_ortholog_evidence_uses_absolute_alt_support(tmp_path: Path) -> None:
    path = tmp_path / "ortholog_evidence_summary.tsv.gz"
    fields = [
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
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            [
                {
                    "strategy": "s1",
                    "target_context": "cds",
                    "taxonomic_scope": "mammalia",
                    "evidence_unit": "species",
                    "site_aligned_count": depth,
                    "alt_support_count": alt,
                    "gnomad_found_count": found,
                    "gnomad_not_found_count": 1,
                    "gnomad_lookup_failed_count": 3,
                }
                for depth, alt, found in [(10, 1, 1), (10, 5, 2), (20, 10, 3)]
            ]
        )

    available, cells, distributions = read_taxonomic_ortholog_evidence(path)

    assert available
    assert set(cells["taxonomic_scope"]) == {"mammalia"}
    assert set(cells["evidence_unit"]) == {"species"}
    assert cells["alt_label"].str.contains("%").sum() == 0
    assert int(cells[cells["quantile_count"].eq(2)]["gnomad_eligible_count"].sum()) == 9
    site = distributions[distributions["metric"].eq("site_aligned")].set_index("value")
    alt = distributions[distributions["metric"].eq("exact_alt")].set_index("value")
    assert int(site.loc[10, "variant_count"]) == 11
    assert int(site.loc[20, "variant_count"]) == 7
    assert int(alt.loc[1, "variant_count"]) == 5
    assert int(distributions[distributions["metric"].eq("exact_alt")]["variant_count"].sum()) == 18


def test_variant_summary_prefers_compact_ortholog_evidence_and_tracks_its_cache(
    tmp_path: Path,
) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    support = tmp_path / "variant_strategy_support.tsv.gz"
    compact = tmp_path / "ortholog_evidence_summary.tsv.gz"
    pd.DataFrame(
        [
            {
                "variant_key": "1:100:A>G",
                "gene_id": "1",
                "event_type": "snv",
                "ref": "A",
                "alt": "G",
                "lookup_status": "ok",
                "strategies": "s1",
            }
        ],
        columns=VARIANT_USECOLS,
    ).fillna("").to_csv(annotations, sep="\t", index=False, compression="gzip")
    pd.DataFrame(
        [
            {
                "variant_key": "1:100:A>G",
                "gene_id": "1",
                "strategy": "s1",
                "alt_support_row_count": 11,
                "alt_support_ortholog_count": 11,
                "site_aligned_ortholog_count": 10,
            }
        ]
    ).to_csv(support, sep="\t", index=False, compression="gzip")

    def write_compact(found: int) -> None:
        pd.DataFrame(
            [
                {
                    "strategy": "s1",
                    "target_context": "cds",
                    "taxonomic_scope": "mammalia",
                    "evidence_unit": "species",
                    "site_aligned_count": 10,
                    "alt_support_count": 5,
                    "gnomad_found_count": found,
                    "gnomad_not_found_count": 1,
                    "gnomad_lookup_failed_count": 0,
                }
            ]
        ).to_csv(compact, sep="\t", index=False, compression="gzip")

    write_compact(2)
    work_dir = tmp_path / "analytics"
    profile = PerformanceProfile(
        tmp_path / "performance.json",
        run_dir=tmp_path,
        report_path=tmp_path / "report.html",
    )
    with profile.stage("Variant summary"):
        summary = build_variant_summary(
            annotations,
            work_dir,
            strategy_label=str,
            variant_strategy_support_path=support,
            ortholog_evidence_summary_path=compact,
            performance_profile=profile,
        )
    profile.finish()

    assert summary.ortholog_evidence_available
    aggregation = next(
        stage
        for stage in profile.stages
        if stage["name"] == "Variant summary aggregation"
    )
    assert aggregation["parent_id"] == profile.stages[0]["id"]
    assert aggregation["metrics"]["duckdb_source_scan_seconds"] >= 0
    assert aggregation["metrics"]["duckdb_compacted_relations_seconds"] >= 0
    assert aggregation["metrics"]["duckdb_allele_gene_row_count"] == 1
    assert aggregation["metrics"]["duckdb_global_allele_row_count"] == 1
    assert aggregation["metrics"]["duckdb_memory_limit"]
    assert aggregation["metrics"]["summary_assembly_seconds"] >= 0
    assert "temporary_sqlite_bytes" not in aggregation["metrics"]
    assert set(summary.ortholog_evidence_cells["taxonomic_scope"]) == {"mammalia"}
    assert int(
        summary.ortholog_evidence_cells[
            summary.ortholog_evidence_cells["quantile_count"].eq(2)
        ]["gnomad_found_count"].sum()
    ) == 2

    cached = build_variant_summary(
        annotations,
        work_dir,
        strategy_label=str,
        variant_strategy_support_path=support,
        ortholog_evidence_summary_path=compact,
    )
    assert cached.cache_hit

    write_compact(3)
    rebuilt = build_variant_summary(
        annotations,
        work_dir,
        strategy_label=str,
        variant_strategy_support_path=support,
        ortholog_evidence_summary_path=compact,
    )
    assert not rebuilt.cache_hit
    assert int(
        rebuilt.ortholog_evidence_cells[
            rebuilt.ortholog_evidence_cells["quantile_count"].eq(2)
        ]["gnomad_found_count"].sum()
    ) == 3


def test_variant_summary_accepts_compact_annotation_schema(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    rows = [
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "event_type": "snv",
            "ref": "A",
            "alt": "G",
            "strategies": "s1,s2",
            "support_row_count": "4",
            "support_ortholog_count": "3",
            "clinvar_id": "VCV1",
            "clinvar_sig": "Benign",
            "clinvar_review_stars": "2",
            "gnomad_af": "0.01",
            "vep_status": "ok",
            "vep_primary_consequence": "synonymous_variant",
        },
        {
            "variant_key": "1:200:C>A",
            "gene_id": "1",
            "event_type": "snv",
            "ref": "C",
            "alt": "A",
            "strategies": "s1",
            "support_row_count": "1",
            "support_ortholog_count": "1",
            "clinvar_id": "VCV2",
            "clinvar_sig": "Pathogenic",
            "clinvar_review_stars": "3",
            "gnomad_af": "",
            "vep_status": "ok",
            "vep_primary_consequence": "missense_variant",
        },
        {
            "variant_key": "1:300:AT>A",
            "gene_id": "1",
            "event_type": "del",
            "ref": "AT",
            "alt": "A",
            "strategies": "s2",
            "support_row_count": "2",
            "support_ortholog_count": "2",
            "clinvar_id": "",
            "clinvar_sig": "",
            "clinvar_review_stars": "",
            "gnomad_af": "",
            "vep_status": "no_target_gene",
            "vep_primary_consequence": "",
        },
    ]
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=VARIANT_USECOLS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in VARIANT_USECOLS})

    summary = build_variant_summary(
        annotations,
        tmp_path / "analytics",
        strategy_label=lambda value: value,
    )

    assert summary.input_row_count == 3
    assert summary.unique_variant_count == 3
    assert summary.all_strategy_variant_count == 1
    assert summary.strategy_record_count == 4
    assert summary.strategies == ["s1", "s2"]
    assert summary.gene_variant_counts.to_dict(orient="records") == [
        {"strategy": "s1", "gene_id": "1", "Variant_Count": 2},
        {"strategy": "s2", "gene_id": "1", "Variant_Count": 2},
    ]
    by_strategy = summary.strategy_stats.set_index("Strategy")
    assert by_strategy.loc["s1", "Ti/Tv"] == 1.0
    assert by_strategy.loc["s2", "Ti/Tv"] == float("inf")
    assert summary.clinvar_found == 2
    assert summary.gnomad_found == 1
    assert summary.consequence_source == "Ensembl VEP"

    cached = build_variant_summary(
        annotations,
        tmp_path / "analytics",
        strategy_label=lambda value: value,
    )
    assert cached.cache_hit
    assert cached.all_strategy_variant_count == 1
    assert cached.gene_variant_counts.equals(summary.gene_variant_counts)
    assert cached.strategy_stats.set_index("Strategy").loc["s2", "Ti/Tv"] == float("inf")
    assert (tmp_path / "analytics" / "variant_summary.json.gz").stat().st_mode & 0o777 == 0o644


def test_variant_summary_uses_vep_consequences_for_all_candidates(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    fields = VARIANT_USECOLS
    rows = [
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "event_type": "snv",
            "ref": "A",
            "alt": "G",
            "strategies": "s1",
            "gnomad_af": "",
            "vep_status": "ok",
            "vep_primary_consequence": "missense_variant",
        },
        {
            "variant_key": "1:200:C>T",
            "gene_id": "1",
            "event_type": "snv",
            "ref": "C",
            "alt": "T",
            "strategies": "s1",
            "gnomad_af": "0.01",
            "vep_status": "no_target_gene",
            "vep_primary_consequence": "",
        },
    ]
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary = build_variant_summary(annotations, tmp_path / "analytics", strategy_label=str)

    assert summary.consequence_source == "Ensembl VEP"
    assert summary.consequence_counts.to_dict(orient="records") == [
        {"strategy": "s1", "value": UNANNOTATED_CONSEQUENCE, "Variant_Count": 1},
        {"strategy": "s1", "value": "missense_variant", "Variant_Count": 1},
    ]
    assert summary.gnomad_consequence_counts.to_dict(orient="records") == [
        {
            "strategy": "s1",
            "gnomad_status": "found",
            "value": UNANNOTATED_CONSEQUENCE,
            "Variant_Count": 1,
        },
        {
            "strategy": "s1",
            "gnomad_status": "not_found",
            "value": "missense_variant",
            "Variant_Count": 1,
        },
    ]


def test_variant_summary_rebuilds_a_corrupt_cache(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=VARIANT_USECOLS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                field: value
                for field, value in {
                    "variant_key": "1:100:A>G",
                    "gene_id": "1",
                    "event_type": "snv",
                    "ref": "A",
                    "alt": "G",
                    "strategies": "s1",
                }.items()
            }
        )

    work_dir = tmp_path / "analytics"
    first = build_variant_summary(annotations, work_dir, strategy_label=str)
    assert not first.cache_hit
    (work_dir / "variant_summary.json.gz").write_bytes(b"not a gzip stream")

    rebuilt = build_variant_summary(annotations, work_dir, strategy_label=str)
    assert not rebuilt.cache_hit
    assert rebuilt.unique_variant_count == 1


def test_variant_summary_adds_per_strategy_pathogenic_support(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    support = tmp_path / "variant_strategy_support.tsv.gz"
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=VARIANT_USECOLS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "variant_key": "1:100:A>G",
                "gene_id": "1",
                "event_type": "snv",
                "ref": "A",
                "alt": "G",
                "lookup_status": "ok",
                "strategies": "s1,s2",
                "clinvar_sig": "Pathogenic",
                "clinvar_review_stars": "2",
                "clinvar_scv_count": "7",
            }
        )
    pd_rows = [
        ["1:100:A>G", "1", "s1", "5", "3"],
        ["1:100:A>G", "1", "s2", "9", "7"],
    ]
    with gzip.open(support, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "variant_key",
                "gene_id",
                "strategy",
                "alt_support_row_count",
                "alt_support_ortholog_count",
            ]
        )
        writer.writerows(pd_rows)

    summary = build_variant_summary(
        annotations,
        tmp_path / "analytics",
        strategy_label=str,
        variant_strategy_support_path=support,
    )

    row = summary.pathogenic_rows.iloc[0]
    assert summary.pathogenic_variant_count == 1
    assert row["support_ortholog_mean"] == 5.0
    assert row["support_ortholog_min"] == 3
    assert row["support_ortholog_max"] == 7
    assert not summary.ortholog_evidence_available


def test_variant_summary_builds_ortholog_evidence_without_failed_lookups(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    support = tmp_path / "variant_strategy_support.tsv.gz"
    features = tmp_path / "target_features.tsv.gz"
    genes = tmp_path / "genes.tsv.gz"
    failures = tmp_path / "failures.tsv.gz"
    variants = [
        ("1:1:A>G", "0.01", 1, 4),
        ("1:2:C>T", "", 2, 4),
        ("1:3:G>A", "", 3, 4),
        ("1:7:T>C", "0.02", 4, 8),
    ]
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=VARIANT_USECOLS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for variant_key, af, _alt, _depth in variants:
            writer.writerow(
                {
                    "variant_key": variant_key,
                    "gene_id": "1",
                    "event_type": "snv",
                    "ref": variant_key.split(":")[-1][0],
                    "alt": variant_key[-1],
                    "lookup_status": "ok",
                    "strategies": "s1",
                    "gnomad_af": af,
                }
            )
    with gzip.open(support, "wt", newline="") as handle:
        fields = [
            "variant_key",
            "gene_id",
            "strategy",
            "alt_support_row_count",
            "alt_support_ortholog_count",
            "site_aligned_ortholog_count",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for variant_key, _af, alt, depth in variants:
            writer.writerow(
                {
                    "variant_key": variant_key,
                    "gene_id": "1",
                    "strategy": "s1",
                    "alt_support_row_count": alt,
                    "alt_support_ortholog_count": alt,
                    "site_aligned_ortholog_count": depth,
                }
            )
    with gzip.open(features, "wt", newline="") as handle:
        fields = ["gene_id", "feature_type", "target_start0", "target_end0"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            [
                {"gene_id": "1", "feature_type": "gene", "target_start0": 0, "target_end0": 10},
                {"gene_id": "1", "feature_type": "cds", "target_start0": 0, "target_end0": 5},
                {"gene_id": "1", "feature_type": "intron", "target_start0": 5, "target_end0": 10},
            ]
        )
    with gzip.open(genes, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["gene_id", "begin"], delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow({"gene_id": "1", "begin": 1})
    with gzip.open(failures, "wt", newline="") as handle:
        fields = ["source", "scope", "chrom", "start", "end", "failure_type", "message"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {"source": "gnomad", "scope": "region", "chrom": "1", "start": 3, "end": 3}
        )

    summary = build_variant_summary(
        annotations,
        tmp_path / "analytics",
        strategy_label=str,
        target_features_path=features,
        genes_path=genes,
        annotation_failures_path=failures,
        variant_strategy_support_path=support,
    )

    assert summary.ortholog_evidence_available
    assert set(summary.ortholog_evidence_cells["quantile_count"]) == {2, 4, 10}
    assert int(
        summary.ortholog_evidence_distributions[
            summary.ortholog_evidence_distributions["metric"].eq("site_aligned")
        ]["variant_count"].sum()
    ) == 4
    median_cells = summary.ortholog_evidence_cells[
        summary.ortholog_evidence_cells["quantile_count"].eq(2)
    ]
    assert int(median_cells["gnomad_eligible_count"].sum()) == 3
    cds = median_cells[median_cells["target_context"].eq("cds")]
    assert int(cds["gnomad_eligible_count"].sum()) == 2
    assert int(cds["gnomad_found_count"].sum()) == 1

    cached = build_variant_summary(
        annotations,
        tmp_path / "analytics",
        strategy_label=str,
        target_features_path=features,
        genes_path=genes,
        annotation_failures_path=failures,
        variant_strategy_support_path=support,
    )
    assert cached.cache_hit
    assert cached.ortholog_evidence_available
    assert len(cached.ortholog_evidence_cells) == len(summary.ortholog_evidence_cells)
    assert len(cached.ortholog_evidence_distributions) == len(
        summary.ortholog_evidence_distributions
    )


def test_variant_summary_aggregates_target_context_and_gnomad_strata(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    features = tmp_path / "target_features.tsv.gz"
    genes = tmp_path / "genes.tsv.gz"
    rows = [
        {"variant_key": "1:1:A>G", "gene_id": "1", "event_type": "snv", "ref": "A", "alt": "G", "lookup_status": "ok", "strategies": "s1", "gnomad_af": "0.01"},
        {"variant_key": "1:7:C>T", "gene_id": "1", "event_type": "snv", "ref": "C", "alt": "T", "lookup_status": "ok", "strategies": "s1", "gnomad_af": ""},
    ]
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=VARIANT_USECOLS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in VARIANT_USECOLS})
    with gzip.open(features, "wt", newline="") as handle:
        fields = ["gene_id", "feature_type", "target_start0", "target_end0"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            [
                {"gene_id": "1", "feature_type": "gene", "target_start0": 0, "target_end0": 10},
                {"gene_id": "1", "feature_type": "cds", "target_start0": 0, "target_end0": 5},
                {"gene_id": "1", "feature_type": "intron", "target_start0": 5, "target_end0": 10},
            ]
        )
    with gzip.open(genes, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene_id", "begin"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({"gene_id": "1", "begin": 1})

    summary = build_variant_summary(
        annotations,
        tmp_path / "analytics",
        strategy_label=str,
        target_features_path=features,
        genes_path=genes,
    )

    contexts = summary.target_context_counts.set_index("target_context")["Variant_Count"].to_dict()
    assert contexts == {"cds": 1, "intron": 1}
    strata = summary.gnomad_context_counts.set_index(["gnomad_status", "target_context"])["Variant_Count"].to_dict()
    assert strata == {("found", "cds"): 1, ("not_found", "intron"): 1}
    assert summary.gnomad_af_summary.iloc[0]["Median"] == -2.0


def test_variant_summary_excludes_failed_gnomad_lookups_from_denominator(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    failures = tmp_path / "failures.tsv.gz"
    rows = [
        {"variant_key": "1:100:A>G", "gene_id": "1", "event_type": "snv", "ref": "A", "alt": "G", "lookup_status": "ok", "strategies": "s1", "gnomad_af": "0.01"},
        {"variant_key": "1:200:C>T", "gene_id": "1", "event_type": "snv", "ref": "C", "alt": "T", "lookup_status": "ok", "strategies": "s1", "gnomad_af": ""},
        {"variant_key": "1:300:G>A", "gene_id": "1", "event_type": "snv", "ref": "G", "alt": "A", "lookup_status": "ok", "strategies": "s1", "gnomad_af": ""},
        {"variant_key": "raw", "gene_id": "1", "event_type": "del", "ref": "A", "alt": "", "lookup_status": "missing_left_anchor", "strategies": "s1", "gnomad_af": ""},
    ]
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=VARIANT_USECOLS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in VARIANT_USECOLS})
    with gzip.open(failures, "wt", newline="") as handle:
        fields = ["source", "scope", "chrom", "start", "end", "failure_type", "message"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {"source": "gnomad", "scope": "region", "chrom": "1", "start": 190, "end": 210}
        )

    summary = build_variant_summary(
        annotations,
        tmp_path / "analytics",
        strategy_label=str,
        annotation_failures_path=failures,
    )

    stats = summary.strategy_stats.iloc[0]
    assert stats["gnomAD Found"] == 1
    assert stats["gnomAD Eligible"] == 2
    assert stats["gnomAD lookup failed"] == 2
    assert stats["gnomAD found %"] == 0.5
    assert summary.gnomad_lookup_failed == 2
    assert set(summary.gnomad_event_counts["gnomad_status"]) == {"found", "not_found"}
    assert set(summary.gnomad_consequence_counts["gnomad_status"]) == {"found", "not_found"}
    assert summary.gnomad_consequence_counts["Variant_Count"].sum() == 2
