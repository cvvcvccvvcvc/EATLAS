from __future__ import annotations

import csv
import gzip
from pathlib import Path

from analytics.core.variant_summary import VARIANT_USECOLS, build_variant_summary


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
            "gnomad_csq": "synonymous_variant",
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
            "gnomad_csq": "missense_variant",
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
            "gnomad_csq": "",
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
    by_strategy = summary.strategy_stats.set_index("Strategy")
    assert by_strategy.loc["s1", "Ti/Tv"] == 1.0
    assert by_strategy.loc["s2", "Ti/Tv"] == float("inf")
    assert summary.clinvar_found == 2
    assert summary.gnomad_found == 1

    cached = build_variant_summary(
        annotations,
        tmp_path / "analytics",
        strategy_label=lambda value: value,
    )
    assert cached.cache_hit
    assert cached.all_strategy_variant_count == 1
    assert cached.strategy_stats.set_index("Strategy").loc["s2", "Ti/Tv"] == float("inf")
    assert (tmp_path / "analytics" / "variant_summary.json.gz").stat().st_mode & 0o777 == 0o644


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
