from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from analytics.analyses import matched_control as controls
from analytics.analyses.matched_control import (
    _build_matched_rows,
    _generate_candidate_controls,
    _matched_ecdf,
    _matched_metric_summary,
    _matched_summary,
)
from analytics.analyses.target_context import read_disjoint_contexts
from analytics.io.performance import PerformanceProfile


def test_build_target_space_null_end_to_end_with_mocked_annotations(
    tmp_path: Path,
    monkeypatch,
) -> None:
    run_dir = tmp_path / "run"
    target_dir = run_dir / "fetch" / "sequences" / "targets"
    target_dir.mkdir(parents=True)
    genes_path = run_dir / "fetch" / "genes.tsv.gz"
    features_path = run_dir / "fetch" / "target_features.tsv.gz"
    annotations_path = run_dir / "annotation" / "variant_annotations.tsv.gz"
    annotations_path.parent.mkdir(parents=True)

    pd.DataFrame(
        [
            {
                "gene_id": "1",
                "genomic_accession": "NC_000001.11",
                "chromosome": "1",
                "begin": 100,
                "end": 119,
                "sequence_length": 20,
            }
        ]
    ).to_csv(genes_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(
        [
            {"gene_id": "1", "feature_type": "gene", "target_start0": 0, "target_end0": 20},
            {"gene_id": "1", "feature_type": "cds", "target_start0": 0, "target_end0": 20},
        ]
    ).to_csv(features_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(
        [
            {
                "variant_key": "1:100:A>G",
                "gene_id": "1",
                "event_type": "snv",
                "ref": "A",
                "alt": "G",
                "strategies": "s1",
                "lookup_status": "ok",
            }
        ]
    ).to_csv(annotations_path, sep="\t", index=False, compression="gzip")
    import gzip

    with gzip.open(target_dir / "1.fa.gz", "wt") as handle:
        handle.write(">1\n" + "A" * 20 + "\n")

    vep_calls = []

    def fake_vep(frame, _cache_path, **kwargs):
        vep_calls.append(kwargs)
        result = frame[["variant_key", "gene_id"]].drop_duplicates().copy()
        result["status"] = "ok"
        result["primary_consequence"] = "missense_variant"
        result["consequence_terms"] = "missense_variant"
        result["transcript_id"] = "NM_1"
        result["mane_select"] = "NM_1"
        result["canonical"] = True
        result["impact"] = "MODERATE"
        result["variant_class"] = "SNV"
        return result, {"status": "complete", "release": "116", "requested": len(result)}

    def fake_conservation(frame, output_path):
        unique = frame[["variant_key"]].drop_duplicates().copy()
        unique["phyloP100way"] = np.arange(len(unique), dtype=float)
        unique.to_csv(output_path, sep="\t", index=False, compression="gzip")
        return unique, {"status": "complete", "track": "phyloP100way"}

    def fake_external_evidence(*, matched, output_path, manifest_path, **_kwargs):
        assert _kwargs["gnomad_cache_dir"] == tmp_path / "gnomad_cache"
        evidence = matched[["variant_key"]].drop_duplicates().copy()
        evidence["clinvar_found"] = False
        evidence["clinvar_classified"] = False
        evidence["clinvar_class"] = ""
        evidence["gnomad_status"] = "ok"
        evidence["gnomad_found"] = False
        evidence["gnomad_af"] = np.nan
        evidence.to_csv(output_path, sep="\t", index=False, compression="gzip")
        manifest = {"complete": True, "gnomad": {"failed_region_count": 0}}
        manifest_path.write_text("{}\n")
        return evidence, manifest

    monkeypatch.setattr(controls, "annotate_vep_consequences", fake_vep)
    monkeypatch.setattr(controls, "_annotate_conservation", fake_conservation)
    monkeypatch.setattr(controls, "build_external_evidence", fake_external_evidence)
    performance_path = run_dir / "analytics" / "performance.json"
    performance = PerformanceProfile(
        performance_path,
        run_dir=run_dir,
        report_path=run_dir / "report.html",
    )

    analysis = controls.build_target_space_null(
        run_dir=run_dir,
        variant_annotations_tsv=annotations_path,
        target_features_tsv=features_path,
        genes_tsv=genes_path,
        target_sequences_dir=target_dir,
        clinvar_vcf=tmp_path / "clinvar.vcf.gz",
        strategies=["s1"],
        sample_size_per_strategy=10,
        resamples=100,
        seed=3,
        gnomad_cache_dir=tmp_path / "gnomad_cache",
        vep_result_cache_dir=tmp_path / "vep_result_cache",
        performance_profile=performance,
    )

    matched = pd.read_csv(analysis.matched_path, sep="\t", compression="gzip")
    assert matched["role"].tolist().count("observed") == 1
    assert matched["role"].tolist().count("control") == 5
    assert set(matched["primary_consequence"]) == {"missense_variant"}
    assert set(matched["alt"]) == {"G"}
    assert analysis.manifest["matched_focal_count"] == 1
    assert vep_calls
    assert all(
        call["vep_result_cache_dir"] == tmp_path / "vep_result_cache"
        for call in vep_calls
    )
    assert analysis.focal_path is not None and analysis.focal_path.exists()
    assert (
        analysis.focal_manifest_path is not None
        and analysis.focal_manifest_path.exists()
    )
    stage_names = {
        stage["name"]
        for stage in json.loads(performance_path.read_text())["stages"]
    }
    assert {
        "Target-null focal sampling",
        "Target-null observed store",
        "Target-null focal VEP",
        "Target-null control VEP",
        "Target-null observed-control exclusion",
        "Target-null matching",
        "Target-null phyloP",
        "Target-null external evidence",
        "Target-null resampling",
    }.issubset(stage_names)

    analysis.manifest_path.unlink()
    analysis.matched_path.unlink()
    analysis.conservation_path.unlink()
    monkeypatch.setattr(
        controls,
        "_sample_focal_snvs",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("focal sample was not cached")
        ),
    )
    controls.build_target_space_null(
        run_dir=run_dir,
        variant_annotations_tsv=annotations_path,
        target_features_tsv=features_path,
        genes_tsv=genes_path,
        target_sequences_dir=target_dir,
        clinvar_vcf=tmp_path / "clinvar.vcf.gz",
        strategies=["s1"],
        sample_size_per_strategy=10,
        resamples=100,
        seed=3,
        gnomad_cache_dir=tmp_path / "gnomad_cache",
        vep_result_cache_dir=tmp_path / "vep_result_cache",
        performance_profile=performance,
    )
    focal_stages = [
        stage
        for stage in json.loads(performance_path.read_text())["stages"]
        if stage["name"] == "Target-null focal sampling"
    ]
    assert focal_stages[-1]["details"] == "cache hit"


def test_target_space_genes_use_accession_chromosome_for_par_loci(
    tmp_path: Path,
) -> None:
    genes_path = tmp_path / "genes.tsv.gz"
    pd.DataFrame(
        [
            {
                "gene_id": "64109",
                "genomic_accession": "NC_000023.11",
                "chromosome": "X,Y",
                "begin": 1_190_490,
                "end": 1_212_649,
                "sequence_length": 22_160,
            }
        ]
    ).to_csv(genes_path, sep="\t", index=False, compression="gzip")

    assert controls._read_genes(genes_path)["64109"]["chrom"] == "X"


def test_contexts_keep_noncoding_exon_separate_from_other_sequence(tmp_path: Path) -> None:
    features = pd.DataFrame(
        [
            {"gene_id": "1", "feature_type": "gene", "target_start0": 0, "target_end0": 20},
            {"gene_id": "1", "feature_type": "exon", "target_start0": 2, "target_end0": 12},
            {"gene_id": "1", "feature_type": "cds", "target_start0": 4, "target_end0": 8},
            {"gene_id": "1", "feature_type": "intron", "target_start0": 12, "target_end0": 18},
        ]
    )
    path = tmp_path / "features.tsv.gz"
    features.to_csv(path, sep="\t", index=False, compression="gzip")

    contexts = read_disjoint_contexts(path, {"1": 20})["1"]

    assert contexts == [
        (0, 2, "other"),
        (2, 4, "other_exon"),
        (4, 8, "cds"),
        (8, 12, "other_exon"),
        (12, 18, "intron"),
        (18, 20, "other"),
    ]


def test_candidate_controls_preserve_gene_context_and_exact_substitution() -> None:
    focal = pd.DataFrame(
        [
            {
                "gene_id": "1",
                "variant_key": "1:101:A>G",
                "context": "cds",
                "ref": "A",
                "alt": "G",
                "primary_consequence": "missense_variant",
            }
        ]
    )

    candidates = _generate_candidate_controls(
        focal,
        {"1": [(0, 10, "cds")]},
        {"1": {"chrom": "1", "begin": 100, "length": 10}},
        {"1": "AAAAAAAAAA"},
        seed=7,
    )

    assert not candidates.empty
    assert set(candidates["gene_id"]) == {"1"}
    assert set(candidates["context"]) == {"cds"}
    assert set(candidates["ref"]) == {"A"}
    assert set(candidates["alt"]) == {"G"}
    assert "1:101:A>G" not in set(candidates["variant_key"])


def test_matched_rows_require_consequence_match_and_exclude_observed_control() -> None:
    focal = pd.DataFrame(
        [
            {
                "focal_id": "f1",
                "strategy": "s1",
                "gene_id": "1",
                "variant_key": "1:101:A>G",
                "context": "cds",
                "primary_consequence": "missense_variant",
                "chrom": "1",
                "pos": 101,
                "target_pos": 1,
                "ref": "A",
                "alt": "G",
                "vep_consequence_terms": "missense_variant",
                "vep_transcript_id": "NM_1",
            }
        ]
    )
    candidates = pd.DataFrame(
        [
            {
                "control_group": "1|1:101:A>G",
                "variant_key": "1:102:A>G",
                "gene_id": "1",
                "context": "cds",
                "primary_consequence": "missense_variant",
                "chrom": "1",
                "pos": 102,
                "target_pos": 2,
                "ref": "A",
                "alt": "G",
                "vep_consequence_terms": "missense_variant",
                "vep_transcript_id": "NM_1",
            },
            {
                "control_group": "1|1:101:A>G",
                "variant_key": "1:103:A>G",
                "gene_id": "1",
                "context": "cds",
                "primary_consequence": "missense_variant",
                "chrom": "1",
                "pos": 103,
                "target_pos": 3,
                "ref": "A",
                "alt": "G",
                "vep_consequence_terms": "missense_variant",
                "vep_transcript_id": "NM_1",
            },
        ]
    )

    result = _build_matched_rows(focal, candidates, {("1:102:A>G", "s1")})

    assert result["variant_key"].tolist() == ["1:101:A>G", "1:103:A>G"]
    assert result["role"].tolist() == ["observed", "control"]


def test_matched_rows_apply_exclusions_per_strategy_and_preserve_control_order() -> None:
    focal = pd.DataFrame(
        [
            {
                "focal_id": f"f{index}",
                "strategy": strategy,
                "gene_id": "1",
                "variant_key": "1:101:A>G",
                "context": "cds",
                "primary_consequence": "missense_variant",
                "chrom": "1",
                "pos": 101,
                "target_pos": 1,
                "ref": "A",
                "alt": "G",
                "vep_consequence_terms": "missense_variant",
                "vep_transcript_id": "NM_1",
            }
            for index, strategy in enumerate(["s1", "s2"])
        ]
    )
    candidates = pd.DataFrame(
        [
            {
                "control_group": "1|1:101:A>G",
                "variant_key": f"1:{position}:A>G",
                "chrom": "1",
                "pos": position,
                "target_pos": position - 100,
                "ref": "A",
                "alt": "G",
                "vep_consequence_terms": "missense_variant",
                "vep_transcript_id": "NM_1",
            }
            for position in range(102, 109)
        ]
    )

    result = _build_matched_rows(
        focal,
        candidates,
        {("1:102:A>G", "s1")},
    )
    by_focal = {
        focal_id: group["variant_key"].tolist()
        for focal_id, group in result.groupby("focal_id", sort=False)
    }

    assert by_focal["f0"] == [
        "1:101:A>G",
        "1:103:A>G",
        "1:104:A>G",
        "1:105:A>G",
        "1:106:A>G",
        "1:107:A>G",
    ]
    assert by_focal["f1"] == [
        "1:101:A>G",
        "1:102:A>G",
        "1:103:A>G",
        "1:104:A>G",
        "1:105:A>G",
        "1:106:A>G",
    ]
    assert result.groupby("focal_id")["option"].apply(list).to_dict() == {
        "f0": [0, 1, 2, 3, 4, 5],
        "f1": [0, 1, 2, 3, 4, 5],
    }


def test_matched_summary_uses_paired_control_options_deterministically() -> None:
    rows = []
    for index, (observed, controls) in enumerate([(2.0, [0.0, 1.0]), (4.0, [1.0, 3.0])]):
        focal_id = f"f{index}"
        rows.append(
            {
                "focal_id": focal_id,
                "strategy": "s1",
                "role": "observed",
                "phyloP100way": observed,
            }
        )
        rows.extend(
            {
                "focal_id": focal_id,
                "strategy": "s1",
                "role": "control",
                "phyloP100way": value,
            }
            for value in controls
        )
    frame = pd.DataFrame(rows)

    first = _matched_summary(frame, ["strategy"], resamples=200, seed=7)
    second = _matched_summary(frame, ["strategy"], resamples=200, seed=7)

    pd.testing.assert_frame_equal(first, second)
    assert first.loc[0, "matched_focals"] == 2
    assert first.loc[0, "observed_median"] == 3.0
    assert first.loc[0, "observed_ci_low"] <= first.loc[0, "observed_median"]
    assert first.loc[0, "observed_ci_high"] >= first.loc[0, "observed_median"]
    assert first.loc[0, "valid_resamples"] == 200
    assert "empirical_p" not in first.columns


def test_matched_summary_preserves_paired_difference_in_bootstrap() -> None:
    rows = []
    for index, observed in enumerate([2.0, 4.0, 8.0, 16.0]):
        focal_id = f"f{index}"
        rows.append(
            {
                "focal_id": focal_id,
                "strategy": "s1",
                "role": "observed",
                "phyloP100way": observed,
            }
        )
        rows.extend(
            {
                "focal_id": focal_id,
                "strategy": "s1",
                "role": "control",
                "phyloP100way": observed - 1.0,
            }
            for _ in range(3)
        )

    summary = _matched_summary(pd.DataFrame(rows), ["strategy"], resamples=200, seed=7)

    assert summary.loc[0, "median_difference"] == 1.0
    assert summary.loc[0, "difference_ci_low"] == 1.0
    assert summary.loc[0, "difference_ci_high"] == 1.0


def test_matched_metric_summary_resamples_one_control_per_focal() -> None:
    frame = pd.DataFrame(
        [
            {"focal_id": "f1", "strategy": "s1", "role": "observed", "found": 1.0},
            {"focal_id": "f1", "strategy": "s1", "role": "control", "found": 0.0},
            {"focal_id": "f1", "strategy": "s1", "role": "control", "found": 1.0},
            {"focal_id": "f2", "strategy": "s1", "role": "observed", "found": 0.0},
            {"focal_id": "f2", "strategy": "s1", "role": "control", "found": 0.0},
        ]
    )

    summary = _matched_metric_summary(frame, ["strategy"], "found", "mean", 200, 7)

    assert summary.loc[0, "matched_focals"] == 2
    assert summary.loc[0, "observed_value"] == 0.5
    assert summary.loc[0, "observed_ci_low"] <= summary.loc[0, "observed_value"]
    assert summary.loc[0, "observed_ci_high"] >= summary.loc[0, "observed_value"]
    assert np.isfinite(summary.loc[0, "difference_ci_low"])
    assert np.isfinite(summary.loc[0, "difference_ci_high"])
    assert summary.loc[0, "valid_resamples"] == 200


def test_matched_ecdf_weights_each_focal_equally() -> None:
    frame = pd.DataFrame(
        [
            {"focal_id": "f1", "strategy": "s1", "role": "observed", "phyloP100way": 0.0},
            {"focal_id": "f1", "strategy": "s1", "role": "control", "phyloP100way": 0.0},
            {"focal_id": "f2", "strategy": "s1", "role": "observed", "phyloP100way": 2.0},
            {"focal_id": "f2", "strategy": "s1", "role": "control", "phyloP100way": 2.0},
            {"focal_id": "f2", "strategy": "s1", "role": "control", "phyloP100way": 2.0},
            {"focal_id": "f2", "strategy": "s1", "role": "control", "phyloP100way": 2.0},
        ]
    )

    ecdf = _matched_ecdf(frame)
    control = ecdf[ecdf["set"] == "Matched target-space null"]
    midpoint = control.loc[(control["phyloP100way"] - 1.0).abs().idxmin()]

    assert midpoint["fraction_leq"] == 0.5
