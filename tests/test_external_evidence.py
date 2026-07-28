from __future__ import annotations

from pathlib import Path

import pandas as pd

from analytics.core import external_evidence
from analytics.core.external_evidence import _annotate_gnomad, categorize_clinvar_sig


def test_clinvar_categories_keep_missing_classification_separate() -> None:
    assert categorize_clinvar_sig("") == ""
    assert categorize_clinvar_sig("Benign") == "B/LB"
    assert categorize_clinvar_sig("Likely_pathogenic") == "P/LP"
    assert categorize_clinvar_sig("Uncertain_significance") == "VUS"
    assert categorize_clinvar_sig("Conflicting_classifications_of_pathogenicity") == "Other"


def test_gnomad_annotation_distinguishes_absence_from_failed_lookup(monkeypatch) -> None:
    variants = pd.DataFrame(
        [
            {"variant_key": "1:100:A>G", "chrom": "1", "pos": 100, "ref": "A", "alt": "G"},
            {"variant_key": "2:200:C>T", "chrom": "2", "pos": 200, "ref": "C", "alt": "T"},
        ]
    )

    def fake_fetch(chrom, _start, _end):
        if chrom == "2":
            raise RuntimeError("network unavailable")
        return [
            {
                "chrom": "1",
                "pos": 100,
                "ref": "A",
                "alt": "G",
                "joint": {"an": 100, "ac": 2},
            }
        ]

    monkeypatch.setattr(external_evidence, "fetch_region_variants_recursive", fake_fetch)

    evidence, summary = _annotate_gnomad(variants)
    evidence = evidence.set_index("variant_key")

    assert evidence.loc["1:100:A>G", "gnomad_status"] == "ok"
    assert bool(evidence.loc["1:100:A>G", "gnomad_found"])
    assert evidence.loc["1:100:A>G", "gnomad_af"] == 0.02
    assert evidence.loc["2:200:C>T", "gnomad_status"] == "error"
    assert not bool(evidence.loc["2:200:C>T", "gnomad_found"])
    assert summary["failed_region_count"] == 1


def test_incomplete_external_evidence_retries_only_failed_alleles(
    tmp_path: Path,
    monkeypatch,
) -> None:
    matched = pd.DataFrame(
        [
            {"variant_key": "1:100:A>G", "chrom": "1", "pos": 100, "ref": "A", "alt": "G"},
            {
                "variant_key": "1:500000:C>T",
                "chrom": "1",
                "pos": 500000,
                "ref": "C",
                "alt": "T",
            },
        ]
    )
    matched_path = tmp_path / "matched.tsv.gz"
    matched.to_csv(matched_path, sep="\t", index=False, compression="gzip")
    clinvar_vcf = tmp_path / "clinvar.vcf.gz"
    clinvar_vcf.write_bytes(b"vcf")
    Path(f"{clinvar_vcf}.tbi").write_bytes(b"index")
    output_path = tmp_path / "evidence.tsv.gz"
    manifest_path = tmp_path / "evidence.manifest.json"

    monkeypatch.setattr(
        external_evidence,
        "_annotate_clinvar",
        lambda variants, _path: pd.DataFrame(
            {
                "variant_key": variants["variant_key"],
                "clinvar_found": False,
                "clinvar_classified": False,
                "clinvar_class": "",
            }
        ),
    )
    calls = []

    def fake_gnomad(variants: pd.DataFrame):
        calls.append(variants["variant_key"].tolist())
        if len(calls) == 1:
            return pd.DataFrame(
                [
                    {
                        "variant_key": "1:100:A>G",
                        "gnomad_status": "ok",
                        "gnomad_found": True,
                        "gnomad_af": 0.1,
                    },
                    {
                        "variant_key": "1:500000:C>T",
                        "gnomad_status": "error",
                        "gnomad_found": False,
                        "gnomad_af": None,
                    },
                ]
            ), {"failed_region_count": 1, "errors": [{"chrom": "1"}]}
        return pd.DataFrame(
            [
                {
                    "variant_key": "1:500000:C>T",
                    "gnomad_status": "ok",
                    "gnomad_found": False,
                    "gnomad_af": None,
                }
            ]
        ), {"failed_region_count": 0, "errors": []}

    monkeypatch.setattr(external_evidence, "_annotate_gnomad", fake_gnomad)

    _first, first_manifest = external_evidence.build_external_evidence(
        matched=matched,
        matched_path=matched_path,
        clinvar_vcf=clinvar_vcf,
        output_path=output_path,
        manifest_path=manifest_path,
    )
    second, second_manifest = external_evidence.build_external_evidence(
        matched=matched,
        matched_path=matched_path,
        clinvar_vcf=clinvar_vcf,
        output_path=output_path,
        manifest_path=manifest_path,
    )

    assert first_manifest["complete"] is False
    assert second_manifest["complete"] is True
    assert calls == [["1:100:A>G", "1:500000:C>T"], ["1:500000:C>T"]]
    assert second_manifest["gnomad"]["cached_ok_allele_count"] == 1
    assert second_manifest["gnomad"]["queried_allele_count"] == 1
    preserved = second.set_index("variant_key").loc["1:100:A>G"]
    assert bool(preserved["gnomad_found"])
    assert preserved["gnomad_af"] == 0.1
