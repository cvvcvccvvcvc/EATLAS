from __future__ import annotations

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
