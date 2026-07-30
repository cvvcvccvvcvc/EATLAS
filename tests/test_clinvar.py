from __future__ import annotations

from genomics.clinvar import record_category, review_stars, significance_class


def test_clinvar_category_distinguishes_absence_from_missing_classification() -> None:
    assert record_category("", found=False) == "Not in ClinVar"
    assert record_category("", found=True) == "Unclassified"
    assert record_category("Benign", found=True) == "B/LB"
    assert record_category("Likely_pathogenic", found=True) == "P/LP"
    assert significance_class("Uncertain_significance") == "VUS"
    assert significance_class("Conflicting_classifications_of_pathogenicity") == "Other"


def test_review_stars_requires_an_exact_unambiguous_mapping() -> None:
    assert review_stars("reviewed_by_expert_panel") == ("3", "mapped")
    assert review_stars("unknown_status") == ("", "unmapped:unknown_status")
    assert review_stars(
        "reviewed_by_expert_panel|criteria_provided,_single_submitter"
    ) == ("", "ambiguous_multiple_review_statuses")
