"""ClinVar review-status and clinical-significance semantics."""

from __future__ import annotations

from pathlib import Path


CLINVAR_REVIEW_STARS = {
    "practice_guideline": "4",
    "reviewed_by_expert_panel": "3",
    "criteria_provided,_multiple_submitters,_no_conflicts": "2",
    "criteria_provided,_multiple_submitters": "2",
    "criteria_provided,_conflicting_classifications": "1",
    "criteria_provided,_conflicting_interpretations": "1",
    "criteria_provided,_single_submitter": "1",
    "no_assertion_criteria_provided": "0",
    "no_assertion_provided": "0",
    "no_classification_provided": "0",
    "no_classification_for_the_individual_variant": "0",
}

CLINVAR_CLASSIFIED_ORDER = ("P/LP", "B/LB", "VUS", "Other")
CLINVAR_CLASS_ORDER = (*CLINVAR_CLASSIFIED_ORDER, "Unclassified", "Not in ClinVar")
PATHOGENIC_SUBTYPE_ORDER = (
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic / likely pathogenic",
)


def parse_vcf_record_fields(
    line: str,
    *,
    source: str | Path,
    line_number: int,
) -> list[str] | None:
    """Parse one tabix VCF output line without hiding malformed records."""

    record = line.rstrip("\r\n")
    if not record or record.startswith("#"):
        return None
    fields = record.split("\t")
    if len(fields) < 8:
        preview = record[:200] + ("..." if len(record) > 200 else "")
        raise ValueError(
            f"Malformed ClinVar VCF record from {source} at tabix output line "
            f"{line_number}: expected at least 8 tab-separated fields, "
            f"observed {len(fields)}; record={preview!r}"
        )
    return fields


def normalize_review_status(value: object) -> str:
    return str(value or "").strip().lower().replace(" ", "_")


def review_stars(review_status: object) -> tuple[str, str]:
    """Map a complete ClinVar review status to an exact star count."""

    if not review_status:
        return "", "missing"
    statuses = [
        normalize_review_status(item)
        for item in str(review_status).split("|")
        if item
    ]
    if not statuses:
        return "", "missing"

    stars = []
    for status in statuses:
        star_value = CLINVAR_REVIEW_STARS.get(status)
        if star_value is None:
            return "", f"unmapped:{status}"
        stars.append(star_value)

    unique_stars = sorted(set(stars))
    if len(unique_stars) != 1:
        return "", "ambiguous_multiple_review_statuses"
    return unique_stars[0], "mapped"


def significance_class(value: object) -> str | None:
    """Return the unambiguous report class for a non-empty CLNSIG value."""

    text = str(value or "").lower()
    if not text:
        return None
    if "conflicting" in text:
        return "Other"
    if "uncertain" in text or "vus" in text:
        return "VUS"
    benign = "benign" in text
    pathogenic = "pathogenic" in text
    if benign and not pathogenic:
        return "B/LB"
    if pathogenic and not benign:
        return "P/LP"
    return "Other"


def pathogenic_subtype(value: object) -> tuple[str | None, bool]:
    """Return the P/LP subtype and whether ClinVar marks low penetrance."""

    if significance_class(value) != "P/LP":
        return None, False
    normalized = str(value or "").strip().lower().replace(" ", "_")
    low_penetrance = "low_penetrance" in normalized
    classifications: set[str] = set()
    for item in normalized.split("|"):
        base = item.replace(",_low_penetrance", "").replace("/low_penetrance", "")
        if "pathogenic/likely_pathogenic" in base:
            classifications.update({"pathogenic", "likely_pathogenic"})
        elif "likely_pathogenic" in base:
            classifications.add("likely_pathogenic")
        elif "pathogenic" in base:
            classifications.add("pathogenic")
    if classifications == {"pathogenic"}:
        return "Pathogenic", low_penetrance
    if classifications == {"likely_pathogenic"}:
        return "Likely pathogenic", low_penetrance
    return "Pathogenic / likely pathogenic", low_penetrance


def record_category(significance: object, *, found: bool) -> str:
    if not found:
        return "Not in ClinVar"
    return significance_class(significance) or "Unclassified"
