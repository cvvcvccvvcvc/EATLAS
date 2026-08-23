"""Sequence Ontology consequence priority and analysis group definitions."""

from __future__ import annotations

from genomics.vep.terms import VEP_CONSEQUENCE_ORDER

VALIDATION_CONSEQUENCE_OPTIONS = (
    ("all", "All consequences"),
    ("missense", "Missense"),
    ("synonymous", "Synonymous"),
    ("protein_truncating", "Protein-truncating / frameshift"),
    ("canonical_splice", "Canonical splice"),
    ("inframe_protein_altering", "In-frame / protein-altering"),
    ("splice_region", "Splice region"),
    ("intronic", "Intronic"),
    ("utr_noncoding", "UTR / noncoding"),
    ("other", "Other"),
)

VALIDATION_CONSEQUENCE_TERMS = {
    "missense": frozenset({"missense_variant"}),
    "synonymous": frozenset({"synonymous_variant", "stop_retained_variant"}),
    "protein_truncating": frozenset(
        {
            "frameshift_variant",
            "nonsense",
            "start_lost",
            "stop_gained",
            "transcript_ablation",
        }
    ),
    "canonical_splice": frozenset({"splice_acceptor_variant", "splice_donor_variant"}),
    "inframe_protein_altering": frozenset(
        {
            "coding_sequence_variant",
            "inframe_deletion",
            "inframe_insertion",
            "protein_altering_variant",
            "stop_lost",
        }
    ),
    "splice_region": frozenset(
        {
            "splice_donor_5th_base_variant",
            "splice_donor_region_variant",
            "splice_polypyrimidine_tract_variant",
            "splice_region_variant",
        }
    ),
    "intronic": frozenset({"intron_variant"}),
    "utr_noncoding": frozenset(
        {
            "3_prime_UTR_variant",
            "5_prime_UTR_variant",
            "downstream_gene_variant",
            "mature_miRNA_variant",
            "non_coding_transcript_exon_variant",
            "non_coding_transcript_variant",
            "non-coding_transcript_variant",
            "regulatory_region_ablation",
            "regulatory_region_amplification",
            "regulatory_region_variant",
            "TF_binding_site_variant",
            "upstream_gene_variant",
        }
    ),
}

VALIDATION_CONSEQUENCE_BITS = {
    key: 1 << index
    for index, (key, _label) in enumerate(VALIDATION_CONSEQUENCE_OPTIONS)
    if key != "all"
}
RECOGNIZED_VALIDATION_TERMS = frozenset().union(*VALIDATION_CONSEQUENCE_TERMS.values())

DISPLAY_CONSEQUENCE_GROUP_ORDER = (
    "LoF/splice",
    "Missense/inframe",
    "Synonymous",
    "Noncoding/UTR/intron",
    "Other",
    "Not annotated",
)

UNANNOTATED_CONSEQUENCE = "__not_annotated__"

DISPLAY_CONSEQUENCE_GROUP_TERMS = {
    "LoF/splice": (
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "start_lost",
        "stop_gained",
        "stop_lost",
    ),
    "Missense/inframe": (
        "inframe_deletion",
        "inframe_insertion",
        "missense_variant",
        "protein_altering_variant",
    ),
    "Synonymous": (
        "stop_retained_variant",
        "synonymous_variant",
    ),
    "Noncoding/UTR/intron": (
        "3_prime_UTR_variant",
        "5_prime_UTR_variant",
        "intron_variant",
        "non_coding_transcript_exon_variant",
        "splice_region_variant",
    ),
}


def validation_consequence_memberships(terms_text: object) -> set[str]:
    terms = {term.strip() for term in str(terms_text or "").split("|") if term.strip()}
    memberships = {
        group
        for group, accepted_terms in VALIDATION_CONSEQUENCE_TERMS.items()
        if terms & accepted_terms
    }
    if not terms or terms - RECOGNIZED_VALIDATION_TERMS:
        memberships.add("other")
    return memberships


def validation_consequence_memberships_text(terms_text: object) -> str:
    memberships = validation_consequence_memberships(terms_text)
    return "|".join(
        key
        for key, _label in VALIDATION_CONSEQUENCE_OPTIONS
        if key != "all" and key in memberships
    )


def validation_consequence_membership_mask(terms_text: object) -> int:
    return sum(
        VALIDATION_CONSEQUENCE_BITS[group]
        for group in validation_consequence_memberships(terms_text)
    )


def display_consequence_group(value: object) -> str:
    consequence = str(value or "")
    if consequence == UNANNOTATED_CONSEQUENCE:
        return "Not annotated"
    for group, terms in DISPLAY_CONSEQUENCE_GROUP_TERMS.items():
        if consequence in terms:
            return group
    return "Other"
