"""Display order, labels, and colors used across the analytics report."""

from genomics.clinvar import CLINVAR_CLASS_ORDER
from analytics.vep.consequences import (
    DISPLAY_CONSEQUENCE_GROUP_ORDER,
    DISPLAY_CONSEQUENCE_GROUP_TERMS,
)

FEATURE_ORDER = ["gene", "exon", "cds", "utr", "intron"]

PROFILE_FEATURE_ORDER = ["cds", "utr", "intron"]

TARGET_CONTEXT_ORDER = ["cds", "utr", "other_exon", "intron", "other"]

TARGET_CONTEXT_LABELS = {
    "cds": "CDS",
    "utr": "UTR",
    "other_exon": "Other exon",
    "intron": "Intron",
    "other": "Other",
}

TARGET_CONTEXT_COLORS = {
    "CDS": "#2166ac",
    "UTR": "#67a9cf",
    "Other exon": "#d1e5f0",
    "Intron": "#ef8a62",
    "Other": "#bdbdbd",
}

CLINVAR_ORDER = list(CLINVAR_CLASS_ORDER)

CLINVAR_COLORS = {
    "B/LB": "#2ca25f",
    "P/LP": "#de2d26",
    "VUS": "#f1c40f",
    "Other": "#8c8c8c",
    "Unclassified": "#c7c7c7",
    "Not in ClinVar": "#e5e7eb",
}

REVIEW_STAR_ORDER = ["4", "3", "2", "1", "0", "Unmapped"]

REVIEW_STAR_COLORS = {
    "4": "#08519c",
    "3": "#3182bd",
    "2": "#6baed6",
    "1": "#9ecae1",
    "0": "#fdbb84",
    "Unmapped": "#bdbdbd",
}

CONSEQUENCE_GROUP_ORDER = list(DISPLAY_CONSEQUENCE_GROUP_ORDER)

CONSEQUENCE_GROUP_COLORS = {
    "LoF/splice": "#de2d26",
    "Missense/inframe": "#fb6a4a",
    "Synonymous": "#74add1",
    "Noncoding/UTR/intron": "#abd9e9",
    "Other": "#9e9e9e",
    "Not annotated": "#d9d9d9",
}

CONSEQUENCE_GROUP_TERMS = DISPLAY_CONSEQUENCE_GROUP_TERMS

STRATEGY_LABELS = {
    "bwa_pseudoreads_150_75": "BWA pseudo 150/75",
    "minimap2_asm10": "minimap2 asm10",
    "minimap2_asm20": "minimap2 asm20",
    "minimap2_map_ont_pseudoreads_30000_15000": "minimap2 map-ont pseudo 30k/15k",
    "nucmer": "nucmer",
    "precomputed_ensembl_92_mammals_epo_extended": "Ensembl EPO",
}

TAXONOMIC_SCOPE_ORDER = [
    "all",
    "eukaryota",
    "metazoa",
    "vertebrata",
    "tetrapoda",
    "amniota",
    "mammalia",
    "primates",
]

TAXONOMIC_SCOPE_LABELS = {
    "all": "All selected",
    "eukaryota": "Eukaryota",
    "metazoa": "Metazoa",
    "vertebrata": "Vertebrata",
    "tetrapoda": "Tetrapoda",
    "amniota": "Amniota",
    "mammalia": "Mammalia",
    "primates": "Primates",
}

EVIDENCE_UNIT_ORDER = ["ortholog", "species", "genus", "family", "order"]

EVIDENCE_UNIT_LABELS = {
    "ortholog": "Ortholog",
    "species": "Species",
    "genus": "Genus",
    "family": "Family",
    "order": "Order",
}
