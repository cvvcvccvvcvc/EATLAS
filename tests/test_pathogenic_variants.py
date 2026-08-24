from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from analytics.analyses.pathogenic_variants import build_pathogenic_variant_analysis
from analytics.reporting.pathogenic_variants import build_pathogenic_variant_sections
from genomics.clinvar import pathogenic_subtype


def test_pathogenic_subtype_preserves_low_penetrance_flag() -> None:
    assert pathogenic_subtype("Pathogenic") == ("Pathogenic", False)
    assert pathogenic_subtype("Likely_pathogenic") == ("Likely pathogenic", False)
    assert pathogenic_subtype("Pathogenic/Likely_pathogenic") == (
        "Pathogenic / likely pathogenic",
        False,
    )
    assert pathogenic_subtype("Pathogenic,_low_penetrance") == (
        "Pathogenic",
        True,
    )


def test_pathogenic_analysis_builds_conditions_and_snv_support_fraction(
    tmp_path: Path,
) -> None:
    common = {
        "lookup_status": "ok",
        "clinvar_category": "P/LP",
        "clinvar_revstat": "criteria_provided,_single_submitter",
        "clinvar_scv_count": "3",
        "clinvar_allele_id": "1",
        "clinvar_hgvs": "NC_000001.11:g.100A>G",
        "clinvar_disease": "",
        "clinvar_variant_type": "single nucleotide variant",
        "gnomad_af": "0.0001",
        "vep_status": "ok",
        "vep_primary_consequence": "missense_variant",
        "vep_consequence_terms": "missense_variant",
        "vep_transcript_id": "NM_1",
        "vep_mane_select": "NM_1",
        "support_ortholog_mean": 3.0,
        "support_ortholog_min": 2,
        "support_ortholog_max": 4,
    }
    pathogenic_rows = pd.DataFrame(
        [
            {
                **common,
                "variant_key": "1:100:A>G",
                "variant_id": "1:100:A>G",
                "gene_id": "1",
                "event_type": "snv",
                "ref": "A",
                "alt": "G",
                "strategies": "s1,s2",
                "clinvar_sig": "Pathogenic,_low_penetrance",
                "clinvar_review_stars": "2",
                "clinvar_id": "VCV1",
            },
            {
                **common,
                "variant_key": "1:200:A>AT",
                "variant_id": "1:200:A>AT",
                "gene_id": "2",
                "event_type": "insertion",
                "ref": "A",
                "alt": "AT",
                "strategies": "s1",
                "clinvar_sig": "Likely_pathogenic",
                "clinvar_review_stars": "1",
                "clinvar_id": "VCV2",
            },
        ]
    )
    support_rows = pd.DataFrame(
        [
            {
                "variant_key": "1:100:A>G",
                "gene_id": "1",
                "strategy": "s1",
                "alt_support_row_count": 4,
                "alt_support_ortholog_count": 4,
                "alt_support_genus_count": 3,
                "site_aligned_ortholog_count": 10,
            },
            {
                "variant_key": "1:200:A>AT",
                "gene_id": "2",
                "strategy": "s1",
                "alt_support_row_count": 2,
                "alt_support_ortholog_count": 2,
                "alt_support_genus_count": 1,
                "site_aligned_ortholog_count": 8,
            },
        ]
    )
    summary = SimpleNamespace(
        pathogenic_rows=pathogenic_rows,
        pathogenic_variant_count=2,
        pathogenic_support_rows=support_rows,
        pathogenic_consequence_counts=pd.DataFrame(
            [{"strategy": "s1", "value": "missense_variant", "Variant_Count": 1}]
        ),
    )
    universe = pd.DataFrame(
        [
            {
                "variant_key": "1:100:A>G",
                "clinvar_disease_names": "Disease_one|not_provided",
                "clinvar_disease_ids": "MONDO:1|.",
            },
            {
                "variant_key": "1:200:A>AT",
                "clinvar_disease_names": "Disease_two",
                "clinvar_disease_ids": "OMIM:2,MedGen:C2",
            },
        ]
    )
    conservation = pd.DataFrame(
        [
            {"variant_key": "1:100:A>G", "phyloP100way": 2.5},
            {"variant_key": "1:200:A>AT", "phyloP100way": 1.0},
        ]
    )

    analysis = build_pathogenic_variant_analysis(
        summary=summary,
        clinvar_universe=universe,
        conservation_cohort=conservation,
        analytics_dir=tmp_path,
    )

    assert len(analysis.variants) == 2
    first = analysis.variants.set_index("variant_key").loc["1:100:A>G"]
    assert first["pathogenic_subtype"] == "Pathogenic"
    assert bool(first["low_penetrance"])
    assert first["conditions"] == "Disease one"
    assert first["condition_ids"] == "MONDO:1"
    assert set(analysis.condition_counts["condition"]) == {"Disease one", "Disease two"}
    assert len(analysis.evolution_rows) == 1
    assert analysis.evolution_rows.iloc[0]["alt_support_fraction"] == 0.4
    exported = pd.read_csv(analysis.variants_path, sep="\t", compression="gzip")
    assert len(exported) == 2
    rendered = "".join(build_pathogenic_variant_sections(analysis))
    assert "Pathogenic ClinVar Hits" in rendered
    assert "Primary sort" in rendered
    assert "Disease one" in rendered
