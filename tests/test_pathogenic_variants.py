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


def test_pathogenic_analysis_builds_condition_backgrounds_and_unique_snv_support(
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
                "alt_support_family_count": 3,
                "site_aligned_ortholog_count": 10,
            },
            {
                "variant_key": "1:200:A>AT",
                "gene_id": "2",
                "strategy": "s1",
                "alt_support_row_count": 2,
                "alt_support_ortholog_count": 2,
                "alt_support_family_count": 1,
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
                "variant_type": "snv",
                "gene_ids": "1",
                "label_class": "pathogenic",
                "clinvar_disease_ids": "MONDO:1|.",
            },
            {
                "variant_key": "1:200:A>AT",
                "clinvar_disease_names": "Disease_two",
                "variant_type": "indel",
                "gene_ids": "2",
                "label_class": "pathogenic",
                "clinvar_disease_ids": "OMIM:2,MedGen:C2",
            },
        ]
    )
    vcf = tmp_path / "clinvar.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNDN=Disease_one;CLNDISDB=MONDO:1\n"
        "1\t200\t2\tA\tAT\t.\t.\tCLNSIG=Likely_pathogenic;CLNDN=Disease_two;CLNDISDB=OMIM:2,MedGen:C2\n"
    )
    support_rows.loc[len(support_rows)] = {
        **support_rows.iloc[0].to_dict(),
        "gene_id": "3",
        "alt_support_ortholog_count": 3,
    }

    analysis = build_pathogenic_variant_analysis(
        summary=summary,
        clinvar_universe=universe,
        clinvar_vcf=vcf,
        condition_cache_dir=tmp_path / "cache",
        eligible_gene_ids_by_strategy={"s1": {"1", "2"}, "s2": {"1"}},
        analytics_dir=tmp_path,
    )

    assert len(analysis.variants) == 2
    first = analysis.variants.set_index("variant_key").loc["1:100:A>G"]
    assert first["pathogenic_subtype"] == "Pathogenic"
    assert bool(first["low_penetrance"])
    assert first["conditions"] == "Disease one"
    assert first["condition_ids"] == "MONDO:1"
    assert set(
        analysis.condition_counts.loc[analysis.condition_counts["variant_count"].gt(0), "condition"]
    ) == {"Disease one", "Disease two"}
    assert len(analysis.support_rows) == 1
    assert analysis.support_rows.iloc[0]["alt_support_ortholog_count"] == 4
    assert analysis.support_rows.iloc[0]["gene_id"] == "1"
    assert "phylop100way" not in analysis.variants
    exported = pd.read_csv(analysis.variants_path, sep="\t", compression="gzip")
    assert len(exported) == 2
    rendered = "".join(build_pathogenic_variant_sections(analysis))
    assert "Pathogenic ClinVar Hits" in rendered
    assert "Primary sort" in rendered
    assert "Disease one" in rendered
    assert "SNV support rows plotted" not in rendered
    assert "phyloP100way" not in rendered
    from analytics.reporting.document import render_html

    (tmp_path / "pathogenic_smoke.html").write_text(
        render_html(
            [
                (
                    "pathogenic-clinvar-hits",
                    "Pathogenic ClinVar Hits",
                    build_pathogenic_variant_sections(analysis),
                ),
            ]
        )
    )


def test_clinvar_conditions_deduplicate_alleles_identifiers_and_preserve_denominators(
    tmp_path: Path,
) -> None:
    from analytics.analyses.clinvar_conditions import (
        global_condition_distribution,
        parse_conditions,
    )

    vcf = tmp_path / "clinvar.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(
            f"1\t{position}\t1\tA\tG\t.\t.\tCLNSIG={label};CLNDN={names};CLNDISDB={ids}\n"
            for position, label, names, ids in [
                (10, "Pathogenic", "Disease|Disease_alias", "MedGen:C1|MedGen:C1"),
                (10, "Pathogenic", "Disease", "MedGen:C1"),
                (20, "Likely_pathogenic", "not_provided", "."),
                (30, "Benign", "Other", "MedGen:C2"),
                (40, "Pathogenic", "Other", "MedGen:C2"),
                (40, "Benign", "Other", "MedGen:C2"),
            ]
        )
    )
    counts = global_condition_distribution(vcf, tmp_path / "cache")
    snv = counts[counts["variant_type"].eq("snv")]
    assert len(snv) == 1
    assert snv.iloc[0]["variant_count"] == 1
    assert snv.iloc[0]["total_variant_count"] == 2
    assert snv.iloc[0]["named_variant_count"] == 1
    pd.testing.assert_frame_equal(counts, global_condition_distribution(vcf, tmp_path / "cache"))
    assert len(parse_conditions("A|B", "MedGen:C1|MedGen:C1")) == 1


def test_pathogenic_violin_displays_counts_on_log_ticks(tmp_path: Path) -> None:
    import numpy as np
    from analytics.reporting.pathogenic_variants import pathogenic_support_figure
    from analytics.reporting.document import render_html
    from analytics.reporting.components import fig_html

    rows = pd.DataFrame(
        [
            {
                "variant_key": f"1:{100 + index}:A>G",
                "gene_id": "1",
                "strategy": strategy,
                "alt_support_ortholog_count": count,
                "alt_support_family_count": 1,
                "site_aligned_ortholog_count": 400,
            }
            for strategy in ["s1", "s2"]
            for index, count in enumerate([1, 1, 1, 2, 2, 3, 4, 24, 330])
        ]
    )
    figure = pathogenic_support_figure(rows)
    assert all(trace.type == "violin" for trace in figure.data)
    assert figure.data[0].marker.color == figure.data[1].marker.color
    assert list(figure.layout.yaxis.ticktext)[:3] == ["1", "2", "5"]
    assert np.isclose(figure.data[0].y[-1], np.log10(330))
    assert figure.data[0].customdata[-1][2] == 330
    (tmp_path / "violin_smoke.html").write_text(
        render_html([("support", "Exact ALT support", [fig_html(figure)])])
    )
