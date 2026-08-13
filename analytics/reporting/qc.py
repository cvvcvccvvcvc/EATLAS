"""Methods, quality-control definitions, and artifact provenance."""

from __future__ import annotations

import html
from pathlib import Path

import pandas as pd

from analytics.analyses.conservation_validation import SPLINE_DF
from analytics.analyses.conservation_analysis import ConservationAnalysis
from analytics.analyses.matched_control import TargetSpaceNullAnalysis
from analytics.analyses.variant_summary import VariantSummary
from analytics.annotation.consequences import (
    VALIDATION_CONSEQUENCE_OPTIONS as CONSEQUENCE_OPTIONS,
    VALIDATION_CONSEQUENCE_TERMS as CONSEQUENCE_TERMS,
)
from analytics.io.run_inputs import RunInputs, file_size_label, read_json
from .components import format_int, metric_cards, table_html
from .config import (
    CONSEQUENCE_GROUP_ORDER,
    CONSEQUENCE_GROUP_TERMS,
    EVIDENCE_UNIT_LABELS,
    TAXONOMIC_SCOPE_LABELS,
)
from .conservation import hidden_clinvar_association_views
from .matched_control import build_target_space_null_qc_sections
from .variant_profile import pathogenic_variant_table


def consequence_grouping_table(source: str) -> pd.DataFrame:
    def definition(group: str) -> str:
        if group == "Not annotated":
            return f"VEP status is not ok or the selected {source} consequence is empty."
        if group == "Other":
            return f"Any non-empty {source} consequence not listed above."
        return ", ".join(CONSEQUENCE_GROUP_TERMS.get(group, []))

    rows = [
        {
            "Group": group,
            f"{source} consequence values": definition(group),
        }
        for group in CONSEQUENCE_GROUP_ORDER
    ]
    return pd.DataFrame(rows)


def clinvar_class_mapping_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Class": "P/LP",
                "Rule": "CLNSIG contains pathogenic and does not contain benign, uncertain/VUS, or conflicting.",
                "Used for": "ClinVar class composition and pathogenic-only evidence plots.",
            },
            {
                "Class": "B/LB",
                "Rule": "CLNSIG contains benign and does not contain pathogenic, uncertain/VUS, or conflicting.",
                "Used for": "ClinVar class composition and validation enrichment counts.",
            },
            {
                "Class": "VUS",
                "Rule": "CLNSIG contains uncertain or VUS, unless the record is marked conflicting.",
                "Used for": "Shown as uncertainty in class composition; excluded from validation denominator.",
            },
            {
                "Class": "Other",
                "Rule": "Conflicting, mixed, ambiguous, or non-empty CLNSIG values outside the clean classes above.",
                "Used for": "Shown separately in class composition; excluded from validation denominator.",
            },
            {
                "Class": "Unclassified",
                "Rule": "A ClinVar record was found, but the allele has no CLNSIG value.",
                "Used for": "Not included in classified-variant plots or validation denominator.",
            },
            {
                "Class": "Not in ClinVar",
                "Rule": "No exact normalized ClinVar allele record was found.",
                "Used for": "Not included in classified-variant plots or validation denominator.",
            },
        ]
    )


def clinvar_review_star_mapping_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Stars": "4",
                "ClinVar review status values": "practice_guideline",
                "Interpretation": "Practice guideline.",
            },
            {
                "Stars": "3",
                "ClinVar review status values": "reviewed_by_expert_panel",
                "Interpretation": "Reviewed by expert panel.",
            },
            {
                "Stars": "2",
                "ClinVar review status values": (
                    "criteria_provided,_multiple_submitters,_no_conflicts; "
                    "criteria_provided,_multiple_submitters"
                ),
                "Interpretation": "Multiple submitters with criteria; no-conflict status when provided by ClinVar.",
            },
            {
                "Stars": "1",
                "ClinVar review status values": (
                    "criteria_provided,_single_submitter; "
                    "criteria_provided,_conflicting_classifications; "
                    "criteria_provided,_conflicting_interpretations"
                ),
                "Interpretation": "Criteria provided, but lower review confidence or conflicting submissions.",
            },
            {
                "Stars": "0",
                "ClinVar review status values": (
                    "no_assertion_criteria_provided; no_assertion_provided; "
                    "no_classification_provided; no_classification_for_the_individual_variant"
                ),
                "Interpretation": "No assertion criteria or no classification.",
            },
            {
                "Stars": "Unmapped",
                "ClinVar review status values": "Missing, empty, or unrecognized review status.",
                "Interpretation": "Kept visible only when such records are present.",
            },
        ]
    )


def validation_method_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Step": "Universe",
                "Definition": (
                    "ClinVar alleles overlapping fetched target loci, normalized to the same target-context "
                    "variant_key representation as GAPH annotations."
                ),
            },
            {
                "Step": "Variant types",
                "Definition": "SNV and INDEL validation are computed separately; complex/MNV/symbolic alleles are excluded.",
            },
            {
                "Step": "Included labels",
                "Definition": "Only clean B/LB and clean P/LP ClinVar labels enter the validation denominator.",
            },
            {
                "Step": "Excluded labels",
                "Definition": "VUS, missing CLNSIG, conflicting/other, mixed B/P labels, and unnormalizable alleles.",
            },
            {
                "Step": "Observed",
                "Definition": "A ClinVar allele is observed for a strategy when its normalized variant_key is present in that strategy's variant_annotations rows.",
            },
            {
                "Step": "2x2 table",
                "Definition": "[B/LB observed, P/LP observed; B/LB not observed, P/LP not observed].",
            },
            {
                "Step": "Statistics",
                "Definition": (
                    "Raw odds ratio, approximate 95% CI on log(OR) with Haldane 0.5 correction for zero cells, "
                    "and two-sided Fisher exact p-value. Benjamini-Hochberg FDR is computed across strategies "
                    "within each variant-type, target-context, and consequence selection."
                ),
            },
        ]
    )


def validation_consequence_grouping_table(source: str = "Ensembl VEP") -> pd.DataFrame:
    rows = []
    for key, label in CONSEQUENCE_OPTIONS:
        if key == "all":
            continue
        terms = CONSEQUENCE_TERMS.get(key)
        rows.append(
            {
                "Consequence subset": label,
                f"{source} terms": ", ".join(sorted(terms))
                if terms
                else f"Missing {source} consequence or any term not assigned to a named subset.",
            }
        )
    return pd.DataFrame(rows)


def vep_qc_tables(
    candidate_manifest: dict | None,
    validation,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidate = candidate_manifest or {}
    clinvar = dict(getattr(validation, "manifest", {}).get("vep", {}))
    datasets = [
        (
            "Candidates",
            int(candidate.get("row_count", 0)),
            dict(candidate.get("status_counts", {})),
            dict(candidate.get("config", {})),
        ),
        (
            "ClinVar universe",
            int(clinvar.get("allele_count", 0)),
            dict(clinvar.get("status_counts", {})),
            dict(clinvar.get("contract", {})),
        ),
    ]
    configuration_rows = []
    status_rows = []
    for dataset, total, status_counts, config in datasets:
        configuration_rows.append(
            {
                "Dataset": dataset,
                "VEP release": str(config.get("release", "")),
                "Backend": str(config.get("backend", "")),
                "Rows": total,
            }
        )
        statuses = ["ok"] if "ok" in status_counts else []
        statuses.extend(sorted(status for status in status_counts if status != "ok"))
        for status in statuses:
            count = int(status_counts[status])
            status_rows.append(
                {
                    "Dataset": dataset,
                    "VEP status": status,
                    "Rows": count,
                    "Fraction": f"{count / total:.3%}" if total else "",
                }
            )
    return pd.DataFrame(configuration_rows), pd.DataFrame(status_rows)


def conservation_validation_method_table(analysis: ConservationAnalysis | None) -> pd.DataFrame:
    versions = analysis.validation.r_versions if analysis is not None else {}
    return pd.DataFrame(
        [
            {
                "Step": "Variant set",
                "Definition": (
                    "Clean B/LB and P/LP ClinVar SNVs and simple INDELs from the normalized target-locus universe. "
                    "For each strategy, the denominator is restricted to genes with an alignment result for that "
                    "strategy. Within those genes, ALT_observed=0 means that the strategy did not report that exact "
                    "normalized ALT; no per-base callability filter is applied."
                ),
            },
            {
                "Step": "Target contexts",
                "Definition": (
                    "Each normalized ClinVar allele is assigned to one target-locus context using the affected "
                    "target position and the precedence CDS > UTR > exon > intron > other. If the same allele "
                    "overlaps multiple target genes, the highest-priority context is used once, without duplicating "
                    "the allele."
                ),
            },
            {
                "Step": "Consequence subsets",
                "Definition": (
                    "Subsets use release-pinned RefSeq VEP Sequence Ontology terms. A record with terms from multiple groups enters "
                    "each matching group once and also enters All consequences; consequence-view counts are not additive."
                ),
            },
            {
                "Step": "Conservation annotation",
                "Definition": (
                    "phyloP100way is read from the hg38 UCSC bigWig. SNVs use the substituted base; deletions use "
                    "the mean across deleted reference bases excluding the VCF padding base; insertions use the mean "
                    "of the two flanking reference bases. All required bases must have a score."
                ),
            },
            {
                "Step": "Fixed bands",
                "Definition": (
                    "<= -1.30103 nominal acceleration band; (-1.30103, 1.30103) central band; "
                    ">= 1.30103 nominal conservation band. The cutoffs equal signed -log10(0.05) "
                    "for a single-base phyloP score and are descriptive, not genome-wide significance claims."
                ),
            },
            {
                "Step": "Fixed-band statistics",
                "Definition": (
                    "Each band receives a B/LB-vs-P/LP ALT-observed 2x2 table, OR, approximate 95% CI, and two-sided "
                    "Fisher test. A Mantel-Haenszel common OR and CMH test summarize across bands. This is a "
                    "sensitivity analysis because residual phyloP differences can remain within a band."
                ),
            },
            {
                "Step": "Continuous analysis",
                "Definition": (
                    f"Firth logistic regression: logit P(B/LB) = intercept + beta_ALT*ALT_observed + natural "
                    f"spline(phyloP100way, df={SPLINE_DF}). exp(beta_ALT) is the adjusted OR. The 95% CI and p-value "
                    "use profile penalized likelihood, not a Wald approximation."
                ),
            },
            {
                "Step": "Continuous estimability",
                "Definition": (
                    "Both clinical classes, both ALT-observed groups, at least four distinct scores, and overlap "
                    "between the groups' observed phyloP ranges are required. Complete outcome separation is handled "
                    "by Firth penalization rather than by discarding the view."
                ),
            },
            {
                "Step": "Multiplicity",
                "Definition": (
                    "For each analysis, variant-type, target-context, and consequence selection, "
                    "Benjamini-Hochberg correction is applied across strategies. Band-specific Fisher tests are "
                    "corrected across strategies within the same band."
                ),
            },
            {
                "Step": "INDEL interpretation",
                "Definition": (
                    "The fixed thresholds have their nominal single-base p-value interpretation only for SNVs. "
                    "INDEL views apply the same bands to an aggregate score for descriptive comparability."
                ),
            },
            {
                "Step": "Software",
                "Definition": (
                    f"R {versions.get('R', 'not recorded')}; logistf {versions.get('logistf', 'not recorded')}. "
                    "The report does not fall back to ordinary maximum-likelihood logistic regression."
                ),
            },
        ]
    )


def negative_control_method_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Step": "Focal sample",
                "Definition": (
                    "Normalized GAPH SNVs are selected by a stable hash independently per strategy, up to the "
                    "configured engineering cap. All eligible SNVs are used below the cap."
                ),
            },
            {
                "Step": "Target-space matching",
                "Definition": (
                    "Each focal SNV is matched to up to five SNVs from the same gene and target context, with the "
                    "same genomic REF>ALT substitution and the same primary RefSeq VEP consequence. A control is "
                    "excluded when the same GAPH strategy observed it."
                ),
            },
            {
                "Step": "VEP consequence",
                "Definition": (
                    "Ensembl VEP uses RefSeq transcripts and pick_allele_gene. The target Entrez Gene ID is "
                    "selected; the most severe Sequence Ontology term on the picked transcript is the matching key."
                ),
            },
            {
                "Step": "Outcomes",
                "Definition": (
                    "The same matched sets compare phyloP100way, exact-allele gnomAD overlap and AF, and exact-allele "
                    "ClinVar overlap and class composition. Conservation and external evidence are outcomes, not "
                    "matching variables."
                ),
            },
            {
                "Step": "Resampling",
                "Definition": (
                    "Each iteration resamples matched sets with replacement and selects one available control from "
                    "each selected set. GAPH, matched-control, and paired-difference statistics use the same draws. "
                    "The report shows descriptive 95% paired matched-set bootstrap intervals and no inferential "
                    "p-value. ClinVar class proportions exclude records with missing CLNSIG; failed gnomAD regions "
                    "remain missing rather than absent."
                ),
            },
        ]
    )


def feature_coverage_formula_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Metric": "Per-feature breadth",
                "Formula": "coverage_breadth = covered_bases / length_bp",
                "Notes": "Read directly from alignment/feature_coverage.tsv.gz.",
            },
            {
                "Metric": "Weighted breadth",
                "Formula": "sum(covered_bases) / sum(length_bp)",
                "Notes": "Length-weighted aggregate used in the target-bases-covered plot.",
            },
            {
                "Metric": "Per-feature mean depth",
                "Formula": "mean_depth = depth_bases / length_bp",
                "Notes": "Uses pipeline-provided depth_bases; the report does not recompute depth from raw aligner output.",
            },
            {
                "Metric": "Weighted mean ortholog depth",
                "Formula": "sum(depth_bases) / sum(length_bp)",
                "Notes": "Main ortholog-depth plot metric.",
            },
            {
                "Metric": "Median feature metrics",
                "Formula": "median(coverage_breadth), median(mean_depth), median(orthologs_covered)",
                "Notes": "Computed over feature rows within each strategy and feature type.",
            },
            {
                "Metric": "Main feature classes",
                "Formula": "CDS, UTR, intron",
                "Notes": "Exon and gene aggregates are omitted because they overlap CDS/UTR/intron and would be redundant.",
            },
        ]
    )


def build_methods_sections(
    inputs: RunInputs,
    out_html: Path,
    variant_summary: VariantSummary,
    cov: pd.DataFrame,
    failures: pd.DataFrame,
    annotation_manifest: dict,
    alignment_manifest: dict,
    validation=None,
    conservation_analysis: ConservationAnalysis | None = None,
    negative_controls: TargetSpaceNullAnalysis | None = None,
    report_timings: list[dict[str, object]] | None = None,
    taxonomy_summary: pd.DataFrame | None = None,
    report_profile_path: Path | None = None,
    candidate_vep_manifest: dict | None = None,
) -> list[str]:
    files = [
        ("Analysis Root", inputs.run_dir),
        ("Fetch Manifest", inputs.fetch_manifest_json),
        ("Variant Annotations", inputs.variant_annotations_tsv),
        ("Variant Strategy Support", inputs.variant_strategy_support_tsv),
        ("Ortholog Evidence Summary", inputs.ortholog_evidence_summary_tsv),
        ("Target Features", inputs.target_features_tsv),
        ("Target Sequences", inputs.target_sequences_dir),
        ("Feature Coverage", inputs.feature_coverage_tsv),
        ("Alignment Segments", inputs.alignment_segments_tsv),
        ("Strategy Summary", inputs.strategy_summary_tsv),
        ("Taxonomy Summary", inputs.taxonomy_summary_tsv),
        ("Annotation Manifest", inputs.annotation_manifest_json),
        ("Alignment Manifest", inputs.alignment_manifest_json),
        (
            "Bulk VEP Manifest",
            inputs.run_dir / "analytics" / "vep_consequences" / "manifest.json",
        ),
        ("Output HTML", out_html),
    ]
    if inputs.cohort_manifest_json is not None:
        files.insert(1, ("Resolved Cohort Manifest", inputs.cohort_manifest_json))
    if validation is not None:
        files.extend(
            [
                ("ClinVar VEP Universe", validation.universe_path),
                ("ClinVar VEP Manifest", validation.manifest_path),
            ]
        )
    if conservation_analysis is not None:
        files.extend(
            [
                ("Candidate phyloP distributions", conservation_analysis.candidate.distributions_path),
                ("Candidate phyloP histograms", conservation_analysis.candidate.histograms_path),
                ("Candidate phyloP manifest", conservation_analysis.candidate.manifest_path),
            ]
        )
    file_rows = [
        {"Key": label, "Path": str(path), "Exists": path.exists(), "Size": file_size_label(path)}
        for label, path in files
    ]

    ok_events = int(annotation_manifest.get("event_key_status_counts", {}).get("ok", 0))
    missing_left_anchor = int(annotation_manifest.get("event_key_status_counts", {}).get("missing_left_anchor", 0))
    vep_config, vep_status = vep_qc_tables(candidate_vep_manifest, validation)
    candidate_vep_ok = int(
        (candidate_vep_manifest or {}).get("status_counts", {}).get("ok", 0)
    )
    clinvar_vep = dict(getattr(validation, "manifest", {}).get("vep", {}))
    clinvar_vep_ok = int(clinvar_vep.get("status_counts", {}).get("ok", 0))
    sections = [
        "<h2>QC</h2>",
        metric_cards(
            [
                (
                    "VEP release",
                    str((candidate_vep_manifest or {}).get("config", {}).get("release", "")),
                ),
                ("Candidate rows VEP ok", format_int(candidate_vep_ok)),
                ("ClinVar alleles VEP ok", format_int(clinvar_vep_ok)),
                ("Event keys normalized", format_int(ok_events)),
                ("Missing left anchor", format_int(missing_left_anchor)),
                ("gnomAD regions failed", format_int(annotation_manifest.get("gnomad_region_failure_count", 0))),
                ("Candidate contexts excluded from gnomAD", format_int(variant_summary.gnomad_lookup_failed)),
                (
                    "ClinVar cached variants",
                    "Not pooled"
                    if inputs.is_cohort
                    else format_int(annotation_manifest.get("clinvar_cached_variant_count", 0)),
                ),
                (
                    "gnomAD cached variants",
                    "Not pooled"
                    if inputs.is_cohort
                    else format_int(annotation_manifest.get("gnomad_cached_variant_count", 0)),
                ),
                ("Feature coverage rows", format_int(len(cov))),
            ]
        ),
    ]
    if inputs.is_cohort and inputs.cohort_manifest_json is not None:
        resolved_cohort = read_json(inputs.cohort_manifest_json)
        provenance_rows = [
            {
                "Run": str(member.get("label", "")),
                "Run directory": str(member.get("run_dir", "")),
                "Scientific fingerprint": str(member.get("fingerprint", "")),
                "Accepted genes": int(member.get("requested_gene_count", 0)),
                "Target genes": int(member.get("target_gene_count", 0)),
            }
            for member in resolved_cohort.get("members", [])
            if isinstance(member, dict)
        ]
        sections.extend(
            [
                "<details open><summary>Cohort provenance</summary>",
                "<p class=\"lead\">All listed completed runs were validated before "
                "pooling. Accepted gene IDs are disjoint; scientific statistics are "
                "recomputed over their union.</p>",
                table_html(
                    pd.DataFrame(provenance_rows),
                    classes="table table-sm table-striped",
                ),
                "</details>",
            ]
        )
    sections.append("<details><summary>VEP annotation coverage</summary>")
    sections.append(
        "<p class=\"lead\">Candidate and ClinVar consequence analyses use release-pinned "
        "Ensembl VEP with RefSeq transcripts and the target-gene consequence selection rule. "
        "A finalized artifact covers every input row; non-<code>ok</code> statuses remain "
        "explicit in the report QC.</p>"
    )
    sections.append(table_html(vep_config, classes="table table-sm table-striped"))
    sections.append(table_html(vep_status, classes="table table-sm table-striped"))
    sections.append("</details>")
    if not failures.empty:
        sections.append("<h3>Annotation Failures</h3>")
        sections.append(table_html(failures, classes="table table-sm table-striped", max_rows=50))
    if taxonomy_summary is not None and not taxonomy_summary.empty:
        shown = taxonomy_summary.copy()
        shown["Taxonomic scope"] = shown["taxonomic_scope"].map(
            lambda value: TAXONOMIC_SCOPE_LABELS.get(str(value), str(value))
        )
        shown["Evidence unit"] = shown["evidence_unit"].map(
            lambda value: EVIDENCE_UNIT_LABELS.get(str(value), str(value))
        )
        shown = shown.rename(
            columns={
                "gene_count": "Genes",
                "ortholog_count": "Selected ortholog rows",
                "taxon_count": "Distinct taxa",
                "unit_count": "Distinct units",
                "orthologs_per_gene_median": "Median orthologs/gene",
                "units_per_gene_median": "Median units/gene",
            }
        )[
            [
                "Taxonomic scope",
                "Evidence unit",
                "Genes",
                "Selected ortholog rows",
                "Distinct taxa",
                "Distinct units",
                "Median orthologs/gene",
                "Median units/gene",
            ]
        ]
        sections.append("<details><summary>Taxonomic evidence scope and grouping</summary>")
        sections.append(table_html(shown, classes="table table-sm table-striped"))
        sections.append("</details>")
    pathogenic_table = pathogenic_variant_table(variant_summary.pathogenic_rows)
    if not pathogenic_table.empty:
        shown = min(len(pathogenic_table), 100)
        sections.append(
            "<details><summary>Top "
            f"{format_int(shown)} of {format_int(variant_summary.pathogenic_variant_count)} "
            "unique P/LP variants</summary>"
        )
        sections.append(
            "<p>Sorted by ClinVar review stars, then supporting SCV count.</p>"
        )
        sections.append(
            table_html(pathogenic_table, classes="table table-sm table-striped", max_rows=100)
        )
        sections.append("</details>")
    if report_timings:
        sections.append("<details><summary>Report computation timing</summary>")
        sections.append(
            "<p>Wall and CPU durations describe this invocation. Process peak RSS is the "
            "high-water mark reached by the report process at the end of each stage, not "
            "memory allocated exclusively by that stage.</p>"
        )
        if report_profile_path is not None:
            sections.append(
                "<p>Detailed nested profile: <code>"
                f"{html.escape(str(report_profile_path))}</code>.</p>"
            )
        sections.append(
            table_html(
                pd.DataFrame(report_timings),
                classes="table table-sm table-striped",
            )
        )
        sections.append("</details>")
    sections.append("<details><summary>Input files and loaded row counts</summary>")
    sections.append(table_html(pd.DataFrame(file_rows), classes="table table-sm table-striped"))
    sections.append(
        table_html(
            pd.DataFrame(
                [
                    {"Metric": "Variant-context rows read", "Value": variant_summary.input_row_count},
                    {"Metric": "Unique candidate variants", "Value": variant_summary.unique_variant_count},
                    {"Metric": "Strategy-supported variant records", "Value": variant_summary.strategy_record_count},
                    {"Metric": "Feature coverage rows loaded", "Value": len(cov)},
                    {"Metric": "Annotation failure rows", "Value": len(failures)},
                    {"Metric": "Alignment event mode", "Value": alignment_manifest.get("alignment_event_mode", "")},
                ]
            ),
            classes="table table-sm table-striped",
        )
    )
    sections.append("</details>")
    sections.append("<details><summary>ClinVar class mapping</summary>")
    sections.append(
        "<p class=\"lead\">The report collapses raw <code>clinvar_sig</code> values into conservative plotting classes.</p>"
    )
    sections.append(table_html(clinvar_class_mapping_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    sections.append("<details><summary>ClinVar review stars</summary>")
    sections.append(
        "<p class=\"lead\">Review-star plots use the normalized star value written during annotation from ClinVar review status.</p>"
    )
    sections.append(table_html(clinvar_review_star_mapping_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    sections.append(f"<details><summary>{variant_summary.consequence_source} consequence grouping</summary>")
    sections.append(
        f"<p class=\"lead\">The Candidate Profile consequence plots use "
        f"{variant_summary.consequence_source} annotations and group them as follows. "
        "Only the release-pinned VEP annotations define these groups.</p>"
    )
    sections.append(
        table_html(
            consequence_grouping_table(variant_summary.consequence_source),
            classes="table table-sm table-striped",
        )
    )
    sections.append("</details>")
    sections.append("<details><summary>Target-context assignment</summary>")
    sections.append(
        "<p class=\"lead\">Candidate and ClinVar validation positions are assigned to one exclusive target context. "
        "Overlapping transcript features use the precedence CDS &gt; UTR &gt; exon &gt; intron; exon sequence "
        "outside CDS/UTR is labelled Other exon, and remaining target sequence is labelled Other. "
        "The ClinVar Association selector exposes All, CDS, UTR, and Intron; Other exon and Other remain included "
        "in the All denominator.</p>"
    )
    sections.append("</details>")
    sections.append("<details><summary>ClinVar validation denominator and statistics</summary>")
    if validation is not None:
        manifest = validation.manifest
        cohort_flow = pd.DataFrame(
            [
                {"Cohort step": "Raw SNV/INDEL alleles", "Alleles": manifest.get("raw_allele_count", 0)},
                {"Cohort step": "Excluded VUS", "Alleles": manifest.get("excluded_vus_count", 0)},
                {"Cohort step": "Excluded missing CLNSIG", "Alleles": manifest.get("excluded_missing_count", 0)},
                {"Cohort step": "Excluded other/conflicting", "Alleles": manifest.get("excluded_other_count", 0)},
                {"Cohort step": "Included B/LB", "Alleles": manifest.get("benign_count", 0)},
                {"Cohort step": "Included P/LP", "Alleles": manifest.get("pathogenic_count", 0)},
                {"Cohort step": "Final validation cohort", "Alleles": manifest.get("usable_allele_count", 0)},
            ]
        )
        sections.append(table_html(cohort_flow, classes="table table-sm table-striped"))
    sections.append(table_html(validation_method_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    sections.append("<details><summary>ClinVar association modes and consequence subsets</summary>")
    sections.append(
        "<p class=\"lead\">All three modes use the same normalized ClinVar allele cohort. They differ only "
        "in whether and how phyloP100way is included.</p>"
    )
    sections.append(table_html(conservation_validation_method_table(conservation_analysis), classes="table table-sm table-striped"))
    consequence_source = validation.consequence_source if validation is not None else "Ensembl VEP"
    sections.append(f"<h4>{consequence_source} consequence subsets</h4>")
    sections.append(
        table_html(
            validation_consequence_grouping_table(consequence_source),
            classes="table table-sm table-striped",
        )
    )
    if conservation_analysis is not None:
        visibility_summary, hidden_views = hidden_clinvar_association_views(
            conservation_analysis.validation
        )
        sections.append("<h4>Adaptive selector visibility</h4>")
        sections.append(
            "<p>Consequence options are hidden from the interactive view only when no strategy has an "
            "estimable result for the selected analysis, variant type, and target context. Hidden combinations "
            "remain listed here; no minimum sample-size threshold is applied.</p>"
        )
        sections.append(
            table_html(visibility_summary, classes="table table-sm table-striped")
        )
        if not hidden_views.empty:
            sections.append(
                table_html(hidden_views, classes="table table-sm table-striped")
            )
    sections.append("</details>")
    sections.append("<details><summary>Candidate-wide phyloP stratification</summary>")
    sections.append(
        table_html(
            pd.DataFrame(
                [
                    {
                        "Step": "Eligible candidate alleles",
                        "Definition": "Normalized lookup_status=ok SNVs and indels with a defined phyloP score basis.",
                    },
                    {
                        "Step": "Allele score",
                        "Definition": (
                            "SNV substituted base; deletion mean across deleted reference bases without VCF padding; "
                            "insertion mean across the two flanking bases. All required bases must have a score."
                        ),
                    },
                    {
                        "Step": "Unit and strata",
                        "Definition": (
                            "Unique variant_key x strategy records with a completed gnomAD lookup, split by presence "
                            "of an exact gnomAD AF annotation. Failed lookups are excluded from both strata."
                        ),
                    },
                    {
                        "Step": "Distribution",
                        "Definition": (
                            "Exact percentiles from 0 through 100, plus relative-frequency histograms using shared "
                            "Found/Not-found bins selected by the Freedman-Diaconis rule and capped at 80 bins. "
                            "Box plots use Tukey 1.5-IQR whiskers."
                        ),
                    },
                    {
                        "Step": "Shared read",
                        "Definition": (
                            "A cold report reads the union of candidate and ClinVar-required positions from bigWig once; "
                            "both analyses reuse that positional score map."
                        ),
                    },
                ]
            ),
            classes="table table-sm table-striped",
        )
    )
    sections.append("</details>")
    sections.append("<details><summary>Negative-control construction</summary>")
    sections.append(
        "<p class=\"lead\">The target-space null uses normalized SNVs and preserves gene, target context, "
        "substitution, and allele-specific functional consequence.</p>"
    )
    sections.append(table_html(negative_control_method_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    sections.append("<details><summary>Feature coverage formulas</summary>")
    sections.append(
        "<p class=\"lead\">Candidate Profile coverage plots use the normalized feature-level table emitted by the alignment stage.</p>"
    )
    sections.append(table_html(feature_coverage_formula_table(), classes="table table-sm table-striped"))
    sections.append("</details>")
    if validation is not None:
        validation_files = [
            ("ClinVar universe", validation.universe_path),
            ("ClinVar universe manifest", validation.manifest_path),
        ]
        if validation.observed_memberships_path is not None:
            validation_files.append(
                ("ClinVar observed memberships", validation.observed_memberships_path)
            )
        if validation.observed_memberships_manifest_path is not None:
            validation_files.append(
                (
                    "ClinVar observed-membership manifest",
                    validation.observed_memberships_manifest_path,
                )
            )
        regions_bed = validation.manifest.get("regions_bed", "")
        if regions_bed:
            validation_files.append(("ClinVar target regions", Path(regions_bed)))
        sections.append("<details><summary>Validation cache files</summary>")
        sections.append(
            table_html(
                pd.DataFrame(
                    [
                        {"Key": label, "Path": str(path), "Exists": path.exists(), "Size": file_size_label(path)}
                        for label, path in validation_files
                    ]
                ),
                classes="table table-sm table-striped",
            )
        )
        sections.append("</details>")
    if conservation_analysis is not None:
        conservation_files = [
            ("Conservation allele annotations", conservation_analysis.annotations_path),
            ("Conservation annotation manifest", conservation_analysis.manifest_path),
        ]
        sections.append("<details><summary>Conservation cache files</summary>")
        sections.append(
            table_html(
                pd.DataFrame(
                    [
                        {"Key": label, "Path": str(path), "Exists": path.exists(), "Size": file_size_label(path)}
                        for label, path in conservation_files
                    ]
                ),
                classes="table table-sm table-striped",
            )
        )
        track_rows = [
            {
                "Track": item.get("track", ""),
                "Status": item.get("status", ""),
                "Annotated positions": item.get("annotated_positions", ""),
                "Unique positions": item.get("unique_positions", ""),
                "Annotated alleles": item.get("annotated_variants", ""),
                "Missing alleles": item.get("missing_variants", ""),
                "Blocks": item.get("block_count", ""),
                "Failed blocks": item.get("failed_block_count", ""),
                "Open seconds": item.get("open_seconds", ""),
                "Read seconds": item.get("read_seconds", ""),
                "Error": item.get("error", ""),
                "URL": item.get("url", ""),
            }
            for item in conservation_analysis.manifest.get("tracks", [])
        ]
        if track_rows:
            sections.append(table_html(pd.DataFrame(track_rows), classes="table table-sm table-striped"))
        sections.append("</details>")
    if negative_controls is not None:
        sections.extend(build_target_space_null_qc_sections(negative_controls))
        control_files = [
            ("Negative-control manifest", negative_controls.manifest_path),
            ("Target-space-null rows", negative_controls.matched_path),
            ("Target-space-null phyloP annotations", negative_controls.conservation_path),
            ("VEP consequence cache", negative_controls.vep_cache_path),
            ("Target-space-null external evidence", negative_controls.external_evidence_path),
            ("External-evidence manifest", negative_controls.external_evidence_manifest_path),
        ]
        if negative_controls.focal_path is not None:
            control_files.append(
                ("Target-space-null focal sample", negative_controls.focal_path)
            )
        if negative_controls.focal_manifest_path is not None:
            control_files.append(
                (
                    "Target-space-null focal manifest",
                    negative_controls.focal_manifest_path,
                )
            )
        sections.append("<details><summary>Negative-control cache files</summary>")
        sections.append(
            table_html(
                pd.DataFrame(
                    [
                        {"Key": label, "Path": str(path), "Exists": path.exists(), "Size": file_size_label(path)}
                        for label, path in control_files
                    ]
                ),
                classes="table table-sm table-striped",
            )
        )
        sections.append("</details>")
    return sections
