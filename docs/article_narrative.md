# GAPH Article Narrative

This document defines the main claim, validation logic, and unresolved design
decisions for the GAPH study.

## Central Hypothesis

> Observation of the exact human alternative allele in orthologous sequence is
> an allele-specific evolutionary signal that may identify tolerated human
> variants beyond the information summarized by site-level conservation.

The study asks whether this signal is associated with benign ClinVar
classification and whether it adds information beyond site-level conservation.
Ortholog observation alone does not prove clinical benignity or causality.

## Main Results Narrative

### 1. Define the signal

GAPH aligns orthologous gene sequences and records whether an alignment reports
the exact human ALT allele:

```text
ALT_observed = at least one ortholog alignment reports the exact human ALT
```

This is a positive-evidence feature. `ALT_observed=0` means only that no exact
ALT was reported. It does not imply REF support, positional constraint, or ALT
intolerance. Alignment completeness and evidence depth are reported as QC, not
as biological outcome classes.

**Main evidence:** workflow diagram and a compact dataset summary.

### 2. Characterize where GAPH produces signal

Before testing clinical enrichment, describe:

- consequence composition;
- conservation scores of GAPH candidates and of a target-space null matched by
  gene, target context, exact substitution, and allele-specific consequence;
- taxonomic proximity and the number of independent taxonomic evidence units
  carrying the ALT.

Define a high-confidence tier from the taxonomic evidence before examining its
ClinVar association. Keep the complete signal in descriptive results and state
the tier rule explicitly.

An excess of synonymous variants or a shift toward less constrained positions
would show where the method tends to find candidates. Neither result would
validate the method by itself.

**Main evidence:** sorted consequence composition and conservation-distribution
figures with denominators shown.

### 3. Test association with ClinVar interpretation

Test whether ClinVar B/LB alleles are more likely than P/LP alleles to have the
exact ALT observed in orthologs.

Report the broad eligible-SNV analysis first. Then stratify independently by
exclusive target-locus contexts from the fetched NCBI features and by
release-pinned RefSeq VEP consequence groups. Keep ClinVar MC and gnomAD
consequence fields as provenance, but do not mix them into VEP-defined strata.
Compute an odds ratio only where both B/LB and P/LP are represented. Missense
SNVs are the primary consequence class unless the final cohort supports another
class chosen in advance.

For each strategy, restrict the ClinVar denominator to genes for which that
strategy produced an alignment result. Within those genes, absence of a
reported exact ALT is `ALT_observed=0` with the detection-only interpretation
defined above. Report the full cohort, then compare high- and lower-confidence
ALT observations separately with the common no-ALT-observed group.

For each analysis, show the four cell counts, odds ratio, 95% confidence
interval, and Fisher exact p-value. Interpret the effect size and its
uncertainty, not significance alone.

**Main evidence:** an interactive forest plot containing the broad estimate and
estimable target-context and consequence-specific estimates.

### 4. Test whether the signal exceeds site-level conservation

The key alternative explanation is:

```text
low constraint -> more ortholog substitutions and more benign human variants
```

Use phyloP as the primary site-level score and model it continuously. Estimate
the association of `ALT_observed` with B/LB after adjusting for nonlinear
phyloP using Firth logistic regression with a natural spline. Report a profile
penalized-likelihood confidence interval and test.

Fixed phyloP bands at signed `-log10(0.05)` thresholds and their
Mantel-Haenszel summary are a descriptive sensitivity analysis, not the primary
proof. The same thresholds in INDEL views apply to aggregate allele scores and
do not retain their nominal single-base p-value interpretation.

**Main evidence:** adjusted OR with 95% confidence interval, accompanied by a
plot showing conservation overlap between `ALT_observed` groups.

The claim is supported only within the conservation range where the groups
overlap. If the adjusted effect disappears, the correct conclusion is that the
binary feature largely reproduces site-level conservation.

## Manuscript Assembly

The Introduction should distinguish candidate-allele evidence from site-level
conservation. Methods define the pipeline, cohort, confidence tier, and tests.
Results follow the four steps above. The Discussion should note that ClinVar
concordance is not independent proof of safety because submitted
classifications may themselves use evolutionary evidence.

## Main Figures And Tables

1. **Figure 1:** GAPH workflow and exact-allele feature definition.
2. **Figure 2:** taxonomic support, confidence tiers, consequence, and
   conservation.
3. **Figure 3:** full-cohort and confidence-tier ClinVar enrichment ORs with
   target-context and consequence sensitivities.
4. **Figure 4:** conservation overlap and adjusted ALT-observed association.
5. **Table 1:** cohort flow and final B/LB and P/LP counts.

## Supplementary Evidence

The supplement should contain:

- normalization exclusions, ClinVar review status, and cohort flow details;
- per-strategy alignment completeness and ortholog-evidence-depth distributions;
- descriptive SNV heatmaps relating site-aligned evidence-unit count and
  absolute exact-ALT support to gnomAD overlap, stratified by CDS, UTR, and
  intron context and selectable taxonomic scope/grouping, accompanied by
  pooled empirical distributions of both evidence counts;
- gene-level contribution and full consequence-class counts; if a few genes
  dominate the estimate, repeat it after removing each gene in turn and use a
  gene-cluster bootstrap only if a cluster-aware confidence interval is needed;
- fixed-band conservation sensitivity results, selector-level cohort counts,
  and target-context assignment rules;
- target-space-null construction, matching yield, consequence-specific counts,
  and descriptive phyloP, gnomAD, and ClinVar comparisons with paired
  matched-set bootstrap intervals for GAPH, controls, and their difference;
- gnomAD found-versus-not-found summaries with failed region lookups excluded
  from both the numerator and denominator;
- strategy overlap, runtime, storage, and alignment QC;
- additional matching or robustness analyses only if overlap diagnostics make
  them necessary.

## Strategy Development

The current alignment strategies are method-development candidates, not final
algorithms. Report them all while comparing behavior, but do not choose a
winner from p-values alone. The final method should use the smallest strategy
set that provides stable alignment completeness, nonredundant allele evidence,
acceptable runtime, and reproducible validation.

If ClinVar results influence strategy selection, the selected strategy must be
confirmed on data not used for that selection.

## Data Resource Direction

A searchable database is a later dissemination product, not evidence for the
current validation. Its core should expose normalized human alleles, target
context, supporting taxa and taxonomic evidence units, strategy provenance,
confidence tier, and versioned annotations. Implementation should follow the
final confidence and strategy contracts so the public schema does not freeze
provisional definitions.

## Open Decisions

- final validation gene panel and its prespecified inclusion rule;
- primary strategy or minimal complementary strategy set;
- final ClinVar review-status threshold after cohort-size assessment;
- prespecified high-confidence signal rule based on taxonomic evidence;
- scope and release model for the downstream searchable database;
- whether later work tests added value over CADD or trains a predictor.

These decisions may change as evidence accumulates. Changes should simplify or
strengthen the central argument, not add parallel hypotheses.
