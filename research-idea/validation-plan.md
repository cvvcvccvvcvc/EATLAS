# GAPH Allele-Specific Validation Plan

## Purpose

This document turns the research conversation into a concrete validation plan
for the next stages of the GAPH project.

The central claim is not that GAPH reimplements conservation. The useful claim
is narrower:

> Exact orthologous observation of the human alternative allele contains
> allele-specific tolerance information beyond ordinary site-level
> conservation.

Conservation scores such as phyloP, phastCons, and GERP ask:

```text
Is this genomic position constrained?
```

GAPH allele features ask:

```text
Has this exact alternative allele been observed at the orthologous position?
```

These are related but not identical questions. The first validation goal is to
show that this distinction matters in ClinVar benign/pathogenic variants.

## Main Dataset

The initial validation dataset should be deliberately simple:

- GRCh38 SNVs only.
- ClinVar benign / likely benign versus pathogenic / likely pathogenic.
- VUS, conflicting interpretations, drug response, association-only, risk
  factor, and unsupported labels excluded.
- Variants restricted to genes already fetched and aligned by GAPH.
- Prefer variants with at least some ClinVar review support.

The initial label convention:

```text
benign      = tolerated class
pathogenic  = deleterious class
```

The first analysis does not need CADD. CADD becomes relevant later if we want
to claim improvement over a strong existing pathogenicity score.

## Required Pipeline Outputs

The current pipeline already emits most of what is needed:

- `alignment_segments.tsv.gz`: which orthologs cover target positions.
- `alignment_events.tsv.gz`: which orthologs carry non-reference events.
- `ortholog_alignment_summary.tsv.gz`: ortholog universe per gene/strategy.
- `taxonomy_presets.tsv.gz`: broad taxonomy groups.
- `target_features.tsv.gz`: target-local feature context.
- `feature_coverage.tsv.gz`: coverage summaries by target feature.

The validation feature builder should produce one row per variant and strategy,
including:

- `gaph_all_depth`
- `gaph_all_alt_count`
- `gaph_all_alt_fraction`
- `gaph_all_ref_count`
- `gaph_all_ref_fraction`
- `gaph_all_other_count`
- `gaph_all_indel_count`
- `gaph_all_no_call_count`
- `gaph_all_callable_fraction`
- `gaph_all_entropy`
- the same fields for taxonomy groups when available
- target feature context, for example CDS, UTR, intron, exon, gene

The first binary feature is derived from existing columns:

```text
ALT_observed = gaph_all_alt_count > 0
```

It may still be useful to materialize this column in reports, because it is the
main biological variable in the first experiments.

## Experiment 1: Global ALT-Observed Enrichment

### Hypothesis

ClinVar benign SNVs are more likely than ClinVar pathogenic SNVs to have their
ALT allele observed in at least one ortholog.

### Test

Build a 2 by 2 table:

```text
                         ClinVar benign    ClinVar pathogenic
ALT_observed = 1                a                    b
ALT_observed = 0                c                    d
```

Compute:

- odds ratio
- Fisher exact test p-value
- 95% confidence interval for the odds ratio

### Expected Result

The desired direction is:

```text
OR > 1
```

Interpretation:

> Variants whose ALT allele has already been observed in orthologous sequence
> are enriched among ClinVar benign variants.

This is only the first sanity check. It shows that GAPH contains biological
signal, but it does not yet show that the signal is independent of
conservation.

## Experiment 2: Conservation-Stratified Enrichment

### Rationale

The global enrichment test could be explained by ordinary conservation:

```text
low-conservation positions -> more ortholog variation -> more benign variants
```

If that is the whole story, GAPH is mostly rediscovering site-level
conservation. The key scientific question is whether ALT observation remains
informative among variants with comparable conservation.

### Hypothesis

Within the same conservation stratum, ClinVar benign SNVs are still more likely
than pathogenic SNVs to have `ALT_observed = 1`.

### Inputs

At least one conservation score is needed:

- phyloP
- phastCons
- GERP

The first version can use whichever one is easiest to obtain reliably. Later,
the same analysis can be repeated for multiple scores.

### Binning

Split variants into conservation bins, for example:

```text
low conservation
medium conservation
high conservation
very high conservation
```

The bins can be quantile-based for balanced sample sizes, or biologically
thresholded if a score-specific convention is used. Quantile bins are safer for
the first implementation because they avoid sparse tables.

### Test

For each conservation bin, repeat the same 2 by 2 table:

```text
                         ClinVar benign    ClinVar pathogenic
ALT_observed = 1                a                    b
ALT_observed = 0                c                    d
```

Compute within each bin:

- odds ratio
- Fisher exact p-value
- 95% confidence interval
- number of variants in each cell

Optionally compute a pooled stratified odds ratio with a
Cochran-Mantel-Haenszel test. This is conceptually aligned with the 2 by 2
analysis and is easier to explain than jumping directly to a logistic model.

### Expected Result

The strongest outcome:

```text
OR > 1 in multiple conservation bins
```

especially if the effect persists in high-conservation bins.

Interpretation:

> Even among variants with similar site-level conservation, orthologous
> observation of the exact ALT allele is enriched in benign variants.

This supports the main claim:

> GAPH captures allele-specific evolutionary tolerance, not only generic
> conservation.

## Coverage And Confidence Controls

Coverage analysis is not a separate proof of the main hypothesis. It protects
the enrichment result from a technical artifact.

### Problem

`ALT_observed = 0` can mean two different things:

```text
The ALT allele was not observed despite good ortholog coverage.
```

or:

```text
The position had too few callable orthologs to observe much of anything.
```

Likewise, `ALT_observed = 1` becomes easier to obtain when many orthologs are
covered. With more observations, the chance of seeing some alternative state is
higher.

### Required Summaries

Every enrichment report should show, separately for benign and pathogenic
variants:

- distribution of `gaph_all_depth`
- distribution of `gaph_all_callable_fraction`
- number of zero-depth or low-depth variants
- fraction of variants included after each quality filter

### High-Confidence Subset

Repeat Experiment 1 and Experiment 2 on a high-confidence subset, for example:

```text
gaph_all_depth >= 20
gaph_all_callable_fraction >= 0.10
```

The exact thresholds should be reported and can be tuned after inspecting the
depth distribution.

### Interpretation

Good outcome:

```text
The OR remains > 1 after coverage filtering.
```

Stronger outcome:

```text
The OR increases in the high-confidence subset.
```

Concerning outcome:

```text
The effect exists only in low-depth or low-callable regions.
```

In that case, the signal may be an alignment or ascertainment artifact rather
than a robust biological effect.

## Dose-Response Feature Development

After the binary `ALT_observed` test, the next question is whether stronger
ortholog support gives stronger evidence of benignity.

### Candidate Variables

Instead of only:

```text
ALT_observed = gaph_all_alt_count > 0
```

test:

- `gaph_all_alt_count`
- `gaph_all_alt_fraction`
- clade-specific ALT support, for example primates or other mammals
- entropy over REF / ALT / other / indel states

### Simple Dose-Response Bins

Example count bins:

```text
0 ALT-supporting orthologs
1 ALT-supporting ortholog
2-4 ALT-supporting orthologs
5+ ALT-supporting orthologs
```

Example fraction bins:

```text
0
low
medium
high
```

For each bin, report:

- benign fraction
- pathogenic fraction
- odds ratio versus the zero-support group
- confidence interval

### Expected Result

The desired pattern:

```text
higher ALT support -> higher benign enrichment
```

This would support using richer GAPH features rather than only a binary
indicator.

## Region-Specific Analysis

Region-specific analysis is not the first proof of the hypothesis. It answers a
different question:

> Where does the method work best, and where is it weak?

Useful strata:

- coding sequence
- UTR
- intronic
- splice-adjacent regions, if available
- whole gene background

Expected pattern:

- coding and splice-adjacent variants may show clearer signal
- UTR and intronic variants may be noisier
- some noncoding regions may suffer from weaker orthology and turnover

This analysis should be presented after the global and conservation-stratified
tests.

## Later Extensions

These are valuable but should not block the first validation report.

### Phylogenetic Weighting

`ALT_observed` treats all orthologs equally. This is biologically crude:

```text
20 closely related rodents are not 20 independent evolutionary observations.
```

Later features should account for:

- broad clade support
- nearest species carrying ALT
- independent ALT emergence events
- phylogenetic depth of ALT support

The first approximation can use existing taxonomy groups. A more rigorous
version would use a species tree.

### Alignment Quality Weighting

ALT support from a weak or ambiguous alignment should not count the same as ALT
support from a high-confidence alignment.

Later features can include:

- identity-weighted support
- MAPQ-weighted support
- exclusion or downweighting of QC-flagged alignments
- strategy agreement across minimap2, nucmer, and BWA-derived evidence

### CADD And Predictive Benchmarking

CADD is not required for the first two enrichment tests.

It becomes important if the claim changes from:

```text
GAPH carries allele-specific biological signal.
```

to:

```text
GAPH improves pathogenicity prediction beyond existing predictors.
```

That later benchmark should compare:

```text
baseline features only
GAPH features only
baseline + GAPH features
baseline + shuffled GAPH features
```

The baseline can include CADD and conservation scores. The main metrics should
be AUROC, AUPRC, and top-k enrichment, preferably with gene-held-out or
chromosome-held-out splits.

### Adjusted Logistic Regression

Adjusted logistic regression is a later optional summary model, not the first
central argument.

It can estimate an adjusted odds ratio for `ALT_observed` while controlling for
continuous variables:

```text
coverage depth
callable fraction
conservation score
region type
```

For the immediate plan, conservation-stratified 2 by 2 enrichment is more
transparent and should be prioritized.

## First Report To Build

The first report script should consume the variant-level GAPH feature table and
produce an HTML or Markdown report with:

1. Dataset summary.
2. Global `ALT_observed` enrichment table.
3. OR, Fisher exact p-value, and 95% CI.
4. Conservation-stratified enrichment tables, if conservation columns are
   present.
5. Coverage and callable-depth distributions.
6. High-confidence subset enrichment.
7. Dose-response summaries for ALT count and ALT fraction.
8. Region-specific summaries, if feature context is present.

The report should also write machine-readable TSV/JSON outputs so results can
be reused in the thesis and later plots.

## Decision Gates

### Gate 1

Proceed if:

```text
global ALT_observed OR > 1
```

and the confidence interval is reasonably above or not strongly crossing 1.

### Gate 2

Proceed to stronger claims if:

```text
ALT_observed OR remains > 1 within conservation bins
```

or a pooled stratified test supports an effect after accounting for
conservation strata.

### Gate 3

Proceed to feature engineering if:

```text
dose-response suggests stronger ALT support gives stronger benign enrichment
```

### Gate 4

Proceed to CADD/incremental benchmark only after:

```text
the enrichment signal survives conservation and coverage checks
```

This keeps the project staged: first prove the biological signal, then refine
the feature, then test predictive improvement over existing methods.
