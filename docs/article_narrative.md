# GAPH Article Narrative

This is the evolving scientific spine of the GAPH study. It defines what the
main paper must establish and keeps analytics work aligned with one coherent
argument. It is the canonical source for the study's current hypotheses,
validation logic, and unresolved scientific decisions.

## Central Hypothesis

> Observation of the exact human alternative allele in orthologous sequence is
> an allele-specific evolutionary signal that may identify tolerated human
> variants beyond the information summarized by site-level conservation.

The study tests association and incremental information. It does not claim
that ortholog observation proves clinical benignity or constitutes a causal
effect.

## Main Results Narrative

### 1. Define the signal

Introduce GAPH as a reproducible pipeline that aligns orthologous gene
sequences and records whether the exact human REF or ALT allele is observed.
The first feature is deliberately simple:

```text
ALT_observed = at least one callable ortholog carries the exact human ALT
```

The unit of analysis is one human allele, strategy, and gene context.
The annotation layer also preserves the number of distinct ALT-supporting
orthologs for each such unit so that a later dose-response analysis does not
depend on raw alignment-event publication.

**Main evidence:** workflow diagram and a compact dataset summary.

### 2. Characterize where GAPH produces signal

Describe the candidate variants before testing clinical enrichment:

- consequence-class composition;
- conservation-score distribution relative to the callable target background;
- fraction of callable alleles with `ALT_observed = 1` across conservation.

An excess of synonymous variants and a shift toward less constrained or
accelerated positions are biologically expected and scientifically informative.
They define the method's operating domain; they do not by themselves validate
the method or invalidate its usefulness.

**Main evidence:** sorted consequence composition and conservation-distribution
figures with denominators shown.

### 3. Test association with ClinVar interpretation

Ask whether ClinVar B/LB alleles are more likely than P/LP alleles to have the
exact ALT observed in orthologs.

First report the broad eligible-SNV analysis. It summarizes the overall output
of the method but may reflect consequence composition. Then show results by
consequence class. Compute an odds ratio only where both B/LB and P/LP are
represented; missense SNVs are the primary comparable class unless the final
cohort supports another prespecified class.

For each reported analysis show the four cell counts, odds ratio, 95% confidence
interval, and Fisher exact p-value. The scientific emphasis is the effect size
and its uncertainty, not significance alone.

**Main evidence:** one forest plot containing the broad estimate and the
estimable consequence-specific estimates.

### 4. Test whether the signal exceeds site-level conservation

The key alternative explanation is:

```text
low constraint -> more ortholog substitutions and more benign human variants
```

Use phyloP as the prespecified primary site-level score and model it
continuously. Estimate the association of `ALT_observed` with B/LB after
adjusting for nonlinear phyloP and the number of callable orthologs. The final
model must account for gene-level dependence when the validation panel is large
enough to support it. Sparse or separated data require a bias-reduced or
regularized logistic estimate with an appropriate confidence interval.

Quantile-bin and Mantel-Haenszel results are descriptive sensitivity analyses,
not the primary proof. GERP and phastCons are secondary conservation
sensitivities.

**Main evidence:** adjusted OR with 95% confidence interval, accompanied by a
plot showing conservation overlap between `ALT_observed` groups.

The claim is supported only within the conservation range where the groups
overlap. If the adjusted effect disappears, the correct conclusion is that the
binary feature largely reproduces site-level conservation.

## Manuscript Assembly

The Introduction moves from the need for evolutionary tolerance evidence to
the gap between site-level scores and candidate-allele evidence. Methods define
the pipeline, cohort, and tests. Results follow the four steps above without
adding parallel hypotheses. The Discussion interprets the observed operating
domain and states that ClinVar concordance is not independent proof of safety,
because submitted classifications may themselves use evolutionary evidence.

## Main Figures And Tables

1. **Figure 1:** GAPH workflow and exact-allele feature definition.
2. **Figure 2:** consequence and conservation landscape of GAPH candidates.
3. **Figure 3:** broad and consequence-specific ClinVar enrichment ORs.
4. **Figure 4:** conservation overlap and adjusted ALT-observed association.
5. **Table 1:** cohort flow and final B/LB and P/LP counts.

## Supplementary Evidence

Supplementary material protects the main inference without creating competing
stories:

- normalization exclusions, ClinVar review status, and cohort flow details;
- per-strategy callability and ortholog-depth distributions;
- gene-level and full consequence-class counts;
- GERP, phastCons, and quantile-stratified sensitivity results;
- strategy overlap, runtime, storage, and alignment QC;
- additional matching or robustness analyses only if overlap diagnostics make
  them necessary.

## Strategy Development

The current alignment strategies are method-development candidates, not final
algorithms. Report them all while comparing behavior, but do not choose a
winner from p-values alone. The final method should use the smallest strategy
set that provides stable callability, nonredundant allele evidence, acceptable
runtime, and reproducible validation.

If ClinVar results influence strategy selection, the selected strategy must be
confirmed on data not used for that selection.

## Open Decisions

- final validation gene panel and its prespecified inclusion rule;
- primary strategy or minimal complementary strategy set;
- final ClinVar review-status threshold after cohort-size assessment;
- callable-ortholog eligibility threshold;
- whether later work extends from association analysis to a trained predictor.

These decisions may change as evidence accumulates. Changes should simplify or
strengthen the central argument, not add parallel hypotheses.
