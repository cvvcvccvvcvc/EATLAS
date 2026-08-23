# Validation Design

## Scientific Hypothesis

GAPH is not trying to rediscover conservation alone. The useful question is
whether ortholog-alignment evidence adds variant-level predictive signal after
existing pathogenicity predictors and conservation annotations are known.

Null hypothesis:

```text
GAPH-derived evolutionary features do not improve held-out prediction metrics
after accounting for CADD and conservation features.
```

Alternative hypothesis:

```text
GAPH-derived evolutionary features provide independent signal, measurable as
improved held-out AUROC/AUPRC and enrichment metrics in ablation experiments.
```

## Baseline Context

CADD is the strongest practical first baseline because it already aggregates
many annotations and is available as public scores and offline scripts:

- official CADD information: https://cadd.gs.washington.edu/info
- official CADD downloads: https://cadd.gs.washington.edu/download
- CADD offline scripts: https://github.com/kircherlab/CADD-scripts

The validation should initially use CADD as a frozen external score rather than
trying to reproduce the full CADD training stack. Full CADD reproduction is a
larger infrastructure project because the offline annotations and prescored
files are large and version-specific.

## Dataset

Recommended first target:

- GRCh38 SNVs
- ClinVar high-confidence pathogenic/likely pathogenic versus benign/likely
  benign labels
- variants restricted to genes successfully fetched and aligned by GAPH
- VUS and conflicting classifications excluded from supervised evaluation

The first pass should avoid mixing train/test variants from the same gene when
possible. A random variant split can overstate performance because nearby
variants share annotations, gene context, and sometimes label ascertainment.

## GAPH Feature Rationale

The current pipeline preserves three complementary Stage 2 evidence relations
inside each partition:

- `alignment_segments.tsv.gz`: which orthologs cover target intervals
- `alignment_events.tsv.gz`: compact observed non-reference events in
  target-local and GRCh38 coordinates
- `event_ortholog_support.tsv.gz`: the exact positive supporters of each event

For a candidate variant, this lets us distinguish:

```text
ortholog covers position and has exact ALT event -> ortholog supports ALT
ortholog covers position with no event           -> ortholog supports human REF
ortholog covers position with another event      -> ortholog supports other
ortholog does not cover position                 -> no-call
```

Canonical variant identity is supplied by the Stage 3 event map; taxonomy is
joined from the canonical Stage 1 table. The archived feature builder predates
this contract and must be migrated before the experiment is rerun.

This is different from a pure conservation score because the features can encode
allele direction, clade distribution, and alignment support. Without explicit
ancestral-state reconstruction, "direction" should be interpreted as phylogeny-
stratified support rather than a definitive ancestral/derived claim.

## Feature Groups

Initial features should be per variant and per alignment strategy:

- callable ortholog depth
- REF support count and fraction
- ALT support count and fraction
- other-SNV support count and fraction
- indel-overlap support count and fraction
- Shannon entropy over REF/ALT/other/indel support classes
- the same counts/fractions stratified by taxonomy groups:
  primates, other mammals, non-mammal vertebrates, and other/unknown
- optional full ortholog universe and no-call count when summaries are supplied
- target feature context and feature-level coverage/depth when
  `target_features.tsv.gz` and `feature_coverage.tsv.gz` are supplied

Later extensions:

- ancestral-state inference using explicit outgroups
- strategy agreement features across minimap2/nucmer/BWA comparators
- alignment quality weighting using identity, MAPQ, and QC flags
- protein or transcript-aware consequences

BWA pseudoreads should be interpreted separately in early experiments because
the current comparator can emit event rows without segment rows. Event-only
evidence can support an observed ALT/other event, but cannot establish REF
support at positions without events.

## Ablation

Run the same model and split on these feature sets:

| Name | Features |
| --- | --- |
| `baseline` | CADD/conservation/other existing annotations |
| `gaph` | only `gaph_*` columns |
| `baseline_plus_gaph` | baseline plus `gaph_*` |
| `baseline_plus_shuffled_gaph` | baseline plus row-shuffled `gaph_*` negative control |

Primary metrics:

- AUROC
- AUPRC
- fold improvement in AUPRC over baseline

Secondary checks:

- calibration
- top-k precision or pathogenic enrichment
- per-consequence and per-gene-family slices when labels are large enough

## Leakage Controls

Important failure modes:

- variant-level random splits that place many same-gene variants in both train
  and test
- using labels or annotations from a newer source in a baseline score trained on
  overlapping data without a date-aware interpretation
- comparing public CADD scores against labels that were likely available during
  CADD model development
- letting variant consequence or gene identity dominate the result without
  checking gene-held-out performance
- treating missing GAPH coverage as evidence of reference support

The minimum defensible report should include random CV only as a smoke check and
gene- or chromosome-aware CV as the main result.
