# Analytics Contract

Read this document when changing analytics input compatibility, cache identity,
scientific derivations, or report interpretation. Use
`report_generation.md` only for the operator launch procedure.

## Ownership And Input Boundary

`analytics/` is a read-only consumer of one or more completed GAPH runs. It
does not repair pipeline evidence, fetch missing pipeline-owned annotations, or
write below a source run.

Each source run must contain:

- canonical Stage 1 target, selected-ortholog, feature, and taxonomy evidence;
- partitioned Stage 2 summaries, segments, compact events, and exact support;
- the Stage 3 partitioned ClinVar/gnomAD/VEP variant dataset and
  event-to-variant map;
- a complete root run manifest and its bound evidence inventory.

Multi-run analyses require the same current pipeline contracts, target and
reference identities, alignment strategy configuration, ClinVar and gnomAD
contracts, VEP backend/release, and variant columns. Accepted Gene IDs must be
disjoint. Missing evidence, incompatible runs, overlap, obsolete schemas, or an
incomplete shard set fail explicitly.

There is no fallback to old pipeline aggregates, a synthetic cohort, or a
separate bulk-VEP artifact.

## Immutability And Verification

A completed source run is immutable. Before using a source cache, analytics:

1. validates the root manifest and its binding to `evidence_inventory.json`;
2. verifies every inventoried evidence file by content;
3. records a verification marker only under the external analytics root.

Later runs may reuse that marker while the exact file set and filesystem
metadata remain unchanged. A metadata change triggers full verification; a
content or membership mismatch fails.

Large variant shards stay in source runs and are queried as one virtual table.
Analytics materializes only reusable derivations and genuine run-set results.
It does not concatenate source evidence, create target-FASTA symlink trees, or
modify source metadata.

## Workspace And Identity

Every command requires an external `--analytics-root`:

```text
<analytics-root>/
  .strategy_report.lock
  cache/<source-id>/
    evidence_inventory.verified.json
    calculations/<calculation-id>/
      alignment_aggregates/
      annotation_support/
      taxonomy_summary/
  cache/calculations/<calculation-id>/clinvar_conditions/
  cache/reference_identity/
  analyses/<analysis-id>/
    manifest.json
    derived/
    reports/
    performance/
  slurm/
```

- `source-id` depends only on completed pipeline provenance.
- `calculation-id` depends on explicit algorithm/cache versions and the
  Python, R, and computation-library versions actually loaded.
- `analysis-id` depends on the unordered source-ID set, scientific options,
  and calculation identity.
- A report name selects an HTML file within one analysis; it does not create a
  second scientific cache.

Every scientific calculation change must bump its owning cache version, or the
top-level analytics semantics version when no narrower owner exists.
Documentation and HTML-only changes do not require a scientific cache bump.

One process owns an analytics root at a time. A concurrent writer exits rather
than racing shared caches. The lock file is only a rendezvous point; process
exit releases ownership even if the file remains.

A local phyloP BigWig is identified by streamed SHA-256 content. Its cached
filesystem identity avoids repeated hashing while the file is unchanged.
Moving identical content preserves analysis identity; replacing it invalidates
affected artifacts.

## Reusable Derivations

Analytics derives and fingerprints:

- strategy summary and feature coverage from Stage 2 summaries and segments;
- taxonomy summary from selected orthologs and canonical Stage 1 taxonomy;
- variant-strategy support and taxonomic depth/support from exact Stage 2
  evidence plus Stage 3 event lineage;
- `analyses/<analysis-id>/derived/pathogenic_clinvar_hits.tsv.gz`, the complete
  pathogenic/likely-pathogenic allele detail and support table;
- run-set scientific tables used by report sections.

Completed continuous Firth result and distribution tables are cached by their
exact cohort, observed memberships, eligible genes, strategy order, model
definition, and calculation runtime. Rebuilding presentation or retrying a
later failed report stage therefore does not fit identical models again.

Annotation-support preparation is partition-local. Each worker reads only the
variant shards for its evidence partition, uses one DuckDB thread, and
atomically checkpoints its small derived tables. An identical rerun reuses
complete checkpoints and rebuilds missing ones. Partial checkpoint directories
are invalid and completed checkpoints are removed after the final per-source
cache is assembled.

Pipeline manifests may expose integrity counts, but analytics derives scientific
counts from row-level evidence. Report-specific thresholds, bins, taxonomic
scopes, ranks, and plotting tables never become pipeline source data.

## External Evidence Semantics

ClinVar and gnomAD record fields are allele-level. Repeated gene-context rows
are reconciled by canonical `variant_key`: a successful non-empty value is
shared, while conflicting non-empty values or a non-numeric gnomAD allele
frequency fail. VEP consequences, target context, and lookup outcome remain
gene-context evidence.

The indexed ClinVar VCF used for validation must match the portable content
identity recorded by every source run. It may be located at a different path
from the pipeline launch.

Rows with non-`ok` pipeline VEP status remain explicit and appear as
`Not annotated` where a complete denominator is required. The optional
target-space null may issue new VEP and gnomAD requests for generated controls;
it uses the pipeline-pinned VEP release, and its backend/release are part of the
cache identity. Local VEP is probed before an analysis workspace is created.

## Report Semantics

### Firth Model Diagnostics

Each continuous model retains its R warnings in the result reason. Models with
nonconvergence warnings or invalid effects, confidence limits, or p-values are
`not_estimable`; their inferential values are missing and do not enter BH
correction. Other warnings remain visible as `estimated_warning` alongside the
usable estimate. Missing or duplicate model results are an error rather than
silently leaving a pending result.

### Pathogenic ClinVar Hits

Pathogenic and likely-pathogenic subtypes remain distinct in detail tables but
are combined in plots. The tab owns review-strength, molecular-effect,
condition, exact-support, and complete paginated detail views; other tabs do
not repeat them.

Exact-support distributions count one SNV per strategy. When an allele has
several target-gene contexts, the row with maximum ALT-supporting ortholog count
is selected and Gene ID breaks ties. Violin densities use log10 counts while
axes and hover show original counts; sparse distributions use points and a box.
phyloP availability does not restrict this tab.

Conditions compare unique pathogenic/likely-pathogenic alleles with either
ClinVar alleles in eligible target genes or the complete pinned GRCh38 VCF.
Both arms use the selected variant type. Denominators include alleles without a
named condition, and one allele may contribute several conditions, so displayed
fractions need not sum to one. Disease identity prefers MedGen, then MONDO, then
OMIM; names are used only when identifiers are absent.

The whole-VCF background is streamed locally and cached per calculation
identity. It counts distinct VCF alleles, merges repeated records, and excludes
mixed benign/likely-benign and pathogenic/likely-pathogenic classifications.
It describes variants and associated conditions in the VCF, not all clinical
assertion records.

### Basic Filtering

Variant type and filter selection apply to the whole tab. Retention and gnomAD
share a strategy selector; ClinVar has independent strategy, adjustment,
target-context, and consequence selectors.

`Compare all strategies` uses each strategy's own candidate denominator.
`Union (any strategy)` counts each normalized allele once and accepts it when
at least one calling strategy/context passes. Union odds ratios come from a new
allele-level contingency table over the union of eligible genes; they are never
averages of strategy odds ratios.

Exact ALT support and strategy support use minimum thresholds. Supporting
families count distinct known family IDs; an unresolved family is not a new
family. Site-aligned ortholog filters support both maximum and minimum
thresholds. Family and site-depth filters are SNV-only. Within an
allele/strategy, maximum support is used for minimum thresholds and minimum site
depth for maximum thresholds. Non-calls have missing scores and cannot pass a
maximum threshold; zero known families is a valid called-SNV score. Invalid
required support fails validation rather than silently removing rows.

ClinVar thresholds are derived after cohort selection and, for adjusted
analyses, finite phyloP selection. Minimum-threshold membership changes at
score + 1; maximum-threshold membership changes at the observed score.
Unestimable results retain their reasons and break the odds-ratio line. Boundary
triangles denote zero or infinite odds ratios. Odds ratio above one means
relative enrichment of benign/likely-benign over pathogenic/likely-pathogenic
among retained calls. BH q-values cover thresholds and strategies within each
filter, variant type, context, consequence, and adjustment mode; confidence
intervals remain pointwise.

Retention includes failed gnomAD lookups; the overlap denominator excludes
them. Missing lookup and unestimable odds ratio are not negative results.

## Code Owners And Checks

- Input resolution, evidence verification, workspace identity, and artifact
  contracts: `analytics/io/`
- Evidence-derived reusable tables: `analytics/derivations/`,
  `analytics/io/alignment_aggregates.py`, and
  `analytics/io/annotation_support.py`
- Scientific calculations: `analytics/analyses/`
- Presentation only: `analytics/reporting/`
- Orchestration and CLI: `analytics/strategy_report.py`

Run the test matching the changed owner. Cross-boundary changes normally need
`tests/test_analysis_inputs.py`, the relevant cache/analysis test, and
`tests/test_strategy_report.py`. The full environment command is documented
in `run_validation.md`.
