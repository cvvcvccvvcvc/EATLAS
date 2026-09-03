# Report Generation

Use this runbook to build a report from one or more compatible completed GAPH
runs. Candidate VEP annotation is already part of the Nextflow pipeline; no
separate candidate-annotation command is required.

## Data Boundary

A completed pipeline run is immutable source data. Analytics reads it but must
never create, modify, or delete files below its run directory. If pipeline
evidence must change, create a new pipeline run instead of patching the old one.

Every report command therefore requires one external `--analytics-root`. There
is no environment-variable default and no output-path override. The workspace
contains reusable source caches and analysis-specific outputs:

```text
<analytics-root>/
  .strategy_report.lock
  cache/<source-id>/
    alignment_aggregates/
    annotation_support/
    taxonomy_summary/
  analyses/<analysis-id>/
    manifest.json
    derived/
    reports/
    performance/
  slurm/
```

`source-id` depends on completed pipeline provenance. `analysis-id` depends on
the unordered set of source IDs and scientific report options. Repeating the
same analysis reuses the same workspace. Report names only select an HTML file
inside `reports/`; they do not create a second scientific cache.

Large variant shards remain in source runs and are queried as one virtual
dataset. Small gene, feature, failure, coverage, and support tables are read
directly from each source. Analytics materializes only reusable per-source
derivations and genuine run-set results; it does not build a synthetic combined
run, copy source TSV files, or create target-FASTA symlink trees.

One report process owns an analytics root at a time. A concurrent report for
the same root exits immediately instead of writing shared caches concurrently;
reports using different roots remain independent. The lock file itself is only
a stable rendezvous point. Process exit, including a crash, releases ownership,
so the file may safely remain in place.

Annotation-support preparation is partition-local. Each worker reads only the
variant-annotation shards owned by its evidence partition, uses one DuckDB
thread, and atomically checkpoints the two small derived partition tables.
After interruption, a repeated identical command reuses completed checkpoints
and continues with the missing partitions; incomplete partition directories are
never accepted. Checkpoints are removed after the final per-source cache has
been assembled, so they do not become a second durable dataset.

## Submit A Report

First pass the local-to-cluster revision gate in `docs/pipeline_launch.md`; do
not use a merely clean but unfetched cluster checkout. On the ITMO controller:

```bash
cd /nfs/home/$USER/gaph_v2
git pull --ff-only
source "$HOME/.gaph_v2_cluster_env.sh"

RUN="$GAPH_ROOT/results/<run-name>"
ANALYTICS_ROOT="$GAPH_ROOT/analytics"
INTENDED_COMMIT=<paste-the-full-commit-captured-locally>

bash analytics/slurm/submit_strategy_report.sh \
  --analytics-root "$ANALYTICS_ROOT" \
  --run-dir "$RUN" \
  --report-name strategy_compare \
  --expected-commit "$INTENDED_COMMIT"
```

Tracked and staged code changes must be committed before submission. The Slurm
worker verifies both the submitted commit and the working tree again when the
job starts; changing the checkout while a report is queued causes an explicit
failure instead of running different code.

The current launcher requires external `--analytics-root` and accepts repeated
`--run-dir` arguments. If the cluster launcher rejects those arguments or asks
for a historical `--cohort-manifest`/`--cohort-root` workflow, stop: the cluster
checkout is obsolete. Fetch and fast-forward to the intended commit instead of
adapting the report request to the old interface. Once the worker starts,
confirm that the `Git commit:` line in Slurm stdout matches the intended
revision before allowing a long report to continue.

For several runs, repeat `--run-dir`; there is no separate cohort mode or
cohort manifest:

```bash
bash analytics/slurm/submit_strategy_report.sh \
  --analytics-root "$ANALYTICS_ROOT" \
  --run-dir "$GAPH_ROOT/results/panel_a" \
  --run-dir "$GAPH_ROOT/results/panel_b" \
  --report-name combined_panel \
  --expected-commit "$INTENDED_COMMIT"
```

Source runs must use the same current pipeline contracts, target/reference
inputs, strategy configuration, ClinVar and gnomAD contracts, VEP
backend/release, and variant columns. Accepted Gene IDs must be disjoint.
Overlap and incompatibility fail explicitly; analytics never silently
deduplicates them.

ClinVar and gnomAD record fields are allele-level evidence. Analytics reconciles
their repeated values by canonical `variant_key`: a successful non-empty value
is shared across its gene contexts, while conflicting non-empty values or a
non-numeric gnomAD allele frequency fail explicitly. VEP consequences, target
context, and lookup outcome remain gene-context evidence.

`--clinvar-vcf` defaults to `CLINVAR_VCF`, then to
`assets/reference/clinvar/clinvar.vcf.gz`. Its indexed contents must match the
portable ClinVar identity recorded by every source run; the current file does
not need to be at the path used when the pipeline ran.

With the local VEP backend, report startup probes the configured executable and
release-specific cache before creating the analysis workspace or derived data.
An unavailable executable, cache, or release fails the report immediately.

Pass scientific report options unchanged after `--`. For example:

```bash
bash analytics/slurm/submit_strategy_report.sh \
  --analytics-root "$ANALYTICS_ROOT" \
  --run-dir "$RUN" \
  --report-name strategy_compare_with_target_null \
  --expected-commit "$INTENDED_COMMIT" \
  -- \
  --target-space-null \
  --target-space-null-sample-size 5000 \
  --target-space-null-resamples 500
```

`--annotation-support-workers` is operational and does not change the analysis
identity or scientific result. It defaults to `1` outside Slurm and to the
smaller of `4` and `SLURM_CPUS_PER_TASK` in a Slurm job. Override it after `--`,
or set `GAPH_ANNOTATION_SUPPORT_WORKERS`. Keep it within the CPUs requested from
the launcher; each worker may use up to the configured
`GAPH_ANALYTICS_DUCKDB_MEMORY_LIMIT` (default `2GB`).

The report launcher requires absolute paths. `--report-name` may contain only
letters, digits, `.`, `_`, and `-`. The CLI is the source of truth for
scientific options:

```bash
micromamba run -p "$GAPH_ROOT/envs/analytics" \
  python -m analytics.strategy_report --help
```

## Derived Data

The report derives fingerprinted caches from durable evidence:

- strategy summary and feature coverage from Stage 2 summaries and segments;
- taxonomy summary from selected orthologs and canonical Stage 1 taxonomy;
- variant-strategy support and taxonomic depth/support from exact Stage 2
  evidence plus Stage 3 event-to-variant lineage;
- `derived/pathogenic_clinvar_hits.tsv.gz`, the complete unique P/LP allele
  table used by the second report tab. It includes ClinVar review strength and
  conditions, VEP effects, gnomAD AF, strategy membership, and
  exact ortholog-support summaries;
- scientific tables used by report sections.

The `Pathogenic ClinVar Hits` tab owns the P/LP-only review-star and molecular-
effect plots, the condition and exact-support views, and the paginated
detail table. The table sorts the complete in-browser dataset by configurable
primary and secondary columns; it is not a top-N extract. Candidate Profile and
QC do not repeat these P/LP-only views.

### Pathogenic support and condition backgrounds

P/LP subtypes remain in the detail table but are combined in the plots.
Exact-support distributions count one SNV per strategy, selecting the complete
target-gene row with maximum ALT-supporting ortholog count (gene ID breaks
ties). Violin densities use log10 counts with original-count axis labels and
hover values; sparse distributions use points and a box. phyloP availability
does not restrict this tab. Conservation-adjusted analyses still use phyloP.

Conditions compare unique P/LP alleles against either ClinVar in the strategy's
eligible target genes or the whole pinned GRCh38 VCF. Both arms use the selected
variant type. The denominator includes alleles without named conditions;
multiple conditions per allele mean the displayed fractions need not sum to
one. Named-condition coverage and counts are available on hover. Disease IDs
(MedGen, then MONDO, then OMIM) identify conditions when supplied; names are
used otherwise, without inferred disease groups. The plot selects the top 15
conditions across both arms and supports name search.

The whole-VCF background is streamed locally and cached once per source identity
under `cache/clinvar_conditions/`; it does not trigger VEP or network requests.
It counts distinct VCF alleles, merging repeated records and excluding mixed
B/LB–P/LP classifications. This describes the variants and associated conditions
represented in the VCF, not all ClinVar record types or condition-specific
clinical assertions.

### Basic Filtering

Filter and variant type apply to the whole tab. Retention and gnomAD share a
strategy selector; ClinVar has independent strategy, adjustment, target-context,
and consequence selectors. `Compare all strategies` draws separate curves with
each strategy's own candidate denominator. `Union (any strategy)` counts each
normalized allele once and requires it to pass in at least one calling
strategy/context. Union OR is calculated from a new allele-level contingency
table over the union of eligible target genes, never by averaging strategy ORs.

Exact ALT support and strategy support use minimum thresholds. Supporting
families count distinct known family IDs; an unresolved family is not a new
family. Site-aligned ortholog counts support both maximum and minimum thresholds.
Family and site-depth filters are SNV-only. Within an allele/strategy, maximum
support is used for minimum thresholds and minimum site depth for maximum
thresholds. Non-calls have missing scores and cannot pass a maximum threshold;
zero known families is a valid score for a called SNV. Invalid or missing
required SNV support fails validation rather than removing rows.

ClinVar thresholds are derived after selecting the cohort (and finite phyloP
scores for adjustment). Minimum-threshold membership changes at score + 1;
maximum-threshold membership changes at the observed score. There is no fixed
20-point limit. Cumulative counts avoid repeated scans at each threshold.
Unestimable results retain their reasons and break the OR line. Triangles at
the plot boundary denote zero or infinite OR, not finite point estimates.
OR > 1 indicates relative enrichment of B/LB over P/LP among retained calls.
BH q-values cover thresholds and strategies within each filter, variant type,
context, consequence, and adjustment mode; intervals remain pointwise.

Retention includes failed gnomAD lookups; the overlap denominator excludes them.
Neither a missing lookup nor an unestimable OR is treated as a negative result.

These files never replace pipeline evidence. Missing source partitions,
mismatched schemas, changed files, and invalid empty outputs are errors. There
are no legacy aggregate, cohort-workspace, or separate bulk-VEP fallbacks.

Rows with non-`ok` pipeline VEP status remain explicit and appear as
`Not annotated` where a complete denominator is required. The optional
target-space null may make additional VEP and gnomAD requests for generated
control alleles. It always uses the VEP release pinned by the pipeline
artifacts; both that release and the selected VEP backend are part of its cache
identity.

## Monitoring And Completion

The launcher writes Slurm logs below the external analytics root:

```bash
JOB_ID=<printed-job-id>
squeue -j "$JOB_ID"
tail -n 40 "$ANALYTICS_ROOT/slurm/<report-name>.$JOB_ID.out"
tail -n 40 "$ANALYTICS_ROOT/slurm/<report-name>.$JOB_ID.err"
sacct -j "$JOB_ID" --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

The launcher prints the analytics workspace. A successful analysis contains:

```text
<analytics-root>/analyses/<analysis-id>/reports/<report-name>.html
<analytics-root>/analyses/<analysis-id>/performance/<report-name>.json
```

Slurm must report `COMPLETED` and exit code `0:0`. The performance JSON is
created before cache preparation, then written progressively. During annotation
support it records total, completed, built, and reused partitions, worker count,
elapsed time, and the last completed partition. It is the first artifact to
inspect during a long preparation or after failure. Repeating an identical
command reuses valid caches and completed partition checkpoints.

The launcher defaults to 8 CPUs, 128 GB, six hours, and the `main` partition.
Change these only after measuring a representative report.

## Common Failures

- `Rscript was not found`: submit through the launcher so the declared
  analytics environment is active.
- `Repository commit changed after submission`: the cluster checkout changed
  while the job was queued; resubmit from the intended commit.
- Source compatibility error: use runs produced with one scientific contract
  and disjoint accepted Gene IDs.
- DuckDB disk or memory error: inspect the performance profile and Slurm
  `MaxRSS`; keep the analytics root on shared scratch.
