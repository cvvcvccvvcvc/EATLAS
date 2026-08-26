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

bash analytics/slurm/submit_strategy_report.sh \
  --analytics-root "$ANALYTICS_ROOT" \
  --run-dir "$RUN" \
  --report-name strategy_compare
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
  --report-name combined_panel
```

Source runs must use the same current pipeline contracts, target/reference
inputs, strategy configuration, ClinVar and gnomAD contracts, VEP
backend/release, and variant columns. Accepted Gene IDs must be disjoint.
Overlap and incompatibility fail explicitly; analytics never silently
deduplicates them.

`--clinvar-vcf` defaults to `CLINVAR_VCF`, then to
`assets/reference/clinvar/clinvar.vcf.gz`. Its indexed contents must match the
portable ClinVar identity recorded by every source run; the current file does
not need to be at the path used when the pipeline ran.

Pass scientific report options unchanged after `--`. For example:

```bash
bash analytics/slurm/submit_strategy_report.sh \
  --analytics-root "$ANALYTICS_ROOT" \
  --run-dir "$RUN" \
  --report-name strategy_compare_with_target_null \
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

## Pipeline And Report In One Command

For a new run that should submit a report after successful Nextflow completion:

```bash
IDS="$GAPH_ROOT/inputs/panel.ids"
RUN="$GAPH_ROOT/results/run_name"
WORK="$GAPH_ROOT/work/run_name"
ANALYTICS_ROOT="$GAPH_ROOT/analytics"

bash scripts/slurm/run_and_report.sh \
  --ids-file "$IDS" \
  --run-dir "$RUN" \
  --work-dir "$WORK" \
  --analytics-root "$ANALYTICS_ROOT" \
  --report-name strategy_compare \
  -- \
  --target-space-null \
  --target-space-null-sample-size 5000
```

Arguments before `--` belong to the launcher; arguments after it go unchanged
to the report. The launcher uses `-profile slurm` and `-resume`, preserves the
registry alignment default unless `--alignment-strategies` is supplied, and
does not submit a report after a failed pipeline run. It validates the clean
checkout and analytics path before starting Nextflow.

Optional report-job resources are `--slurm-cpus`, `--slurm-memory`,
`--slurm-time`, and `--slurm-partition`.

## Derived Data

The report derives fingerprinted caches from durable evidence:

- strategy summary and feature coverage from Stage 2 summaries and segments;
- taxonomy summary from selected orthologs and canonical Stage 1 taxonomy;
- variant-strategy support and taxonomic depth/support from exact Stage 2
  evidence plus Stage 3 event-to-variant lineage;
- `derived/pathogenic_clinvar_hits.tsv.gz`, the complete unique P/LP allele
  table used by the second report tab. It includes ClinVar review strength and
  conditions, VEP effects, phyloP100way, gnomAD AF, strategy membership, and
  exact ortholog-support summaries;
- scientific tables used by report sections.

The `Pathogenic ClinVar Hits` tab owns the P/LP-only review-star and molecular-
effect plots, the condition and evolutionary-support views, and the paginated
detail table. The table sorts the complete in-browser dataset by configurable
primary and secondary columns; it is not a top-N extract. Candidate Profile and
QC do not repeat these P/LP-only views.

These files never replace pipeline evidence. Missing source partitions,
mismatched schemas, changed files, and invalid empty outputs are errors. There
are no legacy aggregate, cohort-workspace, or separate bulk-VEP fallbacks.

Rows with non-`ok` pipeline VEP status remain explicit and appear as
`Not annotated` where a complete denominator is required. The optional
target-space null may make additional VEP and gnomAD requests for generated
control alleles.

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
