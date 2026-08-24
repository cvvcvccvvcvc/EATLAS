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

## Submit A Report

On the ITMO controller:

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
does not submit a report after a failed pipeline run.

Optional report-job resources are `--slurm-cpus`, `--slurm-memory`,
`--slurm-time`, and `--slurm-partition`.

## Derived Data

The report derives fingerprinted caches from durable evidence:

- strategy summary and feature coverage from Stage 2 summaries and segments;
- taxonomy summary from selected orthologs and canonical Stage 1 taxonomy;
- variant-strategy support and taxonomic depth/support from exact Stage 2
  evidence plus Stage 3 event-to-variant lineage;
- scientific tables used by report sections.

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
written progressively and is the first artifact to inspect after failure.
Repeating an identical command reuses valid caches.

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
