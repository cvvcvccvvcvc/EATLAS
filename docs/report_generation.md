# Report Generation

Use this runbook to create an analytics report for one completed current run or
a compatible cohort. VEP annotation of pipeline candidates is already part of
the end-to-end Nextflow run; there is no separate candidate-VEP preparation
step.

## One Run

On the ITMO controller, update the checkout and load the cluster environment:

```bash
cd /nfs/home/$USER/gaph_v2
git pull --ff-only
source "$HOME/.gaph_v2_cluster_env.sh"

RUN="$GAPH_ROOT/results/<run-name>"

bash analytics/slurm/submit_strategy_report.sh \
  --run-dir "$RUN" \
  --report-name strategy_compare
```

The launcher requires an absolute run path and a report name containing only
letters, digits, `.`, `_`, or `-`. It rejects a run unless the current
`annotation/manifest.json` and partitioned
`annotation/variant_annotations/manifest.json` contracts are complete.

Pass report CLI arguments unchanged after `--`. For example, the optional
consequence-matched target-space null is enabled with:

```bash
bash analytics/slurm/submit_strategy_report.sh \
  --run-dir "$RUN" \
  --report-name strategy_compare_with_target_null \
  -- \
  --target-space-null \
  --target-space-null-sample-size 5000 \
  --target-space-null-resamples 500
```

Options not supplied by the user keep the `analytics.strategy_report` default.
Use the CLI as the source of truth:

```bash
micromamba run -p "$GAPH_ROOT/envs/analytics" \
  python -m analytics.strategy_report --help
```

## Pipeline And Report In One Command

For a new run that should immediately submit a report after successful
Nextflow completion, use the combined launcher inside `tmux`:

```bash
IDS="$GAPH_ROOT/inputs/panel.ids"
RUN="$GAPH_ROOT/results/run_name"
WORK="$GAPH_ROOT/work/run_name"

bash scripts/slurm/run_and_report.sh \
  --ids-file "$IDS" \
  --run-dir "$RUN" \
  --work-dir "$WORK" \
  --report-name strategy_compare \
  -- \
  --target-space-null \
  --target-space-null-sample-size 5000
```

Arguments before `--` belong to the operational launcher. Arguments after `--`
go unchanged to the report. Omit `--` when no report overrides are needed. The
launcher always uses `-profile slurm` and `-resume`, preserves the registry
alignment default unless `--alignment-strategies` is supplied, and never submits
the report after a failed pipeline run.

Optional report-job resource flags are `--slurm-cpus`, `--slurm-memory`,
`--slurm-time`, and `--slurm-partition`. They use the same defaults as the
report-only launcher when omitted.

## Multi-run Cohort

A cohort manifest contains one or more completed runs. Paths may be absolute or
relative to the manifest:

```json
{
  "schema_version": 1,
  "runs": [
    {"label": "panel-a", "run_dir": "/archive/gaph/run-a"},
    {"label": "panel-b", "run_dir": "/archive/gaph/run-b"}
  ]
}
```

```bash
COHORT_MANIFEST="$GAPH_ROOT/cohorts/panel.json"
COHORT_ROOT="$GAPH_ROOT/cohort_reports"

bash analytics/slurm/submit_strategy_report.sh \
  --cohort-manifest "$COHORT_MANIFEST" \
  --cohort-root "$COHORT_ROOT" \
  --report-name strategy_compare
```

Members must use the same current pipeline contracts, target/reference inputs,
strategy configuration, ClinVar and gnomAD contracts, VEP backend/release and
variant-table columns. Accepted Gene IDs must be disjoint. Overlap and
incompatible evidence fail explicitly; rows are never silently deduplicated.

The cohort keeps large variant shards in the source runs and exposes them as one
virtual DuckDB input. Compact tables, target-sequence links, derived caches, and
the resolved cohort manifest live under the stable cohort output directory.
Taxonomy distinct counts and medians are recomputed from the union of member
Stage 1 evidence rather than added across reports.

## Report Data Boundary

The report reads pipeline evidence and builds fingerprinted derived caches under
`analytics/`:

- strategy summary and feature coverage from Stage 2 summaries/segments;
- taxonomy summary from selected orthologs and canonical Stage 1 taxonomy;
- variant-strategy support and taxonomic depth/support from exact Stage 2
  evidence plus the Stage 3 event-to-variant map;
- scientific analysis tables used by the HTML sections.

These caches do not replace pipeline evidence. Missing source partitions,
mismatched schemas, or changed files are errors. There are no old aggregate or
separate bulk-VEP fallback paths.

The ClinVar validation universe publishes `clinvar_disease_ids` as its final
column. It preserves each non-empty source `CLNDISDB` value: comma-separated
identifiers and pipe-separated disease concepts retain their ClinVar meaning,
while distinct source values merged into one normalized variant are separated
by semicolons.

Rows with a non-`ok` pipeline VEP status remain explicit and appear as
`Not annotated` where consequence plots need a complete denominator. The
target-space null can make additional VEP and gnomAD requests for generated
control alleles; this is why it remains opt-in.

## Monitoring And Completion

The report launcher prints the Slurm job ID and log paths:

```bash
JOB_ID=<printed-job-id>
squeue -j "$JOB_ID"
tail -n 40 "$RUN/reports/slurm/<report-name>.$JOB_ID.out"
tail -n 40 "$RUN/reports/slurm/<report-name>.$JOB_ID.err"
sacct -j "$JOB_ID" --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

A successful single-run report has:

```bash
test -s "$RUN/reports/<report-name>.html"
test -s "$RUN/analytics/performance/<report-name>.json"
```

Slurm must report `COMPLETED` and exit code `0:0`. The performance JSON is
written progressively and is the first artifact to inspect after a report
failure. Repeating an identical report reuses valid fingerprinted caches.

The report launcher defaults to 8 CPUs, 128 GB, six hours, and the `main`
partition. The large memory default reflects a measured 590-gene aggregation;
change it only after measuring a representative report.

## Common Failures

- `Rscript was not found`: submit through the launcher so the declared
  analytics environment is active.
- `Repository commit changed after submission`: the cluster checkout changed
  while the job was queued; resubmit from the intended commit.
- DuckDB disk or memory errors: inspect the performance profile and Slurm
  `MaxRSS`; keep spill and results on the shared scratch allocation.
