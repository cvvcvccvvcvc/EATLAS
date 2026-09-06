# Report Generation

Use this runbook to build an analytics report from one or more compatible,
completed GAPH runs. Candidate VEP annotation is already part of the pipeline;
there is no separate candidate-annotation step.

For source compatibility, cache identity, and scientific interpretation, read
`analytics_contract.md`.

## Required Inputs

Identify:

- one or more completed, immutable source-run directories;
- one external analytics root on shared storage;
- a report name;
- the intended full Git commit captured from the authoritative local checkout;
- any explicitly requested scientific or resource options.

The analytics root must be outside every source run. Source runs must be
complete, compatible, and non-overlapping as defined by
`analytics_contract.md`.

The report still needs an indexed ClinVar VCF for validation. `--clinvar-vcf`
defaults to `CLINVAR_VCF`, then
`assets/reference/clinvar/clinvar.vcf.gz`. Its content identity must match
every source run; the current file path may differ from the pipeline path.

Production reports use the cluster's pinned local VEP and local phyloP
configuration from `itmo_cluster.md`. Startup validates the VEP executable,
release, cache, and source evidence before creating analysis output.

## Pass The Revision Gate

Before every report submission, complete the local-to-cluster revision gate in
`pipeline_launch.md`. This applies even when the source pipeline run already
finished.

On the cluster, the clean checkout's `HEAD` and freshly fetched
`origin/main` must both equal `INTENDED_COMMIT`. If the launcher does not
accept the arguments documented below, stop and synchronize the checkout.
Never adapt a current report request to a historical cohort interface.

## Submit One Run

On the ITMO controller:

```bash
cd "$GAPH_CODE"
source "$HOME/.gaph_v2_cluster_env.sh"

RUN="$GAPH_ROOT/results/<run-name>"
ANALYTICS_ROOT="$GAPH_ROOT/analytics"

bash analytics/slurm/submit_strategy_report.sh \
  --analytics-root "$ANALYTICS_ROOT" \
  --run-dir "$RUN" \
  --report-name strategy_compare \
  --expected-commit "$INTENDED_COMMIT"
```

The launcher requires absolute paths, a report name containing only letters,
digits, dot, underscore, or hyphen, and a full 40-character commit. It verifies
the submitted revision and clean tree. The Slurm worker repeats the revision
check when the job starts, so a queued job cannot silently run changed code.

Use the launcher help for current resource flags and defaults:

```bash
bash analytics/slurm/submit_strategy_report.sh --help
```

## Submit Several Runs

Repeat `--run-dir`; there is no cohort manifest or combined source run:

```bash
bash analytics/slurm/submit_strategy_report.sh \
  --analytics-root "$ANALYTICS_ROOT" \
  --run-dir "$GAPH_ROOT/results/panel_a" \
  --run-dir "$GAPH_ROOT/results/panel_b" \
  --report-name combined_panel \
  --expected-commit "$INTENDED_COMMIT"
```

Compatibility and accepted-Gene-ID overlap are checked before scientific work.
Analytics never silently deduplicates source runs.

## Pass Report Options

Arguments after `--` are passed unchanged to
`python -m analytics.strategy_report`:

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

The target-space null is opt-in because it may perform additional VEP and
gnomAD work. The report CLI is the source of truth for scientific and
calculation-worker options:

```bash
micromamba run -p "$GAPH_ROOT/envs/analytics" \
  python -m analytics.strategy_report --help
```

`--annotation-support-workers` and `--firth-workers` are operational
parallelism controls; they do not change scientific meaning. Keep both within
the CPUs requested from the launcher. Each annotation-support worker may use up
to `GAPH_ANALYTICS_DUCKDB_MEMORY_LIMIT` (default 2 GB per worker).
The main-process DuckDB calculations use `GAPH_DUCKDB_MEMORY_LIMIT` when set;
otherwise they reserve half the Slurm memory allocation, or half DuckDB's
default limit outside Slurm. These are DuckDB budgets, not total process limits:
leave room for Python/R objects, external tools, and concurrent support workers.

## Monitor And Verify

The launcher prints the Slurm job ID, analytics workspace, and log paths:

```bash
JOB_ID=<printed-job-id>
squeue -j "$JOB_ID"
tail -n 40 "$ANALYTICS_ROOT/slurm/<report-name>.$JOB_ID.out"
tail -n 40 "$ANALYTICS_ROOT/slurm/<report-name>.$JOB_ID.err"
sacct -j "$JOB_ID" --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

As soon as the worker starts, confirm that the `Git commit:` line in stdout
equals `INTENDED_COMMIT`.

A successful analysis contains:

```text
<analytics-root>/analyses/<analysis-id>/manifest.json
<analytics-root>/analyses/<analysis-id>/reports/<report-name>.html
<analytics-root>/analyses/<analysis-id>/performance/<report-name>.json
```

Slurm must report `COMPLETED` with exit code `0:0`. The performance profile
is written progressively and is the first artifact to inspect during a long
derivation or after failure. Repeating the identical command reuses valid
source caches, run-set artifacts, and complete partition checkpoints.

The completed source runs must have the same file membership and metadata
before and after report generation.

## Common Failures

- `Repository commit changed after submission`: the checkout changed while
  the job waited; restore the intended clean revision and resubmit.
- Launcher rejects `--analytics-root` or repeated `--run-dir`: the cluster
  checkout is obsolete; fetch and fast-forward it.
- Source compatibility or overlap error: use runs with one current scientific
  contract and disjoint accepted Gene IDs.
- Missing or changed evidence: restore the exact completed run; do not patch it
  or substitute an old aggregate.
- Local VEP preflight error: correct the release-pinned executable/cache
  configuration before rerunning.
- DuckDB disk or memory error: inspect the performance profile and Slurm
  `MaxRSS`; keep the analytics root on shared scratch and adjust measured
  resources rather than changing scientific options.
- `Rscript was not found`: run through the launcher so the declared analytics
  environment is active.
