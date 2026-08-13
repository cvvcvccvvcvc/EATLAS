# Report Generation

Use this runbook when the user asks to generate an analytics report for one
completed GAPH run or a compatible cohort of completed runs. It is the
complete operational path; do not read the wider cluster or analytics
documentation unless this path fails with a concrete error.

## Required User Inputs

The launcher requires exactly one input selector:

- the absolute completed-run directory; or
- an absolute cohort-manifest path;
- a new report name containing only letters, digits, `.`, `_`, or `-`.

Pass every report option stated by the user unchanged. When the user does not
specify an option, retain the Python CLI default. Do not interpret words such as
"full" as permission to enable optional analyses; ask which optional analyses
the user wants.

## Standard Cluster Launch

Publish completed local changes through GitHub before updating the cluster. Do
not transfer tracked code with `rsync`, Git bundles, or ad hoc copies.

```bash
git push origin main

ssh -i ~/.ssh/itmo ilunegov@ctlab.itmo.ru
ssh sphinx

cd /nfs/home/$USER/gaph_v2
git pull --ff-only
source "$HOME/.gaph_v2_cluster_env.sh"

RUN="$GAPH_ROOT/results/slurm_panel590_all_20260720_221912"

bash analytics/slurm/submit_vep_annotation.sh \
  --run-dir "$RUN"

# After bulk VEP has finalized:
bash analytics/slurm/submit_strategy_report.sh \
  --run-dir "$RUN" \
  --report-name strategy_compare_new
```

Bulk VEP is a required, resumable precompute for every report. Its launcher
submits `prepare -> bounded annotation array -> finalize`, pins the repository
commit, and writes logs under
`<run>/analytics/vep_consequences/slurm/`. The defaults are 250,000 rows per
partition, 4 CPUs and 8 GB per annotation task, and at most 4 simultaneous
tasks. Override them only after measuring a representative partition:

```bash
bash analytics/slurm/submit_vep_annotation.sh \
  --run-dir "$RUN" \
  --max-parallel 8 \
  --slurm-memory 12G
```

Repeating the same command resumes at completed partition boundaries. The
finalizer runs after the array stops and succeeds only when every declared
partition is complete and matches one VEP configuration. Do not submit the
HTML report until both files exist and are non-empty:

```bash
test -s "$RUN/analytics/vep_consequences/manifest.json"
test -s "$RUN/analytics/vep_consequences/variant_annotations.vep.tsv.gz"
```

The default command keeps `--target-space-null` disabled because that is the
`analytics.strategy_report` default. Enable it only when requested:

```bash
bash analytics/slurm/submit_strategy_report.sh \
  --run-dir "$RUN" \
  --report-name strategy_compare_with_target_null \
  -- \
  --target-space-null \
  --target-space-null-sample-size 5000 \
  --target-space-null-resamples 500
```

Arguments after `--` are passed unchanged to `analytics.strategy_report`.

## Multi-run Cohort

A cohort manifest makes any positive number of completed runs one analytical
cohort. A one-member cohort is valid. Paths may be absolute or relative to the
manifest:

```json
{
  "schema_version": 1,
  "runs": [
    {"label": "panel-a", "run_dir": "/archive/gaph/run-a"},
    {"label": "panel-b", "run_dir": "/archive/gaph/run-b"}
  ]
}
```

Every member must be a successful, clean, current `--stage all` run with a
finalized bulk-VEP artifact. Before creating analytics artifacts, the report
requires the same pipeline revision, target assembly and GFF contract,
strategy set and parameters, event/support contracts, ClinVar contents,
gnomAD dataset/API contract, and VEP release/backend/columns. Accepted gene IDs
must be disjoint. An overlap is reported with the conflicting genes and run
labels; it is not silently deduplicated because current durable run summaries
cannot subtract one gene correctly.

The large VEP partitions remain in their source runs and are scanned as one
virtual DuckDB input. Compact tables and target-sequence links are assembled
under `<cohort-root>/<cohort-id>/inputs/`; every scientific statistic is then
recomputed over the union. Run-level taxonomy medians and distinct counts are
not pooled because they are not additive. The resolved manifest records every
member fingerprint and appears in report QC.

```bash
COHORT_MANIFEST="$GAPH_ROOT/cohorts/panel590.json"
COHORT_ROOT="$GAPH_ROOT/cohort_reports"

bash analytics/slurm/submit_strategy_report.sh \
  --cohort-manifest "$COHORT_MANIFEST" \
  --cohort-root "$COHORT_ROOT" \
  --report-name strategy_compare
```

The default root is `<cohort-manifest-dir>/cohorts`. Reports, cohort caches,
performance profiles, and temporary spill data stay under the stable cohort ID;
source runs are read-only.

GAPH Browser cumulative releases are dissemination catalogs, not scientific
report inputs. A cohort manifest may describe the same membership, but every
analysis source must resolve to the completed `gaph_v2` run artifacts above.

## Launcher Arguments

### Bulk-VEP launcher

| Argument | Default | Meaning |
|---|---:|---|
| `--run-dir PATH` | required | Absolute completed-run directory. |
| `--partition-size N` | `250000` | Deterministic VEP input rows per array task. A resumed artifact must use the same value. |
| `--max-parallel N` | `4` | Maximum simultaneous VEP array tasks. |
| `--slurm-cpus N` | `4` | CPUs reserved per VEP array task. |
| `--slurm-memory SIZE` | `8G` | Memory per VEP array task. |
| `--slurm-time D-HH:MM:SS` | `01:00:00` | Time limit per VEP array task. |
| `--slurm-partition NAME` | `main` | Slurm partition for every step. |

### Report launcher

| Argument | Default | Meaning |
|---|---:|---|
| `--run-dir PATH` | mutually exclusive | Absolute completed-run directory visible on compute nodes. |
| `--cohort-manifest PATH` | mutually exclusive | Absolute JSON manifest containing the cohort runs. |
| `--cohort-root PATH` | `<manifest-dir>/cohorts` | Root containing stable cohort-ID output directories. Valid only with `--cohort-manifest`. |
| `--report-name NAME` | required | HTML stem and Slurm log stem. The HTML is written under the selected analysis root. |
| `--slurm-cpus N` | `8` | CPUs reserved for the report job. |
| `--slurm-memory SIZE` | `128G` | Slurm memory reservation. The 590-gene aggregation exceeded the effective DuckDB budget of a 32 GB job and used about 65 GB RSS with the larger allocation. |
| `--slurm-time D-HH:MM:SS` | `06:00:00` | Slurm wall-time limit. |
| `--slurm-partition NAME` | `main` | Slurm partition. |
| `--help` | — | Show launcher help without submitting. |
| `--` | — | End launcher arguments; pass all remaining arguments to the report CLI. |

The launcher does not choose scientific settings. It establishes the analytics
environment, records the Git commit, assigns Slurm resources, and submits the
batch job.

## Report CLI Arguments

The input selector and `--report-name` are supplied by the launcher. Put every other
report argument after `--`.

| Argument | Default and allowed values | Meaning |
|---|---|---|
| `--clinvar-vcf PATH` | `assets/reference/clinvar/clinvar.vcf.gz` | Indexed ClinVar VCF used for validation. |
| `--out-html PATH` | unset | Explicit output path. When set, it takes precedence over the launcher-provided report name. |
| `--target-space-null` / `--no-target-space-null` | disabled | Enable or disable the consequence-matched target-space null analysis. It can invoke VEP and gnomAD and can take hours on a cold cache. |
| `--target-space-null-sample-size N` | `25000`, minimum `1` | Maximum deterministic focal-SNV sample per strategy. It is a cap, not necessarily the final cohort size. |
| `--target-space-null-resamples N` | `1000`, minimum `100` | Matched-set bootstrap resampling iterations. |
| `--target-space-null-seed N` | `20260721` | Deterministic target-space-null sampling seed. |
| `--gnomad-cache-dir PATH` | `$GAPH_GNOMAD_CACHE_DIR`, otherwise unset | Shared resumable gnomAD regional cache. |
| `--phylop-bigwig PATH` | `$GAPH_PHYLOP_BIGWIG`; otherwise `$GAPH_ROOT/reference/ucsc/hg38.phyloP100way.bw` when present | Local hg38 phyloP100way BigWig. A supplied path must exist. |
| `--vep-backend {rest,local}` | `$GAPH_VEP_BACKEND`, otherwise `rest` | VEP execution backend. `local` requires a release and cache directory. |
| `--vep-release RELEASE` | `$GAPH_VEP_RELEASE`, otherwise detected from a matching bulk-VEP artifact or REST | Pinned Ensembl VEP release. Required for local VEP and reports using bulk-VEP consequences. |
| `--vep-executable PATH` | `$GAPH_VEP_EXECUTABLE`, otherwise `vep` | Local VEP executable or container wrapper. Used only by the local backend. |
| `--vep-cache-dir PATH` | `$GAPH_VEP_CACHE_DIR`, otherwise unset | Indexed local VEP cache root. Required by the local backend. |
| `--vep-result-cache-dir PATH` | `$GAPH_VEP_RESULT_CACHE_DIR`; otherwise `$GAPH_ROOT/cache/vep_results` | Sparse cross-run cache of completed VEP results. |
| `--vep-result-cache-tile-size-bp N` | `$GAPH_VEP_RESULT_CACHE_TILE_SIZE_BP`, otherwise `1000000`; minimum `1` | Genomic tile size for the shared VEP result cache. |
| `--vep-forks N` | `$GAPH_VEP_FORKS`, otherwise `4`; minimum `1` | Local VEP worker processes for the ClinVar universe and optional target-space null. |
| `--firth-workers N` | `$GAPH_FIRTH_WORKERS`; otherwise allocated Slurm CPUs capped at `8`; minimum `1` | Parallel independent Firth models. |

The report requires a finalized bulk-VEP artifact under
`<run>/analytics/vep_consequences` that matches `<run>/annotation`. It verifies
the manifest before expensive analysis and fails with the VEP preparation and
finalization commands when the artifact is missing. Non-`ok` per-row VEP
statuses are valid in a finalized artifact and are reported in QC. The Slurm
launcher checks that the final manifest and joined VEP table exist before it
submits a report job. Consequence plots retain those rows as a light-grey
`Not annotated` group instead of excluding them from the plotted denominator.
The gnomAD Stratification consequence view uses these same VEP groups and only
completed gnomAD lookups (`found` or `not_found`).

## Monitoring And Completion

The bulk-VEP preparation log prints the annotation-array and finalizer job IDs:

```bash
squeue -u "$USER"
tail -n 30 "$RUN/analytics/vep_consequences/slurm/prepare.<job-id>.out"
tail -n 30 "$RUN/analytics/vep_consequences/slurm/annotate.<array-id>_<task-id>.err"
tail -n 30 "$RUN/analytics/vep_consequences/slurm/finalize.<job-id>.err"
```

The preparation job, every array task, and the finalizer must finish with
Slurm state `COMPLETED` and exit code `0:0`. A failed array causes finalization
to fail on the missing partition rather than publishing a partial artifact;
rerun the launcher to retry only incomplete partitions.

After bulk VEP is finalized, confirm that the report itself starts and passes
its initial preflight:

The report launcher prints the job ID and exact log paths:

```bash
squeue -j "$JOB_ID"
tail -n 40 "$RUN/reports/slurm/<report-name>.$JOB_ID.out"
tail -n 40 "$RUN/reports/slurm/<report-name>.$JOB_ID.err"
```

A successful job must satisfy all of the following:

```bash
sacct -j "$JOB_ID" --format=JobID,State,Elapsed,MaxRSS,ExitCode
test -s "$RUN/reports/<report-name>.html"
test -s "$RUN/analytics/performance/<report-name>.json"
```

The Slurm state must be `COMPLETED` with exit code `0:0`. The performance JSON
is written progressively and is the first place to inspect after a failure.
Repeating the same submission reuses valid run-local and shared caches.

## Known Launch Failures

- `Rscript was not found`: the analytics environment was not added to `PATH`.
  Use the launcher; do not invoke its Python interpreter directly without the
  environment `bin` directory on `PATH`.
- DuckDB `Disk quota exceeded` during Variant Summary on the 590-gene artifact:
  a 32 GB job gave DuckDB only about 16 GB and forced excessive spill. The
  launcher therefore defaults to 128 GB. Increase `--slurm-memory` only when a
  measured run still requires it.
- `Repository commit changed after submission`: the cluster checkout changed
  while the job was queued. Resubmit from the intended checked-out commit.
