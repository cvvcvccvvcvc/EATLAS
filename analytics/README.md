# GAPH Analytics

`analytics/` turns one or more completed evidence-first pipeline runs into
reproducible derived tables and an HTML report. It does not fetch or repair
missing pipeline evidence, and pipeline commands do not import analytics
modules.

For the supported cluster launch sequence, arguments, and monitoring commands,
use `docs/report_generation.md`.

## Package Boundary

```text
analytics/
  strategy_report.py     report CLI and orchestration
  vep/                   report consequence vocabularies and grouping rules
  derivations/           deterministic tables built from pipeline evidence
  io/                    input validation and fingerprinted artifact contracts
  analyses/              scientific calculations
  reporting/             HTML sections and composition
  slurm/                 cluster launchers and workers
```

Cross-boundary variant, ClinVar, gnomAD, and taxonomy semantics live in
`genomics/`. `reporting/` renders already-computed results; it does not own
scientific calculations or external data access.

## Input Contract

The report accepts one or more repeated `--run-dir` arguments. Runs must be
compatible and non-overlapping. Each run must provide:

- canonical Stage 1 target, selected-ortholog, feature, and taxonomy evidence;
- partitioned Stage 2 summaries, segments, compact events, and exact support;
- the Stage 3 partitioned ClinVar/gnomAD/VEP variant dataset and
  event-to-variant map.

Missing evidence, an incomplete shard set, or an obsolete schema is a contract
error. There is no fallback to old pipeline aggregates or a separate bulk-VEP
artifact.

## Derived Data

The pipeline does not publish report-specific coverage, site-depth, support, or
taxonomic counters. A completed source run is immutable: analytics never writes
below it. Derived files live under the one explicitly supplied analytics root:

```text
<analytics-root>/
  cache/<source-id>/      reusable per-source evidence derivations
  analyses/<analysis-id>/
    derived/              run-set scientific artifacts
    reports/              HTML reports
    performance/          timing, memory, I/O, and disk profiles
  slurm/                  report job logs and revision markers
```

Every reusable artifact has a manifest that identifies its inputs and schema.
A cache hit avoids recomputation; a changed input invalidates the artifact.
Large candidate data is streamed or queried with bounded DuckDB operations
instead of being copied into parallel permanent relations.

Taxonomic scope/rank choices and the corresponding counters are analytics
derivations from `fetch/taxonomy.tsv.gz` and exact Stage 2 evidence. They are not
pipeline source data.

## Local Environment

Create the declared environment on first setup:

```bash
micromamba create --yes -f envs/analytics.yml
```

Synchronize an existing named environment after the YAML changes:

```bash
micromamba install --yes -n gaph-v2-analytics -f envs/analytics.yml
```

After the end-to-end pipeline has completed:

```bash
micromamba run -n gaph-v2-analytics python -m analytics.strategy_report \
  --analytics-root /absolute/path/to/analytics \
  --run-dir /absolute/path/to/completed-run \
  --report-name strategy_compare
```

The consequence-matched target-space null remains opt-in because it can perform
additional VEP and gnomAD work:

```bash
micromamba run -n gaph-v2-analytics python -m analytics.strategy_report \
  --analytics-root /absolute/path/to/analytics \
  --run-dir /absolute/path/to/completed-run \
  --report-name strategy_compare_with_target_null \
  --target-space-null
```

Use `python -m analytics.strategy_report --help` for the current CLI contract.

The same environment runs the Python suite:

```bash
micromamba run -n gaph-v2-analytics python -m pytest -q
```
