# GAPH Analytics

`analytics/` turns a completed evidence-first pipeline run into reproducible
derived tables and an HTML report. It does not fetch or repair missing pipeline
evidence, and pipeline commands do not import analytics modules.

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

The report accepts one current completed run or a manifest of compatible,
non-overlapping runs. Each run must provide:

- canonical Stage 1 target, selected-ortholog, feature, and taxonomy evidence;
- partitioned Stage 2 summaries, segments, compact events, and exact support;
- the Stage 3 partitioned ClinVar/gnomAD/VEP variant dataset and
  event-to-variant map.

Missing evidence, an incomplete shard set, or an obsolete schema is a contract
error. There is no fallback to old pipeline aggregates or a separate bulk-VEP
artifact.

## Derived Data

The pipeline does not publish report-specific coverage, site-depth, support, or
taxonomic counters. Analytics derives them under the run's `analytics/`
directory:

```text
analytics/
  alignment_aggregates/  strategy summary and feature coverage
  taxonomy_summary/      selected-ortholog taxonomy summary
  annotation_support/    variant support and ortholog-evidence views
  performance/           progressive timing, memory, I/O, and disk profiles
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
  --run-dir /absolute/path/to/completed-run
```

The consequence-matched target-space null remains opt-in because it can perform
additional VEP and gnomAD work:

```bash
micromamba run -n gaph-v2-analytics python -m analytics.strategy_report \
  --run-dir /absolute/path/to/completed-run \
  --target-space-null
```

Use `python -m analytics.strategy_report --help` for the current CLI contract.

The same environment runs the Python suite:

```bash
micromamba run -n gaph-v2-analytics python -m pytest -q
```
