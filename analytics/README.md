# GAPH Analytics

`analytics/` turns one or more completed evidence-first pipeline runs into
reproducible derived tables and an HTML report. It never imports into pipeline
commands, repairs source evidence, or writes below a completed run.

Read `docs/report_generation.md` for the supported cluster launch and
`docs/analytics_contract.md` for source compatibility, cache identity, and
scientific-result semantics.

## Package Boundary

```text
analytics/
  strategy_report.py     CLI and orchestration
  io/                    source validation, identity, and artifact contracts
  derivations/           deterministic evidence-derived relations
  analyses/              scientific calculations
  reporting/             HTML presentation of computed results
  vep/                   report consequence vocabulary and grouping
  slurm/                 cluster launcher and worker
```

Cross-boundary variant, ClinVar, gnomAD, taxonomy, and VEP-provider semantics
belong to `genomics/`. Presentation modules do not fetch source data or own
scientific calculations.

## Input And Output

The CLI accepts one or more repeated `--run-dir` arguments. Each run must
provide the current Stage 1–3 evidence and a valid root evidence inventory.
Runs must be compatible and have disjoint accepted Gene IDs.

All writable state lives below one explicit external `--analytics-root`:

```text
<analytics-root>/
  cache/                  verified per-source and calculation caches
  analyses/<analysis-id>/
    derived/
    reports/
    performance/
  slurm/
```

Large candidate shards remain in source runs and are scanned in place. Missing,
changed, incomplete, or obsolete evidence is a contract error; there is no
legacy aggregate or bulk-VEP fallback.

## Local Use

Create or synchronize the declared environment:

```bash
micromamba create --yes -f envs/analytics.yml
micromamba install --yes -n gaph-v2-analytics -f envs/analytics.yml
```

Build a report after the source pipeline run completes:

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

Use `python -m analytics.strategy_report --help` for the current CLI. The same
environment runs the suite:

```bash
micromamba run -n gaph-v2-analytics python -m pytest -q
```
