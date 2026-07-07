# GAPH analytics

This package contains analysis and reporting entrypoints for completed GAPH
runs.

Primary report:

```bash
RUN="results/run_all_strategies_20260703_135905"

.venv/bin/python -m analytics.strategy_report \
  --run-dir "$RUN"
```

When an analysis needs durable intermediate tables, write them under
`<run-dir>/analytics/`. The source tree does not keep a default scratch/work
directory.

The strategy report currently writes an SNV-only ClinVar validation cache under:

```text
<run-dir>/analytics/clinvar_universe.snv.tsv.gz
<run-dir>/analytics/clinvar_universe.snv.manifest.json
```
