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

The strategy report writes its ClinVar validation universe under:

```text
<run-dir>/analytics/clinvar_universe.snv_indel.tsv.gz
<run-dir>/analytics/clinvar_universe.snv_indel.manifest.json
<run-dir>/analytics/clinvar_target_regions.bed
```

Validation statistics are computed separately for SNV and INDEL rows.
