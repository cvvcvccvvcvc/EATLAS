scripts/compare_strategies.py

```
RUN="results/run_all_strategies_20260703_135905"

.venv/bin/python scripts/compare_strategies.py \
  --events-tsv "$RUN/alignment/alignment_events.tsv.gz" \
  --feature-coverage-tsv "$RUN/alignment/feature_coverage.tsv.gz" \
  --out-html "$RUN/reports/strategy_compare.html"
```
