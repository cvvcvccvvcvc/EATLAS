Standalone utility scripts live here. The main run-level analytics report lives
under `analytics/`.

```
RUN="results/run_all_strategies_20260703_135905"

.venv/bin/python -m analytics.strategy_report \
  --run-dir "$RUN"
```
