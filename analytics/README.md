# GAPH analytics

This package contains analysis and reporting entrypoints for completed GAPH
runs.

Create the reproducible analytics environment once with:

```bash
micromamba create -f envs/analytics.yml
```

Primary report:

```bash
RUN="results/run_all_strategies_20260703_135905"

micromamba run -n gaph-v2-analytics python -m analytics.strategy_report \
  --run-dir "$RUN"
```

The report computes ClinVar enrichment and conservation-stratified validation
by default. Conservation tracks and quantile bins can be adjusted with:

```bash
micromamba run -n gaph-v2-analytics python -m analytics.strategy_report \
  --run-dir "$RUN" \
  --conservation-tracks phyloP100way,phastCons100way,GERP_RS_92mammals \
  --conservation-bins 4
```

When an analysis needs durable intermediate tables, write them under
`<run-dir>/analytics/`. The source tree does not keep a default scratch/work
directory.

The report reads `annotation/variant_annotations.tsv.gz` in chunks. It uses a
temporary SQLite file under `<run-dir>/analytics/` to deduplicate
variant-strategy records without loading the full annotation table into memory;
the file is removed when aggregation finishes or fails.
The report requires the canonical `alignment/strategy_summary.tsv.gz`; it does
not reconstruct that aggregate from a raw per-ortholog table.

The strategy report writes its ClinVar validation universe under:

```text
<run-dir>/analytics/clinvar_universe.snv_indel.tsv.gz
<run-dir>/analytics/clinvar_universe.snv_indel.manifest.json
<run-dir>/analytics/clinvar_target_regions.bed
```

Validation statistics are computed separately for SNV and INDEL rows. The
conservation-stratified block also writes:

```text
<run-dir>/analytics/clinvar_universe.snv.conservation.tsv.gz
<run-dir>/analytics/clinvar_universe.snv.conservation.manifest.json
```

The conservation cache is SNV-only and is reused on later report runs when the
ClinVar universe and requested tracks are unchanged. Successful track columns
are retained when another remote track fails; later runs retry only failed or
partial tracks until the cache is complete.

Raw p-values remain visible in the report. ClinVar Fisher tests use
Benjamini-Hochberg correction separately within the SNV and INDEL families.
Conservation bin-level Fisher tests and score-by-strategy CMH tests are corrected
as two separate analysis families. Mantel-Haenszel confidence intervals use the
Robins-Breslow-Greenland variance implemented by `statsmodels.StratifiedTable`.
