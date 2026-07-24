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

The report computes ClinVar enrichment, fixed-band and continuous
phyloP100way-adjusted validation, and two sampled SNV background comparators.
The conservation analyses support strategy, variant-type, and ClinVar MC
consequence selectors. The primary view is SNV / Missense.

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
conservation-adjusted blocks also write:

```text
<run-dir>/analytics/clinvar_universe.snv_indel.conservation.tsv.gz
<run-dir>/analytics/clinvar_universe.snv_indel.conservation.manifest.json
```

The conservation cache is allele-level and is reused when the ClinVar universe
and phyloP100way source are unchanged. SNVs use the substituted base, deletions
use the mean across deleted reference bases without the VCF padding base, and
insertions use the mean of their two flanking bases. All required bases must be
scored.

Background comparators are shown in separate report tabs:

- `Matched Callable Background` compares GAPH SNVs with unobserved SNVs matched by
  gene, target context, REF, and callable-species depth bin;
- `Same-Position ALT: Raw` compares the exact GAPH ALT with other unobserved
  SNV ALTs at the same position. This view is descriptive because it does not
  yet adjust for transition/transversion class or context-specific mutability.

Their bounded, deterministic samples and annotations are cached under:

```text
<run-dir>/analytics/negative_control/
```

The default limit is 25,000 focal SNVs per strategy with 1,000 matched-background
resamples. For a faster exploratory run, use
`--negative-control-sample-size` and `--negative-control-permutations`; changing
only the number of resamples reuses the prepared control tables.

The sample size is an engineering cap, not a fixed scientific cohort size. A
run uses every eligible focal SNV when fewer than the cap are available. The
report shows full focal-weighted phyloP ECDFs and descriptive background
resampling intervals; it does not report an inferential p-value for these
comparators.

Raw p-values remain visible for the formal validation analyses. ClinVar Fisher
tests use Benjamini-Hochberg correction separately within the SNV and INDEL families.
Within conservation-adjusted validation, fixed-band Fisher tests, pooled CMH
tests, and continuous-model profile-likelihood-ratio tests are corrected as
three separate families across displayed selector combinations.
Mantel-Haenszel confidence intervals use the Robins-Breslow-Greenland variance
implemented by `statsmodels.StratifiedTable`. Continuous models use Firth
logistic regression (`logistf`) with a three-degree-of-freedom natural spline
for phyloP100way and profile penalized-likelihood confidence intervals.
