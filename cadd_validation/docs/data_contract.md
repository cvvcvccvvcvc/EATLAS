# Data Contract

## Input Variants

`build_features` consumes a TSV with one row per candidate variant. The variant
universe should come from an independent source such as ClinVar or a benchmark
VCF, not from `alignment_events.tsv.gz`.

Required columns:

- `ref`
- `alt`

Coordinate columns, in priority order:

1. `gene_id` plus `target_start0`
2. `gene_id` plus `genomic_start1` or `pos`, using the gene row from
   `target_features.tsv.gz`

Recommended columns:

- `variant_id`
- `gene_id`
- `genomic_accession`
- `genomic_start1`
- `label`

If `variant_id` is absent, it is derived from available coordinate and allele
columns.

Only rows with resolvable target-local coordinates are emitted. Unresolved rows
are counted in the JSON summary and skipped.

## GAPH Stage 2 Inputs

`build_features` reads these published GAPH outputs:

- `alignment_segments.tsv.gz`
- `alignment_events.tsv.gz`
- optional `ortholog_alignment_summary.tsv.gz`
- optional `taxonomy_presets.tsv.gz`
- optional `target_features.tsv.gz`

The feature builder accepts gzip or plain TSV files.

## Feature Output

The output has one row per input variant and alignment strategy.

Core columns:

- `variant_id`
- `gene_id`
- `genomic_accession`
- `genomic_start1`
- `ref`
- `alt`
- `target_start0`
- `target_end0`
- `event_type`
- `strategy`

All model-ready GAPH columns start with `gaph_`. Examples:

- `gaph_all_depth`
- `gaph_all_alt_count`
- `gaph_all_alt_fraction`
- `gaph_primates_ref_fraction`
- `gaph_all_entropy`
- `gaph_all_no_call_count`
- `gaph_feature_context_mean_depth`

Fractions use callable depth as the denominator. If depth is zero, fractions are
`0.0` and entropy is `0.0`.

## Baseline Join Input

`join_baseline` joins feature rows to external annotations by `variant_id` by
default. Use `--left-key` and `--right-key` for another key.

Baseline tables should include labels and CADD/conservation columns such as:

- `label`
- `CADD_RAW`
- `CADD_PHRED`
- `phyloP`
- `phastCons`
- `GERP`

## Evaluation Dataset

`evaluate_ablation` expects:

- a binary label column
- zero or more baseline feature columns
- one or more `gaph_*` feature columns

Accepted positive labels include `1`, `pathogenic`, `likely_pathogenic`,
`p/lp`, and `lp/p`. Accepted negative labels include `0`, `benign`,
`likely_benign`, `b/lb`, and `lb/b`.

It writes:

- `metrics.tsv`
- `predictions.tsv`
- `split_assignments.tsv`
- `summary.json`
