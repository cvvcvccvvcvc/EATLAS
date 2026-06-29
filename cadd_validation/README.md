# CADD Validation

Isolated validation workspace for testing whether GAPH ortholog-alignment
features add predictive value beyond CADD and standard conservation features.

This directory is intentionally independent from the production Nextflow
pipeline. It consumes published GAPH outputs and external annotation tables, but
does not change pipeline code.

## Question

Do variant-level features derived from GAPH ortholog alignments improve
pathogenicity prediction after accounting for existing predictors and
conservation scores?

The primary test is an ablation:

1. baseline features only, for example CADD and conservation scores
2. GAPH ortholog features only
3. baseline plus GAPH features
4. baseline plus shuffled GAPH features as a negative control

The evidence target is a reproducible improvement in held-out AUROC/AUPRC,
especially AUPRC, under gene- or chromosome-aware splits.

## Layout

```text
cadd_validation/
  docs/
    validation_design.md
    data_contract.md
  src/cadd_validation/
    build_features.py
    join_baseline.py
    evaluate_ablation.py
  tests/
    fixtures/
```

## Minimal Workflow

Create a local environment from the repository root:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r cadd_validation/requirements.txt
```

Build an independent ClinVar variant universe for GAPH target loci:

```bash
PYTHONPATH=cadd_validation/src .venv/bin/python -m cadd_validation.build_variant_universe \
  --clinvar-vcf clinvar.vcf.gz \
  --target-features-tsv results/run_001/fetch/target_features.tsv.gz \
  --out-tsv /tmp/gaph_clinvar_variants.tsv
```

Build variant-level GAPH features from existing Stage 2 outputs:

```bash
PYTHONPATH=cadd_validation/src .venv/bin/python -m cadd_validation.build_features \
  --variants-tsv /tmp/gaph_clinvar_variants.tsv \
  --segments-tsv results/run_001/alignment/alignment_segments.tsv.gz \
  --events-tsv results/run_001/alignment/alignment_events.tsv.gz \
  --summaries-tsv results/run_001/alignment/ortholog_alignment_summary.tsv.gz \
  --taxonomy-presets-tsv results/run_001/alignment/taxonomy_presets.tsv.gz \
  --target-features-tsv results/run_001/fetch/target_features.tsv.gz \
  --feature-coverage-tsv results/run_001/alignment/feature_coverage.tsv.gz \
  --out-tsv /tmp/gaph_variant_features.tsv
```

Join CADD/conservation annotations and labels:

```bash
PYTHONPATH=cadd_validation/src .venv/bin/python -m cadd_validation.join_baseline \
  --gaph-features-tsv /tmp/gaph_variant_features.tsv \
  --baseline-tsv baseline_annotations.tsv \
  --out-tsv /tmp/gaph_cadd_validation_dataset.tsv
```

Run the ablation:

```bash
PYTHONPATH=cadd_validation/src .venv/bin/python -m cadd_validation.evaluate_ablation \
  --dataset-tsv /tmp/gaph_cadd_validation_dataset.tsv \
  --label-column label \
  --baseline-features CADD_RAW,CADD_PHRED,phyloP,phastCons,GERP \
  --group-column gene_id \
  --outdir /tmp/gaph_cadd_ablation
```

Use `--strategies minimap2_asm20` in `build_features` when validating one
alignment strategy at a time.

## Scope

The first implementation targets SNVs and exact indel event matches. SNVs are
the recommended first validation target because CADD/conservation baselines and
ClinVar labels are easiest to normalize there. Indel interpretation should be
treated as exploratory until left-normalization and representation checks are
added upstream.
