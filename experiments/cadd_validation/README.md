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
experiments/cadd_validation/
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
.venv/bin/python -m pip install -r experiments/cadd_validation/requirements.txt
mkdir -p experiments/cadd_validation/outputs/pilot
```

Build an independent ClinVar variant universe for GAPH target loci:

```bash
PYTHONPATH=experiments/cadd_validation/src .venv/bin/python -m cadd_validation.build_variant_universe \
  --clinvar-vcf clinvar.vcf.gz \
  --target-features-tsv results/run_001/fetch/target_features.tsv.gz \
  --out-tsv experiments/cadd_validation/outputs/pilot/clinvar_variants.tsv
```

For a pilot where the target genes have not been chosen yet, first select
ClinVar-rich genes with both pathogenic and benign SNVs:

```bash
PYTHONPATH=experiments/cadd_validation/src .venv/bin/python -m cadd_validation.select_clinvar_genes \
  --clinvar-vcf clinvar.vcf.gz \
  --out-genes experiments/cadd_validation/outputs/pilot/gene_ids.txt \
  --out-gene-summary experiments/cadd_validation/outputs/pilot/selected_genes.tsv \
  --out-variants experiments/cadd_validation/outputs/pilot/selected_clinvar_variants.tsv \
  --min-per-label 3 \
  --max-genes 10
```

Build variant-level GAPH features from existing Stage 2 outputs:

```bash
PYTHONPATH=experiments/cadd_validation/src .venv/bin/python -m cadd_validation.sample_variants \
  --variants-tsv experiments/cadd_validation/outputs/pilot/clinvar_variants.tsv \
  --out-tsv experiments/cadd_validation/outputs/pilot/clinvar_variants.sampled.tsv \
  --group-columns gene_id,label \
  --max-per-group 20

PYTHONPATH=experiments/cadd_validation/src .venv/bin/python -m cadd_validation.build_features \
  --variants-tsv experiments/cadd_validation/outputs/pilot/clinvar_variants.sampled.tsv \
  --segments-tsv results/run_001/alignment/alignment_segments.tsv.gz \
  --events-tsv results/run_001/alignment/alignment_events.tsv.gz \
  --summaries-tsv results/run_001/alignment/ortholog_alignment_summary.tsv.gz \
  --taxonomy-tsv results/run_001/fetch/taxonomy.tsv.gz \
  --target-features-tsv results/run_001/fetch/target_features.tsv.gz \
  --feature-coverage-tsv results/run_001/alignment/feature_coverage.tsv.gz \
  --out-tsv experiments/cadd_validation/outputs/pilot/variant_features.tsv
```

Join CADD/conservation annotations and labels:

```bash
PYTHONPATH=experiments/cadd_validation/src .venv/bin/python -m cadd_validation.fetch_cadd_scores \
  --variants-tsv experiments/cadd_validation/outputs/pilot/variant_features.tsv \
  --out-tsv experiments/cadd_validation/outputs/pilot/cadd_scores.tsv

PYTHONPATH=experiments/cadd_validation/src .venv/bin/python -m cadd_validation.join_baseline \
  --gaph-features-tsv experiments/cadd_validation/outputs/pilot/variant_features.tsv \
  --baseline-tsv experiments/cadd_validation/outputs/pilot/cadd_scores.tsv \
  --out-tsv experiments/cadd_validation/outputs/pilot/validation_dataset.tsv
```

Run the ablation:

```bash
PYTHONPATH=experiments/cadd_validation/src .venv/bin/python -m cadd_validation.evaluate_ablation \
  --dataset-tsv experiments/cadd_validation/outputs/pilot/validation_dataset.tsv \
  --label-column label \
  --baseline-features CADD_RAW,CADD_PHRED,phyloP,phastCons,GERP \
  --group-column gene_id \
  --outdir experiments/cadd_validation/outputs/pilot/ablation
```

Use `--strategies minimap2_asm20` in `build_features` when validating one
alignment strategy at a time.

## Scope

The first implementation targets SNVs and exact indel event matches. SNVs are
the recommended first validation target because CADD/conservation baselines and
ClinVar labels are easiest to normalize there. Indel interpretation should be
treated as exploratory until left-normalization and representation checks are
added upstream.
