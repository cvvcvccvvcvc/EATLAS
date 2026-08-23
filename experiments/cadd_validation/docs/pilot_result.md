# Pilot Result

Date: 2026-06-29

This is a first real-data sanity check, not the final CADD-additive validation.
It was produced with the pre-partition pipeline and is retained as a historical
result, not as an executable current runbook.

## Data

GAPH was run on five ClinVar-rich genes:

| Gene ID | Symbol | ClinVar P/LP SNV | ClinVar B/LB SNV |
| --- | --- | ---: | ---: |
| 1080 | CFTR | 743 | 1433 |
| 3043 | HBB | 185 | 763 |
| 3949 | LDLR | 1073 | 949 |
| 5053 | PAH | 666 | 399 |
| 7157 | TP53 | 403 | 1016 |

Pipeline outputs:

- fetch: `/tmp/gaph_cadd_real/gaph_pilot5_fetch`
- alignment: `/tmp/gaph_cadd_real/gaph_pilot5_align_asm20`
- strategy: `minimap2_asm20`
- selected orthologs: 3010
- alignment events: 530810
- alignment failures: 0

Validation sample:

- 200 ClinVar SNV rows
- balanced by `gene_id,label`: 20 benign and 20 pathogenic per gene
- split: `GroupKFold` by `gene_id`

## Runtime

Fetch trace:

| Stage | Runtime |
| --- | ---: |
| LDLR fetch | 47s |
| CFTR fetch | 4m10s |
| HBB fetch | 5s |
| TP53 fetch | 1m03s |
| PAH fetch | 1m30s |
| fetch merge / target features | 9m16s |

Alignment trace:

| Stage | Runtime |
| --- | ---: |
| taxonomy metadata | <1s |
| alignment task build | 1.6s |
| PAH minimap2 | 11s |
| CFTR minimap2 | 49s |
| TP53 minimap2 | 1.5s |
| HBB minimap2 | <1s |
| LDLR minimap2 | 5.3s |
| alignment merge | 6.5s |
| event annotation | 3m20s |

## Result Without CADD

CADD API access was blocked by the execution environment usage limit, and no
local CADD score table was found. Therefore this pilot only tests whether GAPH
features carry signal compared with a constant baseline.

Logistic regression, gene-held-out:

| Feature set | AUROC | AUPRC |
| --- | ---: | ---: |
| constant baseline | 0.500 | 0.500 |
| GAPH-only | 0.648 | 0.642 |
| baseline + shuffled GAPH | 0.357 | 0.411 |

Bootstrap 95% intervals for GAPH-only:

- AUROC: 0.589 to 0.765
- AUPRC: 0.588 to 0.717

Linear SVM showed the same direction:

| Feature set | AUROC | AUPRC |
| --- | ---: | ---: |
| constant baseline | 0.500 | 0.500 |
| GAPH-only | 0.634 | 0.619 |
| baseline + shuffled GAPH | 0.364 | 0.417 |

Held-out gene metrics for logistic GAPH-only:

| Gene ID | AUROC | AUPRC |
| --- | ---: | ---: |
| 1080 | 0.715 | 0.727 |
| 3043 | 0.675 | 0.606 |
| 3949 | 0.627 | 0.633 |
| 5053 | 0.639 | 0.616 |
| 7157 | 0.853 | 0.827 |

## Biological Read

The direction is plausible: benign variants had much more ortholog ALT support
and higher entropy than pathogenic variants.

Mean feature values in the 200-row sample:

| Label | Mean callable depth | Mean ALT fraction | Mean callable fraction | Mean entropy |
| --- | ---: | ---: | ---: | ---: |
| benign | 25.28 | 0.0457 | 0.1562 | 0.2414 |
| pathogenic | 31.94 | 0.0002 | 0.1640 | 0.0510 |

Rows with nonzero ortholog ALT support:

- benign: 31 / 100
- pathogenic: 1 / 100

## Interpretation

This supports a weak but real preliminary claim:

> GAPH ortholog-alignment features contain pathogenicity-relevant signal on this
> small gene-held-out ClinVar pilot.

It does not yet support the stronger claim:

> GAPH features add predictive value beyond CADD/conservation.

That claim still requires a CADD/conservation baseline table for the same
variants.

## Next Required Step

Provide or fetch CADD annotations for
`/tmp/gaph_cadd_real/gaph_features_pilot5_sample200.tsv`, then run:

```bash
mkdir -p experiments/cadd_validation/outputs/pilot_real

PYTHONPATH=experiments/cadd_validation/src .venv/bin/python -m cadd_validation.join_baseline \
  --gaph-features-tsv /tmp/gaph_cadd_real/gaph_features_pilot5_sample200.tsv \
  --baseline-tsv /path/to/cadd_scores.tsv \
  --out-tsv experiments/cadd_validation/outputs/pilot_real/dataset_with_cadd.tsv

PYTHONPATH=experiments/cadd_validation/src .venv/bin/python -m cadd_validation.evaluate_ablation \
  --dataset-tsv experiments/cadd_validation/outputs/pilot_real/dataset_with_cadd.tsv \
  --label-column label \
  --baseline-features CADD_RAW,CADD_PHRED \
  --group-column gene_id \
  --folds 5 \
  --outdir experiments/cadd_validation/outputs/pilot_real/eval_gaph_plus_cadd
```
