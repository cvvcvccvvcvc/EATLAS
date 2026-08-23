# CADD Validation

Historical prototype for testing whether GAPH ortholog-alignment evidence adds
predictive value beyond CADD and standard conservation features.

## Status

The scientific question and the 2026-06-29 pilot are retained, but the current
feature-builder CLI is not compatible with the evidence-first pipeline
contract. It expects old global Stage 2 segment, event, summary, and feature-
coverage files. Current runs publish partitioned evidence plus exact support,
and analytics owns feature coverage and taxonomic counters.

Do not use the old workflow commands as instructions for a current run. The
experiment is isolated from production and has not been migrated merely to keep
an unused path alive.

## Question

Do variant-level features derived from GAPH ortholog evidence improve held-out
pathogenicity prediction after accounting for CADD and conservation scores?

The intended ablation remains:

1. baseline features only;
2. GAPH evidence features only;
3. baseline plus GAPH features;
4. baseline plus shuffled GAPH features as a negative control.

Gene- or chromosome-aware splits are required; AUPRC is the primary metric.

## Preserved Material

- `experiments/cadd_validation/docs/validation_design.md` — scientific
  hypothesis and leakage controls.
- `experiments/cadd_validation/docs/data_contract.md` — archived CLI/output contract and the required
  migration boundary.
- `experiments/cadd_validation/docs/pilot_result.md` — dated result from the
  pre-partition pipeline.
- `src/cadd_validation/` and `tests/` — historical implementation and fixtures.

## Requirements Before Reactivation

A new implementation must:

- consume `alignment/evidence/partitions/*` instead of nonexistent global
  alignment tables;
- use `event_ortholog_support.tsv.gz` for exact supporters and the Stage 3 event
  map for canonical variant identity;
- derive coverage and taxonomy views from the same current analytics domain
  logic rather than accepting pipeline-owned aggregates;
- validate parity on a small current end-to-end run before any new biological
  claim.

Keep all experiment data, scratch files, and reports inside this directory or a
documented external scratch path.
