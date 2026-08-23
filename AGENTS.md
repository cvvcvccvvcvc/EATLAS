# AGENTS.md

Minimal always-loaded guidance for agents working in this repository.

## Purpose

GAPH v2 is one end-to-end Nextflow pipeline with three internal boundaries:

1. fetch normalized target, selected-ortholog, feature, and taxonomy data;
2. align ortholog sequences and preserve normalized row-level evidence;
3. annotate alignment events with ClinVar, gnomAD, and Ensembl VEP evidence.

Completed-run analytics is a separate consumer. It derives scientific counts,
taxonomic views, coverage, support summaries, and reports from durable pipeline
evidence.

## Architecture Invariants

1. Nextflow owns orchestration, resources, retry, resume, and process wiring.
   Python owns deterministic parsing, normalization, merging, and reporting.
2. Expose one end-to-end pipeline execution path. Fetch, alignment, and
   annotation remain internal workflow boundaries; recovery uses `-resume`, not
   standalone stage modes.
3. Keep pipeline outputs evidence-first. Preserve normalized biological
   observations, non-reconstructable outcomes, stable join identities, and
   provenance needed after `work/` is deleted.
4. Report-specific scopes, ranks, counters, thresholds, bins, histograms, and
   plotting tables belong to `analytics/`. Manifest counts are integrity or
   operational-QC snapshots, not substitutes for source evidence.
5. Fetch external metadata once in its owning boundary, publish one stable
   reusable representation, and forbid hidden downstream network lookups for
   the same data.
6. Prefer one current pipeline and analytics implementation path. Reject
   obsolete inputs explicitly; do not add compatibility fallbacks, silent
   schema autodetection, parallel old/new paths, or speculative feature
   branches.
7. Share domain logic through a real package boundary. Stage internal packages
   explicitly and run them with `python -m`; do not use `PYTHONPATH`, path
   bridges, or import shims to compensate for misplaced code.
8. Treat `results/` as durable output and Nextflow `work/` as a disposable
   resume cache. Raw NCBI packages, unpacked `gene.fna`, and native aligner files
   remain temporary task data.
9. Preserve the fixed target and fetch policy unless the user changes the
   design: GRCh38.p14 (`GCF_000001405.40`) and NCBI `--ortholog all`.
10. Alignment defaults come from registry metadata. Do not duplicate the
    default strategy list or reintroduce a special `all` selection.
11. Large annotation runs reuse the shared gnomAD regional cache. On ITMO use
    `$GAPH_ROOT/cache/gnomad` through `--gnomad_cache_dir` or
    `GAPH_GNOMAD_CACHE_DIR`.
12. Pipeline candidate VEP is part of annotation and publishes one partitioned
    durable dataset. Do not recreate a separate bulk-VEP/report precompute.

## Engineering Method

1. Read the current code and only the relevant document before proposing a
   change. Distinguish measured current behavior from historical notes and
   assumptions.
2. Choose the smallest complete change that satisfies the current use case.
   Prefer deletion and a direct data flow over adapters or new abstractions.
3. Do not refactor a working path for hypothetical scale. First reproduce or
   measure the problem on a representative workload; keep known but unmeasured
   risks as explicit follow-up observations.
4. Keep cohesive logic together. Split modules when ownership or dependency
   boundaries improve, not merely because a file is long.
5. Make user-facing parameters and table contracts explicit in schema, tests,
   and the one owning document. Do not duplicate generated CLI defaults across
   several runbooks.
6. Fail with a concrete message on missing inputs, mismatched headers,
   incomplete partitions, or invalid empty output. Never convert a wiring error
   into a smaller successful dataset.
7. Validate in proportion to risk: run the narrowest realistic behavior check,
   relevant tests, and a small end-to-end smoke for workflow changes. Compare
   scientific parity, disk, memory, or runtime when those properties are part
   of the change.
8. Review the final diff for unnecessary branching, duplicated data, stale
   documentation, and unrelated workspace changes. Preserve user-owned dirty
   files.
9. Keep the repository root clean. Use `/tmp`, `/private/tmp`, cluster scratch,
   `work/`, or documented fixtures for transient files.

## Commits

- Commit a finished, narrow, verified change directly with one coherent
  purpose. Stage only its files and name the user-facing outcome.
- Do not auto-commit broad exploration or partially validated cross-cutting
  work; report its state instead.
- This is a single-maintainer repository. After a verified commit, push the
  current branch directly to origin, including existing unpushed commits. Do
  not create a pull request or invent a collaboration workflow.

## Scope-Specific Placement

- Standalone research belongs under `experiments/<name>/`, including its data,
  scratch files, reports, and generated outputs.
- First-article materials belong under `scientific_work/article_1/`. Treat
  `article_narrative.md` as optional context; a structured supervisor document
  is the primary article guide when present.
- Native run-archiving logic stays isolated under `run_archiving/`.

## Common Commands

Local end-to-end run:

```bash
nextflow run . \
  --ids_file assets/inputs/gene_ids/smoke_5_genes.txt \
  --outdir results/run_001 \
  -resume
```

Cluster run:

```bash
nextflow run . -profile slurm \
  --ids_file /path/gene_ids.txt \
  --outdir /mnt/tank/scratch/$USER/gaph_v2/results/run_001 \
  -work-dir /mnt/tank/scratch/$USER/gaph_v2/work/run_001 \
  -resume
```

Task environments from `envs/*.yml` are mandatory. Local runs need no profile;
ITMO runs use `-profile slurm`.

## Progressive Disclosure

- `README.md` — human overview, quick start, outputs, and document index.
- `docs/pipeline_launch.md` — ordinary ITMO launch or resume.
- `docs/report_generation.md` — ordinary report or combined pipeline/report launch.
- `docs/run_validation.md` — smoke tests and failure investigation.
- `docs/project_map.md` — repository structure and ownership.
- `docs/stage1_fetch_contract.md` — fetch selection and durable data contract.
- `docs/stage2_alignment_contract.md` — alignment evidence contract.
- `docs/stage3_annotation_contract.md` — external annotation and variant-shard contract.
- `docs/storage_model.md` — durable data, resume cache, and disk policy.
- `docs/itmo_cluster.md` — first-time ITMO setup and verified infrastructure.
- `run_archiving/README.md` — verified archive, restore, and removal operations.
