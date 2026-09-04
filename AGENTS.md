# AGENTS.md

Always-loaded guidance for agents working in this repository. Read only the
owning document and code named for the current task.

## Purpose

GAPH v2 is one end-to-end Nextflow pipeline with three internal boundaries:

1. fetch normalized target, selected-ortholog, feature, and taxonomy evidence;
2. align ortholog sequences and preserve normalized row-level evidence;
3. annotate alignment events with ClinVar, gnomAD, and Ensembl VEP evidence.

Completed-run analytics is a separate, read-only consumer of durable pipeline
evidence. Run archiving is separate operational tooling.

## Architecture Invariants

1. Nextflow owns orchestration, resources, retry, resume, and process wiring.
   Python owns deterministic parsing, normalization, merging, and reporting.
2. The only production pipeline path is end to end. Fetch, alignment, and
   annotation are internal workflow boundaries; recovery uses `-resume`, not
   standalone stage modes.
3. Keep pipeline output evidence-first: retain normalized observations,
   non-reconstructable outcomes, stable join identities, and provenance needed
   after `work/` is deleted. Report-specific aggregates belong to `analytics/`.
4. Fetch external metadata once in its owning boundary and publish one stable
   representation. Do not add hidden downstream lookups for the same data.
5. Reject obsolete schemas and inputs explicitly. Do not add compatibility
   fallbacks, silent schema detection, parallel old/new paths, or speculative
   feature branches.
6. Shared domain logic belongs in an importable package and is executed with
   `python -m`. Do not use `PYTHONPATH`, path bridges, or import shims.
7. Treat `results/` as immutable durable evidence and Nextflow `work/` as a
   disposable resume cache. Raw provider packages and native aligner files are
   temporary task data.
8. Preserve the fixed target and fetch policy unless the design is explicitly
   changed: GRCh38.p14 (`GCF_000001405.40`) and NCBI `--ortholog all`.
9. Alignment defaults come only from `ALIGNMENT_STRATEGY_REGISTRY` in `main.nf`.
   Do not duplicate a default list or reintroduce an `all` selector.
10. Pipeline candidate VEP is part of annotation and publishes one partitioned
    durable dataset. Do not recreate bulk-VEP or report-precompute paths.
11. Large annotation runs reuse the shared gnomAD cache. On ITMO use
    `$GAPH_ROOT/cache/gnomad` through `--gnomad_cache_dir` or
    `GAPH_GNOMAD_CACHE_DIR`.

## Task Routing

| Task | Owning code | Read | Focused check |
| --- | --- | --- | --- |
| Pipeline parameters or top-level wiring | `nextflow_schema.json`, `main.nf`, `nextflow.config`, `modules/local/` | relevant stage contract | `tests/test_pipeline_launch_validation.py`, `tests/test_canonical_partition_wiring.py` |
| Fetch or taxonomy | `bin/normalize_ids.py`, `bin/fetch_parse_chunk.py`, `bin/build_fetch_dataset.py`, `bin/fetch_taxonomy.py` | `docs/stage1_fetch_contract.md` | matching normalize, fallback, dataset, and taxonomy tests |
| Alignment strategy or evidence | `ALIGNMENT_STRATEGY_REGISTRY`, alignment modules, `bin/alignment_*`, `bin/run_*alignment.py`, `bin/merge_alignment_results.py` | `docs/stage2_alignment_contract.md` | matching alignment runner/schema/merge tests |
| Annotation or external variant evidence | annotation modules, `bin/annotate_events.py`, `bin/annotate_vep_partition.py`, `bin/finalize_annotation_partitions.py`, `genomics/` | `docs/stage3_annotation_contract.md` | matching annotation, ClinVar, gnomAD, and VEP tests |
| Run provenance or evidence integrity | `lib/RunManifest.groovy`, `provenance/`, `modules/local/build_evidence_inventory.nf` | `docs/storage_model.md` | `tests/test_run_manifest.py`, `tests/test_evidence_inventory.py` |
| Analytics inputs, caches, or scientific results | `analytics/io/`, `analytics/derivations/`, `analytics/analyses/` | `docs/analytics_contract.md` | matching analysis/cache tests |
| HTML presentation | `analytics/reporting/`, `analytics/strategy_report.py` | `docs/analytics_contract.md` | `tests/test_strategy_report.py` and section-specific tests |
| Pipeline or report launch | the two shell launchers under `scripts/slurm/` and `analytics/slurm/` | `docs/pipeline_launch.md` or `docs/report_generation.md` | matching launcher test |
| Cluster bootstrap or infrastructure diagnosis | `envs/`, `nextflow.config`, `bin/gaph-vep116` | `docs/itmo_cluster.md` | `docs/run_validation.md` |
| Run archive, restore, or removal | `run_archiving/` | `run_archiving/README.md` | `tests/test_run_archiving.py` |

`README.md` is the human entry point. `docs/project_map.md` explains ownership
when a task crosses boundaries. Current CLI values belong to
`nextflow_schema.json`, `--help`, and launcher `--help`, not prose copies.

## Engineering Method

1. Read the current owning code and document before changing either. Separate
   measured behavior, durable contract, operator procedure, and historical
   evidence.
2. Make the smallest complete change. Prefer deletion and direct data flow over
   adapters or new abstractions.
3. Do not refactor for hypothetical scale. Reproduce or measure a problem on a
   representative workload first.
4. Keep cohesive logic together. Split modules only when ownership or
   dependency boundaries improve.
5. Make user parameters and table contracts explicit in schema, tests, and one
   owning document. Link elsewhere instead of copying them.
6. Fail with a concrete message on missing inputs, mismatched headers,
   incomplete partitions, or invalid empty output.
7. Validate in proportion to risk: use the focused check above, then a small
   end-to-end smoke for workflow changes. Compare scientific parity or resource
   use when those properties are part of the change.
8. Review the final diff for duplicated facts, stale links, unnecessary
   branches, and unrelated workspace changes. Preserve user-owned dirty files.
9. Keep the repository root clean. Use `/tmp`, `/private/tmp`, cluster scratch,
   `work/`, or documented fixtures for transient files.

## Cluster And Release Safety

Task environments from `envs/*.yml` are mandatory. Local runs need no profile;
ordinary ITMO pipeline and report work must use the documented launchers.

Before any cluster submission, pass the revision gate in
`docs/pipeline_launch.md`: fetch authoritative `origin/main`, require the
intended local commit, fetched `origin/main`, and cluster `HEAD` to match, and
require a clean cluster tree. If a documented launcher interface is absent,
resynchronize the checkout; never adapt a launch to obsolete code.

## Commits

- Commit a finished, narrow, verified change with one coherent purpose and
  stage only its files.
- Do not commit broad exploration or partially validated cross-cutting work.
- After a verified commit, push the current branch directly to origin,
  including existing unpushed commits. Do not create a pull request.

## Placement

- Standalone research and dated benchmarks belong under `experiments/<name>/`.
- First-article material belongs under `scientific_work/article_1/`; a
  structured supervisor document outranks optional narrative context.
- Native archive logic stays under `run_archiving/`.
