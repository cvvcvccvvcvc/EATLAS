# AGENTS.md

Minimal always-loaded onboarding for AI agents working in this repository.

## WHY

GAPH v2 is a structured rewrite of a gene-level variant discovery pipeline.
The current implemented scope is:
- stage 1: fetch normalized input data for Entrez Gene IDs
- stage 2: align selected ortholog gene sequences to fixed human target loci
- stage 3: annotate alignment events with ClinVar and gnomAD evidence

Stage 1 produces:
- human target gene sequences on GRCh38.p14 (`GCF_000001405.40`)
- selected non-human ortholog gene sequences from NCBI Datasets
- compact metadata tables describing target genes, selected orthologs,
  rejected ortholog candidates, failures, and run constants

Stage 2 produces:
- compact taxonomy metadata for selected orthologs
- normalized alignment segments
- compact alignment events with exact ortholog support
- per-ortholog alignment summaries for downstream variant-support logic

Stage 3 produces:
- annotated event tables in a separate annotation output layer

## WHAT

Core files:
- `main.nf` - Nextflow DSL2 workflow wiring.
- `nextflow.config` - local/slurm profiles and process resource policy.
- `lib/RunManifest.groovy` - root run provenance and completion manifest.
- `bin/normalize_ids.py` - input ID normalization and chunking.
- `bin/fetch_parse_chunk.py` - NCBI Datasets fetch + package parsing.
- `bin/build_fetch_dataset.py` - final fetch dataset assembly.
- `bin/fetch_taxonomy.py` - compact ortholog taxonomy table.
- `bin/prepare_alignment_tasks.py` - per-gene alignment task preparation.
- `bin/run_minimap2_alignment.py` - minimap2 execution and PAF parsing.
- `bin/run_nucmer_alignment.py` - nucmer execution and comparator parsing.
- `bin/merge_alignment_results.py` - final alignment evidence merge.
- `bin/annotate_events.py` - event key normalization and ClinVar/gnomAD annotation.
- `genomics/` - shared variant, ClinVar, and gnomAD domain logic.
- `analytics/strategy_report.py` - completed-run analytics and HTML report entrypoint.
- `analytics/analyses/` - bounded-memory scientific analyses.
- `analytics/reporting/` - report sections and HTML composition.
- `run_archiving/` - isolated CLI and Slurm wrapper for verified run archival.
- `envs/fetch.yml` - minimal conda environment for stage 1.
- `envs/alignment.yml` - CLI dependencies for stage 2.

## HOW

Default local run for all stages:

```bash
nextflow run . --ids_file assets/inputs/gene_ids/smoke_5_genes.txt --outdir results/run_001 -resume
```

Cluster run:

```bash
nextflow run . -profile slurm --ids_file /path/gene_ids.txt \
  --outdir /scratch/$USER/gaph_v2/results/run_001 \
  -work-dir /scratch/$USER/gaph_v2/work/run_001 -resume
```

Agent workflow rules:
1. Keep this file short and universally relevant.
2. Before changing code, read only relevant files from `docs/`.
3. Treat `results/` as durable output and Nextflow `work/` as disposable resume cache.
4. Do not publish raw NCBI zip files or unpacked `gene.fna` as final outputs.
5. Preserve fixed stage-1 constants unless the user explicitly changes the design:
   GRCh38.p14 (`GCF_000001405.40`) and NCBI `--ortholog all`.
6. Stage 2 native aligner outputs are debug artifacts; do not publish them by
   default.
7. Prefer small, focused changes and validate with a small local Nextflow smoke run.
   Task environments from `envs/*.yml` are mandatory for every run; local runs
   need no profile, while cluster runs use `-profile slurm`.
8. Keep the repository root clean: no ad hoc scripts, reports, downloaded tools,
   or smoke outputs in production paths. Use `/tmp`, `/private/tmp`, `work/`,
   or documented test fixtures.
9. Preserve modular boundaries: Nextflow owns orchestration and process wiring;
   Python owns deterministic parsing/merging/report generation.
10. Prefer registry/config-driven feature selection over scattered booleans.
    Defaults such as alignment strategy selection should mean "all registered"
    rather than a duplicated literal list.
11. Make contracts explicit in schema/docs/tests when changing user-facing
    parameters or output table shapes.
12. Do not silently accept missing inputs, mismatched table headers, or empty
    outputs caused by wiring bugs; fail with a concrete message.
13. Keep commits atomic: each commit should contain one coherent behavior,
    contract, or documentation change and exclude unrelated workspace noise.
14. Commit finished, narrow, verified changes directly once their scope is
    closed. Stage only the relevant files and use a message that names the
    user-facing purpose.
15. Before committing a fix, manually run the narrowest realistic check that
    exercises the changed behavior and verify it behaves as intended, not only
    that it exits successfully. Report any skipped check in the final status.
16. Do not auto-commit broad, exploratory, cross-cutting, or partially
    validated work. Leave it uncommitted with a clear status summary, or ask
    before committing when the boundary is unclear.
17. Keep alignment and annotation as separate stages; `--stage align` must not
    trigger annotation.
18. Put standalone experiments that build on the pipeline or its data under
    `experiments/<experiment_name>/` (create `experiments/` when needed). Keep
    each experiment's code, data, scratch files, reports, and generated outputs
    isolated inside that experiment directory unless an external scratch path is
    explicitly documented.
19. Keep first-article materials under `scientific_work/article_1/`. Treat
   `article_narrative.md` as an optional source of ideas, not a binding
   scientific contract. Once a structured manuscript document from the
   scientific supervisor is present there, use it as the primary article guide
   and refine the manuscript through maintainer feedback.
20. Large annotation runs must reuse the shared gnomAD regional cache. On the
   ITMO cluster use `$GAPH_ROOT/cache/gnomad`; pass it explicitly with
   `--gnomad_cache_dir` or export `GAPH_GNOMAD_CACHE_DIR` before launch.
21. This is a single-maintainer repository: after completing and committing
   verified changes, always push the current branch directly to origin,
   including any existing unpushed commits; do not create pull requests or
   invent collaboration/conflict workflows.
22. Keep durable pipeline outputs evidence-first. Preserve normalized row-level
   biological observations, non-reconstructable outcomes, and stable identifiers
   needed to reproduce downstream analyses after `work/` is deleted. Selected
   taxonomic scopes or ranks, scientific counters, thresholds, histograms, and
   plotting tables belong to `analytics/` and must never replace their source
   evidence. Small manifest counts are allowed only as reproducible integrity or
   operational-QC snapshots, not as scientific data products.
23. Prefer one current, explicit, straightforward implementation path. Do not
   introduce legacy compatibility fallbacks, silent schema or version
   autodetection, parallel old/new code paths, or speculative branching. Reject
   obsolete inputs with a concrete error. Add branches or abstractions only
   when required by a demonstrated current use case.

## Progressive Disclosure

Read these only when relevant:

- `docs/report_generation.md` - read only this operational runbook for an
  ordinary completed-run analytics report launch.
- `docs/pipeline_launch.md` - read only this operational runbook for an ordinary
  ITMO Nextflow launch or resume.
- `docs/project_map.md` - repository structure and ownership.
- `docs/stage1_fetch_contract.md` - stage-1 data model, selection rules, and outputs.
- `docs/stage2_alignment_contract.md` - stage-2 alignment model, outputs, and rationale.
- `docs/run_validation.md` - smoke tests, validation, and pipeline failure
  investigation; not required for an ordinary launch.
- `docs/storage_model.md` - Nextflow cache, result layout, and disk-space policy.
- `docs/itmo_cluster.md` - first-time ITMO setup, verified infrastructure facts,
  transfer, and smoke tests; not required for an ordinary launch.
- `scientific_work/article_1/` - first-article workspace.
  `article_narrative.md` is an optional idea source; the scientific
  supervisor's structured manuscript document, once added, is the primary
  article guide.
