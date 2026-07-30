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
- taxonomy-driven minimap2 preset metadata
- normalized alignment segments
- normalized raw alignment events
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
- `bin/fetch_taxonomy_presets.py` - compact taxonomy-to-preset table.
- `bin/prepare_alignment_tasks.py` - per-gene alignment task preparation.
- `bin/run_minimap2_alignment.py` - minimap2 execution and PAF parsing.
- `bin/run_nucmer_alignment.py` - nucmer execution and comparator parsing.
- `bin/merge_alignment_results.py` - final alignment evidence merge.
- `bin/annotate_events.py` - event key normalization and ClinVar/gnomAD annotation.
- `genomics/` - shared variant, ClinVar, and gnomAD domain logic.
- `analytics/strategy_report.py` - completed-run analytics and HTML report entrypoint.
- `analytics/analyses/` - bounded-memory scientific analyses.
- `analytics/reporting/` - report sections and HTML composition.
- `envs/fetch.yml` - minimal conda environment for stage 1.
- `envs/alignment.yml` - CLI dependencies for stage 2.

## HOW

Default local run for all stages:

```bash
nextflow run . -profile local,conda --ids_file assets/inputs/gene_ids/smoke_5_genes.txt --outdir results/run_001 -resume
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
   Use `-profile local,conda` for normal local pipeline runs so tasks use the
   declared `envs/*.yml` environments instead of the caller's active shell env.
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
19. Before changing scientific hypotheses, validation statistics, or analytics
   report content, read `docs/article_narrative.md`. Keep the main argument,
   supporting analyses, and unresolved design decisions consistent with it;
   update the document when the research design changes materially.

## Progressive Disclosure

Read these only when relevant:

- `docs/project_map.md` - repository structure and ownership.
- `docs/stage1_fetch_contract.md` - stage-1 data model, selection rules, and outputs.
- `docs/stage2_alignment_contract.md` - stage-2 alignment model, outputs, and rationale.
- `docs/run_validation.md` - commands for local/slurm runs and verification.
- `docs/storage_model.md` - Nextflow cache, result layout, and disk-space policy.
- `docs/itmo_cluster.md` - verified ITMO Slurm setup, paths, transfer, and smoke tests.
- `docs/article_narrative.md` - mandatory scientific narrative for validation
  analyses and analytics report development.
