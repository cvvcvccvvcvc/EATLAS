# Project Map

Use this file to navigate the repository before making changes.

## Production Boundary

Production logic is in:
- `main.nf`
- `nextflow.config`
- `lib/*.groovy`
- `modules/local/*.nf`
- `bin/*.py`
- `genomics/`
- `analytics/`
- `envs/*.yml`

`genomics/` is the shared domain library used by both pipeline commands and
completed-run analytics. `analytics/` owns reproducible analyses and report
generation. Standalone research experiments belong under `experiments/`; they may
consume production outputs but keep their scratch data and generated reports
inside their own package directories.

Do not put experiments, ad hoc downloaded data, or smoke-test outputs in the
repository root. Use `/private/tmp`, `/tmp`, or cluster scratch for temporary
runs.

## Entrypoints

- Cluster pipeline launch or resume: `scripts/slurm/run_pipelines.sh`
- Pipeline workflow entrypoint: `main.nf`
- Runtime configuration: `nextflow.config`
- Analytics report: `python -m analytics.strategy_report --analytics-root <root> --run-dir <run> --report-name <name>`
- Cluster report submission: `analytics/slurm/submit_strategy_report.sh`
- Run archive: `python -m run_archiving`
- Local execution: no profile required
- Cluster run profile: `-profile slurm`

The pipeline and report shell launchers above are the two user-facing cluster
commands. `analytics/slurm/strategy_report.sbatch` is the internal Slurm report
worker, not a third launch mode.

Runtime environments:
- `envs/controller.yml` - Nextflow, Java, and Micromamba for the login/controller host.
- `envs/fetch.yml` - Stage 1 task dependencies.
- `envs/alignment.yml` - alignment and ClinVar/gnomAD annotation dependencies.
- `envs/vep.yml` - bounded pipeline VEP task dependencies.
- `envs/analytics.yml` - completed-run analytics and report dependencies.

## Core Modules

- `bin/normalize_ids.py`
  - reads an input ID file
  - accepts whitespace/comma-separated Entrez Gene IDs
  - removes duplicates while preserving first occurrence order
  - writes `input.ids.tsv`, `chunks.tsv`, and `chunks/*.ids.txt`

- `bin/fetch_parse_chunk.py`
  - runs `datasets download gene gene-id --ortholog all --include gene`
  - unpacks the package inside the task work directory
  - parses `data_report.jsonl` and `gene.fna`
  - selects the GRCh38 human target sequence
  - selects one non-human sequence per ortholog GeneID
  - writes compressed per-chunk FASTA and TSV outputs

- `bin/build_fetch_dataset.py`
  - merges per-chunk TSV files
  - copies target and ortholog FASTA files into the internal fetch handoff
  - derives collapsed target structural features from the configured local target assembly GFF3
  - writes final `manifest.json`

- `bin/fetch_taxonomy.py`
  - reads unique ortholog `tax_id` values
  - fetches ordered lineage and direct domain-through-species metadata from NCBI Datasets
  - writes the Stage 1 `taxonomy.tsv.gz` handoff once per fetch dataset

- `bin/prepare_alignment_tasks.py`
  - validates Stage 1 outputs for alignment
  - rewrites per-gene FASTA inputs with stable short sequence IDs
  - writes per-gene task directories and `alignment_tasks.tsv.gz`

- `bin/alignment_table_schema.py` and `bin/alignment_task_io.py`
  - define exact aligner output schemas and shared task I/O validation

- `bin/run_minimap2_alignment.py`
  - runs fixed-preset minimap2 on complete orthologs or deterministic long pseudo-reads
  - reduces long-read placements to a dominant-strand monotonic backbone
  - parses PAF `cs` evidence
  - writes alignment segments, events, summaries, and failures

- `bin/run_nucmer_alignment.py`
  - runs multi-query nucmer comparator without global one-to-one filtering
  - parses SAM-long alignment records and emits CIGAR-normalized events
  - writes the same normalized alignment evidence schema

- `bin/run_bwa_pseudoreads.py` and `bin/bwa_pseudoread_filter.py`
  - generate the fixed BWA pseudoread comparator and retain its monotonic
    target-order alignment backbone

- `bin/merge_alignment_results.py`
  - merges per-gene/per-strategy evidence into bounded genomic partitions
  - publishes the canonical partition contract
  - writes compact events and their `event_group_id`-keyed positive ortholog
    handoff in one index-ordered pass
  - copies final partitions without global evidence recompression

- `bin/annotate_events.py`
  - normalizes alignment events to VCF-style keys using target context
  - annotates events with ClinVar when a VCF is configured
  - emits one event-to-canonical-variant lineage row per compact event
  - streams event rows and fetches gnomAD regions within one bounded partition

- `bin/annotate_vep_partition.py`
  - adds one declared VEP consequence record to each valid variant/gene row
  - preserves source order and writes one independently resumable enriched shard

- `bin/prepare_annotation_contexts.py`
  - validates partition membership and materializes only each partition's
    target metadata and FASTA context

- `bin/finalize_annotation_partitions.py`
  - validates the exact base/VEP shard relation and semantic VEP configuration
  - copies enriched shards and event maps without a global variant-table rewrite

- `analytics/io/alignment_aggregates.py`
  - derives strategy summary and feature coverage from partition evidence

- `analytics/io/annotation_support.py`
  - derives variant-strategy and taxonomic ortholog-evidence report tables from
    Stage 1 taxonomy, Stage 2 evidence, and Stage 3 event lineage

- `analytics/derivations/`
  - deterministic strategy, coverage, taxonomy, and exact-support derivations
    shared by analytics cache builders

## Shared Domain Library

- `genomics/variants.py`
  - canonical variant keys, normalization, and target-context lookup
- `genomics/clinvar.py`
  - ClinVar review-star and significance semantics
- `genomics/gnomad.py`
  - gnomAD API requests and response normalization
- `genomics/gnomad_cache.py`
  - reusable regional gnomAD response cache
- `genomics/gnomad_index.py`
  - derived immutable Parquet fragments for exact-allele analytics lookups
  - preserves complete-tile coverage so absence remains distinct from failure
- `genomics/taxonomy.py`
  - canonical taxonomy schema, normalization, and lineage/rank access
- `genomics/vep/`
  - shared VEP provider integration, release-pinned semantics, and immutable
    cross-run result cache

## Analytics Package

- `analytics/strategy_report.py`
  - command-line contract and report orchestration
- `analytics/vep/`
  - report-facing consequence vocabularies and grouping rules
- `analytics/analyses/`
  - scientific calculations and bounded-memory aggregation
- `analytics/derivations/`
  - reusable deterministic builders for source-evidence-derived tables
- `analytics/io/`
  - run input resolution and atomic artifact contracts
- `analytics/reporting/`
  - report sections, Plotly components, and final HTML document

Analytics artifacts belong under the explicitly selected external analytics
root. `analytics/io/run_inputs.py` resolves one or more immutable source runs,
reuses per-source caches, and assigns an analysis workspace. Presentation
modules do not fetch data or own scientific calculations.

## Run Archiving

`run_archiving/` is isolated operational tooling for copying complete run
directories to an rclone remote, verifying their content, restoring them, and
removing a local copy only after a fresh remote verification. It does not
participate in the Nextflow workflow or analytics package. Its environment,
Slurm wrapper, and usage contract are documented in
`run_archiving/README.md`.

## Output Boundary

Durable final outputs are published under `params.outdir`.

Default end-to-end layout:

```text
results/run_001/
  run_manifest.json
  fetch/
  alignment/
  annotation/
  reports/
    nextflow/
```

Completed-run analytics is stored separately under
`<analytics-root>/cache/` and `<analytics-root>/analyses/`. Source run
directories are read-only after successful completion.

Temporary files that must not be treated as final outputs:
- NCBI raw zip packages
- unpacked `ncbi_dataset/`
- unpacked `gene.fna`
- Nextflow `work/`
- `.nextflow*` local execution metadata

## Assets

Reusable local inputs and reference files live under `assets/`, grouped by role,
provider, and resource family:

```text
assets/inputs/gene_ids/
assets/reference/clinvar/
assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz
```

`assets/inputs/gene_ids/` contains reusable input lists for local runs.
`assets/reference/` contains operational cache/reference inputs. Large
reference files are not Git-tracked source files. Small required lookup tables
under `assets/reference/` can be Git-tracked through `.gitignore` exceptions.
The workflow uses these paths as defaults when matching explicit parameters or
environment variables are not set.

## Documentation Ownership

Each operational question has one primary document:

- ordinary pipeline launch or resume: `docs/pipeline_launch.md`
- ordinary report launch: `docs/report_generation.md`
- smoke tests and failure investigation: `docs/run_validation.md`
- first-time ITMO setup: `docs/itmo_cluster.md`
- durable-versus-temporary data: `docs/storage_model.md`
- fetch, alignment, and annotation semantics: the three stage contract documents

README files provide navigation and package boundaries; they do not duplicate
complete runbooks. Historical experiment notes describe the run that produced
them and are not current production contracts.
