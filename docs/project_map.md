# Project Map

Use this file to navigate the repository before making changes.

## Production Boundary

Production logic is in:
- `main.nf`
- `nextflow.config`
- `lib/*.groovy`
- `bin/*.py`
- `genomics/*.py`
- `analytics/`
- `envs/*.yml`

`genomics/` is the shared domain library used by both pipeline commands and
completed-run analytics. `analytics/` owns reproducible analyses and report
generation. Standalone research experiments live under `experiments/`; they may
consume production outputs but keep their scratch data and generated reports
inside their own package directories.

Do not put experiments, ad hoc downloaded data, or smoke-test outputs in the
repository root. Use `/private/tmp`, `/tmp`, or cluster scratch for temporary
runs.

## Entrypoints

- Workflow entrypoint: `main.nf`
- Runtime configuration: `nextflow.config`
- Analytics report: `python -m analytics.strategy_report --run-dir <run-dir>`
- Bulk VEP annotation: `python -m analytics.vep_annotation`
- Run archive: `python -m run_archiving`
- Local run profile: `-profile local`
- Cluster run profile: `-profile slurm`

Runtime environments:
- `envs/controller.yml` - Nextflow, Java, and Micromamba for the login/controller host.
- `envs/fetch.yml` - Stage 1 task dependencies.
- `envs/alignment.yml` - alignment and annotation task dependencies.

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
  - copies target and ortholog FASTA files into the final layout
  - derives collapsed target structural features from the configured local target assembly GFF3
  - writes final `manifest.json`

- `bin/fetch_taxonomy_presets.py`
  - reads unique ortholog `tax_id` values
  - maps them through `assets/reference/ncbi/taxonomy/taxonomy_classes.json.gz`
  - writes `taxonomy_presets.tsv.gz`

- `bin/prepare_alignment_tasks.py`
  - validates Stage 1 outputs for alignment
  - rewrites per-gene FASTA inputs with stable short sequence IDs
  - writes per-gene task directories and `alignment_tasks.tsv.gz`

- `bin/run_minimap2_alignment.py`
  - runs fixed or taxonomy-adaptive minimap2
  - parses PAF `cs` evidence
  - writes alignment segments, events, summaries, and failures

- `bin/run_nucmer_alignment.py`
  - runs multi-query nucmer comparator without global one-to-one filtering
  - parses SAM-long alignment records and emits CIGAR-normalized events
  - writes the same normalized alignment evidence schema

- `bin/build_ensembl_compara_maf_manifest.py`
  - builds a small run-specific manifest of Ensembl Compara MAF chunks for the
    human chromosomes present in `genes.tsv.gz`
  - reads MAF directory listings and first human rows, not whole MAF files

- `bin/run_ensembl_compara_maf_alignment.py`
  - streams selected Ensembl Compara MAF chunks for one target gene
  - clips MSA evidence to the target gene interval
  - writes the same normalized alignment evidence schema with species as the
    support unit

- `bin/merge_ensembl_compara_maf_gene.py`
  - consolidates all source-chunk fragments for one gene
  - recomputes union-based MAF summaries and gene-local feature coverage

- `bin/merge_alignment_results.py`
  - merges per-gene/per-strategy evidence into bounded genomic partitions
  - keeps only annotation inputs in end-to-end partitions
  - writes report-ready Stage 2 summaries without globally rewriting raw events
  - intersects target features with alignment segments for coverage/depth summaries
  - writes a canonical small per-strategy summary for downstream reports
  - can write compact event support rows when `--compact_alignment_events true`
  - copies optional native outputs only when enabled

- `bin/annotate_events.py`
  - normalizes alignment events to VCF-style keys using target context
  - annotates events with ClinVar when a VCF is configured
  - preserves distinct per-strategy ortholog support in a compact table
  - streams event rows and fetches gnomAD regions within one bounded partition

- `bin/finalize_annotation_partitions.py`
  - streams partition annotations and strategy-support rows into canonical Stage 3 outputs
  - aggregates partition manifests without loading variant rows into memory

## Shared Domain Library

- `genomics/variants.py`
  - canonical variant keys, normalization, and target-context lookup
- `genomics/clinvar.py`
  - ClinVar review-star and significance semantics
- `genomics/gnomad.py`
  - gnomAD API requests and response normalization
- `genomics/gnomad_cache.py`
  - reusable regional gnomAD response cache

## Analytics Package

- `analytics/strategy_report.py`
  - command-line contract and report orchestration
- `analytics/annotation/`
  - VEP integration and consequence vocabularies
- `analytics/analyses/`
  - scientific calculations and bounded-memory aggregation
- `analytics/io/`
  - run input resolution and atomic artifact contracts
- `analytics/reporting/`
  - report sections, Plotly components, and final HTML document

Durable analytics artifacts belong under `<run-dir>/analytics/`. Presentation
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
```

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
assets/reference/ncbi/taxonomy/taxonomy_classes.json.gz
assets/reference/ensembl/compara/release-116/92_mammals.epo_extended/
```

`assets/inputs/gene_ids/` contains reusable input lists for local runs.
`assets/reference/` contains operational cache/reference inputs. Large
reference files are not Git-tracked source files. Small required lookup tables
under `assets/reference/` can be Git-tracked through `.gitignore` exceptions.
The workflow uses these paths as defaults when matching explicit parameters or
environment variables are not set.

## Local Tools

`tools/bin/` is an ignored local directory for symlinks or copies of external
CLI binaries that should not be committed. The workflow resolves the NCBI
Datasets CLI as `DATASETS_BIN`, then `tools/bin/datasets` when present, then
`datasets` on `PATH`.

## Design Direction

This repository uses Nextflow for orchestration and Python for small parsing
utilities. Keep that separation: Nextflow owns task graph, resources, retry, and
resume; Python owns deterministic file parsing and table generation.

See `docs/pipeline_scaling_notes.md` for known scaling risks and future
large-run refactoring directions.
