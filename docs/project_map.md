# Project Map

Use this file to navigate the repository before making changes.

## Production Boundary

Production logic is in:
- `main.nf`
- `nextflow.config`
- `bin/*.py`
- `envs/*.yml`

Do not put experiments, ad hoc downloaded data, or smoke-test outputs in the
repository root. Use `/private/tmp`, `/tmp`, or cluster scratch for temporary
runs.

## Entrypoints

- Workflow entrypoint: `main.nf`
- Runtime configuration: `nextflow.config`
- Local run profile: `-profile local`
- Cluster run profile: `-profile slurm`

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

- `bin/merge_fetch_results.py`
  - merges per-chunk TSV files
  - copies target and ortholog FASTA files into the final layout
  - derives collapsed target structural features from the target assembly GFF3
  - writes final `manifest.json`

- `bin/fetch_taxonomy_presets.py`
  - reads unique ortholog `tax_id` values
  - queries compact NCBI taxonomy summaries
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
  - parses `show-coords` and `show-snps`
  - writes the same normalized alignment evidence schema

- `bin/merge_alignment_results.py`
  - merges per-gene/per-strategy alignment evidence tables
  - intersects target features with alignment segments for coverage/depth summaries
  - copies optional native outputs only when enabled

## Output Boundary

Durable final outputs are published under `params.outdir`.

Default end-to-end layout:

```text
results/run_001/
  fetch/
  alignment/
```

Temporary files that must not be treated as final outputs:
- NCBI raw zip packages
- unpacked `ncbi_dataset/`
- unpacked `gene.fna`
- Nextflow `work/`
- `.nextflow*` local execution metadata

## Design Direction

This repository uses Nextflow for orchestration and Python for small parsing
utilities. Keep that separation: Nextflow owns task graph, resources, retry, and
resume; Python owns deterministic file parsing and table generation.
