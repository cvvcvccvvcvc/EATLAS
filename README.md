# GAPH v2

GAPH v2 currently implements two stages of the gene-level comparative variant
pipeline:

1. fetch human target gene loci and NCBI ortholog gene sequences for Entrez Gene IDs
2. align selected ortholog gene sequences against the fixed human target loci

The pipeline is implemented with Nextflow DSL2. It keeps raw NCBI packages and
native aligner outputs in temporary task work directories by default and
publishes normalized compressed FASTA/TSV outputs.

## Run

Default local execution runs fetch and alignment in one command:

```bash
nextflow run . \
  -profile local \
  --ids_file gene_ids.txt \
  --outdir results/run_test \
  -resume
```

For Slurm, use the same workflow with the `slurm` profile and put `work/` on
scratch storage:

```bash
nextflow run . \
  -profile slurm \
  --ids_file /path/to/gene_ids.txt \
  --outdir /scratch/$USER/gaph_v2/results/run_001 \
  -work-dir /scratch/$USER/gaph_v2/work/run_001 \
  -resume
```

Alignment-only debug mode can reuse an existing fetch result:

```bash
nextflow run . \
  -profile local \
  --stage align \
  --fetch_dir results/run_test/fetch \
  --outdir results/align_debug \
  -resume
```

By default, alignment runs every strategy registered in the workflow. Use a
comma-separated list to run a subset:

```bash
nextflow run . \
  -profile local \
  --stage align \
  --fetch_dir results/run_test/fetch \
  --outdir results/align_minimap2_asm20 \
  --alignment_strategies minimap2_asm20 \
  -resume
```

## Outputs

Default output layout:

```text
results/run_test/
  fetch/
  alignment/
```

Fetch outputs:

- `fetch/manifest.json` - run constants and tool versions.
- `fetch/input.ids.tsv` - normalized input IDs.
- `fetch/chunks.tsv` - deterministic chunk list used for NCBI requests.
- `fetch/genes.tsv.gz` - target gene metadata.
- `fetch/target_features.tsv.gz` - collapsed target structural intervals.
- `fetch/orthologs.selected.tsv.gz` - selected ortholog sequence metadata.
- `fetch/orthologs.candidates.tsv.gz` - candidate records and deterministic selection decisions.
- `fetch/failures.tsv.gz` - gene-level failures.
- `fetch/sequences/targets/*.fa.gz` - GRCh38 target gene sequences.
- `fetch/sequences/orthologs/*.fa.gz` - selected non-human ortholog gene sequences.

Alignment outputs:

- `alignment/manifest.json` - strategies and row counts.
- `alignment/alignment_tasks.tsv.gz` - per-gene alignment task manifest.
- `alignment/taxonomy_presets.tsv.gz` - tax_id to minimap2 preset metadata.
- `alignment/taxonomy_failures.tsv.gz` - taxonomy lookup warnings/failures.
- `alignment/ortholog_alignment_summary.tsv.gz` - one row per gene/ortholog/strategy.
- `alignment/alignment_segments.tsv.gz` - normalized alignment intervals.
- `alignment/feature_coverage.tsv.gz` - coverage/depth by target feature and strategy.
- `alignment/alignment_events.tsv.gz` - raw mismatch/indel evidence.
- `alignment/failures.tsv.gz` - alignment-stage failures.

The target assembly is fixed to GRCh38.p14 (`GCF_000001405.40`). Ortholog
retrieval always uses the complete NCBI ortholog set (`--ortholog all`).

## Storage Model

Nextflow `work/` is the execution cache used by `-resume`. Published `results/`
is the durable data layer for downstream stages.

Raw NCBI zip files, unpacked `gene.fna`, minimap2 PAF, and MUMmer delta/coords
files are not published by default. Set `--keep_native_alignments true` only for
targeted debug or benchmark runs.

For cluster runs, keep `-work-dir` on scratch storage, not in the project
directory or home quota.

## Stage Steps

| Step | Process | What Happens | Durable Output |
| --- | --- | --- | --- |
| 1 | `VALIDATE_IDS` | Read Entrez IDs, remove duplicates, split accepted IDs into chunks. | `fetch/input.ids.tsv`, `fetch/chunks.tsv` |
| 2 | `FETCH_PARSE_CHUNK` | Download one NCBI gene package with `--ortholog all --include gene`; parse `data_report.jsonl` and `gene.fna`; select GRCh38 human target and one sequence per ortholog GeneID. | Per-chunk compressed FASTA/TSV files in `work/` |
| 3 | `MERGE_FETCH_RESULTS` | Merge chunk tables, copy selected per-gene FASTA files, and derive target structural features. | `fetch/` |
| 4 | `FETCH_TAXONOMY_PRESETS` | Build compact tax_id to minimap2 preset metadata. | `alignment/taxonomy_presets.tsv.gz` |
| 5 | `BUILD_ALIGNMENT_TASKS` | Validate fetch outputs and create per-gene alignment inputs with stable sequence IDs. | `alignment/alignment_tasks.tsv.gz` |
| 6 | `ALIGN_MINIMAP2_ASM10` | Fixed minimap2 baseline. | Per-gene normalized evidence in `work/` |
| 7 | `ALIGN_MINIMAP2_ASM20` | More permissive fixed minimap2 baseline. | Per-gene normalized evidence in `work/` |
| 8 | `ALIGN_MINIMAP2_TAXONOMY_ADAPTIVE` | Taxonomy-driven minimap2 presets. | Per-gene normalized evidence in `work/` |
| 9 | `ALIGN_NUCMER_COMPARATOR` | Multi-query nucmer comparator without global one-to-one delta filtering. | Per-gene normalized evidence in `work/` |
| 10 | `ALIGN_BWA_PSEUDOREADS` | Pseudoread comparator evidence. | Per-gene normalized evidence in `work/` |
| 11 | `MERGE_ALIGNMENT_EVIDENCE` | Merge normalized alignment evidence and summarize feature coverage. | `alignment/` |
