# GAPH v2

GAPH v2 currently implements three stages of the gene-level comparative variant
pipeline:

1. fetch human target gene loci and NCBI ortholog gene sequences for Entrez Gene IDs
2. align selected ortholog gene sequences against the fixed human target loci
3. annotate emitted alignment events with ClinVar and gnomAD evidence

The pipeline is implemented with Nextflow DSL2. It keeps raw NCBI packages and
native aligner outputs in temporary task work directories by default and
publishes normalized compressed FASTA/TSV outputs.

## Run

Default local execution runs every stage in one command:

```bash
nextflow run . \
  -profile local,conda \
  --ids_file assets/inputs/gene_ids/smoke_5_genes.txt \
  --outdir results/run_test \
  -resume
```

Use the `conda` profile for normal runs so tasks use `envs/*.yml` through
micromamba/conda instead of whatever Python environment is active in the shell.
By default the workflow resolves the NCBI Datasets CLI as `DATASETS_BIN`, then
`tools/bin/datasets` when present, then `datasets` on `PATH`. Aligner binaries
are expected on `PATH`.
Set persistent local paths once through environment variables when needed:

```bash
export DATASETS_BIN=/path/to/datasets
export ENTREZ_EMAIL=you@example.org
export ENTREZ_API_KEY=your_ncbi_api_key
export GAPH_TARGET_ANNOTATION_GFF3=/path/to/genomic.gff.gz
export ENSEMBL_COMPARA_MAF_MANIFEST=/path/to/ensembl_compara_maf_manifest.tsv.gz
export CLINVAR_VCF=/path/to/clinvar.vcf.gz
```

For local runs, these values can also be stored in an ignored `.env` file using
the format shown in `.env.example`. `FETCH_PARSE_CHUNK` loads `.env` before
calling NCBI Datasets. The API key is passed to `datasets download` as
`--api-key`; the email is recorded as configured contact metadata for NCBI/API
auditability.

If `--target_annotation_gff3` and `GAPH_TARGET_ANNOTATION_GFF3` are unset, the
fetch stage uses
`assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz`.
If `--ensembl_compara_maf_manifest` and `ENSEMBL_COMPARA_MAF_MANIFEST` are
unset, the precomputed Ensembl strategy uses the matching manifest under
`assets/reference/ensembl/compara/` when present, otherwise it builds one during
the run.
If `--clinvar_vcf` and `CLINVAR_VCF` are unset, annotation uses
`assets/reference/clinvar/clinvar.vcf.gz` when present. If no ClinVar VCF is
configured, ClinVar evidence is skipped instead of using a placeholder file.

For Slurm, use the same workflow with the `slurm` profile and put `work/` on
scratch storage:

```bash
nextflow run . \
  -profile slurm,conda \
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

Annotation-only debug mode can reuse an existing alignment event table:

```bash
nextflow run . \
  -profile local \
  --stage annotate \
  --events_tsv results/run_test/alignment/alignment_events.tsv.gz \
  --fetch_dir results/run_test/fetch \
  --outdir results/annotate_debug \
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
  annotation/
```

Fetch outputs:

- `fetch/manifest.json` - run constants and tool versions.
- `fetch/input.ids.tsv` - normalized input IDs.
- `fetch/chunks.tsv` - deterministic chunk list used for NCBI requests.
- `fetch/chunk_metrics.tsv.gz` - per-chunk fetch timing and package-size metrics.
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
- `alignment/alignment_events.tsv.gz` - raw mismatch/indel evidence by default;
  unique event support rows with `--compact_alignment_events true`.
- `alignment/failures.tsv.gz` - alignment-stage failures.

Annotation outputs:

- `annotation/alignment_events_annotated.tsv.gz` - alignment events plus ClinVar and gnomAD annotation columns.

The target assembly is fixed to GRCh38.p14 (`GCF_000001405.40`). Ortholog
retrieval always uses the complete NCBI ortholog set (`--ortholog all`).

## Storage Model

Nextflow `work/` is the execution cache used by `-resume`. Published `results/`
is the durable data layer for downstream stages.

Raw NCBI zip files, unpacked `gene.fna`, minimap2 PAF, MUMmer delta/coords
files, and external Ensembl Compara MAF chunks are not published by default.
Set `--keep_native_alignments true` only for targeted debug or benchmark runs.

For cluster runs, keep `-work-dir` on scratch storage, not in the project
directory or home quota.

## Stage Steps

| Step | Process | What Happens | Durable Output |
| --- | --- | --- | --- |
| 1 | `VALIDATE_IDS` | Read Entrez IDs, remove duplicates, split accepted IDs into chunks. | `fetch/input.ids.tsv`, `fetch/chunks.tsv` |
| 2 | `FETCH_PARSE_CHUNK` | Download one NCBI gene package with `--ortholog all --include gene`; parse `data_report.jsonl` and `gene.fna`; select GRCh38 human target and one sequence per ortholog GeneID. Concurrent downloads are staggered by `--fetch_request_stagger_seconds`. | Per-chunk compressed FASTA/TSV files in `work/`; durable metrics in `fetch/chunk_metrics.tsv.gz` |
| 3 | `BUILD_FETCH_DATASET` | Assemble chunk tables, selected per-gene FASTA files, and target structural features into the final fetch dataset. | `fetch/` |
| 4 | `FETCH_TAXONOMY_PRESETS` | Build compact tax_id to minimap2 preset metadata. | `alignment/taxonomy_presets.tsv.gz` |
| 5 | `BUILD_ALIGNMENT_TASKS` | Validate fetch outputs and create per-gene alignment inputs with stable sequence IDs. | `alignment/alignment_tasks.tsv.gz` |
| 6 | `ALIGN_MINIMAP2_ASM10` | Fixed minimap2 baseline. | Per-gene normalized evidence in `work/` |
| 7 | `ALIGN_MINIMAP2_ASM20` | More permissive fixed minimap2 baseline. | Per-gene normalized evidence in `work/` |
| 8 | `ALIGN_MINIMAP2_TAXONOMY_ADAPTIVE` | Taxonomy-driven minimap2 presets. | Per-gene normalized evidence in `work/` |
| 9 | `ALIGN_NUCMER_COMPARATOR` | Multi-query nucmer comparator without global one-to-one delta filtering. | Per-gene normalized evidence in `work/` |
| 10 | `ALIGN_BWA_PSEUDOREADS` | Pseudoread comparator evidence. | Per-gene normalized evidence in `work/` |
| 11 | `BUILD_ENSEMBL_COMPARA_MAF_MANIFEST` | When selected, build a run-specific manifest for Ensembl Compara MAF chunks covering target chromosomes. | Per-run manifest in `work/` |
| 12 | `ALIGN_ENSEMBL_COMPARA_MAF` | When selected, stream precomputed Ensembl Compara MAF blocks and normalize species-level evidence. | Per-gene normalized evidence in `work/` |
| 13 | `MERGE_ALIGNMENT_EVIDENCE` | Merge normalized alignment evidence and summarize feature coverage. | `alignment/` |
| 14 | `ANNOTATE_EVENTS` | Normalize event keys and annotate with ClinVar/gnomAD evidence. | `annotation/` |
