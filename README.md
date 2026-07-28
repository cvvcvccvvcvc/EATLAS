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
RUN="results/run_all_strategies_$(date +%Y%m%d_%H%M%S)"

nextflow run . \
  -profile local,conda,low_storage \
  --stage all \
  --ids_file assets/inputs/gene_ids/panel_10_genes.txt \
  --outdir "$RUN" \
  --alignment_strategies all
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
`assets/reference/clinvar/clinvar.vcf.gz` when present. Annotation requires a
ClinVar VCF and matching `.tbi`; the workflow fails early when neither the
parameter, environment variable, nor default asset is available.

For Slurm, combine the `slurm` and `low_storage` profiles and put `work/`,
results, and environment caches under the assigned shared scratch allocation.
The ITMO-specific bootstrap and validation procedure is documented in
`docs/itmo_cluster.md`.

```bash
export GAPH_ROOT="/mnt/tank/scratch/$USER/gaph_v2"
export GAPH_WORK_DIR="$GAPH_ROOT/work"
export NXF_CONDA_CACHEDIR="$GAPH_ROOT/conda/envs"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export NXF_HOME="$GAPH_ROOT/nextflow"

nextflow run . \
  -profile slurm,low_storage \
  --ids_file /path/to/gene_ids.txt \
  --outdir "$GAPH_ROOT/results/run_001"
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
  --events_tsv results/align_debug/alignment_events.tsv.gz \
  --segments_tsv results/align_debug/alignment_segments.tsv.gz \
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

The default end-to-end `--stage all` output is intentionally compact.

Fetch outputs:

- `fetch/manifest.json` - run constants and tool versions.
- `fetch/input.ids.tsv` - normalized input IDs.
- `fetch/genes.tsv.gz` - target gene metadata.
- `fetch/target_features.tsv.gz` - compact gene/exon/CDS/UTR/intron intervals.
- `fetch/orthologs.selected.tsv.gz` - selected ortholog sequence metadata.
- `fetch/failures.tsv.gz` - gene-level failures.
- `fetch/sequences/targets/*.fa.gz` - GRCh38 target gene sequences.

Alignment outputs:

- `alignment/manifest.json` - strategies and row counts.
- `alignment/strategy_summary.tsv.gz` - compact per-strategy aggregate.
- `alignment/feature_coverage.tsv.gz` - coverage/depth by target feature and strategy.
- `alignment/failures.tsv.gz` - alignment-stage failures.

Annotation outputs:

- `annotation/variant_annotations.tsv.gz` - compact unique variant-context rows
  with report-relevant ClinVar classification/review evidence and selected
  gnomAD AF/consequence fields.
- `annotation/variant_strategy_support.tsv.gz` - compact per-strategy ALT-support
  counts and, for SNVs, the distinct orthologs aligned at the variant site.
- `annotation/manifest.json` - annotation input, source, row-count, cache, and diagnostic counters.
- `annotation/failures.tsv.gz` - non-fatal external annotation lookup failures, such as gnomAD region fetch errors.

Annotation accepts only concrete A/C/G/T event alleles for external variant
lookup. Non-concrete legacy or external event rows are excluded before lookup
regions are built and counted in `annotation/manifest.json`.

Standalone `--stage fetch` and `--stage align` runs publish the full handoff
tables and ortholog FASTA required to start the following stage separately.

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
