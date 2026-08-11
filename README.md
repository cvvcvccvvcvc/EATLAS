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
RUN="results/run_default_strategies_$(date +%Y%m%d_%H%M%S)"

nextflow run . \
  --stage all \
  --ids_file assets/inputs/gene_ids/panel_10_genes.txt \
  --outdir "$RUN" \
  --alignment_strategies all
```

Every run uses the declared `envs/*.yml` task environments through Micromamba;
local execution needs no profile and never depends on the active shell's Python
or command-line tools. Set machine-specific data paths through environment
variables when needed:

```bash
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
unset, an explicitly selected precomputed Ensembl strategy uses the matching
manifest under `assets/reference/ensembl/compara/` when present, otherwise it
builds one during the run.
If `--clinvar_vcf` and `CLINVAR_VCF` are unset, annotation uses
`assets/reference/clinvar/clinvar.vcf.gz` when present. Annotation requires a
ClinVar VCF and matching `.tbi`; the workflow fails early when neither the
parameter, environment variable, nor default asset is available.

For Slurm, use the `slurm` profile and put `work/`, results, and environment
caches under the assigned shared scratch allocation.
The ITMO-specific bootstrap and validation procedure is documented in
`docs/itmo_cluster.md`.

```bash
export GAPH_ROOT="/mnt/tank/scratch/$USER/gaph_v2"
export GAPH_WORK_DIR="$GAPH_ROOT/work"
export NXF_CONDA_CACHEDIR="$GAPH_ROOT/conda/envs"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export NXF_HOME="$GAPH_ROOT/nextflow"

nextflow run . \
  -profile slurm \
  --ids_file /path/to/gene_ids.txt \
  --outdir "$GAPH_ROOT/results/run_001"
```

Alignment-only debug mode can reuse an existing fetch result:

```bash
nextflow run . \
  --stage align \
  --fetch_dir results/run_test/fetch \
  --outdir results/align_debug \
  -resume
```

Annotation-only debug mode can reuse an existing alignment event table:

```bash
nextflow run . \
  --stage annotate \
  --events_tsv results/align_debug/alignment_events.tsv.gz \
  --event_ortholog_support_tsv results/align_debug/event_ortholog_support.tsv.gz \
  --segments_tsv results/align_debug/alignment_segments.tsv.gz \
  --fetch_dir results/run_test/fetch \
  --outdir results/annotate_debug \
  -resume
```

`alignment_events.tsv.gz` and `event_ortholog_support.tsv.gz` form one handoff.
The compact event table contains one row per target event and strategy; the
sidecar retains its exact supporting orthologs and joins through
`event_group_id`. Standalone annotation requires both files.

By default, `--alignment_strategies all` runs `minimap2_asm10`,
`minimap2_asm20`, `nucmer`, and `bwa_pseudoreads`. The precomputed Ensembl
strategy remains available only when named explicitly. Use a comma-separated
list to select a different set:

```bash
nextflow run . \
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
- `annotation/variant_ortholog_support/*.parquet` - one row per normalized
  variant, strategy, and supporting ortholog, with tax ID, taxname, and
  observation count.
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
is the durable data layer for downstream stages. Successful execution sessions
clean their task work by default; failed or interrupted work remains available
for recovery. See `docs/storage_model.md` for resumed-run cleanup details.

Raw NCBI zip files, unpacked `gene.fna`, minimap2 PAF, MUMmer delta/coords
files, and external Ensembl Compara MAF chunks are not published by default.
Set `--keep_native_alignments true` only for targeted debug or benchmark runs.

For cluster runs, keep `-work-dir` on scratch storage, not in the project
directory or home quota.

## Stage Steps

| Step | Process | What Happens | Durable Output |
| --- | --- | --- | --- |
| 1 | `VALIDATE_IDS` | Read Entrez IDs, remove duplicates, split accepted IDs into chunks. | `fetch/input.ids.tsv`, `fetch/chunks.tsv` |
| 2 | `FETCH_PARSE_CHUNK` | Download one NCBI gene package with `--ortholog all --include gene`; parse `data_report.jsonl` and `gene.fna`; select GRCh38 human target and one sequence per ortholog GeneID. Concurrent request starts are spaced by a fixed 5 seconds. | Per-chunk compressed FASTA/TSV files in `work/`; durable metrics in `fetch/chunk_metrics.tsv.gz` |
| 3 | `BUILD_FETCH_DATASET` | Assemble chunk tables, selected per-gene FASTA files, and target structural features into the final fetch dataset. | `fetch/` |
| 4 | `FETCH_TAXONOMY` | Build compact taxonomy metadata for downstream taxonomic evidence. | `alignment/taxonomy.tsv.gz` |
| 5 | `BUILD_ALIGNMENT_TASKS` | Validate fetch outputs and create per-gene alignment inputs with stable sequence IDs. | `alignment/alignment_tasks.tsv.gz` |
| 6 | `ALIGN_MINIMAP2` | Run the selected fixed asm10/asm20 minimap2 baselines. | Per-gene normalized evidence in `work/` |
| 7 | `ALIGN_NUCMER_COMPARATOR` | Multi-query nucmer comparator without global one-to-one delta filtering. | Per-gene normalized evidence in `work/` |
| 8 | `ALIGN_BWA_PSEUDOREADS` | Fixed pseudoread comparator evidence. | Per-gene normalized evidence in `work/` |
| 9 | `BUILD_ENSEMBL_COMPARA_MAF_MANIFEST` | When selected, build a release-116 EPO Extended manifest covering target chromosomes. | Per-run manifest in `work/` |
| 10 | `ALIGN_ENSEMBL_COMPARA_MAF_CHUNK` | When selected, stream each required MAF chunk once for all overlapping genes. | Per-gene fragments in `work/` |
| 11 | `MERGE_ENSEMBL_COMPARA_MAF_GENE` | Consolidate all EPO fragments for one target gene. | Per-gene normalized evidence in `work/` |
| 12 | `MERGE_ALIGNMENT_EVIDENCE` | Merge normalized alignment evidence and summarize feature coverage. | `alignment/` |
| 13 | `ANNOTATE_EVENTS` | Normalize event keys and annotate with ClinVar/gnomAD evidence. | `annotation/` |
