# GAPH v2

GAPH v2 implements one end-to-end gene-level comparative variant pipeline with
three internal workflow boundaries:

1. fetch human target gene loci and NCBI ortholog gene sequences for Entrez Gene IDs
2. align selected ortholog gene sequences against the fixed human target loci
3. annotate emitted alignment events with ClinVar, gnomAD, and Ensembl VEP evidence

The pipeline is implemented with Nextflow DSL2. It keeps raw NCBI packages and
native aligner outputs in temporary task work directories by default and
publishes normalized compressed FASTA/TSV outputs.

## Run

The pipeline has one end-to-end execution path:

```bash
RUN="results/run_default_strategies_$(date +%Y%m%d_%H%M%S)"

nextflow run . \
  --ids_file assets/inputs/gene_ids/panel_10_genes.txt \
  --outdir "$RUN" \
  --alignment_strategies default
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
export GAPH_VEP_BACKEND=local
export GAPH_VEP_RELEASE=116
export GAPH_VEP_EXECUTABLE=/path/to/gaph-vep116
export GAPH_VEP_CACHE_DIR=/path/to/vep-cache
export GAPH_VEP_RESULT_CACHE_DIR=/path/to/shared/vep-results
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
Candidate VEP annotation is part of the same workflow. Small local runs default
to the REST backend. Production cluster runs use the release-pinned local VEP
configuration from the cluster environment.

For Slurm, use the `slurm` profile and put `work/`, results, and environment
caches under the assigned shared scratch allocation.
The ITMO-specific bootstrap and validation procedure is documented in
`docs/itmo_cluster.md`.

```bash
export GAPH_ROOT="/mnt/tank/scratch/$USER/gaph_v2"
export GAPH_WORK_DIR="$GAPH_ROOT/work"
export GAPH_GNOMAD_CACHE_DIR="$GAPH_ROOT/cache/gnomad"
export NXF_CONDA_CACHEDIR="$GAPH_ROOT/conda/envs"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export NXF_HOME="$GAPH_ROOT/nextflow"

nextflow run . \
  -profile slurm \
  --ids_file /path/to/gene_ids.txt \
  --gnomad_cache_dir "$GAPH_GNOMAD_CACHE_DIR" \
  --outdir "$GAPH_ROOT/results/run_001" \
  -resume
```

By default, `--alignment_strategies default` runs `minimap2_asm10`,
`minimap2_asm20`, `nucmer`, and `bwa_pseudoreads_150_75`. Both the precomputed
Ensembl strategy and `minimap2_map_ont_pseudoreads_30000_15000` remain available
only when named explicitly. Use a comma-separated list to select a different set:

```bash
nextflow run . \
  --ids_file assets/inputs/gene_ids/panel_10_genes.txt \
  --outdir results/run_minimap2_asm20 \
  --alignment_strategies minimap2_asm20 \
  -resume
```

The long-pseudoread comparator is selected with:

```bash
--alignment_strategies minimap2_map_ont_pseudoreads_30000_15000
```

## Outputs

Default durable output layout:

```text
results/run_test/
  run_manifest.json
  fetch/
  alignment/
  annotation/
  reports/nextflow/
```

The published end-to-end output is intentionally compact.

Fetch outputs:

- `fetch/manifest.json` - run constants and tool versions.
- `fetch/input.ids.tsv` - normalized input IDs.
- `fetch/genes.tsv.gz` - target gene metadata.
- `fetch/target_features.tsv.gz` - compact gene/exon/CDS/UTR/intron intervals.
- `fetch/orthologs.selected.tsv.gz` - selected ortholog sequence metadata.
- `fetch/taxonomy.tsv.gz` - canonical status, ordered lineage, and direct domain-through-species metadata for selected ortholog tax IDs.
- `fetch/taxonomy_failures.tsv.gz` - missing taxonomy responses.
- `fetch/failures.tsv.gz` - gene-level failures.
- `fetch/sequences/targets/*.fa.gz` - GRCh38 target gene sequences.

Alignment outputs:

- `alignment/manifest.json` - strategies and row counts.
- `alignment/evidence/partitions/<partition_id>/` - normalized per-ortholog
  summaries, segments, compact events, exact ortholog support, and the partition
  manifest; gzip evidence files are retained without a global rewrite.
- `alignment/failures.tsv.gz` - alignment-stage failures.

Annotation outputs:

- `annotation/variant_annotations/manifest.json` and
  `annotation/variant_annotations/partitions/<partition_id>/<shard_id>.tsv.gz`
  - compact unique variant-context rows with ClinVar, gnomAD, and VEP evidence;
  headered bounded shards are published once without a duplicate global table.
- `annotation/event_variant_map/partitions/<partition_id>/event_variant_map.tsv.gz`
  - one provenance row per compact alignment event, linking its partition-local
  `event_group_id` to the canonical variant key and normalization status.
- `annotation/manifest.json` - annotation input, source, row-count, cache, and diagnostic counters.
- `annotation/failures.tsv.gz` - non-fatal external annotation lookup failures, such as gnomAD region fetch errors.

Annotation accepts only concrete A/C/G/T event alleles for external variant
lookup. Non-concrete event rows are excluded before lookup
regions are built and counted in `annotation/manifest.json`.

The report derives strategy summaries, feature coverage, site depth, taxonomic
evidence, and variant-strategy support under an explicitly selected external
analytics root. These analytic tables are not pipeline outputs, and completed
source run directories remain unchanged.

Fetch, alignment, and annotation are internal workflow boundaries, not separate
CLI modes. Recovery uses Nextflow `-resume` with the same inputs, result
directory, and work directory. Removed stage-selection parameters are rejected
instead of being ignored.

The target assembly is fixed to GRCh38.p14 (`GCF_000001405.40`). Ortholog
retrieval always uses the complete NCBI ortholog set (`--ortholog all`).

## Storage Model

Nextflow `work/` is the execution cache used by `-resume`. Published `results/`
is the durable data layer for analytics and archival. Successful execution sessions
clean their task work by default; failed or interrupted work remains available
for recovery. See `docs/storage_model.md` for resumed-run cleanup details.

Raw NCBI zip files, unpacked `gene.fna`, minimap2 PAF, MUMmer delta/coords
files, and external Ensembl Compara MAF chunks are not published. Inspect a
retained failed/interrupted task work directory when native debugging is needed.

For cluster runs, keep `-work-dir` on scratch storage, not in the project
directory or home quota.

## Internal Workflow Steps

| Step | Process | What Happens | Output Role |
| --- | --- | --- | --- |
| 1 | `VALIDATE_IDS` | Normalize Entrez IDs and plan deterministic fetch chunks. | Intermediate plan |
| 2 | `FETCH_PARSE_CHUNK` | Fetch NCBI gene packages and normalize target/ortholog records. | Intermediate chunk evidence |
| 3 | `BUILD_FETCH_DATASET`, `FETCH_TAXONOMY` | Assemble sequence/metadata handoffs and fetch canonical taxonomy once. | Intermediate fetch dataset |
| 4 | `FINALIZE_FETCH_OUTPUT` | Validate and publish the compact durable fetch layer. | `fetch/` |
| 5 | `BUILD_ALIGNMENT_TASKS` | Create stable per-gene tasks and genomic partitions. | Intermediate task metadata |
| 6 | alignment strategy processes | Run selected minimap2, nucmer, BWA, and opt-in Ensembl evidence producers. | Per-gene normalized evidence |
| 7 | `MERGE_ALIGNMENT_PARTITION` | Compact raw observations into event, exact-support, segment, and summary relations. | Intermediate bounded partitions |
| 8 | `MERGE_ALIGNMENT` | Validate and copy partitions without global recompression or ID rebasing. | `alignment/` |
| 9 | `PREPARE_ANNOTATION_CONTEXTS` | Materialize only each partition's target context. | Intermediate annotation context |
| 10 | `ANNOTATE_EVENTS_PARTITION` | Normalize event keys and annotate unique variant contexts with ClinVar/gnomAD. | Intermediate source shards |
| 11 | `ANNOTATE_VEP_PARTITION` | Add release-declared VEP evidence to bounded source shards. | Resumable enriched shards |
| 12 | `FINALIZE_ANNOTATION` | Validate and publish the partitioned variant dataset plus event-to-variant lineage. | `annotation/` |

## Documentation

- `docs/pipeline_launch.md` — ordinary ITMO launch or resume.
- `docs/report_generation.md` — report-only and combined pipeline/report launch.
- `docs/run_validation.md` — smoke tests, contract checks, and failure diagnosis.
- `docs/stage1_fetch_contract.md` and `docs/stage2_alignment_contract.md` —
  normalized data contracts and their rationale.
- `docs/stage3_annotation_contract.md` — annotation ownership, VEP, and durable
  variant-shard contract.
- `docs/storage_model.md` — durable evidence, resume cache, and disk policy.
- `docs/itmo_cluster.md` — first-time cluster setup and verified infrastructure.
- `docs/project_map.md` — code ownership and repository navigation.
