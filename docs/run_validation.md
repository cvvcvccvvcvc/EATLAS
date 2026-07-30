# Run And Validation

Use this when running or validating the current end-to-end fetch + alignment + annotation
workflow.

## Local Run

```bash
nextflow run . \
  -profile local,conda \
  --ids_file assets/inputs/gene_ids/smoke_5_genes.txt \
  --outdir results/run_001 \
  -resume
```

The default `--stage all` runs every stage. Use the `conda` profile for normal
local runs so tasks use `envs/*.yml` instead of the active shell environment.
Annotation fetches gnomAD data from the live API for clustered event regions and
uses an in-memory lookup bounded to each genomic partition. Set
`GAPH_GNOMAD_CACHE_DIR` or `--gnomad_cache_dir` to reuse complete 25-kb regional
responses across runs and analytics reports. End-to-end runs may process up to
`--annotation_max_forks` partitions concurrently (default: 2).

If command-line tools are not on `PATH`, pass them explicitly:

```bash
nextflow run . \
  -profile local,conda \
  --ids_file assets/inputs/gene_ids/smoke_5_genes.txt \
  --outdir results/run_001 \
  --datasets_bin /path/to/datasets \
  --minimap2_bin /path/to/minimap2 \
  --nucmer_bin /path/to/nucmer \
  -resume
```

For persistent local paths, environment variables keep machine-specific values
out of the repository:

```bash
export DATASETS_BIN=/path/to/datasets
export GAPH_TARGET_ANNOTATION_GFF3=/path/to/genomic.gff.gz
export CLINVAR_VCF=/path/to/clinvar.vcf.gz
export GAPH_WORK_DIR=/path/to/scratch/gaph_v2_work
export GAPH_GNOMAD_CACHE_DIR=/path/to/scratch/gaph_v2_cache/gnomad
```

When neither `--target_annotation_gff3` nor `GAPH_TARGET_ANNOTATION_GFF3` is
set, fetch uses
`assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz`.
When neither `--clinvar_vcf` nor `CLINVAR_VCF` is set, annotation uses
`assets/reference/clinvar/clinvar.vcf.gz` if present and indexed. ClinVar is
required for annotation, so the workflow fails early when no VCF and matching
`.tbi` are available.

The NCBI Datasets CLI is resolved as `DATASETS_BIN`, then
`tools/bin/datasets` when present, then `datasets` on `PATH`.
`FETCH_PARSE_CHUNK` also loads an ignored project `.env` file when present.
Use `.env.example` as the template for `ENTREZ_EMAIL` and `ENTREZ_API_KEY`.

## Cluster Run

The `slurm` profile enables the declared `envs/*.yml` environments through
Micromamba. Every Python process has an explicit environment; compute nodes do
not depend on their system Python.

For a new account, keep reusable environments and run data in the assigned
scratch area. From the repository, create the controller environment once:

```bash
export GAPH_ROOT="/mnt/tank/scratch/$USER/gaph_v2"
mkdir -p "$GAPH_ROOT"/{envs,work,conda,results,micromamba,nextflow}

export CONDA_PKGS_DIRS="$GAPH_ROOT/conda/controller-pkgs"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export NXF_HOME="$GAPH_ROOT/nextflow"

conda env create --prefix "$GAPH_ROOT/envs/controller" -f envs/controller.yml
conda activate "$GAPH_ROOT/envs/controller"
nextflow -version
java -version
micromamba --version
```

Both Nextflow work and its Conda cache must be on storage shared by the Slurm
controller and compute nodes. For the ITMO CT cluster, use the assigned scratch
directory rather than the home quota:

```bash
export GAPH_ROOT="/mnt/tank/scratch/$USER/gaph_v2"
export GAPH_WORK_DIR="$GAPH_ROOT/work"
export NXF_CONDA_CACHEDIR="$GAPH_ROOT/conda/envs"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export NXF_HOME="$GAPH_ROOT/nextflow"

RUN="$GAPH_ROOT/results/run_001"
nextflow run . \
  -profile slurm,low_storage \
  --ids_file /path/to/gene_ids.txt \
  --outdir "$RUN"
```

`GAPH_WORK_DIR` supplies the default Nextflow work path. An explicit
`-work-dir "$GAPH_ROOT/work/run_001"` is also valid and takes precedence.
`NXF_CONDA_CACHEDIR` is intentionally outside the run directory so environments
are built once and reused across runs. `MAMBA_ROOT_PREFIX` keeps downloaded
package and repodata caches out of the home quota. `NXF_HOME` does the same for
Nextflow runtime and plugin files.

The `slurm` profile disables the process-level `scratch` directive. Tasks run in
the shared Nextflow work directory under `GAPH_WORK_DIR` instead of staging
large inputs and outputs through compute-node `/tmp`. This keeps pipeline data
inside the assigned `/mnt/tank/scratch/$USER` allocation. Local profiles retain
their existing task-scratch behavior.

On the ITMO CT cluster, submit Nextflow from `sphinx`; do not run calculations
there directly. The documented `main` partition is the default, and the cluster
instructions do not require a Slurm account or QOS for ordinary jobs, so the
profile does not invent `account`, `queue`, or `clusterOptions` values. Add a QOS
only when the administrators explicitly grant and request one.

See `docs/itmo_cluster.md` for the verified host layout, environment bootstrap,
reference transfer, preflight checks, and staged smoke-test procedure.

Conservative starting parameters:

```bash
--chunk_size 10 --fetch_max_forks 2 --fetch_request_stagger_seconds 5 --alignment_max_forks 4
```

`fetch_max_forks` controls local fetch concurrency. `fetch_request_stagger_seconds`
spaces out the starts of NCBI Datasets download requests so concurrent tasks do
not all submit at the same instant. `fetch_download_retries` handles transient
NCBI stream resets inside one fetch task before Nextflow sees a failure. Tune
these only after measuring disk, runtime, NCBI behavior, and alignment task
memory on the target cluster.

In the `slurm` profile, `executor.queueSize` limits how many jobs Nextflow keeps
submitted to Slurm at once. It does not affect local runs, task CPU count, or
threads inside an aligner process.

Combine `slurm` with `low_storage` to retain process cache for recovery from a
failed or interrupted run and clean the work directory after the workflow
finishes successfully. Resume requires the original work directory and Nextflow
execution metadata; after successful cleanup there is no task cache to reuse.

## End-To-End Smoke Test

Small multi-chunk smoke test:

```bash
printf '59067\n12\n59067\n355\n' > /tmp/gaph_v2_ids.txt

nextflow run . \
  -profile local \
  --ids_file /tmp/gaph_v2_ids.txt \
  --outdir /tmp/gaph_v2_smoke_run \
  --chunk_size 1 \
  --fetch_max_forks 1 \
  --alignment_max_forks 1 \
  -work-dir /tmp/gaph_v2_smoke_work \
  -resume
```

Run only one alignment strategy when debugging strategy-specific failures:

```bash
nextflow run . \
  -profile local \
  --ids_file /tmp/gaph_v2_ids.txt \
  --outdir /tmp/gaph_v2_smoke_run_asm20 \
  --alignment_strategies minimap2_asm20 \
  --chunk_size 1 \
  --fetch_max_forks 1 \
  --alignment_max_forks 1 \
  -work-dir /tmp/gaph_v2_smoke_work_asm20 \
  -resume
```

Expected layout:

```text
/tmp/gaph_v2_smoke_run/
  fetch/
  alignment/
  annotation/
```

Standalone `--stage fetch` expected properties:
- `fetch/manifest.json` exists.
- `input_record_count` is 4.
- `unique_gene_count` is 3.
- `chunk_count` is 3.
- `chunk_metrics.tsv.gz` exists and includes per-chunk download/parse timings.
- `target_gene_count` is 3.
- `target_feature_count` is greater than `target_gene_count`.
- `failure_count` is 0, unless NCBI data changed or the request failed.
- `orthologs.candidates.tsv.gz` has no rows with `tax_id=9606`.

Standalone `--stage align` expected properties:
- `alignment/manifest.json` exists.
- `alignment/manifest.json` `gene_count` equals the length of `gene_ids`, and
  those IDs equal the union of genes eligible for the selected strategies in
  `alignment_tasks.tsv.gz`. Ensembl requires `target_ready=true`; the other
  strategies require `ortholog_ready=true`.
- `strategy_eligible_gene_counts` reports the corresponding denominator for
  every selected strategy.
- `alignment/taxonomy_presets.tsv.gz` has one row per unique ortholog tax_id.
- `alignment/taxonomy_summary.tsv.gz` records run-level scope and evidence-unit
  counts.
- `alignment/ortholog_alignment_summary.tsv.gz` has rows for each enabled strategy.
- `alignment/alignment_segments.tsv.gz` and `alignment/alignment_events.tsv.gz`
  are gzip TSV files with stable headers.
- `alignment/snv_site_depth.tsv.gz` contains one row per observed concrete SNV
  position and strategy, with positive distinct-ortholog depth.
- `alignment/feature_coverage.tsv.gz` has rows for each aligned gene, enabled
  strategy, and target structural feature.
- `alignment/native/` is absent unless `--keep_native_alignments true` was used.

For default end-to-end `--stage all`, fetch and alignment publish only the
analysis-ready subset documented in `docs/storage_model.md`. In that mode,
verify `orthologs.selected.tsv.gz`, target FASTA, `strategy_summary.tsv.gz`,
`feature_coverage.tsv.gz`, manifests, and failure tables instead of expecting
the standalone handoff files.

Annotation expected properties:
- `annotation/variant_annotations.tsv.gz` exists for end-to-end runs.
- `annotation/variant_strategy_support.tsv.gz` contains per-strategy ALT-support
  counts and site-aligned ortholog depth for SNVs.
- `annotation/ortholog_evidence_summary.tsv.gz` contains bounded SNV evidence
  histograms by strategy, target context, taxonomic scope, and evidence unit.
- `annotation/manifest.json` records event and unique variant-context row counts, source metadata, and annotation counters.
- End-to-end annotation records `partition_count`; partition outputs are merged
  by streaming and are not published as duplicate durable tables.
- `annotation/failures.tsv.gz` records non-fatal external lookup failures.
- `clinvar_*` columns are present and populated when matching ClinVar records exist;
  `clinvar_review_stars` is derived from the raw `clinvar_revstat` value.
- `gnomad_*` columns are populated only for variants found in fetched gnomAD regions,
  including the selected AF, its source dataset, and consequence annotation.

## Quick Checks

Inspect manifests:

```bash
jq . results/run_001/fetch/manifest.json
jq . results/run_001/alignment/manifest.json
```

Count fetch table rows:

```bash
python3 - <<'PY'
import csv, gzip
from pathlib import Path
base = Path("results/run_001/fetch")
for name in [
    "genes.tsv.gz",
    "orthologs.selected.tsv.gz",
    "failures.tsv.gz",
]:
    with gzip.open(base / name, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    print(name, len(rows))
PY
```

Count alignment table rows:

```bash
python3 - <<'PY'
import csv, gzip
from pathlib import Path
base = Path("results/run_001/alignment")
for name in [
    "strategy_summary.tsv.gz",
    "feature_coverage.tsv.gz",
    "failures.tsv.gz",
]:
    with gzip.open(base / name, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    print(name, len(rows))
PY
```

## Alignment-Only Debug

Use this only when Stage 1 output already exists and the alignment stage is being
debugged:

```bash
nextflow run . \
  -profile local \
  --stage align \
  --fetch_dir /path/to/fetch \
  --outdir /tmp/gaph_v2_align_debug \
  -work-dir /tmp/gaph_v2_align_debug_work \
  -resume
```

Strategy selection works in alignment-only mode as well:

```bash
nextflow run . \
  -profile local \
  --stage align \
  --fetch_dir /path/to/fetch \
  --outdir /tmp/gaph_v2_align_debug_asm20 \
  --alignment_strategies minimap2_asm20 \
  -work-dir /tmp/gaph_v2_align_debug_work_asm20 \
  -resume
```

## Annotation-Only Debug

Annotation-only mode reuses an existing alignment events table and the matching
fetch directory. The fetch directory is required because ClinVar and gnomAD
lookup normalizes events to VCF keys using `genes.tsv.gz` and
`sequences/targets/*.fa.gz`.

```bash
nextflow run . \
  -profile local \
  --stage annotate \
  --events_tsv /path/to/alignment_events.tsv.gz \
  --segments_tsv /path/to/alignment_segments.tsv.gz \
  --fetch_dir /path/to/fetch \
  --outdir /tmp/gaph_v2_annotate_debug \
  -work-dir /tmp/gaph_v2_annotate_debug_work \
  -resume
```

## After Validation

After a successful run is validated, clean temporary Nextflow work files if disk
space matters:

```bash
nextflow clean -f -keep-logs
```

If a custom `-work-dir` was used and no resume is needed, it can also be removed
manually after validation.
