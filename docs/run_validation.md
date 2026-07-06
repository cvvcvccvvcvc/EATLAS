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
uses an in-memory cache within that single run.

If command-line tools are not on `PATH`, pass them explicitly:

```bash
nextflow run . \
  -profile local,conda \
  --ids_file assets/inputs/gene_ids/smoke_5_genes.txt \
  --outdir results/run_001 \
  --datasets_bin /path/to/datasets \
  --minimap2_bin /path/to/minimap2 \
  --nucmer_bin /path/to/nucmer \
  --show_coords_bin /path/to/show-coords \
  --show_snps_bin /path/to/show-snps \
  -resume
```

For persistent local paths, environment variables keep machine-specific values
out of the repository:

```bash
export DATASETS_BIN=/path/to/datasets
export GAPH_TARGET_ANNOTATION_GFF3=/path/to/genomic.gff.gz
export CLINVAR_VCF=/path/to/clinvar.vcf.gz
export GAPH_WORK_DIR=/path/to/scratch/gaph_v2_work
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

Use the same workflow and put `work/` on scratch storage:

```bash
nextflow run . \
  -profile slurm,conda \
  --ids_file /path/to/gene_ids.txt \
  --outdir /scratch/$USER/gaph_v2/results/run_001 \
  -work-dir /scratch/$USER/gaph_v2/work/run_001 \
  -resume
```

Conservative starting parameters:

```bash
--chunk_size 10 --fetch_max_forks 2 --fetch_request_stagger_seconds 5 --alignment_max_forks 2
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

Use `-profile low_storage` only when preserving the published outputs matters
more than preserving the execution cache. It disables process caching and cleans
the work directory after a successful run, so do not expect `-resume` to reuse
completed tasks from that run.

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

Fetch expected properties:
- `fetch/manifest.json` exists.
- `input_record_count` is 4.
- `unique_gene_count` is 3.
- `chunk_count` is 3.
- `chunk_metrics.tsv.gz` exists and includes per-chunk download/parse timings.
- `target_gene_count` is 3.
- `target_feature_count` is greater than `target_gene_count`.
- `failure_count` is 0, unless NCBI data changed or the request failed.
- `orthologs.candidates.tsv.gz` has no rows with `tax_id=9606`.

Alignment expected properties:
- `alignment/manifest.json` exists.
- `alignment/taxonomy_presets.tsv.gz` has one row per unique ortholog tax_id.
- `alignment/ortholog_alignment_summary.tsv.gz` has rows for each enabled strategy.
- `alignment/alignment_segments.tsv.gz` and `alignment/alignment_events.tsv.gz`
  are gzip TSV files with stable headers.
- `alignment/feature_coverage.tsv.gz` has rows for each aligned gene, enabled
  strategy, and target structural feature.
- `alignment/native/` is absent unless `--keep_native_alignments true` was used.

Annotation expected properties:
- `annotation/alignment_events_annotated.tsv.gz` exists for end-to-end runs.
- `annotation/manifest.json` records row counts, source metadata, and annotation counters.
- `annotation/failures.tsv.gz` records non-fatal external lookup failures.
- `clinvar_*` columns are present and populated when matching ClinVar records exist.
- `gnomad_*` columns are populated only for events found in fetched gnomAD regions.

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
    "target_features.tsv.gz",
    "orthologs.selected.tsv.gz",
    "orthologs.candidates.tsv.gz",
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
    "taxonomy_presets.tsv.gz",
    "ortholog_alignment_summary.tsv.gz",
    "alignment_segments.tsv.gz",
    "feature_coverage.tsv.gz",
    "alignment_events.tsv.gz",
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
