# Run And Validation

Use this when running or validating the current end-to-end fetch + alignment
workflow.

## Local Run

```bash
nextflow run . \
  -profile local \
  --ids_file gene_ids.txt \
  --outdir results/run_001 \
  -resume
```

If command-line tools are not on `PATH`, pass them explicitly:

```bash
nextflow run . \
  -profile local \
  --ids_file gene_ids.txt \
  --outdir results/run_001 \
  --datasets_bin /path/to/datasets \
  --minimap2_bin /path/to/minimap2 \
  --nucmer_bin /path/to/nucmer \
  --show_coords_bin /path/to/show-coords \
  --show_snps_bin /path/to/show-snps \
  -resume
```

## Cluster Run

Use the same workflow and put `work/` on scratch storage:

```bash
nextflow run . \
  -profile slurm \
  --ids_file /path/to/gene_ids.txt \
  --outdir /scratch/$USER/gaph_v2/results/run_001 \
  -work-dir /scratch/$USER/gaph_v2/work/run_001 \
  -resume
```

Conservative starting parameters:

```bash
--chunk_size 10 --fetch_max_forks 1 --alignment_max_forks 2
```

Increase only after measuring disk, runtime, NCBI behavior, and alignment task
memory on the target cluster.

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

Expected layout:

```text
/tmp/gaph_v2_smoke_run/
  fetch/
  alignment/
```

Fetch expected properties:
- `fetch/manifest.json` exists.
- `input_record_count` is 4.
- `unique_gene_count` is 3.
- `chunk_count` is 3.
- `target_gene_count` is 3.
- `failure_count` is 0, unless NCBI data changed or the request failed.
- `orthologs.candidates.tsv.gz` has no rows with `tax_id=9606`.

Alignment expected properties:
- `alignment/manifest.json` exists.
- `alignment/taxonomy_presets.tsv.gz` has one row per unique ortholog tax_id.
- `alignment/ortholog_alignment_summary.tsv.gz` has rows for each enabled strategy.
- `alignment/alignment_segments.tsv.gz` and `alignment/alignment_events.tsv.gz`
  are gzip TSV files with stable headers.
- `alignment/native/` is absent unless `--keep_native_alignments true` was used.

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
for name in ["genes.tsv.gz", "orthologs.selected.tsv.gz", "orthologs.candidates.tsv.gz", "failures.tsv.gz"]:
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

## After Validation

After a successful run is validated, clean temporary Nextflow work files if disk
space matters:

```bash
nextflow clean -f -keep-logs
```

If a custom `-work-dir` was used and no resume is needed, it can also be removed
manually after validation.
