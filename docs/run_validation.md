# Run And Validation

Use this when running or validating the current end-to-end fetch + alignment + annotation
workflow.

For an ordinary ITMO launch or resume, start with the shorter
`docs/pipeline_launch.md`. Use this document for smoke tests, contract checks,
and failure investigation.

## Local Run

```bash
nextflow run . \
  --ids_file assets/inputs/gene_ids/smoke_5_genes.txt \
  --outdir results/run_001 \
  -resume
```

The pipeline always runs fetch, alignment, and annotation in order. Every run uses the declared
`envs/*.yml` task environments through Micromamba, including local runs without
a profile.
Annotation fetches gnomAD data from the live API for clustered event regions and
uses an in-memory lookup bounded to each genomic partition. Set
`GAPH_GNOMAD_CACHE_DIR` or `--gnomad_cache_dir` to reuse complete 25-kb regional
responses across runs and analytics reports. When `GAPH_ROOT` is set, the cache
defaults to `$GAPH_ROOT/cache/gnomad`. End-to-end runs may process up to
`--annotation_max_forks` partitions concurrently (default: 4). Initial memory
is selected from the partition's exact ortholog-support row count; retries add
32 GB per attempt.

For persistent local paths, environment variables keep machine-specific values
out of the repository:

```bash
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

The NCBI Datasets CLI and aligners come from their declared task environments.
`FETCH_PARSE_CHUNK` also loads an ignored project `.env` file when present. Use
`.env.example` as the template for `ENTREZ_EMAIL` and `ENTREZ_API_KEY`.

## Cluster Run

The `slurm` profile changes the executor and scratch policy; the declared
`envs/*.yml` environments are already mandatory globally. Every Python process
has an explicit environment, so compute nodes do not depend on system Python.

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
export GAPH_GNOMAD_CACHE_DIR="$GAPH_ROOT/cache/gnomad"
export NXF_CONDA_CACHEDIR="$GAPH_ROOT/conda/envs"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export NXF_HOME="$GAPH_ROOT/nextflow"

RUN="$GAPH_ROOT/results/run_001"
nextflow run . \
  -profile slurm \
  --ids_file /path/to/gene_ids.txt \
  --gnomad_cache_dir "$GAPH_GNOMAD_CACHE_DIR" \
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
--chunk_size 10 --fetch_max_forks 2 --alignment_max_forks 4 --annotation_max_forks 4
```

`fetch_max_forks` controls local fetch concurrency. The fetch implementation
always spaces request starts by 5 seconds and gives each download up to 4
retries after its initial attempt, with exponential backoff starting at 30
seconds. Tune concurrency only after measuring disk, runtime, NCBI behavior,
and alignment task memory on the target cluster.

In the `slurm` profile, `executor.queueSize` limits how many jobs Nextflow keeps
submitted to Slurm at once. It does not affect local runs, task CPU count, or
threads inside an aligner process.

The default storage policy retains process cache for recovery from a failed or
interrupted run and cleans task work created by a successful execution session.
Resume requires the original work directory and Nextflow execution metadata. A
fresh successful run has no task cache to reuse; a resumed run can retain task
directories from its earlier failed session, which can be removed after recovery
when the work path is dedicated to that run.

## End-To-End Smoke Test

Small multi-chunk smoke test:

```bash
printf '59067\n12\n59067\n355\n' > /tmp/gaph_v2_ids.txt

nextflow run . \
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

Fetch-boundary expected properties:
- `fetch/manifest.json` exists.
- `input_record_count` is 4.
- `unique_gene_count` is 3.
- `chunk_count` is 3.
- `target_gene_count` is 3.
- `target_feature_count` is greater than `target_gene_count`.
- `failure_count` is 0, unless NCBI data changed or the request failed.
- `taxonomy.tsv.gz` has one row per unique selected ortholog `tax_id`, explicit
  `taxonomy_status`, direct domain-through-species ID/name columns, and ordered
  `lineage_tax_ids`; current files do not contain derived `is_*` flags.
- `taxonomy_failures.tsv.gz` exists; `taxonomy_summary.tsv.gz` does not, because
  it is derived by analytics.

Alignment-boundary expected properties:
- `alignment/manifest.json` exists.
- The manifest declares `stage=alignment` and
  `schema=normalized_alignment_evidence_v1`.
- `alignment/manifest.json` `gene_count` equals the length of `gene_ids`, and
  those IDs equal the union of genes eligible for the selected strategies in
  the Stage 1 selection. Ensembl requires a target sequence; the other
  strategies require selected ortholog sequences.
- The Nextflow trace contains exactly one `FETCH_TAXONOMY` task owned by fetch;
  alignment neither receives nor republishes taxonomy.
- `alignment/evidence/partitions/<partition_id>/` exists for every partition in
  the manifest. Each contains exactly its manifest, per-ortholog summary,
  segments, compact events, and exact event-ortholog support.
- Partition manifests declare
  `schema=normalized_alignment_evidence_partition_v1`; their gene sets are
  non-empty and disjoint.
- `event_group_id` is consecutive and partition-local, and support foreign keys
  agree with compact-event counts.
- Global alignment event, segment, summary, feature-coverage, site-depth,
  taxonomy, task, and native-artifact files are absent.

The end-to-end run must publish this Stage 2 contract. Also verify that
`fetch/orthologs.selected.tsv.gz`, canonical taxonomy, target features, and
target FASTA remain available for analytics and annotation, while
`fetch/sequences/orthologs/` is absent from the completed end-to-end output.

Annotation expected properties:
- `annotation/variant_annotations.tsv.gz` exists for end-to-end runs.
- `annotation/event_variant_map/partitions/<partition_id>/event_variant_map.tsv.gz`
  has one row per compact Stage 2 event. Its `event_group_id` values are
  consecutive within the partition; non-concrete alleles have an empty
  `variant_key` and retain their normalization status.
- `annotation/manifest.json` declares `stage=annotation` and
  `schema=normalized_annotation_evidence_v1`.
- Pipeline-owned `variant_strategy_support`, `variant_ortholog_support`, and
  `ortholog_evidence_summary` outputs are absent.
- `annotation/manifest.json` records event and unique variant-context row counts,
  source metadata, useful provider/cache diagnostics, and per-partition phase
  timings.
- `variant_annotations.tsv.gz` is assembled as concatenated gzip members; the
  event map stays partitioned. Neither path globally rewrites source rows.
- `annotation/failures.tsv.gz` records non-fatal external lookup failures.
- `clinvar_*` columns are present and populated when matching ClinVar records exist;
  `clinvar_review_stars` is derived from the raw `clinvar_revstat` value.
- `gnomad_*` columns are populated only for variants found in fetched gnomAD regions,
  including the selected AF, its source dataset, and consequence annotation.

Analytics expected properties after the first report/preflight build:

- `analytics/alignment_aggregates/` contains strategy summary and feature
  coverage derived from partition summaries and segments.
- `analytics/taxonomy_summary/` contains the taxonomy summary derived from
  `fetch/taxonomy.tsv.gz` and `fetch/orthologs.selected.tsv.gz`.
- `analytics/annotation_support/` contains variant-strategy support and
  ortholog-evidence histograms derived from Stage 2 evidence and the Stage 3
  event map.
- Repeating the build without changing inputs is a cache hit. Removing or
  changing any required canonical input fails or invalidates the cache; no old
  pipeline aggregate is accepted instead.

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

Count durable alignment evidence rows:

```bash
python3 - <<'PY'
import csv, gzip
from pathlib import Path
base = Path("results/run_001/alignment/evidence/partitions")
for partition in sorted(path for path in base.iterdir() if path.is_dir()):
    for name in [
        "ortholog_alignment_summary.tsv.gz",
        "alignment_segments.tsv.gz",
        "alignment_events.tsv.gz",
        "event_ortholog_support.tsv.gz",
    ]:
        with gzip.open(partition / name, "rt", newline="") as handle:
            count = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))
        print(partition.name, name, count)
PY
```

## Resume A Failed Boundary

Re-run the same end-to-end command with the same `-work-dir` and `-resume`.
Nextflow restores every completed process whose inputs and code are unchanged;
for example, an annotation failure resumes without re-running successful fetch
or alignment tasks. There are no standalone stage input directories or
stage-selection parameters.

## After Validation

After a successful run is validated, clean temporary Nextflow work files if disk
space matters:

```bash
nextflow clean -f -keep-logs
```

If a custom `-work-dir` was used and no resume is needed, it can also be removed
manually after validation.
