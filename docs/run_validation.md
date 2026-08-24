# Run And Validation

Use this when validating the current end-to-end fetch + alignment +
ClinVar/gnomAD/VEP annotation workflow.

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
is selected from the partition's compact alignment-event count; retries add
32 GB per attempt.

For persistent local paths, environment variables keep machine-specific values
out of the repository:

```bash
export GAPH_TARGET_ANNOTATION_GFF3=/path/to/genomic.gff.gz
export CLINVAR_VCF=/path/to/clinvar.vcf.gz
export GAPH_WORK_DIR=/path/to/scratch/gaph_v2_work
export GAPH_GNOMAD_CACHE_DIR=/path/to/scratch/gaph_v2_cache/gnomad
export GAPH_VEP_BACKEND=local
export GAPH_VEP_RELEASE=116
export GAPH_VEP_EXECUTABLE=/path/to/gaph-vep116
export GAPH_VEP_CACHE_DIR=/path/to/vep-cache
export GAPH_VEP_RESULT_CACHE_DIR=/path/to/shared/vep-results
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

Use `docs/pipeline_launch.md` for an ordinary launch or resume. Use
`docs/itmo_cluster.md` only for first-time environment/reference setup,
compute-node preflight, or infrastructure diagnosis.

The `slurm` profile changes the executor and disables task-local scratch; task
environments remain mandatory. Keep the work directory, run results, task
environment cache, and shared gnomAD cache under the assigned shared scratch
allocation. Do not run pipeline computation directly on the controller host.

Tune concurrency only after a representative run. `fetch_max_forks` controls
NCBI request concurrency, while each alignment process has its own
`alignment_max_forks` limit and `executor.queueSize` bounds submitted Slurm
jobs. These settings are not one global worker count.

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
  run_manifest.json
  fetch/
  alignment/
  annotation/
  reports/nextflow/
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
  `schema=normalized_alignment_evidence_v2`.
- `alignment/manifest.json` `gene_count` equals the length of `gene_ids`, and
  those IDs equal the union of genes eligible for the selected strategies in
  the Stage 1 selection. Every alignment strategy requires selected ortholog
  sequences.
- The Nextflow trace contains exactly one `FETCH_TAXONOMY` task owned by fetch;
  alignment neither receives nor republishes taxonomy.
- `alignment/evidence/partitions/<partition_id>/` exists for every partition in
  the manifest. Each contains exactly its manifest, per-ortholog summary,
  segments, compact events, and exact event-ortholog support.
- Partition manifests declare
  `schema=normalized_alignment_evidence_partition_v2`; their gene sets are
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
- `annotation/variant_annotations/manifest.json` exists and declares
  `schema=gaph_variant_annotation_dataset_v1`, a complete partitioned headered
  gzip dataset, exact fields/counts, VEP configuration, and VEP status counts.
- Every declared
  `annotation/variant_annotations/partitions/<partition_id>/<shard_id>.tsv.gz`
  exists, has the declared compressed size/header, and contains no more than
  250,000 rows.
- `annotation/event_variant_map/partitions/<partition_id>/event_variant_map.tsv.gz`
  has one row per compact Stage 2 event. Its `event_group_id` values are
  consecutive within the partition; non-concrete alleles have an empty
  `variant_key` and retain their normalization status.
- `annotation/manifest.json` declares `stage=annotation` and
  `schema=normalized_annotation_evidence_v3`; its `variant_annotations`
  descriptor exactly matches the child dataset manifest.
- Pipeline-owned `variant_strategy_support`, `variant_ortholog_support`, and
  `ortholog_evidence_summary` outputs are absent.
- `annotation/manifest.json` records event and unique variant-context row counts,
  source metadata, useful provider/cache diagnostics, and per-partition phase
  timings.
- VEP-enriched shards are copied to the durable dataset byte for byte; no
  duplicate global variant table is built. The event map also stays partitioned.
- `annotation/failures.tsv.gz` records non-fatal external lookup failures.
- `clinvar_*` columns are present and populated when matching ClinVar records exist;
  `clinvar_review_stars` is derived from the raw `clinvar_revstat` value.
- `gnomad_*` columns are populated only for variants found in fetched gnomAD regions,
  including the selected AF, its source dataset, and consequence annotation.
- `vep_*` columns and explicit per-row `vep_status` are present. The dataset
  release/backend match the requested or cluster-configured VEP contract.

Analytics expected properties after the first report/preflight build:

- `<analytics-root>/cache/<source-id>/alignment_aggregates/` contains strategy
  summary and feature coverage derived from partition summaries and segments.
- `<analytics-root>/cache/<source-id>/taxonomy_summary/` contains the taxonomy
  summary derived from `fetch/taxonomy.tsv.gz` and
  `fetch/orthologs.selected.tsv.gz`.
- `<analytics-root>/cache/<source-id>/annotation_support/` contains
  variant-strategy support and ortholog-evidence histograms derived from Stage 2
  evidence and the Stage 3 event map.
- `<analytics-root>/analyses/<analysis-id>/` contains only run-set derivations,
  report HTML, performance data, and its analysis manifest.
- The completed source run has exactly the same files and metadata before and
  after analytics.
- Repeating the build without changing inputs is a cache hit. Removing or
  changing any required canonical input fails or invalidates the cache; no old
  pipeline aggregate is accepted instead.

## Quick Checks

Inspect manifests:

```bash
jq . results/run_001/fetch/manifest.json
jq . results/run_001/alignment/manifest.json
jq . results/run_001/annotation/variant_annotations/manifest.json
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

## Speed And Disk Checks

Use a representative gene panel and compare runs with the same inputs,
strategies, references, cache state, and Slurm resources. A faster warm-cache
run is not comparable to a cold-cache run.

Nextflow already publishes the task-level measurements needed for pipeline
review:

```bash
column -t -s $'\t' results/run_001/reports/nextflow/trace.txt | less -S
du -sh results/run_001/{fetch,alignment,annotation,reports}
du -sh /path/to/dedicated/work/run_001
```

Review `realtime`, `%cpu`, `peak_rss`, `rchar`, and `wchar` in the
trace. Check the largest partitions rather than only totals. Durable alignment
must contain one partition tree, not duplicate global copies of summaries,
segments, events, or support.

Every report writes a progressive profile under
`<analytics-root>/analyses/<analysis-id>/performance/<report-name>.json`.
Compare a cold report with an
immediate identical rerun: the second run should scan the same pipeline-owned
VEP shards and report cache hits for completed derived alignment, taxonomy, support, and
analysis artifacts. Investigate a regression only when it reproduces on the
same workload; do not add a new storage format or cache from an unmeasured
hypothesis.

## Resume A Failed Boundary

Re-run the same end-to-end command with the same `-work-dir` and `-resume`.
Nextflow restores every completed process whose inputs and code are unchanged;
for example, an annotation failure resumes without re-running successful fetch
or alignment tasks. There are no standalone stage input directories or
stage-selection parameters.

## After Validation

Successful sessions use the repository's automatic Nextflow cleanup policy.
Failed or interrupted task work is retained for `-resume`. Remove a dedicated
work directory only after the run is verified and recovery is no longer needed;
see `docs/storage_model.md` for the durable/cache boundary.
