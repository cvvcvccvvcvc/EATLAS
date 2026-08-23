# Storage Model

GAPH v2 separates durable biological outputs from temporary execution cache.

## Evidence-First Invariant

Durable pipeline data is normalized evidence, not necessarily a raw provider
dump. Deterministic parsing, coordinate normalization, lossless compaction,
deduplication, partitioning, and compression are allowed when stable identifiers
and exact biological observations remain available.

Report-specific aggregation is not a storage optimization. Taxonomic scope and
rank selections, scientific counters, thresholds, bins, histograms, and plotting
tables are analytics artifacts and must be reproducible after `work/` is deleted.
They may be cached under `analytics/`, but cannot replace or cause deletion of
their durable row-level inputs.

Manifests may repeat inexpensive counts or checksums to detect incomplete or
corrupted output. These values are integrity snapshots, not the canonical source
for scientific analysis. Non-reconstructable outcomes and provenance, such as a
failed lookup, a valid `no_alignment` result, tool versions, and run parameters,
remain part of the pipeline contract even when they are compact.

## Durable Layer

For the default end-to-end run, `params.outdir` is a run root:

```text
results/run_001/
  run_manifest.json
  fetch/
  alignment/
  annotation/
```

`run_manifest.json` is written when the workflow starts and atomically finalized
when Nextflow terminates. It records the launch-time Git commit and dirty state,
selected profiles, schema-declared resolved parameters with secret-like values
redacted, and the workflow completion status. A manifest left in `running`
state identifies an interrupted run whose completion handler did not execute.
Archival tools must treat this file as the provenance source of truth. For
legacy runs without it, the producing commit is unknown and must not be inferred
from the current checkout.

For an end-to-end run, `fetch/` contains:

- compressed target FASTA files
- input and target-gene metadata
- compact target feature intervals
- the compact selected-ortholog provenance table
- one canonical gzip wide taxonomy row per selected ortholog tax ID, its lookup
  failures; no duplicate summary, long classification table, or raw NCBI
  response is published
- fetch failures
- `manifest.json`

`alignment/` contains:

- partitioned normalized Stage 2 source evidence under
  `evidence/partitions/<partition_id>/`
- alignment failures
- `manifest.json`

The report derives and fingerprint-caches strategy, feature-coverage, and
taxonomy summaries under `analytics/`; alignment does not publish duplicate
global copies.

`annotation/` contains:

- the compressed unique variant-context annotation table
- a partitioned event-to-canonical-variant lineage table
- annotation manifest and diagnostic failure table

The variant-context table intentionally stores one compact interpretation
layer: canonical key, gene/event, raw alleles, normalization status, strategy
membership, ClinVar classification/review fields, and the selected gnomAD
AF/source/consequence. Support metrics remain derivable from the alignment
evidence and event map rather than being copied into this table. Provider fields
not used by the report are not duplicated into durable output.

This layer should be kept.

When `GAPH_GNOMAD_CACHE_DIR` or `--gnomad_cache_dir` is set, complete gnomAD
regional responses are also stored in a shared reusable cache. The cache uses
25-kb gzip-compressed tiles under a dataset/reference/schema namespace. It is
neither a run result nor Nextflow resume state: multiple pipeline and analytics
runs may reuse it, and a run remains valid when the cache is absent.

Alignment compacts raw per-aligner event observations once within each bounded
partition. The normalized evidence is copied byte for byte into
`alignment/evidence/partitions/<partition_id>/`: the partition manifest,
per-ortholog summary, segments, compact events, and exact ortholog support.
Local `event_group_id` values are interpreted together with their directory
`partition_id`; they are not globally rebased.

Annotation reads those same partitions in the end-to-end dataflow. It writes
one canonical variant-context annotation table and preserves the partition-local
event-to-variant mapping needed by analytics. It does not materialize
variant-strategy, variant-ortholog, site-depth, or taxonomic histogram products.
In a fresh successful run, disposable inputs and pre-normalization
intermediates removed with task work include:

- selected ortholog FASTA files
- raw per-aligner event tables
- native aligner files

The durable `variant_annotations.tsv.gz` is assembled from partition gzip
members without row parsing or recompression. Event maps stay partitioned, so
their local IDs need no global rewrite.

Selected ortholog FASTA files pass directly from fetch to alignment through the
Nextflow execution cache and are not copied into durable results. A completed
run keeps selected-ortholog metadata and canonical taxonomy. Human target FASTA
remains because annotation and reproducible target-context normalization need
it. Alignment performs no taxonomy lookup and does not receive taxonomy as an
input.

Analytics owns the reproducible derived layer:

- `analytics/alignment_aggregates/` — strategy summary and feature coverage;
- `analytics/taxonomy_summary/` — scope/rank summary from selected orthologs and
  canonical taxonomy;
- `analytics/annotation_support/` — variant-strategy support and taxonomic
  ortholog-evidence histograms.

Each cache records the identities of its source inputs and is rebuilt when they
change. A missing canonical input is an error; no old aggregate is accepted as
a substitute.

## Execution Cache

Nextflow `work/` is a resume cache. By default this repository puts it under
`$GAPH_WORK_DIR` when set, otherwise under the system temporary directory. It can
contain:

- task scripts and logs
- per-chunk intermediate outputs
- temporary task directories

It is useful while developing or recovering a failed run with `-resume`, but it
is not the data product.

The default storage policy keeps the process cache while a run is active or
failed, enables Nextflow successful-run cleanup, and moves terminal end-to-end
fetch and annotation outputs into their published directories instead of
keeping a second copy in `work/`. A failed or interrupted run can use `-resume`
while its work directory and Nextflow execution metadata remain. Cleanup removes
task work created by a successful execution session, so a fresh completed run is
no longer reusable with `-resume`. A resumed run can retain task directories
from its earlier failed session; remove those explicitly after recovery when the
work path is dedicated to that run.

Alignment task directories are metadata-only. They do not duplicate Stage 1
target or ortholog FASTA files. Sequence-based aligner processes receive the
needed per-gene Stage 1 FASTA files as explicit Nextflow inputs and materialize
uncompressed FASTA files inside their own local task workspace.

## Raw NCBI Data

Raw NCBI zip packages and unpacked `gene.fna` files are intentionally not
published. They are temporary task files inside `FETCH_PARSE_CHUNK`.

Reason:
- raw packages are large
- they duplicate normalized outputs
- they are inconvenient downstream inputs
- they increase disk pressure on large gene lists

## Disk Policy

During a run, peak disk usage can exceed final result size because Nextflow keeps
task outputs in `work/` for resume.

Control peak disk with:
- smaller `--chunk_size`
- lower `--fetch_max_forks`
- `-work-dir` on scratch storage
- `GAPH_WORK_DIR=/path/to/scratch/gaph_v2_work`
- `NXF_CONDA_CACHEDIR=/path/to/shared/scratch/gaph_v2_conda`

On Slurm, both paths must be visible from the controller and every compute node.
The Conda cache is reusable infrastructure and should remain outside individual
run directories; successful-run cleanup applies to task work, not that cache.
The repository's `slurm` profile disables task-local scratch so staged task data
stays under the assigned shared work allocation. Local execution retains
task-local scratch for fetch and alignment processes.

Recommended starting point for large runs:

```bash
--chunk_size 10 --fetch_max_forks 2
```

Lower `--fetch_max_forks` if the target cluster, network, scratch filesystem, or
NCBI behavior becomes unstable. Increase only after measuring on the target
cluster. Request starts are always spaced by 5 seconds inside the fetch
implementation.

## What To Keep

Keep:
- `results/.../run_manifest.json`
- `results/.../manifest.json`
- `results/.../*.tsv.gz`
- `results/.../sequences/targets/*.fa.gz`
- `results/.../target_features.tsv.gz`

Usually remove after validation:
- Nextflow `work/`
- `.nextflow*` local execution files in throwaway run directories
- raw NCBI package files if any were produced manually

Rejected ortholog candidate sequences and candidate metadata remain disposable
fetch intermediates; selected-ortholog provenance is retained in
`fetch/orthologs.selected.tsv.gz`.

Raw aligner outputs are not retained as durable results:

- minimap2 `.paf`
- nucmer `.sam`
- Ensembl Compara MAF chunks used by precomputed alignment strategies

For debugging, inspect the task directory of a retained failed or interrupted
run before removing its Nextflow work cache.
