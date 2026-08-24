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
They live in an external analytics workspace and cannot replace or cause
deletion of their durable row-level inputs.

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
  reports/nextflow/
```

Once `run_manifest.json` is finalized successfully, the whole source run is
immutable. Analytics and operational tools may read it, but must not create,
modify, or delete any path below it. A change to pipeline evidence requires a
new run directory. This keeps archived source data stable and prevents report
caches from being backed up as if they were biological evidence.

The pipeline enforces this boundary: an existing completed manifest makes the
launch fail before stage work begins, even with `-resume`. An existing
`running` or `failed` manifest can be continued only with `-resume`; a fresh
launch must use a new output directory.

`run_manifest.json` is written when the workflow starts and atomically finalized
when Nextflow terminates. It records the launch-time Git commit and dirty state,
selected profiles, schema-declared launch parameters with secret-like values
redacted, and the workflow completion status. A manifest left in `running`
state identifies an interrupted run whose completion handler did not execute.
Archival tools must treat this file as the provenance source of truth.

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

The report derives and fingerprint-caches strategy, feature-coverage, support,
and taxonomy summaries under an external analytics root; the pipeline does not
publish duplicate global copies.

`annotation/` contains:

- the partitioned unique variant-context dataset with ClinVar, gnomAD, and VEP
  evidence
- a partitioned event-to-canonical-variant lineage table
- annotation manifest and diagnostic failure table

The variant-context shards intentionally store one compact interpretation
layer: canonical key, gene/event, alleles, external lookup status, strategy
membership, ClinVar classification/review fields, selected gnomAD fields, and
declared VEP consequence fields. Per-event normalization status stays in the
event map.

Support metrics remain derivable from the alignment evidence and event map
rather than being copied into this table. Provider fields not used by the
report are not duplicated into durable output.

The fetch, alignment, and annotation directories together form the durable
pipeline evidence layer and should be kept as one run.

When `GAPH_GNOMAD_CACHE_DIR` or `--gnomad_cache_dir` is set, complete gnomAD
regional responses are also stored in a shared reusable cache. The cache uses
25-kb gzip-compressed tiles under a dataset/reference/schema namespace. It is
neither a run result nor Nextflow resume state: multiple pipeline and analytics
runs may reuse it, and a run remains valid when the cache is absent.

`GAPH_VEP_RESULT_CACHE_DIR` is a second reusable infrastructure cache. It keeps
only complete release/config-matched variant/gene results and is distinct from
the official indexed VEP reference cache. Neither cache replaces the durable
enriched shards inside a completed run.

Alignment compacts raw per-aligner event observations once within each bounded
partition. The normalized evidence is copied byte for byte into
`alignment/evidence/partitions/<partition_id>/`: the partition manifest,
per-ortholog summary, segments, compact events, and exact ortholog support.
Local `event_group_id` values are interpreted together with their directory
`partition_id`; they are not globally rebased.

Annotation reads those same partitions in the end-to-end dataflow. It writes
one canonical partitioned variant-context dataset and preserves the
partition-local event-to-variant mapping needed by analytics. It does not materialize
variant-strategy, variant-ortholog, site-depth, or taxonomic histogram products.
In a fresh successful run, disposable inputs and pre-normalization
intermediates removed with task work include:

- selected ortholog FASTA files
- raw per-aligner event tables
- native aligner files

Pre-VEP annotation rows are bounded temporary shards. Completed VEP shards are
copied once into `annotation/variant_annotations/partitions/<partition_id>/`
without a global merge or recompression. Event maps stay partitioned, so their
local IDs need no global rewrite.

Selected ortholog FASTA files pass directly from fetch to alignment through the
Nextflow execution cache and are not copied into durable results. A completed
run keeps selected-ortholog metadata and canonical taxonomy. Human target FASTA
remains because annotation and reproducible target-context normalization need
it. Alignment performs no taxonomy lookup and does not receive taxonomy as an
input.

Analytics owns a separate reproducible derived layer:

```text
<analytics-root>/
  cache/<source-id>/
  analyses/<analysis-id>/
  slurm/
```

Per-source strategy, feature-coverage, taxonomy, and support derivations live
under `cache/<source-id>/`. Run-set scientific outputs, report HTML, and
performance profiles live under `analyses/<analysis-id>/`. Large source tables
are queried in place; analytics does not create a synthetic combined run.
Every cache records source identities and is rebuilt when they change. A
missing canonical input is an error; no old aggregate is accepted as a
substitute.

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
- smaller `--chunk_size` when fetch packages dominate
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

Lower `--fetch_max_forks` if the target cluster, network, scratch filesystem, or
NCBI behavior becomes unstable. Increase only after measuring on the target
cluster. Request starts are always spaced by 5 seconds inside the fetch
implementation.

## What To Keep

Keep:
- `<run>/run_manifest.json`
- `<run>/fetch/`, `<run>/alignment/`, and `<run>/annotation/`
- `<run>/reports/nextflow/` for execution trace and resource review
- the separate analytics root when derived caches or reports should be retained

The external analytics root can be deleted and rebuilt from durable pipeline
evidence. It is not part of a source-run archive. Do not delete pipeline
evidence in its place.

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

For debugging, inspect the task directory of a retained failed or interrupted
run before removing its Nextflow work cache.
