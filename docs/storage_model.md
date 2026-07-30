# Storage Model

GAPH v2 separates durable biological outputs from temporary execution cache.

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

For a default end-to-end `--stage all` run, `fetch/` contains:

- compressed target FASTA files
- input and target-gene metadata
- compact target feature intervals
- the compact selected-ortholog provenance table
- fetch failures
- `manifest.json`

`alignment/` contains:

- compact per-strategy and feature-coverage summaries
- compact run-level taxonomy scope/unit summary
- alignment failures
- `manifest.json`

`annotation/` contains:

- the compressed unique variant-context annotation table
- compact per-strategy ALT-support counts for every normalized variant
- compact taxonomic ortholog-evidence histograms for report heatmaps
- annotation manifest and diagnostic failure table

The variant-context table intentionally stores one compact interpretation
layer: canonical key, gene/event, raw alleles, normalization status, aggregate
support, strategy membership, ClinVar classification/review fields, and the
selected gnomAD AF/source/consequence. Provider fields not used by the report
are not duplicated into durable output.

This layer should be kept.

When `GAPH_GNOMAD_CACHE_DIR` or `--gnomad_cache_dir` is set, complete gnomAD
regional responses are also stored in a shared reusable cache. The cache uses
25-kb gzip-compressed tiles under a dataset/reference/schema namespace. It is
neither a run result nor Nextflow resume state: multiple pipeline and analytics
runs may reuse it, and a run remains valid when the cache is absent.

Large handoff artifacts remain inside Nextflow `work/` during `--stage all`.
Alignment partitions reduce raw segments to a compact `snv_site_depth.tsv.gz`
table and temporary taxonomy-aware counts containing only observed concrete SNV
positions before annotation. Annotation reduces the latter to a bounded
histogram before publication. The large artifacts are removed after a successful
`low_storage` run:

- selected ortholog FASTA files
- raw alignment events
- alignment segments
- per-ortholog alignment summaries

Standalone `--stage fetch` and `--stage align` runs still publish their full
handoff datasets because a later invocation needs those files.

## Execution Cache

Nextflow `work/` is a resume cache. By default this repository puts it under
`$GAPH_WORK_DIR` when set, otherwise under the system temporary directory. It can
contain:

- task scripts and logs
- per-chunk intermediate outputs
- temporary task directories

It is useful while developing or recovering a failed run with `-resume`, but it
is not the data product.

The `low_storage` profile keeps the default process cache while a run is active
or failed, enables Nextflow successful-run cleanup, and moves terminal
annotation outputs into the published annotation directory instead of keeping a
second copy in `work/`. A failed or interrupted run can use `-resume` while its
work directory and Nextflow execution metadata remain. After a successful run,
cleanup removes the task work, so that run is no longer reusable with `-resume`.

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
- `-profile low_storage` when resume is required only until successful completion

On Slurm, both paths must be visible from the controller and every compute node.
The Conda cache is reusable infrastructure and should remain outside individual
run directories; `low_storage` cleanup applies to task work, not that cache.
The repository's `slurm` profile disables task-local scratch so staged task data
stays under the assigned shared work allocation. Local execution retains
task-local scratch for fetch and alignment processes.

Recommended starting point for large runs:

```bash
--chunk_size 10 --fetch_max_forks 2 --fetch_request_stagger_seconds 5
```

Lower `--fetch_max_forks` if the target cluster, network, scratch filesystem, or
NCBI behavior becomes unstable. Increase only after measuring on the target
cluster.

## What To Keep

Keep:
- `results/.../run_manifest.json`
- `results/.../manifest.json`
- `results/.../*.tsv.gz`
- `results/.../sequences/targets/*.fa.gz`
- `results/.../target_features.tsv.gz`

Keep standalone Stage 1 output until its Stage 2 consumer has finished:
- `results/.../sequences/orthologs/*.fa.gz`

Usually remove after validation:
- Nextflow `work/`
- `.nextflow*` local execution files in throwaway run directories
- raw NCBI package files if any were produced manually

Rejected ortholog candidate sequences are not retained. Rejected candidates are
represented only by metadata rows in `orthologs.candidates.tsv.gz`.

Raw aligner outputs are also not retained by default:

- minimap2 `.paf`
- nucmer `.sam`
- Ensembl Compara MAF chunks used by precomputed alignment strategies

Set `--keep_native_alignments true` only for targeted debug/benchmark runs.
