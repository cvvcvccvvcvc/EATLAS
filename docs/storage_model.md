# Storage Model

GAPH v2 separates durable biological outputs from temporary execution cache.

## Durable Layer

For the default end-to-end run, `params.outdir` is a run root:

```text
results/run_001/
  fetch/
  alignment/
  annotation/
```

`fetch/` contains the normalized data layer for downstream pipeline stages:

- compressed target FASTA files
- compressed selected ortholog FASTA files
- compressed metadata tables
- `manifest.json`

`alignment/` contains:

- compressed alignment segments
- compressed raw alignment events
- per-ortholog alignment summaries
- compact taxonomy preset metadata
- `manifest.json`

`annotation/` contains:

- annotated event tables

This layer should be kept.

## Execution Cache

Nextflow `work/` is a resume cache. By default this repository puts it under
`$GAPH_WORK_DIR` when set, otherwise under the system temporary directory. It can
contain:

- task scripts and logs
- per-chunk intermediate outputs
- temporary task directories

It is useful while developing or recovering a failed run with `-resume`, but it
is not the data product.

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
- `scratch = true` for fetch/parse tasks

Recommended starting point for large runs:

```bash
--chunk_size 10 --fetch_max_forks 2 --fetch_request_stagger_seconds 5
```

Lower `--fetch_max_forks` if the target cluster, network, scratch filesystem, or
NCBI behavior becomes unstable. Increase only after measuring on the target
cluster.

## What To Keep

Keep:
- `results/.../manifest.json`
- `results/.../*.tsv.gz`
- `results/.../sequences/targets/*.fa.gz`
- `results/.../sequences/orthologs/*.fa.gz`

Usually remove after validation:
- Nextflow `work/`
- `.nextflow*` local execution files in throwaway run directories
- raw NCBI package files if any were produced manually

Rejected ortholog candidate sequences are not retained. Rejected candidates are
represented only by metadata rows in `orthologs.candidates.tsv.gz`.

Raw aligner outputs are also not retained by default:

- minimap2 `.paf`
- nucmer `.delta`
- `show-coords`
- `show-snps`
- Ensembl Compara MAF chunks used by precomputed alignment strategies

Set `--keep_native_alignments true` only for targeted debug/benchmark runs.
