# Pipeline Scaling Notes

This document records optimization decisions that are not required for the
current narrow cleanup, but matter for larger runs.

## Raw Events vs Compact Support

Default alignment output remains raw:

```text
alignment_events.tsv.gz
```

In raw mode, one row means one observed event from one strategy for one support
unit. For minimap2/nucmer/BWA the support unit is an ortholog. For Ensembl MAF
it is a species row.

Large runs can enable:

```bash
--compact_alignment_events true
```

In compact mode, the final `alignment_events.tsv.gz` contains one row per unique
target event key plus support counts. This reduces durable disk usage and makes
annotation cheaper, but it drops per-ortholog/native-record traceability from
the final published table. Raw per-task event files remain recoverable through
Nextflow `work/` while the cache is retained.

Future direction: if compact mode becomes the default for production, publish
raw event tables only behind an explicit debug flag such as
`--keep_raw_alignment_events`.

## Partitioned Alignment Outputs

The current final merge process creates one global set of alignment TSV files.
That is simple and convenient for smoke runs, but it does not scale cleanly to
thousands of genes because one process owns all final raw rows.

Future direction:

- write durable alignment outputs partitioned by gene chunk or gene ID;
- write a small manifest/index that lists partitions and row counts;
- build compact support tables through chunk-level aggregation followed by a
  final aggregation step;
- keep raw event rows partitioned as debug evidence rather than one global file.

This keeps memory bounded and avoids a single huge `alignment_events.tsv.gz`.

## BWA Pseudoreads

`bwa_pseudoreads` is intentionally left out of the current cleanup. It currently
does more temporary I/O than the other aligners:

- copies the whole task directory into the output directory;
- generates pseudoreads;
- builds a BWA index per task;
- writes SAM, BAM, sorted BAM, BAM index, mpileup, and VarScan VCF files;
- emits two event streams: `bwa_pseudoreads_pysam` and
  `bwa_pseudoreads_varscan`.

Future direction:

- make the module write to `--outdir` like the other aligner scripts instead of
  copying the task directory;
- stream `bwa mem` into `samtools view/sort` to avoid durable SAM and unsorted
  BAM intermediates;
- make VarScan optional as a separate strategy or submode;
- write normalized summaries/segments consistently with other strategies.

## Ensembl Compara MAF Manifest

The precomputed MAF strategy can build its chunk manifest during a run, but that
requires probing MAF gzip files. For remote HTTP/FTP sources this can be slow and
fragile.

Production runs should prefer a precomputed manifest through:

```bash
--ensembl_compara_maf_manifest /path/to/ensembl_compara_maf_manifest.tsv.gz
```

or:

```bash
export ENSEMBL_COMPARA_MAF_MANIFEST=/path/to/ensembl_compara_maf_manifest.tsv.gz
```

If no explicit manifest is provided, the workflow also checks:

```text
assets/reference/ensembl/compara/release-<release>/<species_set>/ensembl_compara_maf_manifest.tsv.gz
```

If that file is absent, the workflow falls back to building the manifest during
the run.
