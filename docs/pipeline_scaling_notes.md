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

Alignment merging is internally partitioned by genomic-order target groups. A
partition is released as soon as all expected strategy results for its genes are
ready, and compact-event aggregation is bounded to that partition. The final
merge streams partition files and still publishes one global set of alignment
TSVs for compatibility with current downstream tools.

Remaining future direction:

- publish durable alignment outputs as partitions instead of one global table;
- write a small manifest/index that lists partitions and row counts;
- build compact support tables through chunk-level aggregation followed by a
  final aggregation step;
- keep raw event rows partitioned as debug evidence rather than one global file.

This keeps memory bounded and avoids a single huge `alignment_events.tsv.gz`.

## BWA Pseudoreads

The BWA pseudoread strategies still do more compute work than the other local
aligners because they generate synthetic reads, build a per-gene BWA index, run
BWA, sort/index BAM, and apply the LIS BAM filter.

Current behavior:

- `bwa_pseudoreads` and `bwa_pseudoreads_varscan` are separate selectable
  strategies;
- VarScan is required only when `bwa_pseudoreads_varscan` is selected;
- `bwa mem` is streamed through `samtools view` and `samtools sort`, avoiding
  durable SAM and unsorted BAM intermediates;
- native BAM, mpileup, and VCF files are kept only with
  `--keep_native_alignments true`;
- normalized segments, summaries, events, failures, and manifest files are
  emitted like the other aligner strategies.

Future direction:

- avoid rebuilding the BWA target index when many strategies or repeated runs
  reuse the same target sequence;
- evaluate whether pseudoread generation can be made sparser for long genes
  without losing sensitivity at exon/feature boundaries.

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
