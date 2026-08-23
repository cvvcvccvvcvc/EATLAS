# Pipeline Scaling Notes

This document records optimization decisions that are not required for the
current narrow cleanup, but matter for larger runs.

## Raw Events vs Compact Support

Aligner tasks emit one raw row per observed event, strategy, and support unit.
For minimap2/nucmer/BWA the support unit is an ortholog; for Ensembl MAF it is a
species row. These raw rows are an internal input to the bounded partition
merge, not a public Stage 2 output.

The partition merge writes `alignment_events.tsv.gz` with one row per unique
target event and strategy plus support counts. A separate
`event_ortholog_support.tsv.gz` retains one positive row per supporting ortholog,
keyed by the compact row's `event_group_id`. The partition merge writes both in
one index-ordered pass, so it does not globally group and sort the almost
unreduced event-by-ortholog relation. Native-record traceability is still
dropped. Raw per-task event files remain recoverable through Nextflow `work/`
while the cache is retained.

End-to-end and standalone alignment publish the same partition tree and do not
write a second global event table. Annotation consumes that tree directly. It
publishes the unique `variant_annotations.tsv.gz` and a partitioned
event-to-variant map, but not another event-by-ortholog relation. Analytics
streams the event, segment, exact-support, and event-map relations to derive the
smaller report tables only when requested.

## Partitioned Alignment Outputs

Alignment merging is partitioned by genomic-order target groups. A partition is
released as soon as all expected strategy results for its genes are ready. Each
partition retains per-ortholog summaries, segments, compact events, and exact
support. The final alignment process validates and copies those gzip files
without global sorting, recompression, or event-ID rebasing. Standalone and
end-to-end modes share this exact path.

## BWA Pseudoreads

The BWA pseudoread strategy still does more compute work than the other local
aligners because it generates synthetic reads, builds a per-gene BWA index, runs
BWA, sorts and indexes BAM, and applies the LIS BAM filter.

Current behavior:

- `bwa_pseudoreads_150_75` extracts BAM/CIGAR-supported events with `pysam`;
- `bwa mem` is streamed directly into `samtools sort`, avoiding durable SAM and
  unsorted BAM intermediates;
- native BAM files are kept only with `--keep_native_alignments true`;
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
