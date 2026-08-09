# Pipeline Scaling Notes

This document records optimization decisions that are not required for the
current narrow cleanup, but matter for larger runs.

## Raw Events vs Compact Support

Standalone `--stage align` output remains raw by default:

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

In compact mode, `alignment_events.tsv.gz` contains one row per unique target
event and strategy plus support counts. A separate
`event_ortholog_support.tsv.gz` retains one positive row per supporting ortholog,
keyed by the compact row's `event_group_id`. The partition merge writes both in
one index-ordered pass, so it does not globally group and sort the almost
unreduced event-by-ortholog relation. Native-record traceability is still
dropped. Raw per-task event files remain recoverable through Nextflow `work/`
while the cache is retained.

End-to-end `--stage all` does not publish a global event table. Partition events
remain in `work/` until annotation consumes them, and Stage 3 preserves compact
per-strategy ALT-support counts in `variant_strategy_support.tsv.gz` plus exact
positive supporters in the `variant_ortholog_support/` Parquet dataset. Exact
support is aggregated inside each annotation partition using local integer IDs;
finalization copies the Parquet parts without expanding and recompressing them.

## Partitioned Alignment Outputs

Alignment merging is internally partitioned by genomic-order target groups. A
partition is released as soon as all expected strategy results for its genes are
ready. In `--stage all`, each partition retains only events needed by annotation
plus strategy/coverage summaries and failures. The final alignment merge
combines only those small summaries and records raw row counts from partition
manifests. It does not rewrite a global event, segment, or per-ortholog summary
table.

Standalone `--stage align` preserves the full global handoff for a later
annotation invocation.

## BWA Pseudoreads

The BWA pseudoread strategy still does more compute work than the other local
aligners because it generates synthetic reads, builds a per-gene BWA index, runs
BWA, sorts and indexes BAM, and applies the LIS BAM filter.

Current behavior:

- `bwa_pseudoreads` extracts BAM/CIGAR-supported events with `pysam`;
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
