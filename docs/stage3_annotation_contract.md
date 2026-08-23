# Stage 3 Annotation Contract

Stage 3 turns partitioned alignment events into durable variant-context evidence.
It is an internal boundary of the one end-to-end pipeline, not a standalone CLI
mode.

## Ownership

Stage 3 owns external evidence whose lookup is required for every candidate
variant:

- canonical event-to-variant normalization;
- ClinVar evidence from the configured indexed VCF;
- gnomAD lookup status and selected normalized fields;
- Ensembl VEP consequence evidence.

It does not own report-specific support thresholds, taxonomic scopes, depth
counters, feature coverage, histograms, or plotting tables. Analytics derives
those views from Stage 1 taxonomy/features, Stage 2 segments and exact support,
and the Stage 3 event-to-variant lineage.

## Inputs

Stage 3 consumes:

- `alignment/evidence/partitions/<partition_id>/` with Stage 2 summaries,
  segments, compact events, and exact event-to-ortholog support;
- the matching target Gene table and target FASTA context from Stage 1;
- an indexed ClinVar VCF;
- the configured gnomAD API and optional shared regional cache;
- one declared VEP backend and release contract.

The cluster uses local release-pinned VEP. Small local runs may use Ensembl REST.
Local VEP requires `vep_release`, `vep_executable`, and `vep_cache_dir`. Completed
variant/gene results can be reused through the shared immutable
`vep_result_cache_dir`; the official VEP reference cache and the result cache are
different resources.

## Processing

For each genomic partition:

1. normalize every compact event against target context and write one
   partition-local `event_group_id` to canonical `variant_key` lineage row;
2. collapse repeated events to unique gene/variant contexts without losing
   strategy membership;
3. fetch bounded ClinVar and gnomAD evidence and retain explicit lookup failures;
4. write deterministic headered gzip shards of at most 250,000 variant-context
   rows;
5. run VEP independently for each shard, preserving input row order and adding
   the declared VEP fields;
6. validate the complete shard set and copy each completed enriched shard into
   the durable dataset without a global merge or recompression.

The shard boundary limits task memory and is also the Nextflow `-resume`
boundary for VEP. Missing, extra, reordered-schema, differently configured, or
row-count-mismatched shards make finalization fail; a partial dataset is never
published.

## Durable Outputs

```text
annotation/
  manifest.json
  failures.tsv.gz
  event_variant_map/
    partitions/<partition_id>/event_variant_map.tsv.gz
  variant_annotations/
    manifest.json
    partitions/<partition_id>/<shard_id>.tsv.gz
```

`annotation/manifest.json` uses
`schema=normalized_annotation_evidence_v3`. Its `variant_annotations` object is
identical to `annotation/variant_annotations/manifest.json` and declares:

- `schema=gaph_variant_annotation_dataset_v1`;
- exact fields and row, shard, and partition counts;
- each relative shard path, compressed size, and row count;
- the semantic VEP backend/release/options contract;
- explicit VEP status counts.

Each shard retains the normalized source fields plus:

- `vep_status`;
- `vep_primary_consequence`;
- `vep_consequence_terms`;
- `vep_transcript_id`;
- `vep_mane_select`;
- `vep_canonical`;
- `vep_impact`;
- `vep_variant_class`.

Non-concrete events remain in the event map with their normalization status and
an empty canonical key. Variant-context rows with invalid keys retain an
explicit `vep_status`; individual non-`ok` VEP outcomes do not make a complete
dataset partial.

The event map uses partition-local `event_group_id`. Its durable join key is
`(partition_id, event_group_id)`. Variant shards are headered and directly
scannable by DuckDB; there is no durable global annotation TSV and no separate
analytics-owned copy of the same candidate table.

## Analytics Contract

Completed-run analytics requires the current root annotation manifest and the
exact child variant dataset manifest. It validates all declared shard files and
scans them as one virtual table. Cohorts reference member dataset manifests and
shards rather than concatenating large variant tables.

Analytics then derives strategy/feature coverage, site depth, taxonomic depth,
ALT support, genus support, and report summaries under fingerprinted caches.
Deleting those caches is safe; deleting Stage 1–3 evidence is not.

## Storage And Performance

Stage 3 avoids two large duplicate relations:

- pre-VEP shards remain temporary Nextflow task outputs;
- enriched shards are published once without constructing a second global gzip
  table.

The shared VEP result cache is reusable infrastructure, not part of a run and
not a substitute for durable shard output. Compare performance only on the same
candidate rows, backend/release, cache state, resources, and partition size.
Use the Nextflow trace for per-task runtime, RSS, and I/O.
