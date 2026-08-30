# Stage 2 Alignment Contract

Stage 2 aligns each selected non-human ortholog gene sequence from Stage 1
against the fixed human target locus and writes normalized alignment evidence
for downstream variant-support analysis.

## Position In The Pipeline

Production execution is a single end-to-end workflow:

```text
Entrez IDs
  -> Stage 1 fetch
  -> Stage 2 alignment
  -> Stage 3 annotation
  -> analytics
```

The stage boundary is internal to the single pipeline execution path. Nextflow
reuses completed fetch and alignment processes through `-resume`; there is no
standalone alignment CLI mode.

## Inputs

Stage 2 consumes the normalized Stage 1 handoff:

- `genes.tsv.gz`
- `target_features.tsv.gz`
- `orthologs.selected.tsv.gz`
- `sequences/targets/<gene_id>.fa.gz`
- `sequences/orthologs/<gene_id>.fa.gz`

The ortholog FASTA directory is an internal Nextflow handoff from fetch to
alignment. It is not part of the completed run's durable `fetch/` directory;
selected-ortholog metadata remains durable there.

`target_features.tsv.gz` is fingerprinted in the final alignment manifest for
the later analytics join; it is not copied into aligner tasks. Taxonomy is not a
Stage 2 input.

Target FASTA records are already in plus genomic orientation on GRCh38.p14.
Ortholog records are aligned as query sequences; aligners decide forward/reverse
orientation.

## Strategies

Runnable strategies are registered in the workflow and can be selected with
`--alignment_strategies`. The value `default` selects every strategy marked as
default-enabled in that registry.

| Strategy | Tool | Policy |
| --- | --- | --- |
| `minimap2_asm10` | minimap2 | Fixed baseline preset for every ortholog. |
| `minimap2_asm20` | minimap2 | More permissive fixed minimap2 preset. |
| `minimap2_map_ont_pseudoreads_30000_15000` | minimap2 | Error-free 30 kb long pseudo-reads at a 15 kb step, aligned with `map-ont` and reduced to a dominant-strand monotonic backbone. |
| `nucmer` | MUMmer/nucmer | Independent comparator using multi-query nucmer output. |
| `bwa_pseudoreads_150_75` | BWA/samtools/pysam | Pseudoread comparator using 150-base reads at a 75-base step. |

The `asm20` and long-pseudoread minimap2 strategies, nucmer, and BWA
pseudoreads are default-enabled. The `asm10` strategy is runnable only when
named explicitly.

No LASTZ, consensus calling, or production variant filtering is part of Stage 2.
Conservation scores such as GERP are not part of alignment; they belong to the
later annotation/analysis layer.

Example selections:

```bash
--alignment_strategies default
--alignment_strategies minimap2_asm20
--alignment_strategies minimap2_asm10,nucmer
--alignment_strategies bwa_pseudoreads_150_75
--alignment_strategies minimap2_map_ont_pseudoreads_30000_15000
```

At least one strategy must be selected. Single-strategy runs are valid; compare
or report layers must treat cross-strategy-only sections as not applicable or
empty rather than failing.

`default` selects `minimap2_asm20`,
`minimap2_map_ont_pseudoreads_30000_15000`, `nucmer`, and
`bwa_pseudoreads_150_75`.

Each bounded partition merge reduces raw per-ortholog observations to one
compact row per unique event and strategy. Every compact row receives a
partition-local `event_group_id`; the exact positive-support sidecar refers to
that ID instead of repeating event coordinates. Raw event rows exist only
inside aligner task outputs before the partition merge.

## Taxonomy Metadata

Stage 1 requests the NCBI taxonomy summary once for the selected tax IDs and
publishes status, ordered lineage, and direct named-rank identifiers from domain
through species. Stage 2 neither consumes nor republishes that table. Analytics
joins the Stage 1 taxonomy to Stage 2 exact ortholog support when a report needs
taxonomic scopes or ranks. Report-specific membership flags and counts are
derived from the lineage. Alignment presets belong to the strategy registry,
not taxonomy.

## Alignment Granularity

One Nextflow task processes one gene and one strategy.

For minimap2:

```text
target gene FASTA vs multi-FASTA of all selected orthologs for that gene
```

This is equivalent in intent to running each ortholog separately with the same
reference and options: minimap2 maps each query independently against the target
index. Multi-query execution avoids rebuilding the target index and avoids
creating thousands of scheduler tasks.

All minimap2 modes run through one `ALIGN_MINIMAP2` process. Strategy metadata
selects the opt-in `asm10`, default `asm20`, or default `map-ont`
long-pseudoread mode; there are no duplicate workflow modules. Selected
minimap2 strategies share that process's `alignment_max_forks` budget.

The long-pseudoread strategy cuts each ortholog deterministically into 30,000
base windows at a 15,000 base step. An ortholog of at most 30,000 bases is
aligned once at its full length. Longer sequences always receive a final window
ending at the sequence boundary, so generation cannot leave a terminal gap.
The generated sequence is copied exactly from the ortholog: no Nanopore error
model is injected. `map-ont` supplies minimap2's long-read seeding and chaining
policy; the strategy is therefore an alignment-geometry comparator, not a
simulation of ONT sequencing accuracy.

After mapping, all records for one source ortholog are reduced first to the
dominant strand and then to a longest monotonic read-order backbone in target
coordinate order. This keeps coherent long-range placements while permitting
accepted secondary records. Coordinates are lifted back to the complete source
ortholog, and repeated identical events from overlapping windows count as one
support from that ortholog.

For nucmer:

```text
target gene FASTA vs multi-FASTA of all selected orthologs for that gene
```

Raw multi-query nucmer is used, but the workflow does not run global
`delta-filter -1`. One-to-one delta filtering is appropriate for comparing two
assemblies, but it is wrong for many orthologs that are all expected to align to
the same human target locus. Nucmer emits SAM-long records and the Python parser
normalizes their CIGAR operations per query sequence. A contiguous CIGAR
insertion or deletion is therefore published as one event rather than one row
per affected base. Identical events repeated by overlapping Nucmer records are
collapsed within an ortholog. The parser adds `unfiltered_nucmer` QC flags.
Events containing IUPAC ambiguity symbols are excluded and counted in the task
manifest; the affected ortholog summary receives `ambiguous_event_allele`.

The BWA comparator `bwa_pseudoreads_150_75` uses 150-base pseudo-reads, a
75-base step, and synthetic PHRED 30. These values live in the strategy
registry and are passed explicitly to the shared BWA runner; they are not
runtime tuning parameters or duplicated runner defaults.

Mapped pseudo-reads are reduced by two mandatory rules: retain the dominant
strand for each ortholog, then retain its monotonic LIS/LDS backbone in target
coordinate order. Reads are not deduplicated merely because they share an
alignment position. There are no mapped-fraction, filtered-fraction, retained-
fraction, or whole-ortholog rejection thresholds.

## Durable Outputs

Alignment publishes one canonical partitioned evidence contract:

| Path | Meaning |
| --- | --- |
| `manifest.json` | Dataset schema, exact selected-strategy-eligible `gene_ids`, enabled strategies, source-context fingerprints, evidence layout, and durable output counts. |
| `evidence/partitions/<partition_id>/manifest.json` | Partition identity, exact gene set, strategies, and normalized evidence counts. |
| `evidence/partitions/<partition_id>/ortholog_alignment_summary.tsv.gz` | One row per gene/ortholog/output strategy when that strategy emits summary evidence. |
| `evidence/partitions/<partition_id>/alignment_segments.tsv.gz` | Normalized alignment intervals. |
| `evidence/partitions/<partition_id>/alignment_events.tsv.gz` | One compact row per unique target-coordinate event and strategy, keyed by a partition-local consecutive `event_group_id`. |
| `evidence/partitions/<partition_id>/event_ortholog_support.tsv.gz` | Exact positive ortholog identities keyed by the partition-local `event_group_id`. |
| `failures.tsv.gz` | Alignment-stage failures. |

Stage 1 remains the owner of taxonomy, selected-ortholog metadata, target
features, and sequences. Stage 2 references that source dataset but does not
republish it. Native aligner outputs and analytic aggregates are not durable
Stage 2 outputs.

The alignment directory is an indivisible internal annotation handoff.
Annotation receives it directly from the alignment workflow, requires
`normalized_alignment_evidence_v2`, and validates every declared partition and
file. Individual tables and prior layouts are not independent pipeline inputs.

Every per-gene result and partition uses plural `gene_ids` and `strategies`,
nested `strategy_parameters`, canonical evidence counts, and an explicit
`alignment_event_mode`. Merge rejects any other schema.

The four normalized per-aligner tables have one exact ordered schema, defined
in `bin/alignment_table_schema.py`: summaries, segments, raw events, and
failures. Partition merge validates their full headers before reading evidence
and rejects missing, extra, or reordered fields. In particular,
`ortholog_gene_id`, `tax_id`, and `strategy` cannot disappear from raw event
input and be replaced with empty values. Final merge applies the same exact
check to compact events and their `event_group_id`-linked ortholog-support
sidecar. This validation does not change the valid Stage 2 output format.

The compact event row contains event identity and event-level QC flags. It does
not repeat tool or preset metadata, which is fixed by the strategy manifest, or
support counters, which analytics derives from
`event_ortholog_support.tsv.gz`. The sidecar retains each exact positive
ortholog identity, its `tax_id`, representative alignment metadata, and the
lossless raw-row multiplicity needed for later support metrics. Taxon names and
lineage are joined once from the canonical Stage 1 taxonomy table.

At the final alignment boundary, the pipeline copies normalized partitions
without parsing, recompression, or global event-ID rebasing. An
`event_group_id` is local to one partition; its durable join identity is
`(partition_id, event_group_id)`. Annotation discovers exactly this layout and
materializes only the small target context required by each partition.

The alignment stage does not compute or publish strategy summaries, feature coverage,
site depth, taxonomic depth, or ALT taxonomic counters. Analytics derives and
caches those views from normalized segments, events, exact support, canonical
Stage 1 taxonomy, and target features. Alignment failures remain operational
evidence. Per-partition phase timings remain in the final manifest; they are
task measurements, not a synthetic wall-clock total.

Stage 3 follows the same evidence-first boundary. Its canonical
`normalized_annotation_evidence_v4` output contains the partitioned
ClinVar/gnomAD/VEP variant dataset, partitioned `event_variant_map`, annotation
`failures.tsv.gz`, and `manifest.json`. The event map preserves the exact
`(partition_id, event_group_id)` to canonical `variant_key` lineage needed to
derive variant-strategy, exact-ortholog, taxonomy, and site-depth views in
analytics. Those derived support tables are not Stage 3 outputs. The complete
annotation contract is documented in `docs/stage3_annotation_contract.md`.

For Minimap2 rows, `native_record_id` is derived from the PAF record content
rather than its output line number. `event_id` combines the strategy, that
stable record identifier, and the event ordinal within the PAF `cs` tag. These
identifiers therefore remain stable when Minimap2 emits identical records in a
different order at another thread count.

Minimap2, Nucmer, and BWA retain events from accepted primary and non-primary
alignment records. Raw event rows and exact ortholog support carry nullable
`mapq` and `native_alignment_type`; the compact event aggregate does not copy
these record-level fields. Minimap2 preserves the literal PAF `tp` value (`P`,
`S`, `I`, or `i`). Nucmer and BWA use the SAM flags to report `primary`,
`secondary`, `supplementary`, or `secondary_supplementary`.

`mapq` is the integer reported by the strategy's aligner and is only
interpretable within that strategy. Stage 2 applies no MAPQ cutoff and does not
convert low-MAPQ records to no-calls. If primary and non-primary records from
the same ortholog emit the same event, the existing representative-selection
rule still prefers the primary record while Stage 2 compaction sums
`support_row_count`. If annotation normalizes multiple compact events to one
canonical variant, analytics unions their ortholog identities and sums that
lossless multiplicity. Alignment role is not duplicated in `qc_flags`.

Strategy summaries are analytic views over
`ortholog_alignment_summary.tsv.gz`, not alignment-owned source evidence.

## Why Segments And Events Are Separate

`alignment_events.tsv.gz` records observed differences. It does not record where
an ortholog matches the human target.

`alignment_segments.tsv.gz` records coverage intervals. Together with compact
events and exact supporters, it is sufficient for analytics to distinguish:

```text
event present at position       -> ortholog supports observed alternative
position covered and no event   -> ortholog supports human/reference allele
position not covered            -> no-call
bad/ambiguous alignment         -> filtered/no-call
```

The pipeline keeps segments in this partitioned form. Analytics computes site
depth and taxonomy-aware counts only
when requested, using canonical Stage 1 taxonomy. Different ALT alleles at one
site therefore share the same segment-derived denominator. Multiple overlapping
segments from the same ortholog must be merged before depth is counted so one
ortholog cannot inflate a per-base or per-feature depth statistic.

## Coordinate Convention

- Segment intervals use 0-based half-open target coordinates:
  `target_start0`, `target_end0`.
- Events include target-local coordinates and GRCh38 coordinates:
  `genomic_accession`, `genomic_start1`, `genomic_end1`.
- Target coordinates always refer to the plus genomic target sequence from
  Stage 1.
- Query strand is stored per segment/event.
- VCF-style left normalization is not performed in Stage 2; it belongs to the
  later variant/support stage.

## Storage Policy

Alignment processes use scratch task space. Raw PAF/SAM/BAM files and generated
pseudoreads are temporary task inputs and are never published as durable output.

Per-gene alignment task directories contain only the task manifest and metadata
needed to run a strategy. They do not copy the Stage 1 target or ortholog FASTA
files. Sequence-based strategies receive the needed per-gene FASTA files as
explicit Nextflow inputs and materialize uncompressed aligner input FASTA files
inside their scratch workspace.

Task preparation receives the Stage 1 sequence directory as one staged input,
not one command-line argument per FASTA. It streams the required
`query_gene_id`-grouped ortholog table and keeps only one gene's metadata in
memory. A repeated, non-contiguous `query_gene_id` is a contract error; old
ungrouped fetch outputs must be regenerated.

Task preparation does not stage `target_features.tsv.gz` into aligner jobs.
Feature coverage is an analytics-owned view derived later from durable segments
and the canonical Stage 1 feature table.

Task preparation also assigns a stable `partition_id` after sorting target genes
by chromosome and genomic interval. `--alignment_partition_size` controls the
maximum number of target genes in each partition (default: 10). The identifier
is carried in task metadata for bounded downstream merge and annotation; it does
not change gene-level alignment behavior.

The same preparation pass sums the validated selected-ortholog
`sequence_length` values for each gene into internal `ortholog_sequence_bp`
task metadata. Nextflow uses that input-size measure to choose the initial
memory request for the two pseudoread strategies. The long-pseudoread minimap2
strategy starts at 8 GB below 150 million bases, 24 GB below 600 million bases,
and 40 GB otherwise. BWA starts at 8 GB below 600 million bases and 16 GB
otherwise. Ordinary assembly minimap2 stays at 8 GB and nucmer at 12 GB.
Retries remain a safety margin for unusual tasks; they are not the primary
mechanism for reaching a predictable memory class.

Alignment results are merged in two bounded levels. Each partition merges at
most `--alignment_partition_size` genes and streams one raw event group at a
time from its SQLite key index into the compact event/support pair. The final
merge streams those disjoint compact partitions through one staged directory
without rebasing their local group IDs. This avoids one global event database
and avoids placing every gene/strategy result path on a single command line
while preserving the Stage 2 logical evidence.

Both merge levels fail closed. A partition must contain exactly one result for
every eligible gene/strategy pair, required TSV inputs must exist with valid
headers, and genes cannot occur in multiple partitions. Every strategy requires
fetched ortholog inputs. The final merge compares the union of partition
`gene_ids` with the genes eligible for at least one selected strategy in
`alignment_tasks.tsv.gz`; incomplete output is an error rather than a smaller
successful dataset.

Durable output is limited to compressed normalized TSV files so large runs do
not duplicate sequence data and native aligner output in `results/`.

Nextflow `work/` remains a resume cache. After validating a run, it can be
cleaned according to `docs/storage_model.md`.
