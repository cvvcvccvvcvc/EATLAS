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
  -> later variant/support stages
```

The debug-only `--stage align --fetch_dir <dir>` mode exists to re-run alignment
from an already published Stage 1 result without downloading NCBI data again.

## Inputs

Stage 2 consumes the normalized Stage 1 outputs:

- `genes.tsv.gz`
- `target_features.tsv.gz`
- `orthologs.selected.tsv.gz`
- `taxonomy.tsv.gz`
- `taxonomy_failures.tsv.gz`
- `taxonomy_summary.tsv.gz`
- `sequences/targets/<gene_id>.fa.gz`
- `sequences/orthologs/<gene_id>.fa.gz`

Target FASTA records are already in plus genomic orientation on GRCh38.p14.
Ortholog records are aligned as query sequences; aligners decide forward/reverse
orientation.

## Strategies

Runnable strategies are registered in the workflow and can be selected with
`--alignment_strategies`. The default value is `all`, meaning every strategy
marked as default-enabled in that registry.

| Strategy | Tool | Policy |
| --- | --- | --- |
| `minimap2_asm10` | minimap2 | Fixed baseline preset for every ortholog. |
| `minimap2_asm20` | minimap2 | More permissive fixed minimap2 preset. |
| `minimap2_map_ont_pseudoreads_30000_15000` | minimap2 | Error-free 30 kb long pseudo-reads at a 15 kb step, aligned with `map-ont` and reduced to a dominant-strand monotonic backbone. |
| `nucmer` | MUMmer/nucmer | Independent comparator using multi-query nucmer output. |
| `bwa_pseudoreads_150_75` | BWA/samtools/pysam | Pseudoread comparator using 150-base reads at a 75-base step. |
| `precomputed_ensembl_92_mammals_epo_extended` | Ensembl Compara MAF | Uses release-pinned precomputed `92_mammals.epo_extended` whole-genome MSA blocks overlapping the human target gene interval. |

The two assembly-mode minimap2 strategies, nucmer, and BWA pseudoreads are
default-enabled. The long-pseudoread and precomputed Ensembl strategies are
runnable only when named explicitly.

No LASTZ, consensus calling, or production variant filtering is part of Stage 2.
Conservation scores such as GERP are not part of alignment; they belong to the
later annotation/analysis layer.

Example selections:

```bash
--alignment_strategies all
--alignment_strategies minimap2_asm20
--alignment_strategies minimap2_asm10,nucmer
--alignment_strategies bwa_pseudoreads_150_75
--alignment_strategies minimap2_map_ont_pseudoreads_30000_15000
--alignment_strategies minimap2_asm20,precomputed_ensembl_92_mammals_epo_extended
```

At least one strategy must be selected. Single-strategy runs are valid; compare
or report layers must treat cross-strategy-only sections as not applicable or
empty rather than failing.

`all` selects `minimap2_asm10`, `minimap2_asm20`, `nucmer`, and
`bwa_pseudoreads_150_75`. The explicitly selected Ensembl strategy uses release 116,
the `92_mammals.epo_extended` set, and at most three concurrent remote chunk
tasks. These values are part of the strategy definition rather than separate
user options.

Each bounded partition merge reduces raw per-ortholog observations to one
compact row per unique event and strategy. Every compact row receives a
partition-local `event_group_id`; the exact positive-support sidecar refers to
that ID instead of repeating event coordinates. Raw event rows exist only
inside aligner task outputs before the partition merge.

## Taxonomy Metadata

Stage 1 requests the NCBI taxonomy summary once for the selected tax IDs and
publishes status, ordered lineage, and direct named-rank identifiers from domain
through species. Stage 2 consumes that handoff and never fetches taxonomy.
These source fields support the current report's taxonomic scope and
evidence-unit controls; they do not change aligner selection. Report-specific
membership flags and counts are derived from the lineage rather than stored in
the canonical taxonomy row. Alignment presets belong to the strategy registry,
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
selects `asm10`, `asm20`, or the opt-in `map-ont` long-pseudoread mode; there
are no duplicate workflow modules. Selected minimap2 strategies share that
process's `alignment_max_forks` budget.

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

For Ensembl Compara MAF:

```text
precomputed whole-genome MSA blocks overlapping the human target gene interval
```

The workflow first builds a small run-specific MAF chunk manifest directly from
`genes.tsv.gz`. One alignment task streams each required source chunk once for
all overlapping target genes and routes normalized rows into gene fragments.
All fragments for a gene are then consolidated before the normal alignment
merge; target intervals are unioned and feature coverage is recomputed from the
complete gene evidence. Transient network and truncated-gzip read failures are
retried inside the same process with block-level continuation: already committed
MAF blocks are skipped on the next network attempt. If a source still cannot be
fully read after all attempts, the task records gene-level failure rows.
Unexpected process failures terminate the
workflow rather than silently producing incomplete gene evidence. Full MAF
chunks are not published as durable outputs.

Before coordinates and events are derived, dot placeholders in non-reference
MAF rows are resolved against the human alignment row. A dot opposite a human
base is treated as that matching base; a dot opposite a human gap is treated as
a gap. This prevents placeholders from becoming artificial insertions or
advancing query coordinates. Indels containing other non-ACGT symbols are kept
out of `alignment_events.tsv.gz` and marked in the ortholog summary QC flags.

The MAF chunk manifest can be supplied with `--ensembl_compara_maf_manifest` or
`ENSEMBL_COMPARA_MAF_MANIFEST`. If neither is set, the workflow checks
`assets/reference/ensembl/compara/release-116/92_mammals.epo_extended/ensembl_compara_maf_manifest.tsv.gz`;
if that file is absent, it builds the manifest during the run.

This strategy is not based on NCBI ortholog GeneIDs. Its support unit is the
species row in the precomputed MSA. Consequently, `ortholog_gene_id` contains
the species name for this strategy, and support-count reports should interpret
it as species support.

When one gene spans multiple MAF source chunks, Stage 2 consolidates those
fragments before downstream merge. Target coverage is calculated from the union
of target intervals. A species-level MSA does not provide one meaningful query
length denominator across all blocks, so consolidated MAF rows leave
`query_length` and `query_coverage` empty and add
`maf_query_coverage_not_applicable`; `aligned_query_bp` remains available.

## Durable Outputs

Standalone `--stage align` publishes the full handoff contract:

| Path | Meaning |
| --- | --- |
| `manifest.json` | Exact selected-strategy-eligible `gene_ids`, per-strategy eligibility counts, enabled strategies, and output counts. |
| `alignment_tasks.tsv.gz` | Per-gene task manifest with separate human-target and fetched-ortholog readiness checks. |
| `taxonomy.tsv.gz` | Compact tax_id-to-lineage mapping. |
| `taxonomy_failures.tsv.gz` | Taxonomy lookup warnings/failures. |
| `taxonomy_summary.tsv.gz` | Run-level ortholog and taxonomic-unit counts by scope. |
| `ortholog_alignment_summary.tsv.gz` | One row per gene/ortholog/output strategy when that strategy emits summary evidence. |
| `strategy_summary.tsv.gz` | Small canonical per-strategy aggregate derived from `ortholog_alignment_summary.tsv.gz` for reports and run inspection. |
| `alignment_segments.tsv.gz` | Normalized alignment intervals. |
| `snv_site_depth.tsv.gz` | Distinct aligned-ortholog depth for each observed concrete SNV position and strategy. |
| `snv_taxonomic_depth.tsv.gz` | Taxonomy-aware aligned depth for each observed concrete SNV position and strategy. |
| `snv_alt_taxonomic_support.tsv.gz` | Taxonomy-aware exact-ALT support for each observed concrete SNV and strategy. |
| `feature_coverage.tsv.gz` | Per-gene, per-strategy coverage and depth over target structural intervals. |
| `alignment_events.tsv.gz` | One compact row per unique target-coordinate event and strategy, keyed by consecutive `event_group_id`. |
| `event_ortholog_support.tsv.gz` | Exact positive ortholog identities keyed by `event_group_id`. |
| `failures.tsv.gz` | Alignment-stage failures. |
| `native/` | Optional raw PAF/SAM files when enabled. |

Native outputs are disabled by default.

The standalone directory is an indivisible Stage 3 handoff. Annotation accepts
the directory through `--alignment_dir`, validates its manifest and table
counts, and verifies that its target-context fingerprints match the supplied
`--fetch_dir`. Individual tables are not accepted as independent pipeline
parameters.

Every per-gene result and partition uses the same manifest keys: plural
`gene_ids` and `strategies`, nested `strategy_parameters`, canonical output
counts, and an explicit `alignment_event_mode`. Merge rejects missing or legacy
singular fields. `feature_coverage.tsv.gz` is mandatory for every alignment
result and is never reconstructed by a later fallback.

The four normalized per-aligner tables have one exact ordered schema, defined
in `bin/alignment_table_schema.py`: summaries, segments, raw events, and
failures. Partition merge validates their full headers before reading evidence
and rejects missing, extra, or reordered fields. In particular,
`ortholog_gene_id`, `tax_id`, and `strategy` cannot disappear from raw event
input and be replaced with empty values. Final merge applies the same exact
check to compact events and their `event_group_id`-linked ortholog-support
sidecar. This validation does not change the valid Stage 2 output format.

The compact event row contains event identity, raw/supporting-ortholog counts,
and event-level QC flags. It does not repeat tool or preset metadata, which is
fixed by the strategy manifest, or distinct taxon/name counts, which are
reproducible from `event_ortholog_support.tsv.gz`. Exact taxonomic identities
remain in that sidecar rather than being replaced by those counts.

In an end-to-end `--stage all` run, each alignment partition writes normalized
source evidence before Stage 3 consumes the same partition. The final alignment
merge copies that evidence without parsing, recompression, or global event-ID
rebasing into `alignment/evidence/partitions/<partition_id>/`. Each partition
contains its manifest, `ortholog_alignment_summary.tsv.gz`,
`alignment_segments.tsv.gz`, compact `alignment_events.tsv.gz`, and exact
`event_ortholog_support.tsv.gz`. An `event_group_id` is local to one partition;
the durable join identity is `(partition_id, event_group_id)`.

The current Stage 3 compatibility path still consumes partition-level SNV site
depth and taxonomy-aware counts so existing annotation tables and reports keep
their established results. Those derived tables are not copied into the
durable evidence dataset. The global `strategy_summary.tsv.gz`,
`feature_coverage.tsv.gz`, `taxonomy_summary.tsv.gz`, and `failures.tsv.gz`
remain published for the same reason. The final alignment manifest describes
the partitioned evidence layout in `normalized_evidence` and retains
per-partition phase timings plus summed task-runtime totals. The totals are not
wall-clock time because partitions can run concurrently.

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
`secondary`, `supplementary`, or `secondary_supplementary`. Ensembl Compara MAF
leaves both fields empty because it has no equivalent native record metadata.

`mapq` is the integer reported by the strategy's aligner and is only
interpretable within that strategy. Stage 2 applies no MAPQ cutoff and does not
convert low-MAPQ records to no-calls. If primary and non-primary records from
the same ortholog emit the same event, the existing representative-selection
rule still prefers the primary record. If Stage 3 normalizes multiple raw
representations to the same variant-support edge, it keeps the metadata from
the first support record while summing `support_row_count`. Alignment role is
not duplicated in `qc_flags`.

`strategy_summary.tsv.gz` contains `summary_row_count`, `gene_count`,
`aligned_summary_row_count`, `event_count`, and `aligned_target_bp` for each
enabled strategy. Reports must read this alignment-owned table, not mutable
files from `reports/`.

## Why Segments And Events Are Separate

`alignment_events.tsv.gz` records observed differences. It does not record where
an ortholog matches the human target.

`alignment_segments.tsv.gz` records coverage intervals. Stage 2 reduces those
intervals to site-depth tables so Stage 3 can distinguish:

```text
event present at position       -> ortholog supports observed alternative
position covered and no event   -> ortholog supports human/reference allele
position not covered            -> no-call
bad/ambiguous alignment         -> filtered/no-call
```

Both standalone Stage 2 and the end-to-end durable partition dataset keep
segments as normalized alignment evidence. The current Stage 3 compatibility
path consumes the already reduced `snv_site_depth.tsv.gz`; it does not recompute
depth from segments. Both standalone and end-to-end runs derive that table at
each alignment partition boundary. Rows are keyed by gene, strategy, and target
position, so different ALT alleles at one site share the same denominator.
For taxonomy-aware reporting, the same boundary also counts distinct aligned
ortholog, species, genus, family, and order units within supported taxonomic
scopes. Exact-ALT counts use the same unit definitions, so multiple orthologs
from one selected taxon cannot inflate a collapsed unit.

`feature_coverage.tsv.gz` intersects `alignment_segments.tsv.gz` with
`target_features.tsv.gz`. Depth is summed after merging overlapping segments
within each ortholog, so overlapping records from the same ortholog do not
inflate the per-base depth. The implementation uses `bedtools merge`,
`bedtools coverage`, and `bedtools intersect`; intermediate BED files are
task-local temporary data and are removed after the summary is written.

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

Alignment processes use scratch task space. Raw PAF/delta files are temporary and
are removed unless `--keep_native_alignments true` is set.

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

The same preparation step partitions `target_features.tsv.gz` into one compact
`target_features.tsv.gz` inside each ready gene task. Alignment jobs read only
that gene-local feature file when computing coverage; the global feature table
is not staged into every aligner job.

Task preparation also assigns a stable `partition_id` after sorting target genes
by chromosome and genomic interval. `--alignment_partition_size` controls the
maximum number of target genes in each partition (default: 10). The identifier
is carried in task metadata for bounded downstream merge and annotation; it does
not change gene-level alignment behavior.

Alignment results are merged in two bounded levels. Each partition merges at
most `--alignment_partition_size` genes and streams one raw event group at a
time from its SQLite key index into the compact event/support pair. The final
merge streams those disjoint compact partitions through one staged directory
and rebases their local group IDs when a standalone global handoff is requested.
This avoids one global event database and avoids placing every gene/strategy
result path on a single command line while preserving the Stage 2 logical
evidence.

Both merge levels fail closed. A partition must contain exactly one result for
every eligible gene/strategy pair, required TSV inputs must exist with valid
headers, and genes cannot occur in multiple partitions. Ensembl Compara is
eligible when the human target is ready; minimap2, nucmer, and BWA additionally
require fetched ortholog inputs. The final merge compares the union of partition
`gene_ids` with the genes eligible for at least one selected strategy in
`alignment_tasks.tsv.gz`; incomplete output is an error rather than a smaller
successful dataset.

Durable output is limited to compressed normalized TSV files so large runs do
not duplicate sequence data and native aligner output in `results/`.

For precomputed MAF strategies, the source MAF files are streamed or read from a
configured local directory. They are treated as external inputs/cache, not final
pipeline results.

Nextflow `work/` remains a resume cache. After validating a run, it can be
cleaned according to `docs/storage_model.md`.
