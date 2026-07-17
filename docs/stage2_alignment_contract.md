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
- `sequences/targets/<gene_id>.fa.gz`
- `sequences/orthologs/<gene_id>.fa.gz`

Target FASTA records are already in plus genomic orientation on GRCh38.p14.
Ortholog records are aligned as query sequences; aligners decide forward/reverse
orientation.

## Strategies

Current strategies are registered in the workflow and can be selected with
`--alignment_strategies`. The default value is `all`, meaning every registered
strategy, not a hard-coded command-line list.

| Strategy | Tool | Policy |
| --- | --- | --- |
| `minimap2_asm10` | minimap2 | Fixed baseline preset for every ortholog. |
| `minimap2_asm20` | minimap2 | More permissive fixed minimap2 preset. |
| `minimap2_taxonomy_adaptive` | minimap2 | Preset chosen from NCBI taxonomy summary. |
| `nucmer` | MUMmer/nucmer | Independent comparator using multi-query nucmer output. |
| `bwa_pseudoreads` | BWA/samtools/pysam | Pseudoread comparator that maps generated ortholog pseudo-reads to the target and extracts BAM/CIGAR-supported events. |
| `precomputed_ensembl_92_mammals_epo_extended` | Ensembl Compara MAF | Uses release-pinned precomputed `92_mammals.epo_extended` whole-genome MSA blocks overlapping the human target gene interval. |

No LASTZ, consensus calling, or production variant filtering is part of Stage 2.
Conservation scores such as GERP are not part of alignment; they belong to the
later annotation/analysis layer.

Example selections:

```bash
--alignment_strategies all
--alignment_strategies minimap2_asm20
--alignment_strategies minimap2_asm10,nucmer
--alignment_strategies bwa_pseudoreads
--alignment_strategies minimap2_asm20,precomputed_ensembl_92_mammals_epo_extended
```

At least one strategy must be selected. Single-strategy runs are valid; compare
or report layers must treat cross-strategy-only sections as not applicable or
empty rather than failing.

`all` means every registered strategy, including precomputed alignment
strategies. Remote Ensembl MAF chunk tasks default to
`--ensembl_compara_maf_max_forks 3`; increase it only after checking network
stability.

Large runs can enable `--compact_alignment_events true` to publish one support
row per unique event and strategy instead of raw per-ortholog event rows. Raw
remains the default because it preserves maximum traceability.

## Taxonomy Presets

Taxonomy enrichment uses the offline class dictionary published in
`assets/reference/ncbi/taxonomy/taxonomy_classes.json.gz`. The alignment stage maps every selected
ortholog `tax_id` to a compact preset group without making taxonomy network
requests during alignment.

Preset policy:

| Group | Detection | minimap2 preset |
| --- | --- | --- |
| `primates` | ancestor `9443` | `asm10` |
| `other_mammals` | ancestor `40674`, not `9443` | `asm20` |
| `other_vertebrates` | ancestor `7742`, not `40674` | `asm20` |
| `other_or_unknown` | fallback | `asm20` |

`Hominidae` is retained as a metadata flag but does not currently change the
preset. The `asm5` preset is intentionally not used by the default adaptive
policy because it is too strict for broad inter-species ortholog alignment; it
can be added later as a separate benchmark strategy.

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

For taxonomy-adaptive minimap2, the task temporarily splits the ortholog FASTA
into preset groups and runs one minimap2 command per non-empty preset group.
These split FASTA files are scratch-only and are not published.

For nucmer:

```text
target gene FASTA vs multi-FASTA of all selected orthologs for that gene
```

Raw multi-query nucmer is used, but the workflow does not run global
`delta-filter -1`. One-to-one delta filtering is appropriate for comparing two
assemblies, but it is wrong for many orthologs that are all expected to align to
the same human target locus. The parser separates evidence by query sequence and
adds `unfiltered_nucmer` QC flags. `show-snps` rows are published as variant
events only when both non-gap alleles are concrete A/C/G/T bases. Rows containing
IUPAC ambiguity symbols are excluded and counted in the task manifest; the
affected ortholog summary receives `ambiguous_event_allele`.

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
retried with block-level resume: already committed MAF blocks are skipped on the
next attempt. If a source still cannot be fully read after all attempts, the
task records gene-level failure rows. Unexpected process failures terminate the
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
`assets/reference/ensembl/compara/release-<release>/<species_set>/ensembl_compara_maf_manifest.tsv.gz`;
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
| `manifest.json` | Exact ready `gene_ids`, enabled output strategies, and output counts. |
| `alignment_tasks.tsv.gz` | Per-gene task manifest and readiness checks. |
| `taxonomy_presets.tsv.gz` | Compact tax_id-to-preset mapping. |
| `taxonomy_failures.tsv.gz` | Taxonomy lookup warnings/failures. |
| `ortholog_alignment_summary.tsv.gz` | One row per gene/ortholog/output strategy when that strategy emits summary evidence. |
| `strategy_summary.tsv.gz` | Small canonical per-strategy aggregate derived from `ortholog_alignment_summary.tsv.gz` for reports and run inspection. |
| `alignment_segments.tsv.gz` | Normalized alignment intervals. |
| `feature_coverage.tsv.gz` | Per-gene, per-strategy coverage and depth over target structural intervals. |
| `alignment_events.tsv.gz` | Raw mismatch/indel events normalized to target coordinates by default; unique event support rows when `--compact_alignment_events true`. |
| `failures.tsv.gz` | Alignment-stage failures. |
| `native/` | Optional raw PAF/delta/coords/snps files when enabled. |

Native outputs are disabled by default.

In an end-to-end `--stage all` run, Stage 3 consumes partitioned events directly
from Nextflow `work/`. The durable `alignment/` directory therefore contains
only `manifest.json`, `strategy_summary.tsv.gz`, `feature_coverage.tsv.gz`, and
`failures.tsv.gz`. The manifest retains the raw row counts even though raw
events, segments, per-ortholog summaries, task metadata, and taxonomy tables are
not copied into the final run directory.

For Minimap2 rows, `native_record_id` is derived from the PAF record content
rather than its output line number. `event_id` combines the strategy, that
stable record identifier, and the event ordinal within the PAF `cs` tag. These
identifiers therefore remain stable when Minimap2 emits identical records in a
different order at another thread count.

`strategy_summary.tsv.gz` contains `summary_row_count`, `gene_count`,
`aligned_summary_row_count`, `event_count`, and `aligned_target_bp` for each
enabled strategy. Reports must read this alignment-owned table, not mutable
files from `reports/`.

## Why Segments And Events Are Separate

`alignment_events.tsv.gz` records observed differences. It does not record where
an ortholog matches the human target.

`alignment_segments.tsv.gz` records coverage intervals. It is required for the
next stage to distinguish:

```text
event present at position       -> ortholog supports observed alternative
position covered and no event   -> ortholog supports human/reference allele
position not covered            -> no-call
bad/ambiguous alignment         -> filtered/no-call
```

This is why Stage 2 stores both interval evidence and event evidence.

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
not one command-line argument per FASTA. For current Stage 1 outputs it streams
the `query_gene_id`-grouped ortholog table and keeps only one gene's metadata in
memory. Older ungrouped fetch outputs remain supported through a compatibility
path that loads that legacy table in memory; regenerate large legacy fetch
datasets before cluster-scale alignment.

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
most `--alignment_partition_size` genes and, when compact events are requested,
performs support aggregation only inside that partition. The final merge then
streams those disjoint partitions through one staged directory. This avoids one
global event database and avoids placing every gene/strategy result path on a
single command line while preserving the existing Stage 2 output tables.

Both merge levels fail closed. A partition must contain exactly one result for
every expected gene/strategy pair, required TSV inputs must exist with valid
headers, and genes cannot occur in multiple partitions. The final merge compares
the union of partition `gene_ids` with the ready genes in
`alignment_tasks.tsv.gz`; incomplete output is an error rather than a smaller
successful dataset.

Durable output is limited to compressed normalized TSV files so large runs do
not duplicate sequence data and native aligner output in `results/`.

For precomputed MAF strategies, the source MAF files are streamed or read from a
configured local directory. They are treated as external inputs/cache, not final
pipeline results.

Nextflow `work/` remains a resume cache. After validating a run, it can be
cleaned according to `docs/storage_model.md`.
