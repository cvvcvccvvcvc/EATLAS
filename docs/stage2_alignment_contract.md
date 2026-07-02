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
| `bwa_pseudoreads_varscan` | BWA/samtools/VarScan | Uses the same BWA pseudoread alignment preparation, then calls variants with VarScan. This strategy is selected independently and is the only BWA mode that requires VarScan. |
| `precomputed_ensembl_92_mammals_epo_extended` | Ensembl Compara MAF | Uses release-pinned precomputed `92_mammals.epo_extended` whole-genome MSA blocks overlapping the human target gene interval. |

No LASTZ, consensus calling, or production variant filtering is part of Stage 2.
Conservation scores such as GERP are not part of alignment; they belong to the
later annotation/analysis layer.

Example selections:

```bash
--alignment_strategies all
--alignment_strategies minimap2_asm20
--alignment_strategies minimap2_asm10,nucmer
--alignment_strategies bwa_pseudoreads,bwa_pseudoreads_varscan
--alignment_strategies minimap2_asm20,precomputed_ensembl_92_mammals_epo_extended
```

At least one strategy must be selected. Single-strategy runs are valid; compare
or report layers must treat cross-strategy-only sections as not applicable or
empty rather than failing.

`all` means every registered strategy, including precomputed alignment
strategies. When using remote Ensembl FTP sources, keep
`--ensembl_compara_maf_max_forks` conservative.

Large runs can enable `--compact_alignment_events true` to publish unique event
support rows instead of raw per-support event rows. Raw remains the default
because it preserves maximum traceability.

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
adds `unfiltered_nucmer` QC flags.

For Ensembl Compara MAF:

```text
precomputed whole-genome MSA blocks overlapping the human target gene interval
```

The workflow first builds a small run-specific MAF chunk manifest for the human
chromosomes present in `genes.tsv.gz`. One alignment task then reads the
candidate MAF chunks for one gene and clips evidence to the target interval
before writing normalized Stage 2 tables. Remote MAF chunks are downloaded into
a local cache with atomic writes and gzip validation before parsing. Full MAF
chunks are not published as durable outputs.

The MAF chunk manifest can be supplied with `--ensembl_compara_maf_manifest` or
`ENSEMBL_COMPARA_MAF_MANIFEST`. If neither is set, the workflow checks
`assets/reference/ensembl/compara/release-<release>/<species_set>/ensembl_compara_maf_manifest.tsv.gz`;
if that file is absent, it builds the manifest during the run.

Remote MAF chunks are cached under
`--ensembl_compara_maf_cache_dir`, defaulting to
`assets/reference/ensembl/compara/release-<release>/<species_set>/maf_files`.
The cache prevents repeated network reads of the same chunk, but it can grow by
hundreds of MB per human chromosome region touched by a run.

This strategy is not based on NCBI ortholog GeneIDs. Its support unit is the
species row in the precomputed MSA. Consequently, `ortholog_gene_id` contains
the species name for this strategy, and support-count reports should interpret
it as species support.

## Durable Outputs

Stage 2 publishes:

| Path | Meaning |
| --- | --- |
| `manifest.json` | Alignment run counts, enabled output strategies, and output counts. |
| `alignment_tasks.tsv.gz` | Per-gene task manifest and readiness checks. |
| `taxonomy_presets.tsv.gz` | Compact tax_id-to-preset mapping. |
| `taxonomy_failures.tsv.gz` | Taxonomy lookup warnings/failures. |
| `ortholog_alignment_summary.tsv.gz` | One row per gene/ortholog/output strategy when that strategy emits summary evidence. |
| `alignment_segments.tsv.gz` | Normalized alignment intervals. |
| `feature_coverage.tsv.gz` | Per-gene, per-strategy coverage and depth over target structural intervals. |
| `alignment_events.tsv.gz` | Raw mismatch/indel events normalized to target coordinates by default; unique event support rows when `--compact_alignment_events true`. |
| `failures.tsv.gz` | Alignment-stage failures. |
| `native/` | Optional raw PAF/delta/coords/snps files when enabled. |

Native outputs are disabled by default.

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
inflate the per-base depth.

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

Durable output is limited to compressed normalized TSV files so large runs do
not duplicate sequence data and native aligner output in `results/`.

For precomputed MAF strategies, the source MAF files are cached/read from a
configured local directory. They are treated as external inputs/cache, not final
pipeline results.

Nextflow `work/` remains a resume cache. After validating a run, it can be
cleaned according to `docs/storage_model.md`.
