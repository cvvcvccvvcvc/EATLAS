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
- `orthologs.selected.tsv.gz`
- `sequences/targets/<gene_id>.fa.gz`
- `sequences/orthologs/<gene_id>.fa.gz`

Target FASTA records are already in plus genomic orientation on GRCh38.p14.
Ortholog records are aligned as query sequences; aligners decide forward/reverse
orientation.

## Strategies

Current strategies are intentionally limited to the baseline benchmark set:

| Strategy | Tool | Policy |
| --- | --- | --- |
| `minimap2_asm10` | minimap2 | Fixed baseline preset for every ortholog. |
| `minimap2_taxonomy_adaptive` | minimap2 | Preset chosen from NCBI taxonomy summary. |
| `nucmer` | MUMmer/nucmer | Independent comparator using multi-query nucmer output. |

No LASTZ, consensus calling, or production variant filtering is part of Stage 2.

## Taxonomy Presets

Taxonomy enrichment uses batch NCBI Datasets taxonomy summary calls:

```bash
datasets summary taxonomy taxon <tax_ids...> --as-json-lines
```

The workflow does not use `--parents` because the normal taxonomy summary row
already includes a compact `parents` list and `classification` block for each
requested taxon. `--parents` expands every ancestor into separate heavy JSONL
records and only supports one taxon per request.

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

## Durable Outputs

Stage 2 publishes:

| Path | Meaning |
| --- | --- |
| `manifest.json` | Alignment run counts, strategies, and output counts. |
| `alignment_tasks.tsv.gz` | Per-gene task manifest and readiness checks. |
| `taxonomy_presets.tsv.gz` | Compact tax_id-to-preset mapping. |
| `taxonomy_failures.tsv.gz` | Taxonomy lookup warnings/failures. |
| `ortholog_alignment_summary.tsv.gz` | One row per gene/ortholog/strategy. |
| `alignment_segments.tsv.gz` | Normalized alignment intervals. |
| `alignment_events.tsv.gz` | Raw mismatch/indel events normalized to target coordinates. |
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

Nextflow `work/` remains a resume cache. After validating a run, it can be
cleaned according to `docs/storage_model.md`.
