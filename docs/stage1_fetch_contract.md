# Stage 1 Fetch Contract

Stage 1 turns a list of Entrez Gene IDs into normalized target and ortholog
sequence data for downstream alignment and annotation.

## Fixed Constants

- Target organism: Homo sapiens (`tax_id=9606`)
- Target assembly: GRCh38.p14 (`GCF_000001405.40`)
- Ortholog scope: NCBI Datasets complete ortholog set (`--ortholog all`)
- NCBI Datasets include mode: `--include gene`

These constants are intentionally not user-facing parameters in the current
design.

## Input

Required:
- `--ids_file`: text file with Entrez Gene IDs, separated by whitespace or commas.

Optional operational parameters:
- `--outdir`: final output directory.
- `--chunk_size`: accepted IDs per NCBI package request.
- `--fetch_max_forks`: max concurrent NCBI fetch/parse tasks. Default is 4.
- `ENTREZ_API_KEY` or `NCBI_API_KEY`: optional NCBI API key passed to
  `datasets download` as `--api-key`.
- `ENTREZ_EMAIL` or `NCBI_EMAIL`: optional contact email recorded as configured
  metadata; NCBI Datasets CLI does not expose an email flag for this command.
- `--target_annotation_gff3`: local NCBI RefSeq GFF3 for
  `GCF_000001405.40`; defaults to `GAPH_TARGET_ANNOTATION_GFF3`, then
  `assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz`.

The fetch implementation fixes its NCBI request policy: request starts are at
least 5 seconds apart, each download gets up to 4 retries after its initial
attempt, and retry backoff starts at 30 seconds. The NCBI Datasets executable
comes from the mandatory task environment.

## Processing Steps

1. `VALIDATE_IDS`
   - Rejects non-positive/non-integer IDs.
   - Keeps first occurrence of each Entrez ID.
   - Records duplicate input rows in `input.ids.tsv`.
   - Splits accepted IDs into deterministic chunk files.

2. `FETCH_PARSE_CHUNK`
   - Downloads one NCBI gene package per chunk:

     ```bash
     datasets download gene gene-id --inputfile <chunk> --ortholog all --include gene
     ```

   - Retries transient or invalid package downloads with exponential backoff.
   - If the chunk package still fails, retries every requested gene as a
     sequential singleton download.
   - Records a terminal singleton download/package failure in
     `failures.tsv.gz` and preserves successful genes from the same chunk.
   - Uses `data_report.jsonl` to map every ortholog GeneID back to the requested
     query GeneID via `geneGroups[].id`.
   - Uses `gene.fna` as the source of genomic gene sequences.
   - Writes per-chunk status, download mode, attempt counts, and timings into
     chunk `manifest.json`.

3. Target selection
   - Selects the requested human GeneID.
   - Requires a genomic location on `GCF_000001405.40`.
   - Writes the target FASTA in plus genomic orientation.
   - Stores the biological gene orientation separately in metadata.

4. Ortholog selection
   - Excludes human records from ortholog outputs.
   - Groups candidates by `(query_gene_id, ortholog_gene_id)`.
   - Selects one sequence record per ortholog GeneID using deterministic priority:
     1. earlier matching annotation/location index
     2. `NC_` accession before `NW_`
     3. longer sequence
     4. lexical accession/range tie-break
   - Writes rejected candidates as metadata rows only, without sequences.

5. `BUILD_FETCH_DATASET`
   - Merges chunk tables.
   - Requires every chunk table to exist and have the same explicit schema;
     malformed or headerless tables fail the run.
   - Copies final per-gene FASTA files.
   - Requires every accepted GeneID to have exactly one terminal outcome:
     a target gene row or a gene-level failure.
   - Requires one manifest for every planned chunk and rejects duplicate or
     missing chunk results.
   - Requires at least one target gene and one selected ortholog because all
     supported downstream alignment strategies need both.
   - Builds compact target structural features from the configured local target assembly GFF3.
   - Emits `chunk_metrics.tsv.gz` in task work for fetch diagnostics.
   - Writes final `manifest.json`.

6. `FETCH_TAXONOMY`
   - Reads the unique `tax_id` values from `orthologs.selected.tsv.gz`.
   - Sends one logical batch request to NCBI Datasets taxonomy summary and
     normalizes the response into one row per selected tax ID.
   - Writes missing per-taxon responses explicitly to
     `taxonomy_failures.tsv.gz`; a failed batch request is retried by Nextflow
     and prevents finalization after retries are exhausted.

Taxonomy is acquired only in the fetch boundary. Alignment and annotation do
not consume it or issue taxonomy requests; analytics joins it to durable
alignment evidence when taxonomic views are requested.

## Final Output Files

There is no standalone fetch output contract. Ortholog FASTA, chunk tables, and
the candidate table remain in Nextflow `work/` only as long as downstream
alignment or `-resume` needs them. The durable `fetch/` directory keeps `manifest.json`,
`input.ids.tsv`, `genes.tsv.gz`, `target_features.tsv.gz`,
`orthologs.selected.tsv.gz`, taxonomy metadata, `failures.tsv.gz`, and target
FASTA files.

| Path | Meaning |
| --- | --- |
| `manifest.json` | Run status, constants, counts, download fallback metrics, and Datasets CLI version(s). |
| `input.ids.tsv` | All input rows, accepted status, duplicate mapping. |
| `genes.tsv.gz` | Target human gene metadata and sequence checksum. |
| `target_features.tsv.gz` | Collapsed target-local structural intervals: gene, exon, CDS, UTR, intron. |
| `orthologs.selected.tsv.gz` | Metadata for selected ortholog sequences. |
| `taxonomy.tsv.gz` | Canonical tax_id-to-lineage and named-rank metadata for selected ortholog taxa. |
| `taxonomy_failures.tsv.gz` | Selected tax IDs absent from the NCBI taxonomy response. |
| `failures.tsv.gz` | Gene-level failures. |
| `sequences/targets/<gene_id>.fa.gz` | Target human genomic sequence. |

`orthologs.selected.tsv.gz` is grouped by `query_gene_id`. This ordering is part
of the Stage 1 contract: Stage 2 can prepare one gene at a time without loading
all ortholog metadata into memory.

`taxonomy.tsv.gz` is one gzip-compressed wide table with one row per unique
selected `tax_id`. Its stable columns are:

| Columns | Meaning |
| --- | --- |
| `tax_id`, `taxonomy_status` | Requested taxon and whether NCBI returned it (`resolved` or `not_returned`). |
| `scientific_name`, `rank`, `group_name` | Direct NCBI summary values for the requested taxon. |
| `domain_id`, `domain_name` through `species_id`, `species_name` | Direct NCBI classification at domain, kingdom, phylum, class, order, family, genus, and species ranks. Missing ranks remain empty. |
| `lineage_tax_ids` | Ordered, comma-separated NCBI lineage from root to the requested taxon, including the requested tax ID for resolved rows. It is empty when the response is absent. |

This source contract deliberately does not store precomputed taxonomic
membership flags or a second long classification table. Membership and scope
counts are inexpensive derivations from `lineage_tax_ids`; keeping a single
wide row avoids duplicate taxonomy facts and repeated joins. Consumers require
this canonical schema and fail when a selected tax ID does not have exactly one
taxonomy row.
Analytics derives and caches `taxonomy_summary.tsv.gz` from this table and
`orthologs.selected.tsv.gz`; it is not part of the Stage 1 handoff.

## Strand Convention

Target FASTA is written in plus genomic orientation. For genes annotated on the
minus strand, the NCBI complement record is reverse-complemented back to plus
genomic sequence.

Metadata still records:
- `orientation`: gene orientation from NCBI
- `sequence_orientation`: `plus`

This keeps future VCF coordinates and REF/ALT alleles in the normal genomic
coordinate system while retaining gene-strand information for interpretation.

Ortholog FASTA records are currently written as provided by NCBI Datasets.
Downstream alignment is expected to handle forward/reverse mapping.

## Target Features

`target_features.tsv.gz` uses target-local 0-based half-open coordinates plus
GRCh38 1-based inclusive coordinates. Exon, CDS, and UTR intervals are collapsed
across transcripts to avoid alternative-transcript double counting. Introns are
the gene body minus the collapsed exon union.
