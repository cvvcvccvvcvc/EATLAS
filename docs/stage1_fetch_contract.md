# Stage 1 Fetch Contract

Stage 1 turns a list of Entrez Gene IDs into normalized target and ortholog
sequence data for downstream pseudo-read generation and variant calling.

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
- `--fetch_max_forks`: max concurrent NCBI fetch/parse tasks. Default is 2.
- `--fetch_request_stagger_seconds`: minimum spacing between starts of NCBI
  Datasets download requests across concurrent local fetch tasks. Default is
  5 seconds.
- `--fetch_download_retries`: in-process retries for each NCBI Datasets
  download, including chunk and singleton fallback requests. Default is 4.
- `--fetch_download_retry_base_seconds`: base exponential backoff interval for
  NCBI Datasets download retries. Default is 30 seconds.
- `--datasets_bin`: path/name for the NCBI Datasets CLI. Defaults to
  `DATASETS_BIN`, then `tools/bin/datasets` when present, otherwise `datasets`
  on `PATH`.
- `ENTREZ_API_KEY` or `NCBI_API_KEY`: optional NCBI API key passed to
  `datasets download` as `--api-key`.
- `ENTREZ_EMAIL` or `NCBI_EMAIL`: optional contact email recorded as configured
  metadata; NCBI Datasets CLI does not expose an email flag for this command.
- `--target_annotation_gff3`: local NCBI RefSeq GFF3 for
  `GCF_000001405.40`; defaults to `GAPH_TARGET_ANNOTATION_GFF3`, then
  `assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz`.

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
   - Copies final per-gene FASTA files.
   - Requires every accepted GeneID to have exactly one terminal outcome:
     a target gene row or a gene-level failure.
   - Requires one manifest for every planned chunk and rejects duplicate or
     missing chunk results.
   - Builds compact target structural features from the configured local target assembly GFF3.
   - Writes `chunk_metrics.tsv.gz` with durable per-chunk timing and package-size
     metrics.
   - Writes final `manifest.json`.

## Final Output Files

The table below is the full standalone `--stage fetch` handoff contract. In an
end-to-end `--stage all` run, the ortholog FASTA, chunk tables, and candidate
table remain in Nextflow `work/` only as long as downstream alignment needs
them. The durable `fetch/` directory keeps `manifest.json`,
`input.ids.tsv`, `genes.tsv.gz`, `target_features.tsv.gz`,
`orthologs.selected.tsv.gz`, `failures.tsv.gz`, and target FASTA files.

| Path | Meaning |
| --- | --- |
| `manifest.json` | Run status, constants, counts, download fallback metrics, and Datasets CLI version(s). |
| `input.ids.tsv` | All input rows, accepted status, duplicate mapping. |
| `chunks.tsv` | Chunk IDs and accepted Gene IDs assigned to each chunk. |
| `genes.tsv.gz` | Target human gene metadata and sequence checksum. |
| `target_features.tsv.gz` | Collapsed target-local structural intervals: gene, exon, CDS, UTR, intron. |
| `orthologs.selected.tsv.gz` | Metadata for selected ortholog sequences. |
| `orthologs.candidates.tsv.gz` | Non-human ortholog candidate records and reject reasons. |
| `failures.tsv.gz` | Gene-level failures. |
| `sequences/targets/<gene_id>.fa.gz` | Target human genomic sequence. |
| `sequences/orthologs/<gene_id>.fa.gz` | Selected ortholog genomic sequences. |

`orthologs.selected.tsv.gz` is grouped by `query_gene_id`. This ordering is part
of the Stage 1 contract: Stage 2 can prepare one gene at a time without loading
all ortholog metadata into memory. `manifest.json` records this guarantee as
`orthologs_selected_grouped_by_query_gene_id=true`.

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
