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
- `--fetch_max_forks`: max concurrent NCBI fetch/parse tasks. Default is 5.
- `--fetch_request_stagger_seconds`: minimum spacing between starts of NCBI
  Datasets download requests across concurrent local fetch tasks. Default is
  5 seconds.
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

   - Uses `data_report.jsonl` to map every ortholog GeneID back to the requested
     query GeneID via `geneGroups[].id`.
   - Uses `gene.fna` as the source of genomic gene sequences.
   - Writes per-chunk timings into chunk `manifest.json`.

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
   - Builds compact target structural features from the configured local target assembly GFF3.
   - Writes `chunk_metrics.tsv.gz` with durable per-chunk timing and package-size
     metrics.
   - Writes final `manifest.json`.

## Final Output Files

| Path | Meaning |
| --- | --- |
| `manifest.json` | Run constants, counts, and Datasets CLI version(s). |
| `input.ids.tsv` | All input rows, accepted status, duplicate mapping. |
| `chunks.tsv` | Chunk IDs and accepted Gene IDs assigned to each chunk. |
| `genes.tsv.gz` | Target human gene metadata and sequence checksum. |
| `target_features.tsv.gz` | Collapsed target-local structural intervals: gene, exon, CDS, UTR, intron. |
| `orthologs.selected.tsv.gz` | Metadata for selected ortholog sequences. |
| `orthologs.candidates.tsv.gz` | Non-human ortholog candidate records and reject reasons. |
| `failures.tsv.gz` | Gene-level failures. |
| `sequences/targets/<gene_id>.fa.gz` | Target human genomic sequence. |
| `sequences/orthologs/<gene_id>.fa.gz` | Selected ortholog genomic sequences. |

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
