# Precomputed Cross-Species Alignments

Research package for a future GAPH alignment strategy that consumes already
computed cross-species genome alignments instead of aligning every ortholog gene
sequence locally.

This directory is intentionally separate from the production Nextflow workflow.
It is for source evaluation, format experiments, and small reproducible
prototypes.

Detailed FTP MAF conclusions, size estimates, and production design notes are
in `research_notes.md`.

## Initial Source Choice

The best first target is Ensembl Compara:

- REST genomic alignment endpoint:
  `https://rest.ensembl.org/documentation/info/genomic_alignment_region`
- REST method discovery:
  `https://rest.ensembl.org/documentation/info/compara_methods`
- REST species-set discovery:
  `https://rest.ensembl.org/documentation/info/compara_species_sets`
- FTP MAF dumps:
  `https://ftp.ensembl.org/pub/release-116/maf/ensembl-compara/multiple_alignments/`
- FTP conservation scores:
  `https://ftp.ensembl.org/pub/release-116/compara/conservation_scores/`

As of the checked FTP index, `release-116` has no top-level `current_maf/`
alias, so production code should pin a release explicitly instead of relying on
an implicit moving target. Useful release-116 MAF sets include
`10_primates.epo`, `44_mammals.epo`, and `92_mammals.epo_extended`.
GERP conservation scores are a downstream annotation/QC resource, not an
alignment strategy input.

Observed useful Compara methods:

| Method | Meaning | Notes |
| --- | --- | --- |
| `EPO` | Enredo-Pecan-Ortheus ancestral/tree alignment | Good first choice for mammals/primates. |
| `EPO_EXTENDED` | Extended EPO tree alignment | Larger species sets, potentially noisier. |
| `PECAN` | Pecan multiple alignment | Alternative MSA source. |
| `CACTUS_HAL` / `CACTUS_DB` | Cactus whole-genome alignment | Worth revisiting if HAL dumps or APIs are practical. |
| `GERP_CONSERVATION_SCORE` | Conservation score track | Downstream annotation/QC evidence, not alignment input. |

REST is the right MVP path because it can fetch a single human genomic interval.
FTP MAF is the scalable production path, but files are chromosome/chunk-sized and
need an interval index before they are ergonomic for gene-level queries.

## Prototype Commands

Fetch a small MSA around a human locus:

```bash
python3 precomputed_alignments/scripts/fetch_ensembl_rest_alignment.py \
  --region 4:122612500-122612650 \
  --method EPO \
  --species-set-group mammals \
  --output precomputed_alignments/data/ensembl_epo_4_122612500_122612650.json
```

Convert the REST JSON into GAPH-like alignment tables:

```bash
python3 precomputed_alignments/scripts/rest_alignment_to_gaph_tables.py \
  --input precomputed_alignments/data/ensembl_epo_4_122612500_122612650.json \
  --outdir precomputed_alignments/outputs/epo_region_4 \
  --gene-id demo_region_4 \
  --strategy ensembl_compara_rest_epo \
  --method EPO \
  --species-set-group mammals
```

Run the existing strategy report on the produced events:

```bash
python3 scripts/compare_strategies.py \
  --events-tsv precomputed_alignments/outputs/epo_region_4/alignment_events.tsv.gz \
  --out-html /private/tmp/gaph_precomputed_alignment_demo.html
```

The prototype writes `data/` and `outputs/` as ignored scratch directories.
The report command above writes outside the repository for the same reason.

Extract a region directly from a release-pinned MAF chunk without keeping the
full `.maf.gz` file locally:

```bash
python3 precomputed_alignments/scripts/extract_ensembl_maf_region.py \
  --maf https://ftp.ensembl.org/pub/release-116/maf/ensembl-compara/multiple_alignments/92_mammals.epo_extended/92_mammals.epo_extended.4_14.maf.gz \
  --human-src homo_sapiens.4 \
  --start1 122600000 \
  --end1 122700000 \
  --output precomputed_alignments/data/maf_92_mammals_chr4_region.json
```


## Current Compatibility

The converter emits the same core table names as Stage 2:

- `alignment_events.tsv.gz`
- `alignment_segments.tsv.gz`
- `ortholog_alignment_summary.tsv.gz`
- `failures.tsv.gz`
- `manifest.json`

For REST MSA rows, `ortholog_gene_id` is currently the species name. That is
enough for per-species support counts and compare-report experiments, but not
yet equivalent to the current gene-level ortholog IDs. A production integration
should map Ensembl species/genomic intervals back to selected ortholog metadata
when that mapping is required.

Coordinates are emitted on the human Ensembl `seq_region` by default. Existing
GAPH Stage 2 outputs use GRCh38 RefSeq accessions from NCBI. To compare variants
across strategies exactly, a production implementation must normalize chromosome
names/accessions through the Stage 1 target metadata or an assembly report.
The prototype has `--genomic-accession` for one-off overrides, but the real
pipeline needs a deterministic mapping table.

## Proposed Production Shape

Keep this as a distinct strategy family, not a wrapper inside minimap2:

```text
FETCH_PRECOMPUTED_ALIGNMENT_INDEX
  -> EXTRACT_PRECOMPUTED_GENE_MSA
  -> NORMALIZE_PRECOMPUTED_EVENTS
  -> existing MERGE_ALIGNMENT
  -> existing ANNOTATE_EVENTS / compare report
```

The Stage 2 strategy registry would eventually register something like:

```text
precomputed_ensembl_epo_mammals
precomputed_ensembl_epo_primates
precomputed_ensembl_epo_extended_mammals
```

This keeps local aligner strategies and precomputed-MSA strategies independent
while preserving the same downstream table contract.

## REST Versus FTP Design

REST advantages:

- simplest MVP;
- no large downloads;
- exact gene/region query;
- good for validating normalization logic.

REST limitations:

- network dependency during pipeline execution;
- rate limits and reproducibility concerns;
- less control over release pinning unless the server/release is explicit.

FTP MAF advantages:

- release-pinned, reproducible source;
- suitable for offline production runs;
- same source can feed many genes.

FTP MAF limitations:

- files are large and split by human chromosome chunks;
- direct random access requires an index;
- MAF blocks need careful filtering to the human target interval.

Recommended path:

1. Validate event normalization with REST on small intervals.
2. Add an offline MAF index builder that records human `seq_region`,
   `start0`, `end0`, file path, and byte/compressed offset if practical.
3. Extract per-gene MSA blocks from the indexed MAF source.
4. Reuse the same normalization code for REST and MAF block records.
5. Integrate into Nextflow only after the standalone package can produce
   stable `alignment_events.tsv.gz` and `alignment_segments.tsv.gz`.
