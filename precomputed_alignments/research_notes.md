# Precomputed Alignment Research Notes

Checked on 2026-06-30 against Ensembl release 116 FTP listings and Ensembl
REST documentation.

## Bottom Line

The production direction should be an offline, release-pinned Ensembl Compara
MAF strategy. REST is useful for validating normalization on tiny intervals, but
it should not be the production data source because it creates a live network
dependency and makes repeatability harder.

This can be built now, but "production MAF" is not just changing the fetcher
from REST to FTP. It needs a downloader/cache, an MAF block extractor, coordinate
normalization, and table-contract checks. The disk requirement is manageable;
the main risk is correctness of interval extraction and coordinate conversion.

## Official Sources Checked

- Ensembl REST genomic alignment endpoint:
  `https://rest.ensembl.org/documentation/info/genomic_alignment_region`
- Ensembl REST Compara methods:
  `https://rest.ensembl.org/documentation/info/compara_methods`
- Ensembl REST Compara species sets:
  `https://rest.ensembl.org/documentation/info/compara_species_sets`
- Ensembl release 116 MAF dumps:
  `https://ftp.ensembl.org/pub/release-116/maf/ensembl-compara/multiple_alignments/`
- Ensembl release 116 conservation scores:
  `https://ftp.ensembl.org/pub/release-116/compara/conservation_scores/`

The top-level FTP index had `release-116/` and no top-level `current_maf/`
alias at the time of inspection. Production runs should therefore pin the
release explicitly.

## Size Estimates

These are compressed sizes from the release-116 FTP directory listings, not
full downloads.

| Dataset | File count including small metadata files | Compressed total | Largest listed file | Notes |
| --- | ---: | ---: | ---: | --- |
| `10_primates.epo` MAF | 132 | 15.2 GiB | 257 MiB | Smallest useful validation set. |
| `44_mammals.epo` MAF | 1204 | 57.8 GiB | 183 MiB | Good production baseline: broad mammals without extended set complexity. |
| `92_mammals.epo_extended` MAF | 1204 | 71.5 GiB | 207 MiB | Most detailed mammal MSA option inspected. Best target if disk is acceptable. |
| `92_mammals.gerp_conservation_score` | 95 | 542.0 GiB | 12.0 GiB | Optional conservation layer, not required for event extraction. Human GRCh38 bigWig alone is 8.9 GiB. |

Practical storage recommendation:

- reserve about 90 GiB for `92_mammals.epo_extended` compressed MAF plus
  manifests and indexes;
- reserve about 70 GiB for `44_mammals.epo`;
- do not decompress whole MAF sets to plain text;
- store MAF outside the repository, for example under a configured cache or
  scratch path;
- publish only normalized GAPH TSV outputs by default.

If we later keep per-gene extracted MSA fragments for debugging, that should be
behind an explicit `keep_native`/debug option because those fragments duplicate
source data.

## Why Not REST For Production

REST is the fastest way to learn the API and validate event normalization. It
worked on a 151 bp region and produced GAPH-like events/segments that the
current compare report accepts.

REST is weak for production:

- it depends on live Ensembl availability during a pipeline run;
- it can be rate-limited;
- release pinning is less explicit than a local release-specific FTP cache;
- repeated pipeline runs would re-query data unless we implement a cache;
- large gene sets would make many small network calls.

REST should remain a test/prototype backend or a fallback for exploratory
single-region work.

## Why FTP MAF Is Not Just A Downloader

The MAF files are large gzip text files split by human chromosome/chunk-like
names. A production implementation must avoid re-scanning huge files separately
for every gene.

The correct design is:

1. Download and checksum a selected release/set once.
2. Build a file/block index over human rows.
3. For a run, group target gene intervals by candidate MAF files.
4. Stream each needed MAF file once and emit only blocks overlapping requested
   target loci.
5. Normalize extracted alignment blocks into the existing Stage 2 table
   contract.

The index does not need to store all MAF sequence text. The first useful index
can be a SQLite table with:

```text
release
species_set
maf_file
human_seq_region
human_start0
human_end0
block_index
block_order
```

If direct random access is needed later, we can add gzip seek indexes or convert
source files to a more random-access-friendly local representation. The first
production version can be streaming and grouped-by-file, which is simpler and
less risky.

## Coordinate And Normalization Risks

Ensembl MAF follows the UCSC MAF coordinate convention:

- `s` line starts are 0-based;
- on the `-` strand, the start is relative to the reverse-complemented source
  sequence;
- Ensembl's tree comments use 1-based Ensembl-style coordinates.

GAPH Stage 2 currently emits target-local coordinates plus GRCh38 RefSeq
accessions from NCBI target metadata. Ensembl MAF human rows use Ensembl
sequence names such as chromosome `4`. Exact cross-strategy comparison therefore
needs a deterministic mapping from Ensembl seq_region names to the same genomic
accession namespace used by the current pipeline.

Event normalization must also handle:

- human row on either strand;
- gaps in human versus gaps in non-human rows;
- ambiguous bases;
- ancestral rows, which should be skipped by default;
- EPO_EXTENDED composite sequences;
- duplicate or fragmented rows for a species;
- insertions represented at the same target anchor across adjacent columns;
- species support semantics, since precomputed whole-genome MSA rows are
  species/genomic intervals, not NCBI ortholog gene records.

## Strategy Semantics

This should be a separate strategy family, not a minimap2 mode.

Candidate strategy names:

```text
precomputed_ensembl_44_mammals_epo
precomputed_ensembl_92_mammals_epo_extended
precomputed_ensembl_10_primates_epo
```

These strategies should emit the same downstream files as local aligners:

```text
alignment_events.tsv.gz
alignment_segments.tsv.gz
ortholog_alignment_summary.tsv.gz
failures.tsv.gz
manifest.json
```

For early comparison reports, `ortholog_gene_id` can be the species name. For
production support/depth semantics, we should decide whether support unit means
species, species interval, or a mapped ortholog gene. Species is probably the
right first definition for whole-genome MSA evidence, but reports must label it
clearly.

## Recommended Next Implementation Step

Build the FTP MAF path inside this standalone package before wiring it into
Nextflow:

1. `download_ensembl_maf_set.py`
   - input: release, species set, output cache dir;
   - reads FTP listing and `MD5SUM`;
   - downloads with resume;
   - verifies checksums;
   - writes a manifest.
2. `index_ensembl_maf.py`
   - streams each `.maf.gz`;
   - records human interval spans and block order in SQLite;
   - does not store full block text.
3. `extract_ensembl_maf_region.py`
   - input: index, region(s), MAF cache;
   - streams only candidate files;
   - emits compact JSONL or MAF fragments for overlapping blocks.
4. `maf_alignment_to_gaph_tables.py`
   - shares normalization code with the REST prototype;
   - writes GAPH Stage 2 table shapes.

After that works for a few known regions, wire it into Nextflow as a separate
alignment strategy and let the existing merge/annotation/compare layers consume
the output.

## Region Test On The Largest MAF Set

A streaming test was run against the release-116
`92_mammals.epo_extended.4_13.maf.gz` remote chunk for human chr4:
`122600000-122700000`.

The extractor did not keep the full source MAF locally. It streamed the remote
gzip file, selected overlapping MAF blocks, and wrote only the extracted JSON
under the ignored `precomputed_alignments/data/` scratch directory.

Observed result:

| Item | Value |
| --- | ---: |
| Source chunk compressed size from FTP listing | 136 MiB |
| Requested interval | 100001 bp |
| Extracted MAF blocks | 4 |
| Human block span after overlap extraction | 122579945-122710746 |
| Species rows in GAPH summaries | 87 |
| Extracted JSON scratch size | 55 MiB |
| Normalized output scratch size | 18 MiB |
| `alignment_events.tsv.gz` rows | 817483 |
| `alignment_segments.tsv.gz` rows | 222376 |
| Event type counts | 541083 SNVs, 163047 deletions, 113353 insertions |

This validates the chunked approach: a single 136 MiB compressed source chunk
can be streamed and discarded, leaving only normalized evidence for the target
region. For whole-genome processing, the durable output should be normalized
TSV/SQLite/Parquet evidence, not retained source chunks.

## GERP Region Test

The matching human GERP scores were fetched from the release-116 92-mammals
human GRCh38 bigWig by HTTP range access with `pyBigWig`.

Input:

```text
https://ftp.ensembl.org/pub/release-116/compara/conservation_scores/92_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw
chr4:122600000-122700000
```

Observed result:

| Item | Value |
| --- | ---: |
| Full human bigWig size from FTP listing | 8.9 GiB |
| Requested interval | 100001 bp |
| Values returned | 95728 |
| Missing bases | 4273 |
| Intervals written | 72525 |
| Output TSV gzip size | 472 KiB |
| Min score | -5.36 |
| Max score | 2.68 |
| Mean score | -0.902747 |

This means GERP does not require downloading the 8.9 GiB human bigWig for each
run. The production path can query exactly the target intervals, write compact
regional score tables, and leave the full 542 GiB cross-species conservation
directory untouched.

## Current Recommendation

If disk space is acceptable, start with `92_mammals.epo_extended` because it is
the most detailed mammal alignment set inspected and is only about 71.5 GiB
compressed. Also keep `44_mammals.epo` in mind as a cleaner baseline: it is
smaller and likely easier to reason about when validating the parser.

Do not start by downloading all GERP bigWigs. The full inspected GERP directory
is about 542 GiB. The practical route is region access over the human GRCh38
bigWig, which supports HTTP range reads with `pyBigWig`; the 100 kb chr4 test
above fetched values without downloading the 8.9 GiB file. Add GERP as a
chunked annotation/QC layer tied to the same target intervals as the extracted
MAF blocks.
