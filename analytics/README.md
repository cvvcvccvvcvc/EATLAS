# GAPH analytics

This package contains analysis and reporting entrypoints for completed GAPH
runs.

## Package structure

```text
analytics/
  strategy_report.py  # CLI contract and orchestration
  vep_annotation.py   # resumable bulk-VEP command
  annotation/         # VEP integration and consequence semantics
  analyses/           # scientific calculations and aggregation
  io/                 # run inputs and atomic artifact contracts
  reporting/          # HTML sections and document composition
```

Cross-stage variant, ClinVar, and gnomAD rules live in the top-level
`genomics/` package. Analytics modules must not be imported by pipeline
commands.

Create the reproducible analytics environment once with:

```bash
micromamba create -f envs/analytics.yml
```

The same environment runs the complete Python test suite:

```bash
micromamba run -n gaph-v2-analytics python -m pytest -q
```

Primary report:

```bash
RUN="results/run_all_strategies_20260703_135905"

micromamba run -n gaph-v2-analytics python -m analytics.strategy_report \
  --run-dir "$RUN"
```

The report presents unadjusted, fixed-band, and continuous phyloP100way ClinVar
association modes in one view. The consequence-matched target-space null is
an explicit opt-in because it uses Ensembl VEP and the gnomAD GraphQL API and
can take hours:

```bash
micromamba run -n gaph-v2-analytics python -m analytics.strategy_report \
  --run-dir "$RUN" \
  --target-space-null \
  --target-space-null-sample-size 5000
```

REST remains the default VEP backend for small runs. Large cluster runs can use
a pinned local RefSeq cache without changing the annotation schema:

```bash
python -m analytics.strategy_report \
  --run-dir "$RUN" \
  --target-space-null \
  --vep-backend local \
  --vep-release 116 \
  --vep-executable "$GAPH_VEP_EXECUTABLE" \
  --vep-cache-dir "$GAPH_VEP_CACHE_DIR" \
  --vep-result-cache-dir "$GAPH_VEP_RESULT_CACHE_DIR" \
  --vep-forks 4
```

The executable may be a wrapper around a containerized VEP. It must expose the
ordinary VEP CLI and be able to read the cache and run analytics paths. Local
VEP uses GRCh38, RefSeq transcripts, `pick_allele_gene`, and the uploaded
normalized REF/ALT alleles; basic consequence annotation does not require a
reference FASTA.

`GAPH_VEP_CACHE_DIR` is the official indexed RefSeq reference cache used by the
VEP executable. `GAPH_VEP_RESULT_CACHE_DIR` is a separate sparse cross-run cache
of completed variant/gene results. The latter stores immutable regional Parquet
fragments with Zstandard compression, namespaces them by the complete VEP
configuration, and never retains transient `no_response` results. Matched
control still keeps its per-run SQLite database as its resume checkpoint.
When `GAPH_VEP_RESULT_CACHE_DIR` is unset but `GAPH_ROOT` is available, reports
default to `$GAPH_ROOT/cache/vep_results`.

Full candidate annotation is a separate resumable precompute, not an implicit
part of HTML generation. Prepare deterministic input partitions once:

```bash
python -m analytics.vep_annotation prepare \
  --run-dir "$RUN" \
  --partition-size 250000
```

Run each one-based partition index independently, normally through a bounded
Slurm array:

```bash
python -m analytics.vep_annotation annotate \
  --run-dir "$RUN" \
  --partition-index 1 \
  --vep-backend local \
  --vep-release 116 \
  --vep-executable "$GAPH_VEP_EXECUTABLE" \
  --vep-cache-dir "$GAPH_VEP_CACHE_DIR" \
  --vep-result-cache-dir "$GAPH_VEP_RESULT_CACHE_DIR" \
  --vep-forks 4
```

After every partition is complete, validate and join the gzip members without
recompressing them:

```bash
python -m analytics.vep_annotation finalize --run-dir "$RUN"
```

An existing finalized artifact can seed the cross-run result cache without
running VEP again:

```bash
python -m analytics.vep_annotation seed-cache \
  --run-dir "$RUN" \
  --vep-result-cache-dir "$GAPH_VEP_RESULT_CACHE_DIR"
```

The final table is
`<run-dir>/analytics/vep_consequences/variant_annotations.vep.tsv.gz`. It keeps
the source annotation columns and adds `vep_*` fields. Prepared inputs,
headerless output partitions, and their manifests remain in the same directory
for resume and audit. A changed source, partition contract, or VEP runtime
configuration is rejected instead of being silently mixed with old results.
When this finalized artifact matches the current pipeline annotation file, the
report selects it automatically. Candidate consequence plots then use primary
RefSeq VEP consequences for every successfully annotated candidate; raw
`gnomad_csq` remains available as provenance. Runs without the artifact retain
the legacy gnomAD-CSQ plots.

The same report invocation annotates the much smaller normalized ClinVar
validation universe with the configured VEP release and caches it as
`<run-dir>/analytics/clinvar_universe.snv_indel.vep.tsv.gz`. ClinVar Association
then uses those RefSeq VEP terms for consequence subsets while retaining raw
ClinVar MC fields in the artifact. A matching source and VEP release reuse the
file without another VEP process.

Set `GAPH_GNOMAD_CACHE_DIR` to the same shared path used by pipeline annotation,
or pass `--gnomad-cache-dir`, so new matched-control reports reuse complete
regional responses instead of requesting them again.

The ClinVar association view compares all strategies and supports variant-type
and consequence selectors. With a completed bulk-VEP artifact these selectors
use RefSeq VEP; legacy runs use ClinVar MC. A second selector exposes the 2x2,
fixed-band, or continuous-distribution data for one strategy at a time.

`Ortholog Evidence` shows SNV-only heatmaps of site-aligned evidence-unit count
and absolute exact-ALT support for CDS, UTR, and intron contexts. The controls
select strategy, taxonomic scope, evidence unit (ortholog, species, genus,
family, or order), and median/quartile/decile grouping. Boundaries are computed
independently within each context. Cell color is the gnomAD found fraction after
excluding failed lookups. Two weighted empirical cumulative distributions below
the heatmaps summarize site-aligned and exact-ALT evidence units across the
eligible CDS, UTR, and intron SNVs; these distributions include failed gnomAD
lookups because gnomAD is not an evidence-depth eligibility criterion. Ensembl
Compara MAF is marked unavailable because its species rows do not carry the NCBI
taxonomy identifiers required by this calculation. Older runs remain reportable
through the legacy all-ortholog view when site depth is present; otherwise the
section is marked unavailable.

When an analysis needs durable intermediate tables, write them under
`<run-dir>/analytics/`. The source tree does not keep a default scratch/work
directory.

Each report invocation progressively writes
`<run-dir>/analytics/performance/<report-stem>.json`. The profile records
top-level and nested wall/CPU timing, process and child peak RSS, block-I/O
counters, cache details, known temporary-file sizes, and the net size change of
the run-local analytics directory. It is updated atomically after every
completed or failed stage, so a long or interrupted cluster job still leaves
useful diagnostics.

The report reads `annotation/variant_annotations.tsv.gz` in chunks. It uses a
temporary SQLite file under `<run-dir>/analytics/` to deduplicate
variant-strategy records without loading the full annotation table into memory;
the file is removed when aggregation finishes or fails. The compact final
aggregation is cached as `<run-dir>/analytics/variant_summary.json.gz` and is
reused while the annotation input and summary schema remain unchanged.
Current runs load ortholog-evidence heatmaps from the compact
`annotation/ortholog_evidence_summary.tsv.gz`; the report reconstructs those
aggregates from `variant_strategy_support.tsv.gz` only for legacy runs that do
not publish the compact table.
The report requires the canonical `alignment/strategy_summary.tsv.gz`; it does
not reconstruct that aggregate from a raw per-ortholog table.

The strategy report writes its ClinVar validation universe under:

```text
<run-dir>/analytics/clinvar_universe.snv_indel.tsv.gz
<run-dir>/analytics/clinvar_universe.snv_indel.manifest.json
<run-dir>/analytics/clinvar_target_regions.bed
```

Validation statistics are computed separately for SNV and INDEL rows. The
conservation-adjusted blocks also write:

```text
<run-dir>/analytics/clinvar_universe.snv_indel.conservation.tsv.gz
<run-dir>/analytics/clinvar_universe.snv_indel.conservation.manifest.json
```

The conservation cache is allele-level and is reused when the ClinVar universe
and phyloP100way source are unchanged. SNVs use the substituted base, deletions
use the mean across deleted reference bases without the VCF padding base, and
insertions use the mean of their two flanking bases. All required bases must be
scored.

Candidate-wide phyloP is summarized separately under:

```text
<run-dir>/analytics/candidate_variants.phyloP100way.distributions.tsv.gz
<run-dir>/analytics/candidate_variants.phyloP100way.histograms.tsv.gz
<run-dir>/analytics/candidate_variants.phyloP100way.manifest.json
```

The report stores exact one-percentile distribution curves and compact
relative-frequency histograms per strategy and gnomAD-hit stratum, not millions
of allele-level score rows. Histogram widths use the Freedman-Diaconis rule on
the combined hit/non-hit scores for a strategy, capped at 80 display bins. On a
cold run it reads the union of candidate and ClinVar-required positions from
the remote bigWig once and reuses that positional map for both candidate
stratification and ClinVar conservation validation.

`Target-Space Null` compares GAPH SNVs with unobserved possible SNVs matched by
gene, target context, exact genomic REF>ALT substitution, and primary RefSeq VEP
consequence. GAPH callability and conservation are not matching variables.
VEP annotations use a release-specific SQLite cache, so interrupted REST or
local runs resume without repeating completed variant/gene pairs.

The same matched focal-control sets are also compared descriptively for exact
allele overlap with gnomAD and ClinVar, gnomAD AF among exact hits, and ClinVar
class composition among records with non-empty CLNSIG. These annotations are
stored once per unique matched allele in:

```text
<run-dir>/analytics/negative_control/target_space_null.external_evidence.tsv.gz
<run-dir>/analytics/negative_control/target_space_null.external_evidence.manifest.json
```

Failed gnomAD regions remain missing and are never interpreted as absence.
An incomplete external-evidence table is resumable: a later report run keeps
successful allele lookups and requests only alleles from failed regions.

To repair failed gnomAD regions in an existing pipeline annotation without
changing the original output, run:

```bash
python bin/complete_gnomad_annotation.py \
  --run-dir "$RUN" \
  --gnomad-cache-dir "$GAPH_GNOMAD_CACHE_DIR"
```

The command writes `<run-dir>/annotation_gnomad_complete/`. Repeating it retries
only regions still listed as failed. The cache argument is optional and also
defaults to `GAPH_GNOMAD_CACHE_DIR`. The recovery directory contains the
corrected variant annotations, failure table, and manifest; immutable strategy
support tables remain in `<run-dir>/annotation/` and are reused by the report.
Build a report from the repaired annotation with:

```bash
micromamba run -n gaph-v2-analytics python -m analytics.strategy_report \
  --run-dir "$RUN" \
  --annotation-dir "$RUN/annotation_gnomad_complete"
```

Its bounded, deterministic samples and annotations are cached under:

```text
<run-dir>/analytics/negative_control/
```

When enabled, the default limit is 25,000 focal SNVs per strategy with 1,000
target-space-null resamples. For a faster exploratory run, use
`--target-space-null-sample-size` and `--target-space-null-resamples`; changing
only the number of resamples reuses the prepared control tables.

The sample size is an engineering cap, not a fixed scientific cohort size. A
run uses every eligible focal SNV when fewer than the cap are available. The
report shows full focal-weighted phyloP ECDFs and descriptive target-space
intervals for all target-space-null outcomes. Each bootstrap replicate resamples
whole matched sets with replacement, chooses one available control from every
selected set, and uses the same draws for the GAPH, control, and paired-difference
statistics. These are descriptive 95% paired matched-set bootstrap intervals;
the report does not assign an inferential p-value to these comparators.

Raw p-values remain visible for the formal validation analyses. For each
analysis mode, SNV/INDEL selection, target-context selection, and ClinVar MC
consequence selection, Benjamini-Hochberg correction is applied across
strategies. Band-specific Fisher tests are corrected across strategies within
the same band. Consequence options with no estimable strategy result for the
current selectors are hidden from the interactive view and enumerated in QC.
Mantel-Haenszel confidence intervals use the Robins-Breslow-Greenland variance
implemented by `statsmodels.StratifiedTable`. Continuous models use Firth
logistic regression (`logistf`) with a three-degree-of-freedom natural spline
for phyloP100way and profile penalized-likelihood confidence intervals.
