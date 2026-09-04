# GAPH v2

GAPH v2 is an evidence-first comparative variant pipeline for Entrez Gene IDs.
One Nextflow DSL2 run:

1. fetches the fixed GRCh38.p14 human loci and the complete NCBI ortholog set;
2. aligns selected ortholog gene sequences with one or more registered methods;
3. annotates normalized alignment events with ClinVar, gnomAD, and Ensembl VEP.

The pipeline publishes durable row-level evidence. A separate analytics command
derives scientific tables and an HTML report without modifying completed runs.

## Quick Start

Create the controller environment once:

```bash
micromamba create --yes -f envs/controller.yml
```

Place the required reference files at the repository defaults or pass their
paths explicitly:

```text
assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz
assets/reference/clinvar/clinvar.vcf.gz
assets/reference/clinvar/clinvar.vcf.gz.tbi
```

Run a small local panel:

```bash
RUN="results/run_$(date +%Y%m%d_%H%M%S)"

micromamba run -n gaph-v2-controller nextflow run . \
  --ids_file assets/inputs/gene_ids/panel_10_genes.txt \
  --outdir "$RUN"
```

Every process uses its declared environment from `envs/`. Local execution needs
no Nextflow profile. The target assembly is fixed to GRCh38.p14
(`GCF_000001405.40`), and ortholog retrieval always uses NCBI
`--ortholog all`.

The GFF3 defaults to the path above. ClinVar is required and defaults to the
indexed VCF above. Small local runs use Ensembl REST for VEP by default;
production cluster runs use a release-pinned local VEP configuration.

Use the schema-generated help for current pipeline parameters and defaults:

```bash
micromamba run -n gaph-v2-controller nextflow run . --helpFull
```

Only NCBI credentials are loaded from an ignored project `.env` by the fetch
task. Copy `.env.example` and set `ENTREZ_EMAIL` and `ENTREZ_API_KEY` there if
needed. Machine paths and pipeline configuration must be exported in the
controller environment or passed as Nextflow parameters; they are not loaded
from `.env` into the controller.

## Cluster Runs

Ordinary ITMO runs use the validated launcher, not a direct `nextflow run`:

```bash
bash scripts/slurm/run_pipelines.sh \
  --results-root "$GAPH_ROOT/results/run_group" \
  --expected-commit "$INTENDED_COMMIT" \
  /absolute/path/to/gene_ids.txt
```

Read [Pipeline launch](docs/pipeline_launch.md) before launching or resuming.
Use [ITMO cluster operations](docs/itmo_cluster.md) only for first-time setup or
infrastructure diagnosis.

## Durable Output

```text
<run>/
  run_manifest.json
  evidence_inventory.json
  fetch/
  alignment/
  annotation/
  reports/nextflow/
```

- `fetch/` keeps normalized target, selected-ortholog, feature, taxonomy, and
  failure evidence plus target FASTA.
- `alignment/` keeps partitioned summaries, segments, compact events, exact
  ortholog support, and failures.
- `annotation/` keeps the partitioned ClinVar/gnomAD/VEP dataset, the exact
  event-to-variant lineage, and lookup failures.
- `evidence_inventory.json` binds every durable evidence file by size and
  SHA-256; `run_manifest.json` records completion and launch provenance.

A completed run is immutable. Analytics writes to a separate
`--analytics-root`; Nextflow `work/` is only execution/resume state. See
[Storage model](docs/storage_model.md) for retention and cleanup rules.

## Analytics

Create the analytics environment once, then build a report from one or more
compatible completed runs:

```bash
micromamba create --yes -f envs/analytics.yml

micromamba run -n gaph-v2-analytics python -m analytics.strategy_report \
  --analytics-root /absolute/path/to/analytics \
  --run-dir /absolute/path/to/completed-run \
  --report-name strategy_compare
```

Read [Report generation](docs/report_generation.md) for the supported launch
procedure and [Analytics contract](docs/analytics_contract.md) for input,
cache, and scientific-result semantics.

## Documentation

| Need | Read |
| --- | --- |
| Launch or resume pipelines on ITMO | [Pipeline launch](docs/pipeline_launch.md) |
| Generate a report | [Report generation](docs/report_generation.md) |
| Understand analytics compatibility and scientific derivations | [Analytics contract](docs/analytics_contract.md) |
| Run tests, smoke checks, or diagnose a failure | [Run and validation](docs/run_validation.md) |
| Set up or diagnose ITMO infrastructure | [ITMO cluster operations](docs/itmo_cluster.md) |
| Find owning code and focused tests | [Project map](docs/project_map.md) |
| Inspect pipeline data contracts | [Fetch](docs/stage1_fetch_contract.md), [alignment](docs/stage2_alignment_contract.md), and [annotation](docs/stage3_annotation_contract.md) contracts |
| Decide what is durable, cached, or removable | [Storage model](docs/storage_model.md) |
| Archive, verify, restore, or remove a run | [Run archiving](run_archiving/README.md) |
