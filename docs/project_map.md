# Project Map

Use this map when the owning boundary of a change is not obvious. Current
values and executable cases remain in code, schema, and tests; this document
routes work to them.

## Production Boundaries

```text
Entrez IDs
  -> main.nf / FETCH_STAGE
  -> main.nf / ALIGNMENT_STAGE
  -> main.nf / PARTITIONED_ANNOTATION_STAGE
  -> immutable run evidence
  -> analytics.strategy_report
```

Nextflow owns process wiring and runtime policy. Python commands own
deterministic transformations. `genomics/` owns cross-boundary biological
semantics. `analytics/` reads completed evidence but never participates in the
pipeline. `run_archiving/` is isolated operational tooling.

## Task To Owner

| Change | Start in | Contract or procedure | Focused tests |
| --- | --- | --- | --- |
| User-facing pipeline parameter | `nextflow_schema.json`, then `nextflow.config` or `main.nf` | relevant stage contract | `test_pipeline_launch_validation.py`, launcher tests |
| Workflow channel or process boundary | `main.nf`, `modules/local/` | relevant stage contract | `test_canonical_partition_wiring.py` plus boundary-specific tests |
| Fetch parsing, selection, features, or taxonomy | `bin/normalize_ids.py`, `bin/fetch_parse_chunk.py`, `bin/build_fetch_dataset.py`, `bin/fetch_taxonomy.py` | `stage1_fetch_contract.md` | normalize, fetch fallback, fetch dataset, taxonomy ownership tests |
| Alignment registry or runner | `ALIGNMENT_STRATEGY_REGISTRY` in `main.nf`, alignment modules, `bin/run_*alignment.py` | `stage2_alignment_contract.md` | runner and runtime tests |
| Alignment table or partition layout | `bin/alignment_table_schema.py`, `bin/merge_alignment_results.py` | `stage2_alignment_contract.md` | schema and merge tests |
| Variant normalization or ClinVar/gnomAD evidence | `bin/annotate_events.py`, `genomics/variants.py`, `genomics/clinvar.py`, `genomics/gnomad*.py` | `stage3_annotation_contract.md` | annotation and provider-specific tests |
| VEP execution or result reuse | `bin/annotate_vep_partition.py`, `genomics/vep/`, `bin/gaph-vep116` | `stage3_annotation_contract.md`, cluster setup when relevant | VEP partition, consequence, runtime, and cache tests |
| Run manifest or evidence inventory | `lib/RunManifest.groovy`, `provenance/evidence_inventory.py` | `storage_model.md` | run-manifest and evidence-inventory tests |
| Analytics source compatibility or cache identity | `analytics/io/run_inputs.py`, `analytics/io/calculation_identity.py`, `analytics/io/artifacts.py` | `analytics_contract.md` | analysis-input, calculation-identity, and artifact tests |
| Reusable evidence derivation | `analytics/derivations/`, `analytics/io/alignment_aggregates.py`, `analytics/io/annotation_support.py` | `analytics_contract.md` | matching cache and derivation tests |
| Scientific calculation | `analytics/analyses/` | `analytics_contract.md` | matching analysis test |
| HTML/report composition | `analytics/reporting/`, `analytics/strategy_report.py` | `analytics_contract.md` | `test_strategy_report.py` and section-specific tests |
| Cluster pipeline/report launcher | `scripts/slurm/run_pipelines.sh`, `analytics/slurm/submit_strategy_report.sh` | launch runbook | matching launcher test |
| Archive/restore/removal | `run_archiving/` | `run_archiving/README.md` | `test_run_archiving.py` |

Test paths above are under `tests/`. Use `rg` for the named symbol or behavior
instead of relying on line numbers.

## Entrypoints And Sources Of Truth

- End-to-end workflow: `main.nf`
- Pipeline parameter names, types, and defaults: `nextflow_schema.json`
- Process resources, environments, retry policy, publication, and profiles:
  `nextflow.config`
- Alignment strategy membership and defaults: `ALIGNMENT_STRATEGY_REGISTRY` in
  `main.nf`
- Ordinary cluster pipeline launch: `scripts/slurm/run_pipelines.sh --help`
- Analytics CLI: `python -m analytics.strategy_report --help`
- Ordinary cluster report submission:
  `analytics/slurm/submit_strategy_report.sh --help`
- Run archive CLI: `python -m run_archiving --help`

Do not copy these mutable values into additional documents. Stage documents own
durable biological and table semantics; launch documents own operator order and
safety gates.

## Package Ownership

- `modules/local/` stages files and invokes deterministic commands; it should
  not contain domain parsing logic.
- `bin/` contains pipeline-only command modules and exact pipeline table
  schemas.
- `genomics/` is the shared domain package for variant, ClinVar, gnomAD,
  taxonomy, and VEP semantics.
- `provenance/` contains reusable durable-evidence inventory logic.
- `analytics/io/` validates immutable source runs and owns artifact identities.
- `analytics/derivations/` builds reusable evidence-derived relations.
- `analytics/analyses/` owns scientific calculations.
- `analytics/reporting/` renders already-computed results.
- `envs/` owns task and controller dependencies.

Shared packages are staged explicitly into Nextflow tasks and invoked with
`python -m`. Do not solve package placement with `PYTHONPATH` or import shims.

## Data Placement

Durable pipeline output lives below one `params.outdir`; completed source runs
are immutable. Analytics lives under an external `--analytics-root`. Nextflow
task work, package caches, provider downloads, and native aligner files are not
source evidence. The complete retention contract is in `storage_model.md`.

Reusable small inputs live under `assets/inputs/`. Operational reference files
live under `assets/reference/` and may be intentionally ignored by Git.
Experiments and dated benchmarks belong under `experiments/<name>/`, including
their generated outputs and notes.

## Documentation Ownership

| Subject | Owner |
| --- | --- |
| Human overview and entry points | `README.md` |
| Agent invariants and routing | `AGENTS.md` |
| Ordinary pipeline launch/resume | `pipeline_launch.md` |
| Ordinary report launch | `report_generation.md` |
| Analytics compatibility, cache, and scientific semantics | `analytics_contract.md` |
| Smoke tests and failure investigation | `run_validation.md` |
| First-time ITMO setup and infrastructure diagnosis | `itmo_cluster.md` |
| Durable fetch/alignment/annotation semantics | the three stage contracts |
| Durable versus temporary data | `storage_model.md` |
| Archive operations and safety | `run_archiving/README.md` |

Historical measurements explain past decisions but are not current
configuration. Keep them under `experiments/` and link them only when a tuning
or infrastructure task needs that evidence.
