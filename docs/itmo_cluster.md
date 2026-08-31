# ITMO Cluster Operations

Use this document for GAPH v2 setup and execution on the ITMO CT Slurm
cluster. It records repository-specific decisions and facts verified on the
cluster. Do not store passwords, API keys, or private SSH keys here.

## Verified Cluster Layout

Verified through 2026-07-16:

- external login host: `ctlab.itmo.ru`
- Slurm controller: `sphinx`; do not run pipeline computations directly there
- Slurm account: `users`
- default partition: `main`
- home: `/nfs/home/$USER`; administrator-provided quota is 20 GB
- durable work allocation: `/mnt/tank/scratch/$USER`; quota is 250 GB
- `/mnt/tank` is shared between the controller and compute nodes
- a test allocation ran on `meduza-1`
- compute-node outbound DNS/TLS connectivity was present
- system Java on `sphinx` was OpenJDK 11
- `normal` QOS permits up to 64 CPU, 512 GB RAM, and 2 GPU per user
- `high_mem` is not authorized for this account; ordinary jobs use `normal`
- 1 CPU / 3 GB and 2 CPU / 8 GB test allocations started successfully
- Micromamba 2.8.1 is installed under the assigned work allocation
- the controller environment provides Nextflow 25.10.4 and OpenJDK 17
- the controller environment occupies about 608 MB and its package cache about
  445 MB

The user association did not report explicit `MaxSubmitJobs` or `MaxJobs`
limits. A controlled test submitted five independent 1 CPU / 1 GB jobs and all
five ran concurrently on `meduza-1`; the historical four-job observation is not
a current scheduler limit. The ordinary `main` partition does not require an
explicit account or partition setting in the Nextflow profile.

Before a run, inspect `squeue -u "$USER"` for obsolete jobs in
`launch failed requeued held`; cancel only jobs confirmed to be obsolete. Four
such 8 GB jobs from an older pipeline were found and removed during setup, but
they were not the root cause of the first smoke-test failure.

The `main` partition reports `DefMemPerNode=UNLIMITED`. A task submitted without
an explicit memory request can therefore inherit hundreds of gigabytes and be
held with `QOSMaxMemoryPerUser`. This happened when partial `withName` blocks in
composed profiles replaced the base process resource configuration. Scratch
policy now uses the `task_scratch` label instead, so the `slurm` profile can
disable task-local scratch without removing CPU, memory, time, retry, or Conda
directives. When changing profiles, verify the effective configuration rather
than assuming repeated selector blocks merge field by field.

The connectivity probe `curl -I https://api.ncbi.nlm.nih.gov` returned HTTP
403. That confirms routing, DNS, and TLS to the host, but not successful access
to every NCBI, Ensembl, UCSC, or gnomAD endpoint. The one-gene pipeline smoke
test remains the authoritative network check.

Useful scheduler inspection commands:

```bash
squeue -u "$USER"
scontrol show partition main
sacctmgr show user "$USER" withassoc format=User,Account,Partition,MaxSubmitJobs,MaxJobs
sacctmgr show qos normal format=Name%20,MaxTRESPU%100,MaxJobsPU,MaxSubmitJobsPU
```

## Resource Model

Nextflow submits one Slurm job for each process task:

- `cpus` reserves CPUs for that task; the command must explicitly use threads
  to benefit from them
- `memory` is a task limit/reservation, not a performance accelerator
- `maxForks` limits concurrent instances of one process name, not the whole
  workflow; the same value on several alignment processes can yield more jobs
  in total
- `executor.queueSize` limits submitted jobs, including pending jobs; Slurm
  still decides how many actually run

Minimap2 and Nucmer receive `task.cpus` explicitly. A two-gene benchmark using
AFDN (GeneID 4301) and BRCA1 (GeneID 672) compared two, three, and four CPUs per
task. Three CPUs reduced individual alignment time by about 31-36% relative to
two CPUs and was only 4-6% slower than four CPUs. Across the four measured
tasks, two and three CPUs consumed approximately the same reserved CPU-seconds,
whereas four CPUs consumed about 27% more than three. Minimap2 and Nucmer
therefore request three CPUs by default. BWA and the other strategies were not
changed by this benchmark.

The two- and three-CPU outputs had identical summaries, analytics-derived
feature coverage, and canonical segment/event evidence. Minimap2 initially differed only in
line-order-derived provenance identifiers; those identifiers are now derived
from PAF content and remain stable across thread counts.

A 20-gene default-strategy smoke completed without task retries or scheduler
resource failures with `--alignment_max_forks 4`. Minimap2 used 2.32-2.53 CPUs
per task at the median and Nucmer used 2.04. This validates the task shape but
does not establish four as a throughput optimum: the verified user QOS permits
64 CPUs, while each alignment task reserves three.

Three later production runs exposed a size-dependent memory tail in the
long-pseudoread minimap2 strategy. Across 1,470 logical tasks, 1,405 genes below
150 million selected-ortholog bases stayed below 6.5 GB RSS without retry. Of
the 63 genes from 150 to 600 million bases, the maximum successful RSS was
21.5 GB; the two genes at or above 600 million bases reached 31.8 GB. Starting
every task at 8 GB caused 38 failed attempts across 29 genes before retries
reached a sufficient allocation. The workflow now records the total selected-
ortholog sequence volume during task preparation and starts these classes at
8, 24, and 40 GB respectively. BWA starts at 16 GB only for the largest class;
all other BWA tasks remain at 8 GB. Fixed retry increments remain as protection
against outliers.

In the same run, annotation partitions containing two genes used 0.38-5.6 GB
RSS for 66,921-815,531 unique variant contexts. The current resource function
uses each partition's compact `alignment_event_count`: 8 or 16 GB at up to 1 or
5 million events, then 32, 48, 64, or 96 GB at up to 15, 30, 40, or more than
40 million events. Each retry adds 32 GB. With eight annotation forks, eight
largest first attempts would request 768 GB in aggregate. Slurm enforces the
verified 512 GB per-user memory limit, so some tasks and larger retries may
queue rather than run simultaneously.

Current memory requests are conservative initial bounds. Tune them from
Nextflow trace `peak_rss` after representative cluster runs. Requesting the
account maximum for every task wastes capacity and can increase queue time.
The trace also records `attempt`, requested `cpus`, `memory`, and `time`, so a
retry can be distinguished from a large successful first allocation.

The 490-gene `all_genes_batch_001_20260830_163720_f4_a21_n8` run on commit
`8a6f666` completed all 2,347 tasks in 4 h 51 min 49 s, with no failed tasks or
fetch, alignment, or annotation failures. Its tested concurrency settings are
now the defaults: 4 fetch tasks, 21 tasks per alignment process, and 8 tasks per
annotation process. Slurm still enforces the 64-CPU and 512-GB per-user limits,
and `maxForks` applies separately to minimap2, nucmer, and BWA. This successful
run does not establish a throughput optimum; further tuning should compare
wall time, retries, peak RSS, and external-service failures on the same workload
and cache state.

## Shared gnomAD Cache

Large annotation runs must reuse the shared regional cache rather than fetch
the same gnomAD regions again after retries or in later runs. As verified on
2026-08-09, the existing ITMO cache is located at:

```text
/mnt/tank/scratch/$USER/gaph_v2/cache/gnomad
```

For `ilunegov` it contained 3,306 tile files and occupied about 710 MB at the
time of verification. The cache is incremental and safe to keep between runs.
Configure it in every controller shell and pass it explicitly in production
commands:

```bash
export GAPH_GNOMAD_CACHE_DIR="$GAPH_ROOT/cache/gnomad"
mkdir -p "$GAPH_GNOMAD_CACHE_DIR"
test -w "$GAPH_GNOMAD_CACHE_DIR"

nextflow run . \
  -profile slurm \
  --ids_file /path/to/gene_ids.txt \
  --gnomad_cache_dir "$GAPH_GNOMAD_CACHE_DIR" \
  --outdir "$RUN"
```

When `GAPH_GNOMAD_CACHE_DIR` is absent but `GAPH_ROOT` is set, the pipeline
defaults to `$GAPH_ROOT/cache/gnomad`. Explicitly passing the parameter remains
recommended because the resolved cache path is then visible in the launch
command and run manifest.

## Storage Policy

Keep only source code and small reference data in home:

```text
/nfs/home/$USER/gaph_v2
```

Keep environments, package caches, Nextflow runtime files, task work, and run
outputs in the assigned work allocation:

```text
/mnt/tank/scratch/$USER/gaph_v2/
  conda/
  envs/
  micromamba/
  nextflow/
  results/
  work/
```

On `meduza-1`, `$TMPDIR` resolved to `/tmp`, backed by the node-local root disk.
At the time of inspection it had 67 GB free. Node-local storage can be faster
than `/mnt/tank` because it avoids network filesystem I/O, but its capacity and
per-user policy are not guaranteed for this project. The repository therefore
sets `scratch = false` in the `slurm` profile. Large task data remains under
`GAPH_WORK_DIR`; local profiles retain task scratch.

The default storage policy keeps process cache available while a run is active
or failed, then cleans task work created by a successful execution session. This
supports recovery with `-resume` without retaining fresh successful-run work
indefinitely. After a resumed run succeeds, task directories from its earlier
failed session can remain and may be removed when its dedicated work path is no
longer needed.

## Login

From the local machine:

```bash
ssh -i ~/.ssh/itmo -o IdentitiesOnly=yes ilunegov@ctlab.itmo.ru
ssh ilunegov@sphinx
```

Home and `/mnt/tank/scratch` are shared, so files placed there through the login
host are visible from `sphinx` and compute nodes.

## Repository And References

Clone the repository into home:

```bash
cd /nfs/home/$USER
git clone https://github.com/cvvcvccvvcvc/EATLAS.git gaph_v2
cd gaph_v2
git switch main
```

ClinVar and the RefSeq GFF3 are intentionally ignored by Git. Copy them from the
local repository while preserving their paths. Run this on the local Mac from
the repository root:

```bash
ssh -i ~/.ssh/itmo -o IdentitiesOnly=yes \
  ilunegov@ctlab.itmo.ru \
  'mkdir -p /nfs/home/ilunegov/gaph_v2'

rsync -aP --relative \
  -e 'ssh -i ~/.ssh/itmo -o IdentitiesOnly=yes' \
  ./assets/reference/clinvar/clinvar.vcf.gz \
  ./assets/reference/clinvar/clinvar.vcf.gz.tbi \
  ./assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz \
  ilunegov@ctlab.itmo.ru:/nfs/home/ilunegov/gaph_v2/
```

The three files occupy about 249 MB. Verify their integrity before the first
run. On the local Mac:

```bash
shasum -a 256 \
  assets/reference/clinvar/clinvar.vcf.gz \
  assets/reference/clinvar/clinvar.vcf.gz.tbi \
  assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz
```

On `sphinx`, run the corresponding command from the repository root:

```bash
sha256sum \
  assets/reference/clinvar/clinvar.vcf.gz \
  assets/reference/clinvar/clinvar.vcf.gz.tbi \
  assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz
```

The hashes must match line by line. Do not transfer local `results/`, `work/`,
`.nextflow/` or Conda environments.

Expected SHA-256 values for the files transferred on 2026-07-16:

```text
0ac3f9e8084b43ad09d367a20b54ada4cc9e592846f471a7c2d698e6dbf7b71a  clinvar.vcf.gz
4af776bf0c7ca2cb613fd62a64aae550f0d48fc7fdfb11d95c43e403092020c6  clinvar.vcf.gz.tbi
4920f0eae7e2197c50b67a201e06d657387137b49dd60f474b4f1d5b29334051  genomic.gff.gz
```

## Controller Environment

Bootstrap Micromamba on `sphinx` and keep the controller environment and all
package caches in `/mnt/tank/scratch`:

```bash
export GAPH_ROOT="/mnt/tank/scratch/$USER/gaph_v2"
mkdir -p "$GAPH_ROOT"/{bin,cache/gnomad,cache/vep_results,conda,envs,micromamba,nextflow,reference/ucsc,results,work}

curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xvj -C "$GAPH_ROOT" bin/micromamba

export PATH="$GAPH_ROOT/bin:$PATH"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export NXF_HOME="$GAPH_ROOT/nextflow"

cd /nfs/home/$USER/gaph_v2
micromamba create --yes \
  --prefix "$GAPH_ROOT/envs/controller" \
  --file envs/controller.yml

micromamba create --yes \
  --prefix "$GAPH_ROOT/envs/analytics" \
  --file envs/analytics.yml

micromamba run -p "$GAPH_ROOT/envs/controller" nextflow -version
micromamba run -p "$GAPH_ROOT/envs/controller" java -version
micromamba run -p "$GAPH_ROOT/envs/analytics" python -m analytics.strategy_report --help
micromamba info
```

The active Java must be version 17 from `envs/controller.yml`, not the system
Java 11. In `micromamba info`, both the root prefix and package cache must be
under `$GAPH_ROOT`.

Set the runtime variables in every shell that launches Nextflow:

```bash
export GAPH_CODE="/nfs/home/$USER/gaph_v2"
export GAPH_ROOT="/mnt/tank/scratch/$USER/gaph_v2"
export GAPH_WORK_DIR="$GAPH_ROOT/work"
export GAPH_GNOMAD_CACHE_DIR="$GAPH_ROOT/cache/gnomad"
export GAPH_PHYLOP_BIGWIG="$GAPH_ROOT/reference/ucsc/hg38.phyloP100way.bw"
export GAPH_VEP_BACKEND="local"
export GAPH_VEP_RELEASE="116"
export GAPH_VEP_EXECUTABLE="$GAPH_ROOT/bin/gaph-vep116"
export GAPH_VEP_CACHE_DIR="$GAPH_ROOT/reference/vep"
export GAPH_VEP_RESULT_CACHE_DIR="$GAPH_ROOT/cache/vep_results"
export GAPH_VEP_RESULT_CACHE_TILE_SIZE_BP="1000000"
export GAPH_VEP_FORKS="4"
export NXF_CONDA_CACHEDIR="$GAPH_ROOT/conda/envs"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export CONDA_PKGS_DIRS="$MAMBA_ROOT_PREFIX/pkgs"
export NXF_HOME="$GAPH_ROOT/nextflow"
export PATH="$GAPH_ROOT/bin:$PATH"
```

These exports are stored in `$HOME/.gaph_v2_cluster_env.sh` on the cluster. Run
`source "$HOME/.gaph_v2_cluster_env.sh"` after each login. NCBI credentials are
configured in the ignored `$GAPH_CODE/.env`; the file was verified with mode
`600`. Never print its values into logs or commit the file.

If NCBI credentials are used, place them in the ignored `$GAPH_CODE/.env`, set
permissions to `600`, and never commit the file.

## Local phyloP For Analytics

Keep one shared copy of the official UCSC hg38 phyloP100way BigWig under
`$GAPH_ROOT/reference/ucsc`. The file is about 9.2 GB. Downloading is explicit;
report generation never downloads reference files implicitly:

```bash
mkdir -p "$GAPH_ROOT/reference/ucsc"
PHYLOP_PARTIAL="$GAPH_ROOT/reference/ucsc/hg38.phyloP100way.bw.partial"
rsync --partial --append-verify --info=progress2 \
  rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw \
  "$PHYLOP_PARTIAL"

if printf '%s  %s\n' \
  43858006bdf98145b6fd239490bd0478 \
  "$PHYLOP_PARTIAL" \
  | md5sum --check -; then
  mv "$PHYLOP_PARTIAL" "$GAPH_ROOT/reference/ucsc/hg38.phyloP100way.bw"
fi
```

Set `GAPH_PHYLOP_BIGWIG` in `$HOME/.gaph_v2_cluster_env.sh`. The report also
discovers this exact path automatically when `GAPH_ROOT` is set. Candidate,
ClinVar, and target-null caches record the local file path, size, and
modification time; replacing the BigWig invalidates the affected caches.

Verified on 2026-08-09 against the saved remote target-null annotations from
the 590-gene run: all 172,183 allele rows matched exactly, including 458
positions without a score. For the same 172,043 unique positions and 638 read
blocks, BigWig read time fell from 707.361 seconds over UCSC HTTPS to 6.554
seconds from the shared local file (about 108x); the complete local annotation
pass took 10.184 seconds with about 269 MiB peak RSS. This comparison covers
the phyloP annotation pass, not the rest of report generation.

## Local VEP For Pipeline And Analytics

Verified on 2026-07-29:

- Ensembl VEP image: `ensemblorg/ensembl-vep:release_116.0`
- image path: `$GAPH_ROOT/containers/ensembl-vep_release_116.0.sif` (`281 MB`)
- indexed RefSeq cache: `$GAPH_ROOT/reference/vep/homo_sapiens_refseq/116_GRCh38`
  (`25 GB` as reported by `du`)
- retained download archive: `$GAPH_ROOT/downloads/vep/homo_sapiens_refseq_vep_116_GRCh38.tar.gz`
  (`26,409,563,680` bytes)
- wrapper: `$GAPH_ROOT/bin/gaph-vep116`, installed from
  `$GAPH_CODE/bin/gaph-vep116`

Install or refresh the wrapper without embedding a node-specific runtime path:

```bash
install -m 755 "$GAPH_CODE/bin/gaph-vep116" "$GAPH_ROOT/bin/gaph-vep116"
```

The wrapper discovers `singularity` and then `apptainer` from `PATH`, so the
same `GAPH_VEP_EXECUTABLE` works across compute nodes with different install
locations. `GAPH_CONTAINER_RUNTIME` remains an optional explicit override; it
is not required for the verified cluster setup. Automatic discovery was
verified on 2026-08-09 on `kraken`, where Singularity is installed as
`/usr/bin/singularity`.

The cache archive is the official Ensembl release-116 RefSeq-only GRCh38
archive. Installation validated the complete archive, extracted it outside the
final cache path, and ran VEP `--show_cache_info` before marking it complete.
The installation used 50 GB while retaining both the archive and extracted
cache. The archive is intentionally still present; remove it only as a separate
explicit storage decision.

Basic consequence annotation uses `--offline`, `--cache`, `--refseq`,
`--use_given_ref`, and `--assembly GRCh38`. No FASTA is installed because the
current shared VEP contract does not request HGVS or reference checking. The
`--use_given_ref` flag is required for this BAM-edited RefSeq cache; without it,
VEP automatically requests transcript reference sequences and fails without a
FASTA.

A six-variant APC/BRCA1 smoke covering SNV, insertion, and deletion matched
Ensembl REST release 116 exactly after normalizing VEP's `-` placeholder to an
empty value. Consequence terms, selected RefSeq transcript, canonical/MANE
fields, impact, and variant class agreed for all six variants. Local VEP wrote
the six results in about 3.3 seconds.

The end-to-end Nextflow pipeline runs candidate VEP in bounded Slurm tasks; do
not run it directly on `sphinx`. A 118,549,428-row candidate dataset would
produce 475 shards at the fixed 250,000-row bound. Benchmark a representative
shard before increasing annotation concurrency: each task uses four VEP workers,
and concurrent tasks share the same network-backed 25 GB reference cache.
Completed shard tasks are the `-resume` boundary; the finalizer publishes their
gzip outputs without constructing a global candidate table or full-candidate
SQLite cache. Completed variant/gene results are additionally reused across runs
from `$GAPH_VEP_RESULT_CACHE_DIR`. This sparse Parquet cache is separate from the
official indexed VEP reference cache in `$GAPH_VEP_CACHE_DIR`; incomplete
results are not stored in it.

A 500,000-row benchmark sampled evenly from 40 of the 590-gene artifact's 475
partitions. One-megabase and five-megabase tiles had the same 1.88 MB size and
about 8.0-second lookup time; one-megabase tiles remain the default because
they provide finer pruning for sparse regional requests without a measured
cost.

## Compute-Node Preflight

Before the first pipeline run:

```bash
srun -p main --cpus-per-task=1 --mem=2G --time=00:10:00 --pty bash

hostname
test -r /nfs/home/$USER/gaph_v2/main.nf && echo code-readable
test -w /mnt/tank/scratch/$USER && echo scratch-writable
df -h /mnt/tank/scratch/$USER
curl -I --max-time 20 https://api.ncbi.nlm.nih.gov
exit
```

This is a connectivity and mount check only. It does not replace the pipeline
smoke test because the real workflow uses several different remote services.

## End-to-End Pipeline Validation

Run Nextflow from `sphinx`. Use `tmux` or another cluster-approved persistent
terminal mechanism so an SSH disconnect does not terminate the controller.

Create a one-gene input in the durable work allocation:

```bash
mkdir -p "$GAPH_ROOT/inputs"
printf '57007\n' > "$GAPH_ROOT/inputs/smoke_1_gene.ids"
```

First run one strategy with all concurrency set to one:

```bash
cd "$GAPH_CODE"
SMOKE_ID="slurm_smoke_1gene_asm20_$(date +%Y%m%d_%H%M%S)"
RUN="$GAPH_ROOT/results/$SMOKE_ID"
WORK="$GAPH_ROOT/work/$SMOKE_ID"

micromamba run -p "$GAPH_ROOT/envs/controller" nextflow run . \
  -profile slurm \
  --ids_file "$GAPH_ROOT/inputs/smoke_1_gene.ids" \
  --outdir "$RUN" \
  --gnomad_cache_dir "$GAPH_GNOMAD_CACHE_DIR" \
  --alignment_strategies minimap2_asm20 \
  --fetch_max_forks 1 \
  --alignment_max_forks 1 \
  --annotation_max_forks 1 \
  -work-dir "$WORK"
```

Verify the run before increasing scope:

```bash
test -s "$RUN/fetch/manifest.json"
test -s "$RUN/alignment/manifest.json"
test -s "$RUN/annotation/manifest.json"
test -s "$RUN/annotation/variant_annotations/manifest.json"
du -sh "$RUN" "$WORK" 2>/dev/null || du -sh "$RUN"
```

Then validate in this order:

1. one gene with `--alignment_strategies default`
2. 10-20 representative large genes with all default-enabled strategies
3. a medium scaling run whose trace is used to tune concurrency
4. only then a production-scale panel

Do not start with 20,000 genes. The local tests validate pipeline behavior, but
the Slurm scheduler, cluster network, shared filesystem, and external services
must be measured together before production scaling.

## Monitoring And Analytics

For an ordinary pipeline launch, use `docs/pipeline_launch.md`. For an ordinary
analytics report submission, use `docs/report_generation.md`. The commands
below are infrastructure-level reference material.

Monitor orchestration from `sphinx` without running analysis there:

```bash
squeue -u "$USER"
tail -f "$GAPH_CODE/.nextflow.log"
du -sh "$GAPH_ROOT/work" "$RUN"
```

Run large analytics reports through a Slurm allocation, not directly on
`sphinx`. The launchers use `$GAPH_ROOT/envs/analytics`, created from
`envs/analytics.yml` during setup.
