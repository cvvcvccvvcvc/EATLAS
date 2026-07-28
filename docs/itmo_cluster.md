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

The two- and three-CPU outputs had identical summaries, feature coverage, and
canonical segment/event evidence. Minimap2 initially differed only in
line-order-derived provenance identifiers; those identifiers are now derived
from PAF content and remain stable across thread counts.

A 20-gene all-strategies run completed without task retries or scheduler
resource failures with `--alignment_max_forks 4`. Minimap2 used 2.32-2.53 CPUs
per task at the median, Nucmer used 2.04, and peak RSS stayed below 2.4 GB for
all sequence aligners. Four concurrent tasks per alignment process are
therefore the default on the verified cluster setup; this setting does not
change the separate Ensembl MAF or annotation concurrency limits.

In the same run, annotation partitions containing two genes used 0.38-5.6 GB
RSS for 66,921-815,531 unique variant contexts. Memory tracked context count
closely, while gene counts alone did not predict it. The default partition size
is therefore 10 genes, and partitioned annotation requests 24 GB on its first
attempt, retrying at 48 and 72 GB after memory-related failures. With the
default two annotation forks, the maximum retry reservation remains below the
verified 512 GB per-user Slurm memory limit.

Current memory requests are conservative initial bounds. Tune them from
Nextflow trace `peak_rss` after representative cluster runs. Requesting the
account maximum for every task wastes capacity and can increase queue time.

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

The `low_storage` profile disables process cache and cleans successful task work.
It minimizes retained data but does not provide useful `-resume` behavior after
a successful run.

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
`.nextflow/`, Conda environments, or the ignored `tools/bin/datasets` symlink.

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
mkdir -p "$GAPH_ROOT"/{bin,cache/gnomad,conda,envs,micromamba,nextflow,results,work}

curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xvj -C "$GAPH_ROOT" bin/micromamba

export PATH="$GAPH_ROOT/bin:$PATH"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export NXF_HOME="$GAPH_ROOT/nextflow"

cd /nfs/home/$USER/gaph_v2
micromamba create --yes \
  --prefix "$GAPH_ROOT/envs/controller" \
  --file envs/controller.yml

micromamba run -p "$GAPH_ROOT/envs/controller" nextflow -version
micromamba run -p "$GAPH_ROOT/envs/controller" java -version
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

## Staged Pipeline Validation

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
RUN="$GAPH_ROOT/results/slurm_smoke_1gene_asm20_$(date +%Y%m%d_%H%M%S)"

micromamba run -p "$GAPH_ROOT/envs/controller" nextflow run . \
  -profile slurm,low_storage \
  --stage all \
  --ids_file "$GAPH_ROOT/inputs/smoke_1_gene.ids" \
  --outdir "$RUN" \
  --alignment_strategies minimap2_asm20 \
  --fetch_max_forks 1 \
  --alignment_max_forks 1 \
  --annotation_max_forks 1
```

Verify the run before increasing scope:

```bash
test -s "$RUN/fetch/manifest.json"
test -s "$RUN/alignment/manifest.json"
test -s "$RUN/annotation/manifest.json"
test -s "$RUN/annotation/variant_annotations.tsv.gz"
du -sh "$RUN" "$GAPH_WORK_DIR"
```

Then validate in this order:

1. one gene with `--alignment_strategies all`
2. 10-20 representative large genes with all strategies
3. a medium scaling run whose trace is used to tune concurrency
4. only then a production-scale panel

Do not start with 20,000 genes. The local tests validate pipeline behavior, but
the Slurm scheduler, cluster network, shared filesystem, and external services
must be measured together before production scaling.

A two-gene align-only regression exposed and verified a multi-gene partition
merge fix. Both genes now reach `MERGE_ALIGNMENT_PARTITION` and
`MERGE_ALIGNMENT`; the AFDN/BRCA1 Minimap2 ASM20 plus Nucmer run produced
1,411,399 raw events, 2,750 ortholog-strategy summaries, no failures, and a
24 MB durable alignment dataset.

Annotation of the two-gene, seven-strategy AFDN/BRCA1 dataset collapsed
4,379,013 support rows to 470,359 variant contexts in 2m54s and used 3.3 GB
peak RSS. Annotation tasks therefore request 8 GB and partitioned annotation
retries OOM-like failures with a larger allocation. Representative scaling
runs should keep partitions small until memory scaling has been measured on a
larger panel.

## Monitoring And Analytics

Monitor orchestration from `sphinx` without running analysis there:

```bash
squeue -u "$USER"
tail -f "$GAPH_CODE/.nextflow.log"
du -sh "$GAPH_ROOT/work" "$RUN"
```

Run large analytics reports through a Slurm allocation, not directly on
`sphinx`. The analytics environment is defined separately in
`envs/analytics.yml`.
