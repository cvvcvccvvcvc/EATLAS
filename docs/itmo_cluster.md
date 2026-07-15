# ITMO Cluster Operations

Use this document for GAPH v2 setup and execution on the ITMO CT Slurm
cluster. It records repository-specific decisions and facts verified on the
cluster. Do not store passwords, API keys, or private SSH keys here.

## Verified Cluster Layout

Verified on 2026-07-15:

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

The user association did not report explicit `MaxSubmitJobs` or `MaxJobs`
limits. The ordinary `main` partition did not require an explicit QOS, account,
or partition setting in the Nextflow profile. Do not request a high-resource
QOS unless the cluster administrator explicitly authorizes it.

The connectivity probe `curl -I https://api.ncbi.nlm.nih.gov` returned HTTP
403. That confirms routing, DNS, and TLS to the host, but not successful access
to every NCBI, Ensembl, UCSC, or gnomAD endpoint. The one-gene pipeline smoke
test remains the authoritative network check.

Useful scheduler inspection commands:

```bash
squeue -u "$USER"
scontrol show partition main
sacctmgr show user "$USER" withassoc format=User,Account,Partition,MaxSubmitJobs,MaxJobs
sacctmgr show qos format=Name%20,MaxJobsPU,MaxWall
```

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

Clone the repository into home after its Git remote is available:

```bash
cd /nfs/home/$USER
git clone <repository-url> gaph_v2
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

## Controller Environment

Conda is user-managed on this cluster. From the repository on `sphinx`, keep
both the controller environment and all package caches in `/mnt/tank/scratch`:

```bash
export GAPH_CODE="/nfs/home/$USER/gaph_v2"
export GAPH_ROOT="/mnt/tank/scratch/$USER/gaph_v2"

mkdir -p "$GAPH_ROOT"/{conda,envs,micromamba,nextflow,results,work}

export CONDA_PKGS_DIRS="$GAPH_ROOT/conda/controller-pkgs"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export NXF_HOME="$GAPH_ROOT/nextflow"

cd "$GAPH_CODE"
conda env create --prefix "$GAPH_ROOT/envs/controller" -f envs/controller.yml
conda activate "$GAPH_ROOT/envs/controller"

nextflow -version
java -version
micromamba --version
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
export NXF_CONDA_CACHEDIR="$GAPH_ROOT/conda/envs"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export NXF_HOME="$GAPH_ROOT/nextflow"
```

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

nextflow run . \
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
