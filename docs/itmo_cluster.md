# ITMO Cluster Operations

Use this document only for first-time setup or infrastructure diagnosis on the
ITMO CT Slurm cluster. Ordinary pipeline and report work belongs in
`pipeline_launch.md` and `report_generation.md`.

Dated capacity tests and performance measurements are preserved in
`../experiments/itmo_cluster_baseline/README.md`; they explain past decisions
but are not current scheduler or pipeline configuration.

## Revalidate The Cluster

The last verified layout used:

- external login host `ctlab.itmo.ru`;
- Slurm controller `sphinx`;
- shared home under `/nfs/home/$USER`;
- shared work allocation under `/mnt/tank/scratch/$USER`;
- ordinary partition `main` and QOS `normal`.

These are infrastructure facts, not repository contracts. Recheck them during
first setup or after an administrator change:

```bash
squeue -u "$USER"
scontrol show partition main
sacctmgr show user "$USER" withassoc \
  format=User,Account,Partition,MaxSubmitJobs,MaxJobs
sacctmgr show qos normal \
  format=Name%20,MaxTRESPU%100,MaxJobsPU,MaxSubmitJobsPU
```

Do not infer usable concurrency from a historical run. Current task resources,
fork limits, retry policy, and Slurm queue size belong to `nextflow.config` and
`nextflow_schema.json`. Slurm still decides which submitted jobs run.

The repository's `slurm` profile disables task-local scratch so staged task
data stays on the shared work allocation. This matters because compute-node
`$TMPDIR` capacity and policy are not guaranteed. Keep source code and small
reference files in home; keep environments, caches, task work, analytics, and
results under the assigned work allocation.

## Login

On the local machine:

```bash
ITMO_USER=ilunegov
ssh -i ~/.ssh/itmo -o IdentitiesOnly=yes "$ITMO_USER@ctlab.itmo.ru"
ssh sphinx
```

Home and `/mnt/tank/scratch` must be visible from `sphinx` and compute nodes.
Never run pipeline or report computation directly on `sphinx); it hosts the
Nextflow controller and submission commands.

## Repository And Required References

Clone the authoritative repository into cluster home:

```bash
cd "/nfs/home/$USER"
git clone git@github.com:cvvcvccvvcvc/EATLAS.git gaph_v2
cd gaph_v2
git switch main
```

ClinVar and the RefSeq GFF3 are intentionally ignored by Git. Transfer them
from the local repository while preserving paths:

```bash
ITMO_USER=ilunegov

ssh -i ~/.ssh/itmo -o IdentitiesOnly=yes \
  "$ITMO_USER@ctlab.itmo.ru" \
  "mkdir -p /nfs/home/$ITMO_USER/gaph_v2"

rsync -aP --relative \
  -e 'ssh -i ~/.ssh/itmo -o IdentitiesOnly=yes' \
  ./assets/reference/clinvar/clinvar.vcf.gz \
  ./assets/reference/clinvar/clinvar.vcf.gz.tbi \
  ./assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz \
  "$ITMO_USER@ctlab.itmo.ru:/nfs/home/$ITMO_USER/gaph_v2/"
```

Verify content on both machines. The files used for the verified setup have:

```text
0ac3f9e8084b43ad09d367a20b54ada4cc9e592846f471a7c2d698e6dbf7b71a  clinvar.vcf.gz
4af776bf0c7ca2cb613fd62a64aae550f0d48fc7fdfb11d95c43e403092020c6  clinvar.vcf.gz.tbi
4920f0eae7e2197c50b67a201e06d657387137b49dd60f474b4f1d5b29334051  genomic.gff.gz
```

Generate local hashes with `shasum -a 256` and cluster hashes with
`sha256sum`. A different intentional reference version requires a pipeline
contract change; do not silently accept a mismatch.

Do not transfer tracked source files, `results/`, `work/`, `.nextflow/`, or
Conda environments outside the documented Git/reference flow.

## Controller And Analytics Environments

Keep Micromamba, environments, and package caches on shared work storage:

```bash
export GAPH_ROOT="/mnt/tank/scratch/$USER/gaph_v2"
mkdir -p "$GAPH_ROOT"/{bin,cache/gnomad,cache/vep_results,conda,containers,downloads/vep,envs,micromamba,nextflow,reference/ucsc,reference/vep,results,singularity/cache,singularity/tmp,work}

curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xvj -C "$GAPH_ROOT" bin/micromamba

export PATH="$GAPH_ROOT/bin:$PATH"
export MAMBA_ROOT_PREFIX="$GAPH_ROOT/micromamba"
export NXF_HOME="$GAPH_ROOT/nextflow"

cd "/nfs/home/$USER/gaph_v2"
micromamba create --yes \
  --prefix "$GAPH_ROOT/envs/controller" \
  --file envs/controller.yml
micromamba create --yes \
  --prefix "$GAPH_ROOT/envs/analytics" \
  --file envs/analytics.yml

micromamba run -p "$GAPH_ROOT/envs/controller" nextflow -version
micromamba run -p "$GAPH_ROOT/envs/controller" java -version
micromamba run -p "$GAPH_ROOT/envs/analytics" \
  python -m analytics.strategy_report --help
```

The active Java and Nextflow versions are owned by `envs/controller.yml`.
`micromamba info` must show the root prefix and package cache below
`$GAPH_ROOT`.

## Persistent Runtime Configuration

Store non-secret cluster paths in `$HOME/.gaph_v2_cluster_env.sh`:

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

Source this file after every login. NCBI credentials belong in ignored
`$GAPH_CODE/.env`, which is read only by the fetch task. Set mode `600` and
never print or commit it.

## Shared gnomAD And VEP Result Caches

Create and retain the shared infrastructure caches:

```bash
mkdir -p "$GAPH_GNOMAD_CACHE_DIR" "$GAPH_VEP_RESULT_CACHE_DIR"
test -w "$GAPH_GNOMAD_CACHE_DIR"
test -w "$GAPH_VEP_RESULT_CACHE_DIR"
```

The gnomAD cache stores complete regional responses and is safe to reuse across
pipeline and analytics runs. The VEP result cache stores only complete
release/config-matched variant-gene results. Neither cache is part of a source
run or a substitute for durable annotation shards.

## Local phyloP

Analytics never downloads phyloP implicitly. Keep one shared official UCSC
hg38 phyloP100way BigWig:

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
  mv "$PHYLOP_PARTIAL" "$GAPH_PHYLOP_BIGWIG"
fi
```

Analytics identifies this reference by content, not path. Moving unchanged
content preserves cache identity; replacing it invalidates affected artifacts.

## Local VEP

The current cluster contract is release 116 RefSeq GRCh38 with:

```text
$GAPH_ROOT/containers/ensembl-vep_release_116.0.sif
$GAPH_ROOT/reference/vep/homo_sapiens_refseq/116_GRCh38/
$GAPH_ROOT/bin/gaph-vep116
```

The image and indexed cache are large externally provisioned references; verify
them before a first run. Install the repository wrapper:

```bash
install -m 755 "$GAPH_CODE/bin/gaph-vep116" "$GAPH_ROOT/bin/gaph-vep116"
test -s "$GAPH_ROOT/containers/ensembl-vep_release_116.0.sif"
test -d "$GAPH_ROOT/reference/vep/homo_sapiens_refseq/116_GRCh38"
```

The wrapper discovers `singularity` and then `apptainer`; set
`GAPH_CONTAINER_RUNTIME` only when explicit override is required. It binds
`$GAPH_ROOT` into the container. Validate it on a compute node, not on the
controller:

```bash
srun -p main --cpus-per-task=1 --mem=2G --time=00:10:00 \
  "$GAPH_VEP_EXECUTABLE" \
  --offline --cache --refseq --assembly GRCh38 \
  --dir_cache "$GAPH_VEP_CACHE_DIR" --show_cache_info
```

Pipeline VEP uses the indexed cache with `--offline`, `--cache`, `--refseq`,
`--use_given_ref`, and `--assembly GRCh38`. It does not require a FASTA
because the current annotation contract does not request HGVS or reference
checking.

## Compute-Node Preflight

Before the first pipeline run or after a mount/network change:

```bash
srun -p main --cpus-per-task=1 --mem=2G --time=00:10:00 --pty bash

hostname
test -r "/nfs/home/$USER/gaph_v2/main.nf" && echo code-readable
test -w "/mnt/tank/scratch/$USER" && echo scratch-writable
df -h "/mnt/tank/scratch/$USER"
curl -I --max-time 20 https://api.ncbi.nlm.nih.gov
exit
```

An HTTP response proves routing, DNS, and TLS only. The one-gene pipeline smoke
is the authoritative integration check across scheduler, mounts, environments,
and remote services.

## One-Gene Pipeline Smoke

Run from a persistent `tmux` session on `sphinx`:

```bash
source "$HOME/.gaph_v2_cluster_env.sh"
cd "$GAPH_CODE"
mkdir -p "$GAPH_ROOT/inputs"
printf '57007\n' > "$GAPH_ROOT/inputs/smoke_1_gene.ids"

SMOKE_ID="slurm_smoke_1gene_$(date +%Y%m%d_%H%M%S)"
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

Verify the durable boundary:

```bash
test -s "$RUN/run_manifest.json"
test -s "$RUN/evidence_inventory.json"
test -s "$RUN/fetch/manifest.json"
test -s "$RUN/alignment/manifest.json"
test -s "$RUN/annotation/manifest.json"
test -s "$RUN/annotation/variant_annotations/manifest.json"
```

Then validate the broader sequence in `run_validation.md` before production
scale. Ordinary production runs must switch back to
`scripts/slurm/run_pipelines.sh`.

## Diagnosis

Use bounded inspection from `sphinx`:

```bash
squeue -u "$USER"
tail -n 50 "$RUN/reports/nextflow/nextflow.log"
du -sh "$RUN" "$WORK" 2>/dev/null || du -sh "$RUN"
```

Inspect the Nextflow trace before changing resources. `maxForks` applies per
process and is not a global worker count; requested CPUs must also be used by
the task command. Scheduler memory limits can leave valid large tasks pending.

For ordinary recovery, rerun the documented launcher command. Do not create an
ad hoc launcher, run analysis on `sphinx`, or delete a failed run's dedicated
work directory before recovery is no longer needed.
