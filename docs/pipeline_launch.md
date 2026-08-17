# Pipeline Launch

Use this runbook for an ordinary GAPH Nextflow launch or resume on the ITMO CT
cluster. Use `docs/run_validation.md` only for smoke tests, validation, or a
concrete pipeline failure.

## Required User Inputs

Obtain or identify:

- the input Gene ID file;
- the durable result directory;
- the durable Nextflow work directory;
- the requested stage, if it is not the default `all`;
- any explicitly requested alignment strategies or concurrency overrides.

When an option is not specified, keep the Nextflow/config default. In
particular, the default stage is `all` and the default strategy selection runs
`minimap2_asm10`, `minimap2_asm20`, `nucmer`, and `bwa_pseudoreads_150_75`. The
long-pseudoread `minimap2_map_ont_pseudoreads_30000_15000` and precomputed
Ensembl strategies must be selected explicitly.

## Connect And Update

Completed code changes reach the cluster through GitHub only:

```bash
git push origin main

ssh -i ~/.ssh/itmo ilunegov@ctlab.itmo.ru
ssh sphinx

cd /nfs/home/$USER/gaph_v2
git pull --ff-only
source "$HOME/.gaph_v2_cluster_env.sh"
```

Do not copy tracked source files with `rsync`, Git bundles, or ad hoc archives.

## Start A Persistent Controller Session

Nextflow is a long-running controller process and must run in `tmux` or another
cluster-approved persistent terminal session on `sphinx`:

```bash
tmux new -s gaph_run_name
```

Inside the session:

```bash
cd "$GAPH_CODE"
source "$HOME/.gaph_v2_cluster_env.sh"

IDS="$GAPH_ROOT/inputs/panel.ids"
RUN="$GAPH_ROOT/results/run_name"
WORK="$GAPH_ROOT/work/run_name"

micromamba run -p "$GAPH_ROOT/envs/controller" nextflow run . \
  -profile slurm \
  --stage all \
  --ids_file "$IDS" \
  --outdir "$RUN" \
  --gnomad_cache_dir "$GAPH_GNOMAD_CACHE_DIR" \
  -work-dir "$WORK" \
  -resume
```

Add only options requested by the user or required by a concrete run, for
example:

```bash
--alignment_strategies minimap2_asm20,nucmer \
--fetch_max_forks 2 \
--alignment_max_forks 4 \
--annotation_max_forks 4
```

For the fixed opt-in long-pseudoread strategy, use:

```bash
--alignment_strategies minimap2_map_ont_pseudoreads_30000_15000
```

Detach without stopping Nextflow with `Ctrl-b d`. Reattach with:

```bash
tmux attach -t gaph_run_name
```

## Resume

Resume with the same `RUN`, `WORK`, inputs, and relevant parameters. Use the
same command with `-resume`; do not invent a new work directory for a resume.

```bash
tmux attach -t gaph_run_name
```

If the original controller process ended, start a new persistent session and
rerun the original command with `-resume`.

## Minimal Monitoring

```bash
squeue -u "$USER"
tail -n 50 "$GAPH_CODE/.nextflow.log"
```

Do not run pipeline computation directly on `sphinx`; Nextflow submits compute
tasks to Slurm. Do not perform broad disk, environment, or scheduler audits
before every ordinary launch. Use the deeper checks in `docs/run_validation.md`
or `docs/itmo_cluster.md` only for a first-time setup or a concrete failure.
