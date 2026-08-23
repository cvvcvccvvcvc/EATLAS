# Pipeline Launch

Use this runbook for an ordinary GAPH Nextflow launch or resume on the ITMO CT
cluster. Use `docs/run_validation.md` only for smoke tests, validation, or a
concrete pipeline failure.

## Required User Inputs

Obtain or identify:

- the input Gene ID file;
- the durable result directory;
- a dedicated Nextflow work directory retained while resume may be needed;
- any explicitly requested alignment strategies or concurrency overrides.

The cluster environment must also declare the release-pinned local VEP
executable, indexed cache, and optional shared result cache described in
`docs/itmo_cluster.md`. Candidate VEP is part of annotation, not a later report
precompute.

When an option is not specified, keep the Nextflow/config default. The pipeline
always runs end to end, and the default strategy selection runs
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

## Run And Report Together

When a report should be submitted immediately after a successful pipeline run,
use the combined launcher in the persistent session:

```bash
bash scripts/slurm/run_and_report.sh \
  --ids-file "$IDS" \
  --run-dir "$RUN" \
  --work-dir "$WORK" \
  --report-name strategy_compare \
  -- \
  --target-space-null
```

Everything after `--` is a report argument and is forwarded unchanged. Omit
the target-space-null option unless the user requested it. See
`docs/report_generation.md` for report and cohort details.

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
