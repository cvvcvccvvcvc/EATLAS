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
long-pseudoread `minimap2_map_ont_pseudoreads_30000_15000` must be selected
explicitly.

## Connect And Update

Completed code changes reach the cluster through GitHub only:

```bash
git fetch origin main
git status --short --branch
git merge-base --is-ancestor origin/main HEAD
git push origin main
INTENDED_COMMIT=$(git rev-parse HEAD)
test "$(git rev-parse origin/main)" = "$INTENDED_COMMIT"

ssh -i ~/.ssh/itmo ilunegov@ctlab.itmo.ru
ssh sphinx

cd /nfs/home/$USER/gaph_v2
git fetch origin main
git merge --ff-only origin/main
source "$HOME/.gaph_v2_cluster_env.sh"
test -z "$(git status --porcelain)"
git rev-parse HEAD
```

The final cluster hash must equal the exact `INTENDED_COMMIT` captured in the
authoritative local checkout. Do not infer currency from a clean cluster tree,
from a commit-looking run name, or by comparing cluster `HEAD` with the
cluster's `origin/main` before fetching: remote-tracking refs can be stale.

This revision gate applies again before report submission, even if the source
pipeline run is already complete. After a pipeline creates `run_manifest.json`,
verify that its `git_commit` equals `INTENDED_COMMIT`. After a report worker
starts, verify that the `Git commit:` line in its Slurm stdout equals the same
commit.

Treat any mismatch between the documented command and the cluster launcher as
a stale-checkout failure. Stop and synchronize the checkout. In particular, do
not remove current arguments, substitute a historical cohort workflow, or
otherwise reshape the requested run merely to make an older launcher accept
the command.

Do not copy tracked source files with `rsync`, Git bundles, or ad hoc archives.
A source run intended for analytics must start from a clean Git working tree,
including no untracked files. The run manifest records this state, and analytics
rejects dirty source runs instead of treating unreproducible evidence as valid.

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
  --analytics-root "$GAPH_ROOT/analytics" \
  --report-name strategy_compare \
  -- \
  --target-space-null
```

Everything after `--` is a report argument and is forwarded unchanged. Omit
the target-space-null option unless the user requested it. See
`docs/report_generation.md` for report and multi-run details.

Before starting Nextflow, the combined launcher verifies the clean checkout,
the external analytics path, and its own report arguments. It exposes the
standard pipeline configuration plus optional `--alignment-strategies`. When a
measured run requires concurrency overrides, use the direct Nextflow command
above and submit its report separately after completion.

## Resume

Resume with the same `RUN`, `WORK`, inputs, and relevant parameters. Use the
same command with `-resume`; do not invent a new work directory for a resume.
The pipeline refuses to reuse a successfully completed `RUN`; use a new run
directory when source evidence must be regenerated.

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
