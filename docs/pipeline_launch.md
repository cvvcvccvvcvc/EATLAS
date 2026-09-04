# Pipeline Launch

Use this runbook to launch or resume one or more GAPH pipeline runs on the ITMO
CT cluster. Ordinary cluster runs use `scripts/slurm/run_pipelines.sh`; use
direct Nextflow commands only for the smoke tests and failure investigation in
`docs/run_validation.md`.

## Required User Inputs

Obtain or identify:

- one or more input Gene ID files, in execution order;
- one durable results root for their run directories;
- the intended full Git commit from the authoritative local checkout;
- any explicitly requested alignment strategies or concurrency overrides.

The cluster environment must also declare the release-pinned local VEP
executable, indexed cache, and optional shared result cache described in
`docs/itmo_cluster.md`. Candidate VEP is part of annotation, not a later report
precompute.

When an option is not specified, keep the schema/config default. The pipeline
always runs end to end. `ALIGNMENT_STRATEGY_REGISTRY` in `main.nf` owns the
current default strategy set; `stage2_alignment_contract.md` explains each
strategy's scientific policy.

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
INTENDED_COMMIT=<paste-the-full-commit-captured-locally>
test -z "$(git status --porcelain)"
git rev-parse HEAD
```

The final cluster hash must equal the exact `INTENDED_COMMIT` captured in the
authoritative local checkout. Do not infer currency from a clean cluster tree,
from a commit-looking run name, or by comparing cluster `HEAD` with the
cluster's `origin/main` before fetching: remote-tracking refs can be stale.
The pipeline launcher repeats the fetch and equality checks before doing work.

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

Inside the session, pass one input for a single run or several inputs for a
sequential series:

```bash
cd "$GAPH_CODE"
source "$HOME/.gaph_v2_cluster_env.sh"

RESULTS_ROOT="$GAPH_ROOT/results/all_genes"

bash scripts/slurm/run_pipelines.sh \
  --results-root "$RESULTS_ROOT" \
  --expected-commit "$INTENDED_COMMIT" \
  "$GAPH_ROOT/inputs/all_genes/batch_001.txt" \
  "$GAPH_ROOT/inputs/all_genes/batch_002.txt"
```

The input basename becomes the run name, such as `batch_001`. The launcher
creates `$RESULTS_ROOT/batch_001` and uses one internal work directory below
`$GAPH_WORK_DIR` for the group. It runs only one pipeline at a time and stops at
the first failure.

Add only options requested by the user or required by a concrete run:

```bash
--alignment-strategies minimap2_asm20,nucmer \
--fetch-max-forks 2 \
--alignment-max-forks 4 \
--annotation-max-forks 4
```

Ordinary launches use the configured concurrency defaults. Override
the three fork limits only
when a specific run needs different limits.

Detach without stopping Nextflow with `Ctrl-b d`. Reattach with:

```bash
tmux attach -t gaph_run_name
```

## Resume

Rerun the same launcher command. It checks and skips successfully completed
runs, resumes the first incomplete run with its recorded Nextflow session, and
then continues in the original order. It refuses changed input paths, explicit
launcher settings, result paths, or Git provenance. Use a new results root when
completed evidence must be regenerated.

```bash
tmux attach -t gaph_run_name
```

If the original controller ended, start a new persistent session and rerun the
same command. Submit reports separately after the required pipeline runs have
completed; see `docs/report_generation.md`.

## Minimal Monitoring

```bash
squeue -u "$USER"
tail -n 50 "$RESULTS_ROOT/batch_001/reports/nextflow/nextflow.log"
```

Do not run pipeline computation directly on `sphinx`; Nextflow submits compute
tasks to Slurm. Do not perform broad disk, environment, or scheduler audits
before every ordinary launch. Use the deeper checks in `docs/run_validation.md`
or `docs/itmo_cluster.md` only for a first-time setup or a concrete failure.

Use the launcher help for the current accepted operational overrides:

```bash
bash scripts/slurm/run_pipelines.sh --help
```
