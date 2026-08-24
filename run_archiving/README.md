# Run Archiving

This package archives complete pipeline run directories to an rclone remote.
It is operational tooling: it does not modify or import pipeline or analytics
code.

An archive preserves every regular file below the run directory. Nextflow
`work/` remains outside the run and is never archived. Symlinks and embedded
execution-cache directories are rejected instead of being silently skipped.
New runs must have a root `run_manifest.json` with `status=complete`,
`success=true`, and `exit_status=0`.

A completed source run is immutable. Analytics belongs in its separate external
workspace and is not included in a run archive. The archiver rejects a top-level
`analytics/` directory and rejects report content other than
`reports/nextflow/`, so derived reports cannot be silently archived as source
evidence.

Remote layout:

```text
<remote-root>/
  runs/
    <run-id>/
      data/                 exact run contents
      _archive/
        manifest.json       file sizes and content hashes
        MD5SUMS
        SHA256SUMS
        COMPLETE.json       written only after remote verification
```

## One-time cluster setup

Create the isolated environment and log directory:

```bash
source "$HOME/.gaph_v2_cluster_env.sh"
mkdir -p "$GAPH_ROOT/envs" "$GAPH_ROOT/logs/archive"

micromamba create --yes \
  --prefix "$GAPH_ROOT/envs/run-archiving" \
  --file "$GAPH_CODE/run_archiving/environment.yml"
```

Configure Google Drive interactively:

```bash
"$GAPH_ROOT/envs/run-archiving/bin/rclone" config
```

Use a dedicated remote name such as `gdrive`. On a headless cluster, answer
`n` to automatic browser configuration and run the displayed
`rclone authorize drive` command on a computer with a browser. For the narrowest
access, choose the scope that permits access only to files created by rclone.

Protect the resulting token:

```bash
chmod 700 "$HOME/.config" "$HOME/.config/rclone"
chmod 600 "$HOME/.config/rclone/rclone.conf"
```

Create the archive root and verify account quota:

```bash
"$GAPH_ROOT/envs/run-archiving/bin/rclone" mkdir gdrive:GAPH
"$GAPH_ROOT/envs/run-archiving/bin/rclone" about gdrive:
```

Add this non-secret setting to `$HOME/.gaph_v2_cluster_env.sh`:

```bash
export GAPH_ARCHIVE_REMOTE="gdrive:GAPH"
```

Do not put the OAuth token in Git, command-line arguments, Slurm exports, or
job logs. `rclone.conf` contains the refresh token and must remain mode `600`.

## Archive

Run a local dry run first:

```bash
"$GAPH_ROOT/envs/run-archiving/bin/python" -m run_archiving archive \
  --run-dir "$GAPH_ROOT/results/<run-id>" \
  --dry-run
```

Submit the real transfer as a small Slurm job:

```bash
sbatch run_archiving/slurm/run_archiving.sbatch archive \
  --run-dir "$GAPH_ROOT/results/<run-id>"
```

The command is resumable at file boundaries. A repeated command verifies and
returns `already_archived` when the same run is already complete. It refuses to
replace different data under an existing run ID.

Historical runs created before `run_manifest.json` existed require an explicit
exception:

```bash
sbatch run_archiving/slurm/run_archiving.sbatch archive \
  --run-dir "$GAPH_ROOT/results/<legacy-run-id>" \
  --allow-legacy-run
```

The exception applies only when the root manifest is absent. A manifest that
says `running` or `failed` is always rejected.

## List archives

```bash
"$GAPH_ROOT/envs/run-archiving/bin/python" -m run_archiving list
```

The command reads only the small completion markers and reports archive IDs,
completion times, file counts, and sizes. It does not re-hash archived data.
Incomplete uploads without `COMPLETE.json` are intentionally omitted.

## Verify and restore

```bash
"$GAPH_ROOT/envs/run-archiving/bin/python" -m run_archiving verify \
  --run-id <run-id>

"$GAPH_ROOT/envs/run-archiving/bin/python" -m run_archiving restore \
  --run-id <run-id> \
  --destination "$GAPH_ROOT/results_restored/<run-id>"
```

`verify` checks every remote file against its MD5 and checks total file count
and bytes. `restore` downloads into `<destination>.partial`, can resume that
partial copy, verifies it with SHA-256, and only then renames it to the final
destination.

## Remove the cluster copy

Removal is deliberately separate from upload:

```bash
"$GAPH_ROOT/envs/run-archiving/bin/python" -m run_archiving remove-local \
  --run-dir "$GAPH_ROOT/results/<run-id>" \
  --confirm-run-id <run-id>
```

The command performs a fresh remote checksum check, requires the run to be a
direct child of `$GAPH_ROOT/results`, and requires an exact run-ID confirmation.
It then atomically renames the run to a quarantine path in the same directory,
re-hashes it there, and removes it only when it still matches the archive.

The archiver verifies workflow completion and transfer integrity, not scientific
validity. The operator remains responsible for reviewing stage-level failures
and deciding whether a completed run should be retained.
