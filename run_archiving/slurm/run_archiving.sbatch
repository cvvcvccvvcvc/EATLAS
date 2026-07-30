#!/usr/bin/env bash
#SBATCH --job-name=run-archive
#SBATCH --partition=main
#SBATCH --account=users
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/tank/scratch/%u/gaph_v2/logs/archive/%x-%j.log

set -euo pipefail

source "$HOME/.gaph_v2_cluster_env.sh"

archive_env="${GAPH_ARCHIVE_ENV:-$GAPH_ROOT/envs/run-archiving}"
if [[ ! -x "$archive_env/bin/python" ]]; then
    printf 'Archive environment is missing: %s\n' "$archive_env" >&2
    printf 'Create it from run_archiving/environment.yml first.\n' >&2
    exit 2
fi

cd "$GAPH_CODE"
exec "$archive_env/bin/python" -m run_archiving "$@"
