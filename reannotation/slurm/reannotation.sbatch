#!/usr/bin/env bash
set -euo pipefail
project_root=$1
expected_commit=$2
shift 2
source "$HOME/.gaph_v2_cluster_env.sh"
: "${GAPH_ROOT:?GAPH_ROOT is required}"
[[ $(git -C "$project_root" rev-parse HEAD) == "$expected_commit" && -z $(git -C "$project_root" status --porcelain) ]] || {
  printf 'Checkout changed since submission\n' >&2; exit 2;
}
export PATH="$GAPH_ROOT/envs/reannotation/bin:$PATH"
export GAPH_GNOMAD_CACHE_DIR="$GAPH_ROOT/cache/gnomad"
export RCLONE_TPSLIMIT=2
export RCLONE_TPSLIMIT_BURST=2
export PYTHONUNBUFFERED=1
cd "$project_root"
printf 'Code commit: %s\n' "$expected_commit"
exec python -m reannotation run "$@"
