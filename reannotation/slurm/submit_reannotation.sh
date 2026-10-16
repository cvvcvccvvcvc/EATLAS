#!/usr/bin/env bash
set -euo pipefail

if [[ ${1:-} == --help ]]; then
  printf 'Usage: submit_reannotation.sh --expected-commit SHA -- <reannotation run arguments>\n'
  exit 0
fi
[[ ${1:-} == --expected-commit && ${2:-} =~ ^[0-9a-f]{40}$ && ${3:-} == -- ]] || {
  printf 'Expected --expected-commit FULL_SHA -- <reannotation run arguments>\n' >&2; exit 2;
}
expected_commit=$2
shift 3
project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
source "$HOME/.gaph_v2_cluster_env.sh"
: "${GAPH_ROOT:?GAPH_ROOT is required}"
[[ -z $(git -C "$project_root" status --porcelain) ]] || { printf 'Checkout must be clean\n' >&2; exit 2; }
git -C "$project_root" fetch origin main >/dev/null
[[ $(git -C "$project_root" rev-parse HEAD) == "$expected_commit" &&
   $(git -C "$project_root" rev-parse origin/main) == "$expected_commit" ]] || {
  printf 'Expected commit must equal clean HEAD and freshly fetched origin/main\n' >&2; exit 2;
}
[[ -x "$GAPH_ROOT/envs/reannotation/bin/python" ]] || { printf 'Create envs/reannotation.yml environment first\n' >&2; exit 2; }
mkdir -p "$GAPH_ROOT/reannotation/slurm"
sbatch --parsable --job-name=gaph-gnomad-reannotation --partition=main --cpus-per-task=2 --mem=8G --time=7-00:00:00 \
  --output="$GAPH_ROOT/reannotation/slurm/reannotation.%j.out" \
  --error="$GAPH_ROOT/reannotation/slurm/reannotation.%j.err" \
  "$project_root/reannotation/slurm/reannotation.sbatch" "$project_root" "$expected_commit" "$@"
