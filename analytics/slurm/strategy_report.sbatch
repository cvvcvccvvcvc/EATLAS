#!/usr/bin/env bash

set -euo pipefail

if (( $# < 4 )); then
  printf 'Expected RUN_DIR REPORT_NAME GIT_COMMIT PROJECT_ROOT [report arguments...]\n' >&2
  exit 2
fi

run_dir=$1
report_name=$2
expected_commit=$3
project_root=$4
shift 4
report_args=("$@")

source "$HOME/.gaph_v2_cluster_env.sh"

: "${GAPH_ROOT:?GAPH_ROOT is not set by .gaph_v2_cluster_env.sh}"
analytics_env="$GAPH_ROOT/envs/analytics"
analytics_python="$analytics_env/bin/python"
export PATH="$analytics_env/bin:$PATH"

actual_commit=$(git -C "$project_root" rev-parse HEAD)
if [[ "$actual_commit" != "$expected_commit" ]]; then
  printf 'Repository commit changed after submission: expected %s, found %s\n' \
    "$expected_commit" "$actual_commit" >&2
  exit 1
fi
[[ -x "$analytics_python" ]] || {
  printf 'Analytics Python is missing: %s\n' "$analytics_python" >&2
  exit 1
}

mkdir -p "$run_dir/reports"
printf '%s\n' "$actual_commit" > "$run_dir/reports/$report_name.git_commit"

command=(
  "$analytics_python" -m analytics.strategy_report
  --run-dir "$run_dir"
  --report-name "$report_name"
  "${report_args[@]}"
)

printf 'Git commit: %s\n' "$actual_commit"
printf 'Command:'
printf ' %q' "${command[@]}"
printf '\n'

cd "$project_root"
exec "${command[@]}"
