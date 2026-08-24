#!/usr/bin/env bash

set -euo pipefail

if (( $# < 5 )); then
  printf 'Expected ANALYTICS_ROOT REPORT_NAME GIT_COMMIT PROJECT_ROOT RUN_COUNT RUN_DIR... [report arguments...]\n' >&2
  exit 2
fi

analytics_root=$1
report_name=$2
expected_commit=$3
project_root=$4
run_count=$5
shift 5
[[ "$run_count" =~ ^[1-9][0-9]*$ ]] || {
  printf 'Invalid source-run count: %s\n' "$run_count" >&2
  exit 2
}
(( $# >= run_count )) || {
  printf 'Expected %s source-run arguments\n' "$run_count" >&2
  exit 2
}
run_dirs=("${@:1:run_count}")
shift "$run_count"
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
git -C "$project_root" diff --quiet || {
  printf 'Repository has tracked changes after report submission\n' >&2
  exit 1
}
git -C "$project_root" diff --cached --quiet || {
  printf 'Repository has staged changes after report submission\n' >&2
  exit 1
}
[[ -x "$analytics_python" ]] || {
  printf 'Analytics Python is missing: %s\n' "$analytics_python" >&2
  exit 1
}

mkdir -p "$analytics_root/slurm"
printf '%s\n' "$actual_commit" > "$analytics_root/slurm/$report_name.git_commit"

command=(
  "$analytics_python" -m analytics.strategy_report
  --analytics-root "$analytics_root"
  --report-name "$report_name"
)
for run_dir in "${run_dirs[@]}"; do
  command+=(--run-dir "$run_dir")
done
command+=("${report_args[@]}")

printf 'Git commit: %s\n' "$actual_commit"
printf 'Command:'
printf ' %q' "${command[@]}"
printf '\n'

cd "$project_root"
exec "${command[@]}"
