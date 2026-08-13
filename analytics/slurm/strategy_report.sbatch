#!/usr/bin/env bash

set -euo pipefail

if (( $# < 6 )); then
  printf 'Expected SOURCE_KIND SOURCE_PATH COHORT_ROOT REPORT_NAME GIT_COMMIT PROJECT_ROOT [report arguments...]\n' >&2
  exit 2
fi

source_kind=$1
source_path=$2
cohort_root=$3
report_name=$4
expected_commit=$5
project_root=$6
shift 6
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

if [[ "$source_kind" == run ]]; then
  analysis_root=$source_path
  source_args=(--run-dir "$source_path")
  mkdir -p "$analysis_root/reports"
  commit_marker="$analysis_root/reports/$report_name.git_commit"
elif [[ "$source_kind" == cohort ]]; then
  analysis_root=$cohort_root
  source_args=(--cohort-manifest "$source_path" --cohort-root "$cohort_root")
  mkdir -p "$analysis_root/slurm"
  commit_marker="$analysis_root/slurm/$report_name.git_commit"
else
  printf 'Unknown report source kind: %s\n' "$source_kind" >&2
  exit 2
fi
printf '%s\n' "$actual_commit" > "$commit_marker"

command=(
  "$analytics_python" -m analytics.strategy_report
  "${source_args[@]}"
  --report-name "$report_name"
  "${report_args[@]}"
)

printf 'Git commit: %s\n' "$actual_commit"
printf 'Command:'
printf ' %q' "${command[@]}"
printf '\n'

cd "$project_root"
exec "${command[@]}"
