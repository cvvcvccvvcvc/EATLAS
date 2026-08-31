#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  analytics/slurm/submit_strategy_report.sh \
    --analytics-root /absolute/path/to/analytics \
    --run-dir /absolute/path/to/run [--run-dir /absolute/path/to/another_run ...] \
    --report-name report_name \
    --expected-commit FULL_GIT_COMMIT \
    [--slurm-cpus 8] [--slurm-memory 128G] [--slurm-time 06:00:00] \
    [--slurm-partition main] [-- <analytics.strategy_report arguments>]

Arguments after -- are passed unchanged to analytics.strategy_report.
The analytics root, source runs, and report name are managed by this launcher.
EOF
}

fail() {
  printf 'Error: %s\n' "$*" >&2
  exit 2
}

canonical_destination() {
  local path=$1
  local suffix=""
  while [[ ! -e "$path" ]]; do
    suffix="/$(basename "$path")$suffix"
    path=$(dirname "$path")
  done
  [[ -d "$path" ]] || return 1
  path=$(cd "$path" && pwd -P) || return 1
  printf '%s%s\n' "$path" "$suffix"
}

analytics_root=""
run_dirs=()
report_name=""
expected_commit=""
slurm_cpus=8
slurm_memory=128G
slurm_time=06:00:00
slurm_partition=main
report_args=()

while (( $# > 0 )); do
  case "$1" in
    --analytics-root)
      (( $# >= 2 )) || fail "--analytics-root requires a value"
      analytics_root=$2
      shift 2
      ;;
    --run-dir)
      (( $# >= 2 )) || fail "--run-dir requires a value"
      run_dirs+=("$2")
      shift 2
      ;;
    --report-name)
      (( $# >= 2 )) || fail "--report-name requires a value"
      report_name=$2
      shift 2
      ;;
    --expected-commit)
      (( $# >= 2 )) || fail "--expected-commit requires a value"
      expected_commit=$2
      shift 2
      ;;
    --slurm-cpus)
      (( $# >= 2 )) || fail "--slurm-cpus requires a value"
      slurm_cpus=$2
      shift 2
      ;;
    --slurm-memory)
      (( $# >= 2 )) || fail "--slurm-memory requires a value"
      slurm_memory=$2
      shift 2
      ;;
    --slurm-time)
      (( $# >= 2 )) || fail "--slurm-time requires a value"
      slurm_time=$2
      shift 2
      ;;
    --slurm-partition)
      (( $# >= 2 )) || fail "--slurm-partition requires a value"
      slurm_partition=$2
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    --)
      shift
      report_args=("$@")
      break
      ;;
    *)
      fail "unknown launcher argument '$1'; put strategy_report arguments after --"
      ;;
  esac
done

[[ -n "$analytics_root" ]] || fail "--analytics-root is required"
[[ "$analytics_root" = /* ]] || fail "--analytics-root must be an absolute path"
(( ${#run_dirs[@]} > 0 )) || fail "at least one --run-dir is required"
requested_analytics_root=$analytics_root
analytics_root=$(canonical_destination "$analytics_root") || fail \
  "--analytics-root cannot be resolved: $requested_analytics_root"
resolved_run_dirs=()
for run_dir in "${run_dirs[@]}"; do
  [[ "$run_dir" = /* ]] || fail "--run-dir must be an absolute path: $run_dir"
  [[ -d "$run_dir" ]] || fail "run directory does not exist: $run_dir"
  run_dir=$(cd "$run_dir" && pwd -P)
  case "$analytics_root/" in
    "$run_dir/"*) fail "--analytics-root must be outside source run: $run_dir" ;;
  esac
  [[ -s "$run_dir/run_manifest.json" ]] || fail "missing run manifest: $run_dir"
  [[ -s "$run_dir/annotation/manifest.json" ]] || fail \
    "missing annotation manifest: $run_dir"
  [[ -s "$run_dir/annotation/variant_annotations/manifest.json" ]] || fail \
    "missing finalized variant annotations: $run_dir"
  resolved_run_dirs+=("$run_dir")
done
run_dirs=("${resolved_run_dirs[@]}")
[[ -n "$report_name" ]] || fail "--report-name is required"
[[ "$expected_commit" =~ ^[0-9a-f]{40}$ ]] || fail \
  "--expected-commit must be a full 40-character Git commit"
[[ "$report_name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || fail \
  "--report-name may contain only letters, digits, dot, underscore, and hyphen"
[[ "$slurm_cpus" =~ ^[1-9][0-9]*$ ]] || fail "--slurm-cpus must be a positive integer"
[[ -n "$slurm_memory" ]] || fail "--slurm-memory must not be empty"
[[ -n "$slurm_time" ]] || fail "--slurm-time must not be empty"
[[ -n "$slurm_partition" ]] || fail "--slurm-partition must not be empty"
command -v sbatch >/dev/null || fail "sbatch was not found; run this on the Slurm controller"

for argument in "${report_args[@]}"; do
  case "$argument" in
    --analytics-root|--analytics-root=*|--run-dir|--run-dir=*|--report-name|--report-name=*|--expected-commit|--expected-commit=*)
      fail "$argument is managed by the launcher and must appear before --"
      ;;
  esac
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/../.." && pwd)
batch_script="$script_dir/strategy_report.sbatch"
[[ -f "$batch_script" ]] || fail "missing Slurm batch script: $batch_script"
git_status=$(git -C "$project_root" status --porcelain=v1 --untracked-files=normal) || fail \
  "cannot inspect repository status: $project_root"
[[ -z "$git_status" ]] || fail "report submission requires a clean working tree"
git -C "$project_root" fetch origin main >/dev/null || fail "cannot fetch authoritative origin/main"
git_commit=$(git -C "$project_root" rev-parse HEAD) || fail "cannot resolve repository HEAD"
origin_commit=$(git -C "$project_root" rev-parse origin/main) || fail "cannot resolve origin/main"
[[ "$git_commit" = "$expected_commit" ]] || fail \
  "cluster HEAD $git_commit does not match expected commit $expected_commit"
[[ "$origin_commit" = "$expected_commit" ]] || fail \
  "fetched origin/main $origin_commit does not match expected commit $expected_commit"

log_dir="$analytics_root/slurm"
mkdir -p "$log_dir"
job_tag=${report_name:0:40}
job_id=$(sbatch --parsable \
  --job-name="gaph-report-$job_tag" \
  --partition="$slurm_partition" \
  --cpus-per-task="$slurm_cpus" \
  --mem="$slurm_memory" \
  --time="$slurm_time" \
  --output="$log_dir/$report_name.%j.out" \
  --error="$log_dir/$report_name.%j.err" \
  "$batch_script" \
  "$analytics_root" \
  "$report_name" \
  "$git_commit" \
  "$project_root" \
  "${#run_dirs[@]}" \
  "${run_dirs[@]}" \
  "${report_args[@]}")

printf 'Submitted report job %s\n' "$job_id"
printf 'Analytics workspace: %s\n' "$analytics_root"
printf 'Logs: %s/%s.%s.{out,err}\n' "$log_dir" "$report_name" "$job_id"
