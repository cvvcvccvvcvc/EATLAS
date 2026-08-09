#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  analytics/slurm/submit_strategy_report.sh \
    --run-dir /absolute/path/to/run \
    --report-name report_name \
    [--slurm-cpus 8] [--slurm-memory 128G] [--slurm-time 06:00:00] \
    [--slurm-partition main] [-- <analytics.strategy_report arguments>]

Arguments after -- are passed unchanged to analytics.strategy_report.
--run-dir and --report-name are reserved launcher arguments.
EOF
}

fail() {
  printf 'Error: %s\n' "$*" >&2
  exit 2
}

run_dir=""
report_name=""
slurm_cpus=8
slurm_memory=128G
slurm_time=06:00:00
slurm_partition=main
report_args=()

while (( $# > 0 )); do
  case "$1" in
    --run-dir)
      (( $# >= 2 )) || fail "--run-dir requires a value"
      run_dir=$2
      shift 2
      ;;
    --report-name)
      (( $# >= 2 )) || fail "--report-name requires a value"
      report_name=$2
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

[[ -n "$run_dir" ]] || fail "--run-dir is required"
[[ "$run_dir" = /* ]] || fail "--run-dir must be an absolute path visible to compute nodes"
[[ -d "$run_dir" ]] || fail "run directory does not exist: $run_dir"
[[ -n "$report_name" ]] || fail "--report-name is required"
[[ "$report_name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || fail \
  "--report-name may contain only letters, digits, dot, underscore, and hyphen"
[[ "$slurm_cpus" =~ ^[1-9][0-9]*$ ]] || fail "--slurm-cpus must be a positive integer"
[[ -n "$slurm_memory" ]] || fail "--slurm-memory must not be empty"
[[ -n "$slurm_time" ]] || fail "--slurm-time must not be empty"
[[ -n "$slurm_partition" ]] || fail "--slurm-partition must not be empty"
command -v sbatch >/dev/null || fail "sbatch was not found; run this launcher on the Slurm controller"

for argument in "${report_args[@]}"; do
  case "$argument" in
    --run-dir|--run-dir=*|--report-name|--report-name=*)
      fail "$argument is managed by the launcher and must appear before --"
      ;;
  esac
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/../.." && pwd)
batch_script="$script_dir/strategy_report.sbatch"
[[ -f "$batch_script" ]] || fail "missing Slurm batch script: $batch_script"

git -C "$project_root" diff --quiet || fail "tracked working-tree changes must be committed before submission"
git -C "$project_root" diff --cached --quiet || fail "staged changes must be committed before submission"
git_commit=$(git -C "$project_root" rev-parse HEAD)

log_dir="$run_dir/reports/slurm"
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
  "$run_dir" \
  "$report_name" \
  "$git_commit" \
  "$project_root" \
  "${report_args[@]}")

printf 'Submitted report job %s\n' "$job_id"
printf 'Report: %s/reports/%s.html\n' "$run_dir" "${report_name%.html}"
printf 'Logs: %s/%s.%s.{out,err}\n' "$log_dir" "$report_name" "$job_id"
