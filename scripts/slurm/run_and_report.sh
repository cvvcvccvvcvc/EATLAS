#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/slurm/run_and_report.sh \
    --ids-file /absolute/path/to/gene_ids.txt \
    --run-dir /absolute/path/to/results/run_name \
    --work-dir /absolute/path/to/work/run_name \
    --report-name report_name \
    [--alignment-strategies strategy_a,strategy_b] \
    [--slurm-cpus N] [--slurm-memory SIZE] \
    [--slurm-time D-HH:MM:SS] [--slurm-partition NAME] \
    [-- <analytics.strategy_report arguments>]

Runs the complete pipeline with Slurm and -resume. The report is submitted only
after Nextflow exits successfully. Arguments after -- are forwarded unchanged to
analytics/slurm/submit_strategy_report.sh.
EOF
}

fail() {
  printf 'Error: %s\n' "$*" >&2
  exit 2
}

ids_file=""
run_dir=""
work_dir=""
report_name=""
alignment_strategies=""
report_slurm_cpus=""
report_slurm_memory=""
report_slurm_time=""
report_slurm_partition=""
report_args=()

while (( $# > 0 )); do
  case "$1" in
    --ids-file)
      (( $# >= 2 )) || fail "--ids-file requires a value"
      ids_file=$2
      shift 2
      ;;
    --run-dir)
      (( $# >= 2 )) || fail "--run-dir requires a value"
      run_dir=$2
      shift 2
      ;;
    --work-dir)
      (( $# >= 2 )) || fail "--work-dir requires a value"
      work_dir=$2
      shift 2
      ;;
    --report-name)
      (( $# >= 2 )) || fail "--report-name requires a value"
      report_name=$2
      shift 2
      ;;
    --alignment-strategies)
      (( $# >= 2 )) || fail "--alignment-strategies requires a value"
      alignment_strategies=$2
      shift 2
      ;;
    --slurm-cpus)
      (( $# >= 2 )) || fail "--slurm-cpus requires a value"
      report_slurm_cpus=$2
      shift 2
      ;;
    --slurm-memory)
      (( $# >= 2 )) || fail "--slurm-memory requires a value"
      report_slurm_memory=$2
      shift 2
      ;;
    --slurm-time)
      (( $# >= 2 )) || fail "--slurm-time requires a value"
      report_slurm_time=$2
      shift 2
      ;;
    --slurm-partition)
      (( $# >= 2 )) || fail "--slurm-partition requires a value"
      report_slurm_partition=$2
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

[[ -n "$ids_file" ]] || fail "--ids-file is required"
[[ -n "$run_dir" ]] || fail "--run-dir is required"
[[ -n "$work_dir" ]] || fail "--work-dir is required"
[[ -n "$report_name" ]] || fail "--report-name is required"

[[ "$ids_file" = /* ]] || fail "--ids-file must be an absolute path"
[[ "$run_dir" = /* ]] || fail "--run-dir must be an absolute path"
[[ "$work_dir" = /* ]] || fail "--work-dir must be an absolute path"
[[ -f "$ids_file" ]] || fail "IDs file does not exist: $ids_file"
[[ "$report_name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || fail \
  "--report-name may contain only letters, digits, dot, underscore, and hyphen"
[[ -z "$report_slurm_cpus" || "$report_slurm_cpus" =~ ^[1-9][0-9]*$ ]] || fail \
  "--slurm-cpus must be a positive integer"

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/../.." && pwd)
report_launcher="$project_root/analytics/slurm/submit_strategy_report.sh"
[[ -f "$report_launcher" ]] || fail "missing report launcher: $report_launcher"

cluster_env="$HOME/.gaph_v2_cluster_env.sh"
[[ -f "$cluster_env" ]] || fail "cluster environment file not found: $cluster_env"
source "$cluster_env"
: "${GAPH_ROOT:?GAPH_ROOT is not set by $cluster_env}"
command -v micromamba >/dev/null || fail "micromamba was not found"

pipeline_command=(
  micromamba run -p "$GAPH_ROOT/envs/controller"
  nextflow run "$project_root"
  -profile slurm
  --ids_file "$ids_file"
  --outdir "$run_dir"
  -work-dir "$work_dir"
  -resume
)
if [[ -n "$alignment_strategies" ]]; then
  pipeline_command+=(--alignment_strategies "$alignment_strategies")
fi

report_command=(
  bash "$report_launcher"
  --run-dir "$run_dir"
  --report-name "$report_name"
)
if [[ -n "$report_slurm_cpus" ]]; then
  report_command+=(--slurm-cpus "$report_slurm_cpus")
fi
if [[ -n "$report_slurm_memory" ]]; then
  report_command+=(--slurm-memory "$report_slurm_memory")
fi
if [[ -n "$report_slurm_time" ]]; then
  report_command+=(--slurm-time "$report_slurm_time")
fi
if [[ -n "$report_slurm_partition" ]]; then
  report_command+=(--slurm-partition "$report_slurm_partition")
fi
if (( ${#report_args[@]} > 0 )); then
  report_command+=(-- "${report_args[@]}")
fi

printf 'Pipeline command:'
printf ' %q' "${pipeline_command[@]}"
printf '\n'
cd "$project_root"
"${pipeline_command[@]}"

printf 'Report submission command:'
printf ' %q' "${report_command[@]}"
printf '\n'
"${report_command[@]}"
