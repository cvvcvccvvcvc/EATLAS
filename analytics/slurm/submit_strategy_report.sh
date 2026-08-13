#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  analytics/slurm/submit_strategy_report.sh \
    (--run-dir /absolute/path/to/run | --cohort-manifest /absolute/path/to/cohort.json) \
    --report-name report_name \
    [--cohort-root /absolute/path/to/cohort/output/root] \
    [--slurm-cpus 8] [--slurm-memory 128G] [--slurm-time 06:00:00] \
    [--slurm-partition main] [-- <analytics.strategy_report arguments>]

Arguments after -- are passed unchanged to analytics.strategy_report.
The input selector, --cohort-root, and --report-name are reserved launcher arguments.
EOF
}

fail() {
  printf 'Error: %s\n' "$*" >&2
  exit 2
}

run_dir=""
cohort_manifest=""
cohort_root=""
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
    --cohort-manifest)
      (( $# >= 2 )) || fail "--cohort-manifest requires a value"
      cohort_manifest=$2
      shift 2
      ;;
    --cohort-root)
      (( $# >= 2 )) || fail "--cohort-root requires a value"
      cohort_root=$2
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

if [[ -n "$run_dir" && -n "$cohort_manifest" ]]; then
  fail "--run-dir and --cohort-manifest are mutually exclusive"
fi
if [[ -z "$run_dir" && -z "$cohort_manifest" ]]; then
  fail "one of --run-dir or --cohort-manifest is required"
fi
if [[ -n "$run_dir" ]]; then
  [[ "$run_dir" = /* ]] || fail "--run-dir must be an absolute path visible to compute nodes"
  [[ -d "$run_dir" ]] || fail "run directory does not exist: $run_dir"
  [[ -z "$cohort_root" ]] || fail "--cohort-root can only be used with --cohort-manifest"
  source_kind=run
  source_path=$run_dir
  log_dir="$run_dir/reports/slurm"
else
  [[ "$cohort_manifest" = /* ]] || fail "--cohort-manifest must be an absolute path visible to compute nodes"
  [[ -f "$cohort_manifest" ]] || fail "cohort manifest does not exist: $cohort_manifest"
  if [[ -z "$cohort_root" ]]; then
    cohort_root="$(dirname "$cohort_manifest")/cohorts"
  fi
  [[ "$cohort_root" = /* ]] || fail "--cohort-root must be an absolute path visible to compute nodes"
  source_kind=cohort
  source_path=$cohort_manifest
  log_dir="$cohort_root/slurm"
fi
[[ -n "$report_name" ]] || fail "--report-name is required"
[[ "$report_name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || fail \
  "--report-name may contain only letters, digits, dot, underscore, and hyphen"
[[ "$slurm_cpus" =~ ^[1-9][0-9]*$ ]] || fail "--slurm-cpus must be a positive integer"
[[ -n "$slurm_memory" ]] || fail "--slurm-memory must not be empty"
[[ -n "$slurm_time" ]] || fail "--slurm-time must not be empty"
[[ -n "$slurm_partition" ]] || fail "--slurm-partition must not be empty"
command -v sbatch >/dev/null || fail "sbatch was not found; run this launcher on the Slurm controller"

if [[ "$source_kind" == run ]]; then
  vep_dir="$run_dir/analytics/vep_consequences"
  [[ -s "$vep_dir/manifest.json" && -s "$vep_dir/variant_annotations.vep.tsv.gz" ]] || fail \
    "missing finalized bulk VEP artifact under $vep_dir; run: bash analytics/slurm/submit_vep_annotation.sh --run-dir '$run_dir'"
fi

for argument in "${report_args[@]}"; do
  case "$argument" in
    --run-dir|--run-dir=*|--cohort-manifest|--cohort-manifest=*|--cohort-root|--cohort-root=*|--report-name|--report-name=*)
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
  "$source_kind" \
  "$source_path" \
  "$cohort_root" \
  "$report_name" \
  "$git_commit" \
  "$project_root" \
  "${report_args[@]}")

printf 'Submitted report job %s\n' "$job_id"
if [[ "$source_kind" == run ]]; then
  printf 'Report: %s/reports/%s.html\n' "$run_dir" "${report_name%.html}"
else
  printf 'Cohort report root: %s/<cohort-id>/reports/%s.html\n' "$cohort_root" "${report_name%.html}"
fi
printf 'Logs: %s/%s.%s.{out,err}\n' "$log_dir" "$report_name" "$job_id"
