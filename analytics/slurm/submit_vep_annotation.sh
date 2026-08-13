#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  analytics/slurm/submit_vep_annotation.sh \
    --run-dir /absolute/path/to/completed-run \
    [--partition-size 250000] [--max-parallel 4] \
    [--slurm-cpus 4] [--slurm-memory 8G] [--slurm-time 01:00:00] \
    [--slurm-partition main]

Submits the complete resumable bulk-VEP chain:
  prepare -> bounded annotation array -> finalize
EOF
}

fail() {
  printf 'Error: %s\n' "$*" >&2
  exit 2
}

run_dir=""
partition_size=250000
max_parallel=4
slurm_cpus=4
slurm_memory=8G
slurm_time=01:00:00
slurm_partition=main

while (( $# > 0 )); do
  case "$1" in
    --run-dir)
      (( $# >= 2 )) || fail "--run-dir requires a value"
      run_dir=$2
      shift 2
      ;;
    --partition-size)
      (( $# >= 2 )) || fail "--partition-size requires a value"
      partition_size=$2
      shift 2
      ;;
    --max-parallel)
      (( $# >= 2 )) || fail "--max-parallel requires a value"
      max_parallel=$2
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
    *)
      fail "unknown argument '$1'"
      ;;
  esac
done

[[ -n "$run_dir" ]] || fail "--run-dir is required"
[[ "$run_dir" = /* ]] || fail "--run-dir must be an absolute path visible to compute nodes"
[[ -d "$run_dir" ]] || fail "run directory does not exist: $run_dir"
[[ -s "$run_dir/run_manifest.json" ]] || fail "missing run manifest: $run_dir/run_manifest.json"
[[ -s "$run_dir/annotation/variant_annotations.tsv.gz" ]] || fail \
  "missing annotation input: $run_dir/annotation/variant_annotations.tsv.gz"
[[ "$partition_size" =~ ^[1-9][0-9]*$ ]] || fail "--partition-size must be a positive integer"
[[ "$max_parallel" =~ ^[1-9][0-9]*$ ]] || fail "--max-parallel must be a positive integer"
[[ "$slurm_cpus" =~ ^[1-9][0-9]*$ ]] || fail "--slurm-cpus must be a positive integer"
[[ -n "$slurm_memory" ]] || fail "--slurm-memory must not be empty"
[[ -n "$slurm_time" ]] || fail "--slurm-time must not be empty"
[[ -n "$slurm_partition" ]] || fail "--slurm-partition must not be empty"
command -v sbatch >/dev/null || fail "sbatch was not found; run this launcher on the Slurm controller"

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/../.." && pwd)
batch_script="$script_dir/vep_annotation.sbatch"
[[ -f "$batch_script" ]] || fail "missing Slurm worker: $batch_script"

git -C "$project_root" diff --quiet || fail "tracked working-tree changes must be committed before submission"
git -C "$project_root" diff --cached --quiet || fail "staged changes must be committed before submission"
git_commit=$(git -C "$project_root" rev-parse HEAD)

read -r run_status run_success run_exit_status < <(python3 -c \
  'import json, sys; manifest=json.load(open(sys.argv[1], encoding="utf-8")); print(manifest.get("status", ""), str(manifest.get("success", False)).lower(), manifest.get("exit_status", ""))' \
  "$run_dir/run_manifest.json")
[[ "$run_status" == complete && "$run_success" == true && "$run_exit_status" == 0 ]] || fail \
  "run manifest is not successfully complete: $run_dir/run_manifest.json"

vep_manifest="$run_dir/analytics/vep_consequences/manifest.json"
vep_output="$run_dir/analytics/vep_consequences/variant_annotations.vep.tsv.gz"
if [[ -s "$vep_manifest" && -s "$vep_output" ]]; then
  python3 -c \
    'import json, os, sys; manifest=json.load(open(sys.argv[1], encoding="utf-8")); source=os.stat(sys.argv[2]); output=os.stat(sys.argv[3]); expected_source={"path": os.path.realpath(sys.argv[2]), "size_bytes": source.st_size, "mtime_ns": source.st_mtime_ns}; expected_output={"size_bytes": output.st_size, "mtime_ns": output.st_mtime_ns}; valid=manifest.get("status")=="complete" and manifest.get("source")==expected_source and manifest.get("output")==expected_output; raise SystemExit(0 if valid else 1)' \
    "$vep_manifest" "$run_dir/annotation/variant_annotations.tsv.gz" "$vep_output" || fail \
    "existing finalized bulk-VEP artifact does not match the current annotation input"
  printf 'Bulk VEP is already finalized and matches the current annotation input:\n'
  printf '  %s\n' "$vep_output"
  exit 0
fi

log_dir="$run_dir/analytics/vep_consequences/slurm"
mkdir -p "$log_dir"
run_id=$(basename "$run_dir")
job_tag=$(printf '%s' "$run_id" | tr -c 'A-Za-z0-9._-' '-')
job_tag=${job_tag:0:36}

job_id=$(sbatch --parsable \
  --job-name="gaph-vep-$job_tag" \
  --partition="$slurm_partition" \
  --cpus-per-task=2 \
  --mem=8G \
  --time=02:00:00 \
  --output="$log_dir/prepare.%j.out" \
  --error="$log_dir/prepare.%j.err" \
  "$batch_script" prepare \
  "$run_dir" "$git_commit" "$project_root" \
  "$partition_size" "$max_parallel" \
  "$slurm_partition" "$slurm_cpus" "$slurm_memory" "$slurm_time")

job_id=${job_id%%;*}
printf 'Submitted bulk-VEP preparation job %s\n' "$job_id"
printf 'The annotation-array and finalizer job IDs will appear in:\n'
printf '  %s/prepare.%s.out\n' "$log_dir" "$job_id"
printf 'Final artifact:\n'
printf '  %s\n' "$vep_output"
