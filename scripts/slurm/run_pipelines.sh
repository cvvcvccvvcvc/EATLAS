#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/slurm/run_pipelines.sh \
    --results-root /absolute/path/to/results/group \
    --expected-commit FULL_GIT_COMMIT \
    [--alignment-strategies strategy_a,strategy_b] \
    [--fetch-max-forks N] [--alignment-max-forks N] \
    [--annotation-max-forks N] \
    /absolute/path/to/batch_001.txt [/absolute/path/to/batch_002.txt ...]

Runs one or more complete pipelines sequentially on Slurm. Each input filename
becomes a run name below --results-root. Repeating the same command skips
completed runs and resumes the first incomplete run.
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

inspect_manifest() {
  python3 - "$1" <<'PY'
import hashlib
import json
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
parameters = payload.get("parameters")
if not isinstance(parameters, dict):
    raise ValueError(f"run manifest has no parameter map: {path}")

if payload.get("status") == "complete":
    descriptor = payload.get("evidence_inventory")
    if (
        not isinstance(descriptor, dict)
        or set(descriptor) != {"path", "schema_version", "size_bytes", "sha256"}
        or descriptor.get("path") != "evidence_inventory.json"
        or descriptor.get("schema_version") != 1
        or isinstance(descriptor.get("size_bytes"), bool)
        or not isinstance(descriptor.get("size_bytes"), int)
        or descriptor["size_bytes"] < 0
        or not isinstance(descriptor.get("sha256"), str)
        or re.fullmatch(r"[0-9a-f]{64}", descriptor["sha256"]) is None
    ):
        raise ValueError(f"completed run has an invalid evidence inventory descriptor: {path}")
    inventory_path = path.parent / descriptor["path"]
    if inventory_path.is_symlink() or not inventory_path.is_file():
        raise ValueError(f"completed run has no evidence inventory: {inventory_path}")
    inventory_bytes = inventory_path.read_bytes()
    if (
        len(inventory_bytes) != descriptor["size_bytes"]
        or hashlib.sha256(inventory_bytes).hexdigest() != descriptor["sha256"]
    ):
        raise ValueError(
            f"completed run evidence inventory differs from its descriptor: {inventory_path}"
        )

values = (
    payload.get("schema_version"),
    payload.get("pipeline"),
    payload.get("status"),
    payload.get("success"),
    payload.get("exit_status"),
    payload.get("session_id"),
    payload.get("git_commit"),
    payload.get("git_dirty"),
    parameters.get("ids_file"),
    parameters.get("outdir"),
    parameters.get("alignment_strategies"),
    parameters.get("fetch_max_forks"),
    parameters.get("alignment_max_forks"),
    parameters.get("annotation_max_forks"),
)
print("\x1f".join("" if value is None else str(value).lower() if isinstance(value, bool) else str(value) for value in values))
PY
}

results_root=""
expected_commit=""
alignment_strategies=""
fetch_max_forks=""
alignment_max_forks=""
annotation_max_forks=""
ids_files=()

while (( $# > 0 )); do
  case "$1" in
    --results-root)
      (( $# >= 2 )) || fail "--results-root requires a value"
      results_root=$2
      shift 2
      ;;
    --expected-commit)
      (( $# >= 2 )) || fail "--expected-commit requires a value"
      expected_commit=$2
      shift 2
      ;;
    --alignment-strategies)
      (( $# >= 2 )) || fail "--alignment-strategies requires a value"
      alignment_strategies=$2
      shift 2
      ;;
    --fetch-max-forks)
      (( $# >= 2 )) || fail "--fetch-max-forks requires a value"
      fetch_max_forks=$2
      shift 2
      ;;
    --alignment-max-forks)
      (( $# >= 2 )) || fail "--alignment-max-forks requires a value"
      alignment_max_forks=$2
      shift 2
      ;;
    --annotation-max-forks)
      (( $# >= 2 )) || fail "--annotation-max-forks requires a value"
      annotation_max_forks=$2
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    -*)
      fail "unknown launcher argument: $1"
      ;;
    *)
      ids_files+=("$1")
      shift
      ;;
  esac
done

[[ -n "$results_root" ]] || fail "--results-root is required"
[[ "$results_root" = /* ]] || fail "--results-root must be an absolute path"
[[ "$expected_commit" =~ ^[0-9a-f]{40}$ ]] || fail \
  "--expected-commit must be a full 40-character Git commit"
(( ${#ids_files[@]} > 0 )) || fail "at least one IDs file is required"
[[ -z "$alignment_strategies" || "$alignment_strategies" != *[[:space:]]* ]] || fail \
  "--alignment-strategies must be a comma-separated value without spaces"
for value in "$fetch_max_forks" "$alignment_max_forks" "$annotation_max_forks"; do
  [[ -z "$value" || "$value" =~ ^[1-9][0-9]*$ ]] || fail \
    "fork limits must be positive integers"
done

resolved_results_root=$(canonical_destination "$results_root") || fail \
  "--results-root cannot be resolved as a directory destination: $results_root"
run_group=$(basename "$resolved_results_root")
[[ "$run_group" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || fail \
  "--results-root basename may contain only letters, digits, dot, underscore, and hyphen"

run_names=()
resolved_ids_files=()
for ids_file in "${ids_files[@]}"; do
  [[ "$ids_file" = /* ]] || fail "IDs file must be an absolute path: $ids_file"
  [[ -s "$ids_file" ]] || fail "IDs file does not exist or is empty: $ids_file"
  resolved_ids_file=$(cd "$(dirname "$ids_file")" && printf '%s/%s\n' "$PWD" "$(basename "$ids_file")")
  filename=$(basename "$resolved_ids_file")
  case "$filename" in
    *.txt|*.ids|*.tsv) run_name=${filename%.*} ;;
    *) run_name=$filename ;;
  esac
  [[ "$run_name" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || fail \
    "input filename does not produce a safe run name: $filename"
  if [[ ${run_names[0]+present} ]]; then
    for existing_name in "${run_names[@]}"; do
      [[ "$existing_name" != "$run_name" ]] || fail "duplicate run name: $run_name"
    done
  fi
  run_names+=("$run_name")
  resolved_ids_files+=("$resolved_ids_file")
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/../.." && pwd)
git_status=$(git -C "$project_root" status --porcelain=v1 --untracked-files=normal) || fail \
  "cannot inspect repository status: $project_root"
[[ -z "$git_status" ]] || fail "pipeline launch requires a clean working tree"
git -C "$project_root" fetch origin main >/dev/null || fail "cannot fetch authoritative origin/main"
actual_commit=$(git -C "$project_root" rev-parse HEAD) || fail "cannot resolve repository HEAD"
origin_commit=$(git -C "$project_root" rev-parse origin/main) || fail "cannot resolve origin/main"
[[ "$actual_commit" = "$expected_commit" ]] || fail \
  "cluster HEAD $actual_commit does not match expected commit $expected_commit"
[[ "$origin_commit" = "$expected_commit" ]] || fail \
  "fetched origin/main $origin_commit does not match expected commit $expected_commit"

cluster_env="$HOME/.gaph_v2_cluster_env.sh"
[[ -f "$cluster_env" ]] || fail "cluster environment file not found: $cluster_env"
source "$cluster_env"
: "${GAPH_ROOT:?GAPH_ROOT is not set by $cluster_env}"
work_base=${GAPH_WORK_DIR:-"$GAPH_ROOT/work"}
work_root="$work_base/$run_group"
command -v micromamba >/dev/null || fail "micromamba was not found"
command -v python3 >/dev/null || fail "python3 was not found"
mkdir -p "$resolved_results_root" "$work_root"

for index in "${!resolved_ids_files[@]}"; do
  ids_file=${resolved_ids_files[$index]}
  run_name=${run_names[$index]}
  run_dir="$resolved_results_root/$run_name"
  manifest="$run_dir/run_manifest.json"
  resume_session=""
  run_alignment_strategies=$alignment_strategies
  run_fetch_max_forks=$fetch_max_forks
  run_alignment_max_forks=$alignment_max_forks
  run_annotation_max_forks=$annotation_max_forks

  git_status=$(git -C "$project_root" status --porcelain=v1 --untracked-files=normal) || fail \
    "cannot inspect repository status: $project_root"
  [[ -z "$git_status" ]] || fail "pipeline launch requires a clean working tree"
  actual_commit=$(git -C "$project_root" rev-parse HEAD) || fail "cannot resolve repository HEAD"
  [[ "$actual_commit" = "$expected_commit" ]] || fail \
    "cluster HEAD changed during the pipeline series"

  if [[ -e "$manifest" ]]; then
    [[ -f "$manifest" ]] || fail "run manifest is not a file: $manifest"
    manifest_info=$(inspect_manifest "$manifest") || fail "cannot read run manifest: $manifest"
    IFS=$'\x1f' read -r schema pipeline status success exit_status session_id \
      manifest_commit git_dirty manifest_ids manifest_outdir manifest_strategies \
      manifest_fetch_forks manifest_alignment_forks manifest_annotation_forks \
      <<< "$manifest_info"
    [[ "$schema" = 3 && "$pipeline" = gaph_v2 ]] || fail \
      "unsupported run manifest: $manifest"
    [[ "$manifest_commit" = "$expected_commit" && "$git_dirty" = false ]] || fail \
      "run $run_name was not created from the expected clean commit"
    [[ "$manifest_ids" = "$ids_file" && "$manifest_outdir" = "$run_dir" ]] || fail \
      "run $run_name does not match its current input or result path"

    if [[ "$status" = complete && "$success" = true && "$exit_status" = 0 ]]; then
      printf 'Skipping completed run %s\n' "$run_name"
      continue
    fi
    [[ "$status" = running || "$status" = failed ]] || fail \
      "run $run_name has invalid status: $status"
    [[ -n "$session_id" ]] || fail "run $run_name has no resumable Nextflow session"
    [[ -d "$work_root" ]] || fail "resume work directory is missing: $work_root"

    if [[ -n "$alignment_strategies" && "$alignment_strategies" != "$manifest_strategies" ]]; then
      fail "--alignment-strategies differs from the incomplete run $run_name"
    fi
    if [[ -n "$fetch_max_forks" && "$fetch_max_forks" != "$manifest_fetch_forks" ]]; then
      fail "--fetch-max-forks differs from the incomplete run $run_name"
    fi
    if [[ -n "$alignment_max_forks" && "$alignment_max_forks" != "$manifest_alignment_forks" ]]; then
      fail "--alignment-max-forks differs from the incomplete run $run_name"
    fi
    if [[ -n "$annotation_max_forks" && "$annotation_max_forks" != "$manifest_annotation_forks" ]]; then
      fail "--annotation-max-forks differs from the incomplete run $run_name"
    fi
    run_alignment_strategies=$manifest_strategies
    run_fetch_max_forks=$manifest_fetch_forks
    run_alignment_max_forks=$manifest_alignment_forks
    run_annotation_max_forks=$manifest_annotation_forks
    resume_session=$session_id
  fi

  mkdir -p "$run_dir/reports/nextflow"
  pipeline_command=(
    micromamba run -p "$GAPH_ROOT/envs/controller"
    nextflow -log "$run_dir/reports/nextflow/nextflow.log"
    run "$project_root"
    -profile slurm
    --ids_file "$ids_file"
    --outdir "$run_dir"
    -work-dir "$work_root"
  )
  [[ -z "$resume_session" ]] || pipeline_command+=(-resume "$resume_session")
  [[ -z "$run_alignment_strategies" ]] || pipeline_command+=(--alignment_strategies "$run_alignment_strategies")
  [[ -z "$run_fetch_max_forks" ]] || pipeline_command+=(--fetch_max_forks "$run_fetch_max_forks")
  [[ -z "$run_alignment_max_forks" ]] || pipeline_command+=(--alignment_max_forks "$run_alignment_max_forks")
  [[ -z "$run_annotation_max_forks" ]] || pipeline_command+=(--annotation_max_forks "$run_annotation_max_forks")

  printf '%s run %s (%s/%s)\n' \
    "$([[ -n "$resume_session" ]] && printf 'Resuming' || printf 'Starting')" \
    "$run_name" "$((index + 1))" "${#resolved_ids_files[@]}"
  cd "$project_root"
  "${pipeline_command[@]}"

  [[ -f "$manifest" ]] || fail "pipeline exited successfully without a run manifest: $run_name"
  manifest_info=$(inspect_manifest "$manifest") || fail "cannot read completed run manifest: $manifest"
  IFS=$'\x1f' read -r _ _ status success exit_status _ manifest_commit git_dirty _ _ _ _ _ _ \
    <<< "$manifest_info"
  [[ "$status" = complete && "$success" = true && "$exit_status" = 0 ]] || fail \
    "pipeline exited successfully but run $run_name is not complete"
  [[ "$manifest_commit" = "$expected_commit" && "$git_dirty" = false ]] || fail \
    "completed run $run_name has unexpected Git provenance"
done

printf 'All %s pipeline run(s) completed successfully\n' "${#resolved_ids_files[@]}"
