#!/usr/bin/env bash

set -euo pipefail

if (( $# < 4 )); then
  printf 'Expected MODE RUN_DIR GIT_COMMIT PROJECT_ROOT [mode arguments...]\n' >&2
  exit 2
fi

mode=$1
run_dir=$2
expected_commit=$3
project_root=$4
shift 4

source "$HOME/.gaph_v2_cluster_env.sh"

: "${GAPH_ROOT:?GAPH_ROOT is not set by .gaph_v2_cluster_env.sh}"
analytics_python="$GAPH_ROOT/envs/analytics/bin/python"
[[ -x "$analytics_python" ]] || {
  printf 'Analytics Python is missing: %s\n' "$analytics_python" >&2
  exit 1
}

actual_commit=$(git -C "$project_root" rev-parse HEAD)
if [[ "$actual_commit" != "$expected_commit" ]]; then
  printf 'Repository commit changed after submission: expected %s, found %s\n' \
    "$expected_commit" "$actual_commit" >&2
  exit 1
fi

cd "$project_root"
artifact_dir="$run_dir/analytics/vep_consequences"
log_dir="$artifact_dir/slurm"
batch_script="$project_root/analytics/slurm/vep_annotation.sbatch"
run_id=$(basename "$run_dir")
job_tag=$(printf '%s' "$run_id" | tr -c 'A-Za-z0-9._-' '-')
job_tag=${job_tag:0:30}

case "$mode" in
  prepare)
    if (( $# != 6 )); then
      printf 'Prepare mode expects PARTITION_SIZE MAX_PARALLEL PARTITION CPUS MEMORY TIME\n' >&2
      exit 2
    fi
    partition_size=$1
    max_parallel=$2
    slurm_partition=$3
    slurm_cpus=$4
    slurm_memory=$5
    slurm_time=$6

    "$analytics_python" -m analytics.vep_annotation prepare \
      --run-dir "$run_dir" \
      --partition-size "$partition_size" >/dev/null

    read -r partition_count row_count < <(
      "$analytics_python" -c \
        'import json, sys; plan=json.load(open(sys.argv[1], encoding="utf-8")); print(plan["partition_count"], plan["row_count"])' \
        "$artifact_dir/plan.json"
    )
    [[ "$partition_count" =~ ^[1-9][0-9]*$ ]] || {
      printf 'Invalid VEP partition count: %s\n' "$partition_count" >&2
      exit 1
    }
    printf 'Prepared %s rows in %s partitions.\n' "$row_count" "$partition_count"

    annotation_job=$(sbatch --parsable \
      --job-name="gaph-vep-array-$job_tag" \
      --partition="$slurm_partition" \
      --array="1-${partition_count}%${max_parallel}" \
      --cpus-per-task="$slurm_cpus" \
      --mem="$slurm_memory" \
      --time="$slurm_time" \
      --output="$log_dir/annotate.%A_%a.out" \
      --error="$log_dir/annotate.%A_%a.err" \
      "$batch_script" annotate "$run_dir" "$expected_commit" "$project_root")
    annotation_job=${annotation_job%%;*}

    finalize_job=$(sbatch --parsable \
      --job-name="gaph-vep-finalize-$job_tag" \
      --partition="$slurm_partition" \
      --dependency="afterany:$annotation_job" \
      --cpus-per-task=2 \
      --mem=8G \
      --time=02:00:00 \
      --output="$log_dir/finalize.%j.out" \
      --error="$log_dir/finalize.%j.err" \
      "$batch_script" finalize "$run_dir" "$expected_commit" "$project_root")
    finalize_job=${finalize_job%%;*}

    printf 'Submitted annotation array %s (1-%s, at most %s concurrent).\n' \
      "$annotation_job" "$partition_count" "$max_parallel"
    printf 'Submitted finalizer %s after array completion.\n' "$finalize_job"
    ;;
  annotate)
    (( $# == 0 )) || {
      printf 'Annotate mode does not accept additional arguments.\n' >&2
      exit 2
    }
    : "${SLURM_ARRAY_TASK_ID:?Annotate mode requires a Slurm array task}"
    : "${GAPH_VEP_RELEASE:?GAPH_VEP_RELEASE is not set}"
    : "${GAPH_VEP_EXECUTABLE:?GAPH_VEP_EXECUTABLE is not set}"
    : "${GAPH_VEP_CACHE_DIR:?GAPH_VEP_CACHE_DIR is not set}"
    : "${GAPH_VEP_RESULT_CACHE_DIR:?GAPH_VEP_RESULT_CACHE_DIR is not set}"

    "$analytics_python" -m analytics.vep_annotation annotate \
      --run-dir "$run_dir" \
      --partition-index "$SLURM_ARRAY_TASK_ID" \
      --vep-backend local \
      --vep-release "$GAPH_VEP_RELEASE" \
      --vep-executable "$GAPH_VEP_EXECUTABLE" \
      --vep-cache-dir "$GAPH_VEP_CACHE_DIR" \
      --vep-result-cache-dir "$GAPH_VEP_RESULT_CACHE_DIR" \
      --vep-result-cache-tile-size-bp "${GAPH_VEP_RESULT_CACHE_TILE_SIZE_BP:-1000000}" \
      --vep-forks "${GAPH_VEP_FORKS:-4}" >/dev/null
    printf 'Completed VEP partition %s.\n' "$SLURM_ARRAY_TASK_ID"
    ;;
  finalize)
    (( $# == 0 )) || {
      printf 'Finalize mode does not accept additional arguments.\n' >&2
      exit 2
    }
    "$analytics_python" -m analytics.vep_annotation finalize \
      --run-dir "$run_dir" >/dev/null
    "$analytics_python" -c \
      'import json, sys; manifest=json.load(open(sys.argv[1], encoding="utf-8")); print("Finalized {} rows across {} partitions.".format(manifest["row_count"], manifest["partition_count"]))' \
      "$artifact_dir/manifest.json"
    ;;
  *)
    printf 'Unknown VEP Slurm mode: %s\n' "$mode" >&2
    exit 2
    ;;
esac
