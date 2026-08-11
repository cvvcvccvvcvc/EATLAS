process FETCH_PARSE_CHUNK {
    label 'task_scratch'
    tag { meta.id }

    input:
    tuple val(meta), path(chunk_file)
    path fetch_script
    val request_throttle_dir

    output:
    tuple val(meta), path("fetch_*"), emit: chunk_dirs

    script:
    """
    chunk_name=\$(basename "${chunk_file}" .ids.txt)
    if [[ -f "${projectDir}/.env" ]]; then
        set -a
        source "${projectDir}/.env"
        set +a
    fi

    python3 "${fetch_script}" \\
        --ids-file "${chunk_file}" \\
        --outdir "fetch_\${chunk_name}" \\
        --request-throttle-dir "${request_throttle_dir}"
    """
}
