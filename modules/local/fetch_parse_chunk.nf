process FETCH_PARSE_CHUNK {
    label 'task_scratch'
    tag { meta.id }

    input:
    tuple val(meta), path(chunk_file)
    path fetch_script

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
        --datasets-bin "${params.datasets_bin}" \\
        --target-assembly-accession "${params.target_assembly_accession}" \\
        --target-assembly-name "${params.target_assembly_name}" \\
        --target-tax-id "${params.target_tax_id}" \\
        --request-stagger-seconds "${params.fetch_request_stagger_seconds}" \\
        --request-throttle-dir "${workflow.workDir}/.gaph/ncbi_fetch_throttle" \\
        --download-retries "${params.fetch_download_retries}" \\
        --download-retry-base-seconds "${params.fetch_download_retry_base_seconds}"
    """
}
