process ALIGN_ENSEMBL_COMPARA_MAF_CHUNK {
    tag { meta.id }

    input:
    tuple val(meta), path(chunk_task_dir)
    path align_script

    output:
    tuple val(meta), path("align_ensembl_compara_maf_${meta.id}/gene_results/gene_*"), emit: ensembl_compara_maf_gene_fragments

    script:
    def resultDir = "align_ensembl_compara_maf_${meta.id}"
    """
    export PYTHONPATH="${projectDir}/bin:\${PYTHONPATH:-}"
    python3 "${align_script}" \\
        --chunk-task-dir "${chunk_task_dir}" \\
        --outdir "${resultDir}" \\
        --strategy precomputed_ensembl_92_mammals_epo_extended \\
        --release "${params.ensembl_compara_maf_release}" \\
        --species-set "${params.ensembl_compara_maf_species_set}" \\
        --method "${params.ensembl_compara_maf_method}" \\
        --timeout "${params.ensembl_compara_maf_timeout_seconds}" \\
        --retries "${params.ensembl_compara_maf_retries}" \\
        --retry-base-seconds "${params.ensembl_compara_maf_retry_base_seconds}" \\
        --retry-max-seconds "${params.ensembl_compara_maf_retry_max_seconds}"
    """
}
