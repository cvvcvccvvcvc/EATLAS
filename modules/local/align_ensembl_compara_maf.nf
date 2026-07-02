process ALIGN_ENSEMBL_COMPARA_MAF {
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir)
    path maf_manifest
    path align_script

    output:
    tuple val(meta), path("align_ensembl_compara_maf_${meta.id}"), emit: ensembl_compara_maf_result_dirs

    script:
    def resultDir = "align_ensembl_compara_maf_${meta.id}"
    """
    python3 "${align_script}" \\
        --task-dir "${task_dir}" \\
        --maf-manifest "${maf_manifest}" \\
        --outdir "${resultDir}" \\
        --strategy precomputed_ensembl_92_mammals_epo_extended \\
        --release "${params.ensembl_compara_maf_release}" \\
        --species-set "${params.ensembl_compara_maf_species_set}" \\
        --method "${params.ensembl_compara_maf_method}" \\
        --maf-cache-dir "${params.ensembl_compara_maf_cache_dir}"
    """
}
