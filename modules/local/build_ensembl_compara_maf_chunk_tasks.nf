process BUILD_ENSEMBL_COMPARA_MAF_CHUNK_TASKS {
    tag "ensembl_compara_maf_chunk_tasks"

    input:
    path genes_tsv
    path maf_manifest
    path chunk_task_script, stageAs: 'bin/prepare_ensembl_compara_maf_chunk_tasks.py'
    path bin_sources, stageAs: 'bin/*'

    output:
    path "maf_chunk_tasks/*", emit: chunk_task_dirs
    path "maf_chunk_tasks.tsv.gz", emit: chunk_tasks
    path "maf_chunk_tasks_manifest.json", emit: manifest

    script:
    """
    python3 -m bin.prepare_ensembl_compara_maf_chunk_tasks \\
        --maf-manifest "${maf_manifest}" \\
        --outdir . \\
        --genes-tsv "${genes_tsv}"
    """
}
