process BUILD_ENSEMBL_COMPARA_MAF_CHUNK_TASKS {
    tag "ensembl_compara_maf_chunk_tasks"

    input:
    path task_dirs
    path maf_manifest
    path chunk_task_script

    output:
    path "maf_chunk_tasks/*", emit: chunk_task_dirs
    path "maf_chunk_tasks.tsv.gz", emit: chunk_tasks
    path "maf_chunk_tasks_manifest.json", emit: manifest

    script:
    def taskDirList = task_dirs instanceof List ? task_dirs : [task_dirs]
    def taskDirArgs = taskDirList.collect { "--task-dir \"${it}\"" }.join(' ')
    """
    python3 "${chunk_task_script}" \\
        --maf-manifest "${maf_manifest}" \\
        --outdir . \\
        --strategy precomputed_ensembl_92_mammals_epo_extended \\
        ${taskDirArgs}
    """
}
