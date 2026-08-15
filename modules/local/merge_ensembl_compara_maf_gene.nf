process MERGE_ENSEMBL_COMPARA_MAF_GENE {
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir), path(fragment_dirs, stageAs: 'fragments/*')
    path merge_script
    path alignment_table_schema

    output:
    tuple val(meta), path("merge_ensembl_compara_maf_${meta.id}"), emit: gene_result_dirs

    script:
    def resultDir = "merge_ensembl_compara_maf_${meta.id}"
    """
    export PYTHONPATH="${projectDir}/bin:\${PYTHONPATH:-}"
    python3 "${merge_script}" \\
        --gene-id "${meta.gene_id}" \\
        --task-dir "${task_dir}" \\
        --fragment-root fragments \\
        --outdir "${resultDir}"
    """
}
