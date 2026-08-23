process MERGE_ENSEMBL_COMPARA_MAF_GENE {
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir), path(fragment_dirs, stageAs: 'fragments/*')
    path merge_script, stageAs: 'bin/merge_ensembl_compara_maf_gene.py'
    path bin_sources, stageAs: 'bin/*'

    output:
    tuple val(meta), path("merge_ensembl_compara_maf_${meta.id}"), emit: gene_result_dirs

    script:
    def resultDir = "merge_ensembl_compara_maf_${meta.id}"
    """
    python3 -m bin.merge_ensembl_compara_maf_gene \\
        --gene-id "${meta.gene_id}" \\
        --task-dir "${task_dir}" \\
        --fragment-root fragments \\
        --outdir "${resultDir}"
    """
}
