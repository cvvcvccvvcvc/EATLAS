process ALIGN_ENSEMBL_COMPARA_MAF_CHUNK {
    label 'task_scratch'
    tag { meta.id }

    input:
    tuple val(meta), path(chunk_task_dir)
    path align_script, stageAs: 'bin/run_ensembl_compara_maf_chunk_alignment.py'
    path bin_sources, stageAs: 'bin/*'

    output:
    tuple val(meta), path("align_ensembl_compara_maf_${meta.id}/gene_results/gene_*"), emit: ensembl_compara_maf_gene_fragments

    script:
    def resultDir = "align_ensembl_compara_maf_${meta.id}"
    """
    python3 -m bin.run_ensembl_compara_maf_chunk_alignment \\
        --chunk-task-dir "${chunk_task_dir}" \\
        --outdir "${resultDir}"
    """
}
