process ALIGN_MINIMAP2_ASM10 {
    label 'task_scratch'
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir), path(source_target_fasta, stageAs: 'source_target.fa.gz'), path(source_ortholog_fasta, stageAs: 'source_ortholog.fa.gz')
    path minimap2_script

    output:
    tuple val(meta), path("align_minimap2_asm10_${meta.id}"), emit: asm10_result_dirs

    script:
    def resultDir = "align_minimap2_asm10_${meta.id}"
    """
    python3 "${minimap2_script}" \\
        --task-dir "${task_dir}" \\
        --source-target-fasta "${source_target_fasta}" \\
        --source-ortholog-fasta "${source_ortholog_fasta}" \\
        --outdir "${resultDir}" \\
        --strategy minimap2_asm10 \\
        --mode fixed \\
        --fixed-preset asm10 \\
        --minimap2-bin "${params.minimap2_bin}" \\
        --target-features "${task_dir}/target_features.tsv.gz" \\
        --keep-native "${params.keep_native_alignments}"
    """
}
