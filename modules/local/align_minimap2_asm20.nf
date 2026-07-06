process ALIGN_MINIMAP2_ASM20 {
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir), path(source_target_fasta, stageAs: 'source_target.fa.gz'), path(source_ortholog_fasta, stageAs: 'source_ortholog.fa.gz')
    path minimap2_script
    path target_features

    output:
    tuple val(meta), path("align_minimap2_asm20_${meta.id}"), emit: asm20_result_dirs

    script:
    def resultDir = "align_minimap2_asm20_${meta.id}"
    """
    python3 "${minimap2_script}" \\
        --task-dir "${task_dir}" \\
        --source-target-fasta "${source_target_fasta}" \\
        --source-ortholog-fasta "${source_ortholog_fasta}" \\
        --outdir "${resultDir}" \\
        --strategy minimap2_asm20 \\
        --mode fixed \\
        --fixed-preset asm20 \\
        --minimap2-bin "${params.minimap2_bin}" \\
        --target-features "${target_features}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}
