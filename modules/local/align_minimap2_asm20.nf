process ALIGN_MINIMAP2_ASM20 {
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir)
    path minimap2_script

    output:
    tuple val(meta), path("align_minimap2_asm20_${meta.id}"), emit: asm20_result_dirs

    script:
    def resultDir = "align_minimap2_asm20_${meta.id}"
    """
    python3 "${minimap2_script}" \\
        --task-dir "${task_dir}" \\
        --outdir "${resultDir}" \\
        --strategy minimap2_asm20 \\
        --mode fixed \\
        --fixed-preset asm20 \\
        --minimap2-bin "${params.minimap2_bin}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}
