process ALIGN_MINIMAP2_ASM10 {
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir)
    path minimap2_script

    output:
    tuple val(meta), path("align_minimap2_asm10"), emit: asm10_result_dirs

    script:
    """
    python3 "${minimap2_script}" \\
        --task-dir "${task_dir}" \\
        --outdir "align_minimap2_asm10" \\
        --strategy minimap2_asm10 \\
        --mode fixed \\
        --fixed-preset asm10 \\
        --minimap2-bin "${params.minimap2_bin}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}
