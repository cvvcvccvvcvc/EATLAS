process ALIGN_MINIMAP2_ASM20 {
    tag { task_dir.baseName }

    input:
    path task_dir
    path minimap2_script

    output:
    path "align_minimap2_asm20", emit: asm20_result_dirs

    script:
    """
    python3 "${minimap2_script}" \\
        --task-dir "${task_dir}" \\
        --outdir "align_minimap2_asm20" \\
        --strategy minimap2_asm20 \\
        --mode fixed \\
        --fixed-preset asm20 \\
        --minimap2-bin "${params.minimap2_bin}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}
