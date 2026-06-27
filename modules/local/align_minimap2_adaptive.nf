process ALIGN_MINIMAP2_ADAPTIVE {
    tag { task_dir.baseName }

    input:
    path task_dir
    path minimap2_script

    output:
    path "align_minimap2_taxonomy_adaptive", emit: adaptive_result_dirs

    script:
    """
    python3 "${minimap2_script}" \\
        --task-dir "${task_dir}" \\
        --outdir "align_minimap2_taxonomy_adaptive" \\
        --strategy minimap2_taxonomy_adaptive \\
        --mode adaptive \\
        --minimap2-bin "${params.minimap2_bin}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}
