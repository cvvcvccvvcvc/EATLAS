process ALIGN_NUCMER_COMPARATOR {
    tag { task_dir.baseName }

    input:
    path task_dir
    path nucmer_script

    output:
    path "align_nucmer", emit: nucmer_result_dirs

    script:
    """
    python3 "${nucmer_script}" \\
        --task-dir "${task_dir}" \\
        --outdir "align_nucmer" \\
        --nucmer-bin "${params.nucmer_bin}" \\
        --show-coords-bin "${params.show_coords_bin}" \\
        --show-snps-bin "${params.show_snps_bin}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}
