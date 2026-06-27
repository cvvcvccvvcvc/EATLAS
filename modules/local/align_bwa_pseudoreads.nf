process ALIGN_BWA_PSEUDOREADS {
    tag { task_dir.baseName }

    input:
    path task_dir
    path bwa_script

    output:
    path "align_bwa", emit: bwa_result_dirs

    script:
    """
    cp -r "${task_dir}" align_bwa
    python3 "${bwa_script}" \\
        --task-dir "align_bwa" \\
        --bwa-bin "${params.bwa_bin}" \\
        --samtools-bin "${params.samtools_bin}" \\
        --varscan-jar "${params.varscan_jar}"
    """
}
