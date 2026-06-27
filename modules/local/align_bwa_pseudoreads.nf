process ALIGN_BWA_PSEUDOREADS {
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir)
    path bwa_script

    output:
    tuple val(meta), path("align_bwa"), emit: bwa_result_dirs

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
