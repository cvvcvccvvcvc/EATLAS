process ALIGN_BWA_PSEUDOREADS {
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir)
    path bwa_script
    path bam_filtering_script

    output:
    tuple val(meta), path("align_bwa_${meta.id}"), emit: bwa_result_dirs

    script:
    def resultDir = "align_bwa_${meta.id}"
    """
    cp -r "${task_dir}" "${resultDir}"
    export PYTHONPATH="\$PWD:\${PYTHONPATH:-}"
    python3 "${bwa_script}" \\
        --task-dir "${resultDir}" \\
        --bwa-bin "${params.bwa_bin}" \\
        --samtools-bin "${params.samtools_bin}" \\
        --varscan-jar "${params.varscan_jar}"
    """
}
