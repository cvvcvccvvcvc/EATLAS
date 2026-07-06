process ALIGN_BWA_PSEUDOREADS {
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir), path(source_target_fasta, stageAs: 'source_target.fa.gz'), path(source_ortholog_fasta, stageAs: 'source_ortholog.fa.gz')
    path bwa_script
    path bam_filtering_script
    val selected_strategies

    output:
    tuple val(meta), path("align_bwa_${meta.id}"), emit: bwa_result_dirs

    script:
    def resultDir = "align_bwa_${meta.id}"
    """
    export PYTHONPATH="\$PWD:\${PYTHONPATH:-}"
    python3 "${bwa_script}" \\
        --task-dir "${task_dir}" \\
        --source-target-fasta "${source_target_fasta}" \\
        --source-ortholog-fasta "${source_ortholog_fasta}" \\
        --outdir "${resultDir}" \\
        --strategies "${selected_strategies}" \\
        --bwa-bin "${params.bwa_bin}" \\
        --samtools-bin "${params.samtools_bin}" \\
        --varscan-jar "${params.varscan_jar}" \\
        --varscan-min-coverage "${params.bwa_varscan_min_coverage}" \\
        --varscan-min-reads2 "${params.bwa_varscan_min_reads2}" \\
        --varscan-min-var-freq "${params.bwa_varscan_min_var_freq}" \\
        --pseudoread-len "${params.bwa_pseudoread_len}" \\
        --pseudoread-step "${params.bwa_pseudoread_step}" \\
        --pseudoread-phred "${params.bwa_pseudoread_phred}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}
