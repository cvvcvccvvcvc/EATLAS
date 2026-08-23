process ALIGN_BWA_PSEUDOREADS {
    label 'task_scratch'
    tag { "${meta.id}:${meta.strategy}" }

    input:
    tuple val(meta), path(task_dir), path(source_target_fasta, stageAs: 'source_target.fa.gz'), path(source_ortholog_fasta, stageAs: 'source_ortholog.fa.gz')
    path bwa_script, stageAs: 'bin/run_bwa_pseudoreads.py'
    path bin_sources, stageAs: 'bin/*'

    output:
    tuple val(meta), path("align_${meta.strategy}_${meta.id}"), emit: bwa_result_dirs

    script:
    def resultDir = "align_${meta.strategy}_${meta.id}"
    """
    python3 -m bin.run_bwa_pseudoreads \\
        --task-dir "${task_dir}" \\
        --source-target-fasta "${source_target_fasta}" \\
        --source-ortholog-fasta "${source_ortholog_fasta}" \\
        --outdir "${resultDir}" \\
        --strategy "${meta.strategy}" \\
        --pseudoread-len "${meta.pseudoread_len}" \\
        --pseudoread-step "${meta.pseudoread_step}" \\
        --pseudoread-phred "${meta.pseudoread_phred}" \\
        --threads "${task.cpus}"
    """
}
