process ALIGN_MINIMAP2 {
    label 'task_scratch'
    tag { "${meta.id}:${meta.strategy}" }

    input:
    tuple val(meta), path(task_dir), path(source_target_fasta, stageAs: 'source_target.fa.gz'), path(source_ortholog_fasta, stageAs: 'source_ortholog.fa.gz')
    path minimap2_script
    path alignment_table_schema

    output:
    tuple val(meta), path("align_${meta.strategy}_${meta.id}"), emit: result_dirs

    script:
    def resultDir = "align_${meta.strategy}_${meta.id}"
    """
    python3 "${minimap2_script}" \\
        --task-dir "${task_dir}" \\
        --source-target-fasta "${source_target_fasta}" \\
        --source-ortholog-fasta "${source_ortholog_fasta}" \\
        --outdir "${resultDir}" \\
        --strategy "${meta.strategy}" \\
        --preset "${meta.preset}" \\
        --pseudoread-len "${meta.pseudoread_len}" \\
        --pseudoread-step "${meta.pseudoread_step}" \\
        --threads "${task.cpus}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}
