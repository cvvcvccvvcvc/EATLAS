process ALIGN_NUCMER_COMPARATOR {
    label 'task_scratch'
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir), path(source_target_fasta, stageAs: 'source_target.fa.gz'), path(source_ortholog_fasta, stageAs: 'source_ortholog.fa.gz')
    path nucmer_script

    output:
    tuple val(meta), path("align_nucmer_${meta.id}"), emit: nucmer_result_dirs

    script:
    def resultDir = "align_nucmer_${meta.id}"
    """
    python3 "${nucmer_script}" \\
        --task-dir "${task_dir}" \\
        --source-target-fasta "${source_target_fasta}" \\
        --source-ortholog-fasta "${source_ortholog_fasta}" \\
        --outdir "${resultDir}" \\
        --threads "${task.cpus}" \\
        --target-features "${task_dir}/target_features.tsv.gz" \\
        --keep-native "${params.keep_native_alignments}"
    """
}
