process ALIGN_NUCMER_COMPARATOR {
    label 'task_scratch'
    tag { meta.id }

    input:
    tuple val(meta), path(task_dir), path(source_target_fasta, stageAs: 'source_target.fa.gz'), path(source_ortholog_fasta, stageAs: 'source_ortholog.fa.gz')
    path nucmer_script, stageAs: 'bin/run_nucmer_alignment.py'
    path bin_sources, stageAs: 'bin/*'

    output:
    tuple val(meta), path("align_nucmer_${meta.id}"), emit: nucmer_result_dirs

    script:
    def resultDir = "align_nucmer_${meta.id}"
    """
    python3 -m bin.run_nucmer_alignment \\
        --task-dir "${task_dir}" \\
        --source-target-fasta "${source_target_fasta}" \\
        --source-ortholog-fasta "${source_ortholog_fasta}" \\
        --outdir "${resultDir}" \\
        --threads "${task.cpus}"
    """
}
