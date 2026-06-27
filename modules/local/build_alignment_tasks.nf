process BUILD_ALIGNMENT_TASKS {
    tag "alignment_tasks"

    input:
    path genes_tsv
    path orthologs_tsv
    path "sequences/*", stageAs: 'sequences/*'
    path taxonomy_presets
    path prepare_script

    output:
    path "alignment_tasks.tsv.gz", emit: alignment_tasks
    path "tasks/task_*", emit: task_dirs

    script:
    """
    targetArgs=\$(find sequences/*/targets -name "*.fa.gz" | sed 's/^/--target-fasta /' | tr '\\n' ' ')
    orthologArgs=\$(find sequences/*/orthologs -name "*.fa.gz" | sed 's/^/--ortholog-fasta /' | tr '\\n' ' ')
    python3 "${prepare_script}" \\
        --genes-tsv "${genes_tsv}" \\
        --orthologs-tsv "${orthologs_tsv}" \\
        --taxonomy-presets "${taxonomy_presets}" \\
        --outdir . \\
        \$targetArgs \\
        \$orthologArgs
    """
}
