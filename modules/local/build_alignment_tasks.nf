process BUILD_ALIGNMENT_TASKS {
    tag "alignment_tasks"

    input:
    path genes_tsv
    path orthologs_tsv
    path fetch_manifest
    path sequences_dir, stageAs: 'sequences'
    path taxonomy_presets
    path prepare_script

    output:
    path "alignment_tasks.tsv.gz", emit: alignment_tasks
    path "tasks/task_*", emit: task_dirs

    script:
    """
    python3 "${prepare_script}" \\
        --genes-tsv "${genes_tsv}" \\
        --orthologs-tsv "${orthologs_tsv}" \\
        --fetch-manifest "${fetch_manifest}" \\
        --taxonomy-presets "${taxonomy_presets}" \\
        --outdir . \\
        --sequences-dir sequences
    """
}
