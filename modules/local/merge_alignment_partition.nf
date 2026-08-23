process MERGE_ALIGNMENT_PARTITION {
    tag { meta.partition_id }

    input:
    tuple val(meta), path(result_dirs, stageAs: 'results/*')
    path alignment_tasks
    val expected_strategies
    path merge_script, stageAs: 'bin/merge_alignment_results.py'
    path bin_sources, stageAs: 'bin/*'

    output:
    tuple val(meta), path("${meta.partition_id}"), emit: partition_dirs

    script:
    """
    python3 -m bin.merge_alignment_results \\
        --result-root results \\
        --partition-id "${meta.partition_id}" \\
        --alignment-tasks "${alignment_tasks}" \\
        --expected-gene-ids "${meta.gene_ids.join(',')}" \\
        --expected-strategies "${expected_strategies}" \\
        --outdir "${meta.partition_id}"
    """
}
