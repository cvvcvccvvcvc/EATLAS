process MERGE_ALIGNMENT_PARTITION {
    tag { meta.partition_id }

    input:
    tuple val(meta), path(result_dirs, stageAs: 'results/*')
    path alignment_tasks
    val expected_strategies
    path merge_script
    path alignment_table_schema

    output:
    tuple val(meta), path("${meta.partition_id}"), emit: partition_dirs

    script:
    """
    export PYTHONPATH="\$PWD:\${PYTHONPATH:-}"
    python3 "${merge_script}" \\
        --result-root results \\
        --partition-id "${meta.partition_id}" \\
        --alignment-tasks "${alignment_tasks}" \\
        --expected-gene-ids "${meta.gene_ids.join(',')}" \\
        --expected-strategies "${expected_strategies}" \\
        --outdir "${meta.partition_id}"
    """
}
