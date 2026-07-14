process MERGE_ALIGNMENT_PARTITION {
    tag { meta.partition_id }

    input:
    tuple val(meta), path(result_dirs, stageAs: 'results/*')
    path merge_script

    output:
    tuple val(meta), path("${meta.partition_id}"), emit: partition_dirs

    script:
    def compactEventsArg = params.compact_alignment_events ? "--compact-events" : ""
    """
    python3 "${merge_script}" \\
        --result-root results \\
        --partition-id "${meta.partition_id}" \\
        --outdir "${meta.partition_id}" \\
        ${compactEventsArg}
    """
}
