process MERGE_ALIGNMENT_PARTITION {
    tag { meta.partition_id }

    input:
    tuple val(meta), path(result_dirs, stageAs: 'results/*')
    val expected_strategies
    path taxonomy_presets
    path merge_script

    output:
    tuple val(meta), path("${meta.partition_id}"), emit: partition_dirs

    script:
    def outputProfile = params.stage == 'all' ? "annotation-input" : "full"
    def compactEventsArg = (params.compact_alignment_events || outputProfile == 'annotation-input') ? "--compact-events" : ""
    def taxonomyArg = outputProfile == 'annotation-input' ? "--taxonomy-presets \"${taxonomy_presets}\"" : ""
    """
    python3 "${merge_script}" \\
        --result-root results \\
        --partition-id "${meta.partition_id}" \\
        --expected-gene-ids "${meta.gene_ids.join(',')}" \\
        --expected-strategies "${expected_strategies}" \\
        --output-profile "${outputProfile}" \\
        --outdir "${meta.partition_id}" \\
        ${compactEventsArg} \\
        ${taxonomyArg}
    """
}
