process MERGE_ALIGNMENT_PARTITION {
    tag { meta.partition_id }

    input:
    tuple val(meta), path(result_dirs, stageAs: 'results/*')
    path alignment_tasks
    val expected_strategies
    path taxonomy
    path merge_script
    path feature_coverage_script
    path taxonomic_evidence_script

    output:
    tuple val(meta), path("${meta.partition_id}"), emit: partition_dirs

    script:
    def outputProfile = params.stage == 'all' ? "annotation-input" : "full"
    def taxonomyArg = outputProfile == 'annotation-input' ? "--taxonomy \"${taxonomy}\"" : ""
    """
    export PYTHONPATH="\$PWD:\${PYTHONPATH:-}"
    python3 "${merge_script}" \\
        --result-root results \\
        --partition-id "${meta.partition_id}" \\
        --alignment-tasks "${alignment_tasks}" \\
        --expected-gene-ids "${meta.gene_ids.join(',')}" \\
        --expected-strategies "${expected_strategies}" \\
        --output-profile "${outputProfile}" \\
        --outdir "${meta.partition_id}" \\
        ${taxonomyArg}
    """
}
