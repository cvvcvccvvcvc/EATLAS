process MERGE_ALIGNMENT {
    tag "merge_alignment"

    input:
    path alignment_tasks
    path taxonomy_presets
    path taxonomy_failures
    path target_features
    path result_dirs
    path merge_script

    output:
    path "manifest.json", emit: manifest
    path "alignment_tasks.tsv.gz", emit: alignment_tasks
    path "taxonomy_presets.tsv.gz", emit: taxonomy_presets
    path "taxonomy_failures.tsv.gz", emit: taxonomy_failures
    path "ortholog_alignment_summary.tsv.gz", emit: summaries
    path "strategy_summary.tsv.gz", emit: strategy_summary
    path "alignment_segments.tsv.gz", emit: segments
    path "feature_coverage.tsv.gz", emit: feature_coverage
    path "alignment_events.tsv.gz", emit: events
    path "failures.tsv.gz", emit: failures
    path "native", optional: true, emit: native_outputs

    script:
    def resultDirList = result_dirs instanceof List ? result_dirs : [result_dirs]
    def resultDirArgs = resultDirList.collect { "--result-dir \"${it}\"" }.join(' ')
    def compactEventsArg = params.compact_alignment_events ? "--compact-events" : ""
    """
    python3 "${merge_script}" \\
        --alignment-tasks "${alignment_tasks}" \\
        --taxonomy-presets "${taxonomy_presets}" \\
        --taxonomy-failures "${taxonomy_failures}" \\
        --target-features "${target_features}" \\
        --outdir . \\
        ${compactEventsArg} \\
        ${resultDirArgs}
    """
}
