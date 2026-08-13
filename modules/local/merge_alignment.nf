process MERGE_ALIGNMENT {
    tag "merge_alignment"

    input:
    path alignment_tasks, stageAs: "source/alignment_tasks.tsv.gz"
    path taxonomy, stageAs: "source/taxonomy.tsv.gz"
    path taxonomy_failures, stageAs: "source/taxonomy_failures.tsv.gz"
    path taxonomy_summary, stageAs: "source/taxonomy_summary.tsv.gz"
    path source_genes, stageAs: "source/genes.tsv.gz"
    path source_target_features, stageAs: "source/target_features.tsv.gz"
    path result_dirs, stageAs: 'partitions/*'
    val expected_strategies
    path merge_script

    output:
    path "manifest.json", emit: manifest
    path "alignment_tasks.tsv.gz", optional: true, emit: alignment_tasks
    path "taxonomy.tsv.gz", optional: true, emit: taxonomy
    path "taxonomy_failures.tsv.gz", optional: true, emit: taxonomy_failures
    path "taxonomy_summary.tsv.gz", emit: taxonomy_summary
    path "ortholog_alignment_summary.tsv.gz", optional: true, emit: summaries
    path "strategy_summary.tsv.gz", emit: strategy_summary
    path "alignment_segments.tsv.gz", optional: true, emit: segments
    path "feature_coverage.tsv.gz", emit: feature_coverage
    path "alignment_events.tsv.gz", optional: true, emit: events
    path "event_ortholog_support.tsv.gz", optional: true, emit: event_ortholog_support
    path "snv_site_depth.tsv.gz", optional: true, emit: snv_site_depth
    path "snv_taxonomic_depth.tsv.gz", optional: true, emit: snv_taxonomic_depth
    path "snv_alt_taxonomic_support.tsv.gz", optional: true, emit: snv_alt_taxonomic_support
    path "failures.tsv.gz", emit: failures
    path "native", optional: true, emit: native_outputs

    script:
    def outputProfile = params.stage == 'all' ? "report-input" : "full"
    """
    python3 "${merge_script}" \\
        --alignment-tasks "${alignment_tasks}" \\
        --taxonomy "${taxonomy}" \\
        --taxonomy-failures "${taxonomy_failures}" \\
        --taxonomy-summary "${taxonomy_summary}" \\
        --source-genes "${source_genes}" \\
        --source-target-features "${source_target_features}" \\
        --result-root partitions \\
        --expected-strategies "${expected_strategies}" \\
        --output-profile "${outputProfile}" \\
        --outdir .
    """
}
