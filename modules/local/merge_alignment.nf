process MERGE_ALIGNMENT {
    tag "merge_alignment"

    input:
    path alignment_tasks, stageAs: "source/alignment_tasks.tsv.gz"
    path source_genes, stageAs: "source/genes.tsv.gz"
    path source_target_features, stageAs: "source/target_features.tsv.gz"
    path result_dirs, stageAs: 'partitions/*'
    val expected_strategies
    path merge_script, stageAs: 'bin/merge_alignment_results.py'
    path bin_sources, stageAs: 'bin/*'

    output:
    path "manifest.json", emit: manifest
    path "evidence", emit: evidence
    path "failures.tsv.gz", emit: failures

    script:
    """
    python3 -m bin.merge_alignment_results \\
        --alignment-tasks "${alignment_tasks}" \\
        --source-genes "${source_genes}" \\
        --source-target-features "${source_target_features}" \\
        --result-root partitions \\
        --expected-strategies "${expected_strategies}" \\
        --outdir .
    """
}
