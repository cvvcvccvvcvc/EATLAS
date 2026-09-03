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
    path provenance_sources, stageAs: 'provenance/*'

    output:
    path "manifest.json", emit: manifest
    path "evidence", emit: evidence
    path "failures.tsv.gz", emit: failures
    path "alignment.inventory.json", emit: inventory

    script:
    """
    python3 -m bin.merge_alignment_results \\
        --alignment-tasks "${alignment_tasks}" \\
        --source-genes "${source_genes}" \\
        --source-target-features "${source_target_features}" \\
        --result-root partitions \\
        --expected-strategies "${expected_strategies}" \\
        --outdir .
    python3 -m provenance.evidence_inventory create --scope alignment \\
        --output alignment.inventory.json \\
        --input alignment/manifest.json=manifest.json \\
        --input alignment/evidence=evidence \\
        --input alignment/failures.tsv.gz=failures.tsv.gz
    """
}
