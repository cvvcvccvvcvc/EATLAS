process FINALIZE_ANNOTATION {
    tag "annotation_partitions"

    input:
    path partition_dirs, stageAs: 'partitions/*'
    path finalize_script

    output:
    path "variant_annotations.tsv.gz", emit: variant_annotations
    path "variant_strategy_support.tsv.gz", emit: variant_strategy_support
    path "ortholog_evidence_summary.tsv.gz", emit: ortholog_evidence_summary
    path "manifest.json", emit: manifest
    path "failures.tsv.gz", emit: failures

    script:
    """
    python3 "${finalize_script}" \\
        --partition-root partitions \\
        --outdir .
    """
}
