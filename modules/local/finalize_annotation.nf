process FINALIZE_ANNOTATION {
    tag "annotation_partitions"

    input:
    path partition_dirs, stageAs: 'partitions/*'
    path finalize_script

    output:
    path "variant_annotations.tsv.gz", emit: variant_annotations
    path "manifest.json", emit: manifest
    path "failures.tsv.gz", emit: failures

    script:
    """
    python3 "${finalize_script}" \\
        --partition-root partitions \\
        --outdir .
    """
}
