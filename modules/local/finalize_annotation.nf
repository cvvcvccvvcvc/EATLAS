process FINALIZE_ANNOTATION {
    tag "annotation_partitions"

    input:
    path partition_dirs, stageAs: 'partitions/*'
    path vep_shard_dirs, stageAs: 'vep_partitions/*'
    path clinvar_vcf, stageAs: 'references/clinvar.vcf.gz'
    path clinvar_vcf_tbi, stageAs: 'references/clinvar.vcf.gz.tbi'
    path finalize_script

    output:
    path "variant_annotations", emit: variant_annotations
    path "event_variant_map", emit: event_variant_map
    path "manifest.json", emit: manifest
    path "failures.tsv.gz", emit: failures

    script:
    """
    python3 "${finalize_script}" \\
        --partition-root partitions \\
        --vep-root vep_partitions \\
        --clinvar-vcf references/clinvar.vcf.gz \\
        --clinvar-tbi references/clinvar.vcf.gz.tbi \\
        --outdir .
    """
}
