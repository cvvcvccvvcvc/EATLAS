process FINALIZE_ANNOTATION {
    tag "annotation_partitions"

    input:
    path partition_dirs, stageAs: 'partitions/*'
    path vep_shard_dirs, stageAs: 'vep_partitions/*'
    path clinvar_vcf, stageAs: 'references/clinvar.vcf.gz'
    path clinvar_vcf_tbi, stageAs: 'references/clinvar.vcf.gz.tbi'
    path finalize_script, stageAs: 'bin/finalize_annotation_partitions.py'
    path bin_package_init, stageAs: 'bin/__init__.py'
    path genomics_package_init, stageAs: 'genomics/__init__.py'
    path gnomad_source, stageAs: 'genomics/gnomad.py'
    path variants_source, stageAs: 'genomics/variants.py'

    output:
    path "variant_annotations", emit: variant_annotations
    path "event_variant_map", emit: event_variant_map
    path "manifest.json", emit: manifest
    path "failures.tsv.gz", emit: failures

    script:
    """
    python3 -m bin.finalize_annotation_partitions \\
        --partition-root partitions \\
        --vep-root vep_partitions \\
        --clinvar-vcf references/clinvar.vcf.gz \\
        --clinvar-tbi references/clinvar.vcf.gz.tbi \\
        --outdir .
    """
}
