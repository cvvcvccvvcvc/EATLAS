process ANNOTATE_EVENTS_PARTITION {
    tag { meta.partition_id }

    input:
    tuple val(meta), path(alignment_partition), path(genes_tsv), path(target_fastas, stageAs: 'targets/*')
    path annotate_script
    path clinvar_vcf
    path clinvar_vcf_tbi

    output:
    tuple val(meta), path("annotation_${meta.partition_id}"), emit: partition_dirs

    script:
    def resultDir = "annotation_${meta.partition_id}"
    """
    python3 "${annotate_script}" \\
        --events-tsv "${alignment_partition}/alignment_events.tsv.gz" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir targets \\
        --outdir "${resultDir}" \\
        --partition-id "${meta.partition_id}" \\
        --clinvar-vcf "${clinvar_vcf}"
    """
}
