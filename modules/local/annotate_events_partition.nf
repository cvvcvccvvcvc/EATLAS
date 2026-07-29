process ANNOTATE_EVENTS_PARTITION {
    tag { meta.partition_id }

    input:
    tuple val(meta), path(alignment_partition), path(genes_tsv), path(target_fastas, stageAs: 'targets/*')
    path annotate_script
    path clinvar_vcf
    path clinvar_vcf_tbi
    val gnomad_cache_dir

    output:
    tuple val(meta), path("annotation_${meta.partition_id}"), emit: partition_dirs

    script:
    def resultDir = "annotation_${meta.partition_id}"
    """
    GAPH_GNOMAD_CACHE_DIR="${gnomad_cache_dir}" \\
    python3 "${annotate_script}" \\
        --events-tsv "${alignment_partition}/alignment_events.tsv.gz" \\
        --snv-site-depth-tsv "${alignment_partition}/snv_site_depth.tsv.gz" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir targets \\
        --outdir "${resultDir}" \\
        --partition-id "${meta.partition_id}" \\
        --clinvar-vcf "${clinvar_vcf}"
    """
}
