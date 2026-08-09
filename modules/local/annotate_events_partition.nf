process ANNOTATE_EVENTS_PARTITION {
    tag { "${meta.partition_id} support_rows=${meta.annotation_event_ortholog_support_count} base_memory=${meta.annotation_memory_gb}GB" }

    input:
    tuple val(meta), path(alignment_partition), path(genes_tsv), path(target_fastas, stageAs: 'targets/*'), path(target_features, stageAs: 'target_features/*')
    path annotate_script
    path annotation_helpers
    path genomics_sources, stageAs: 'genomics/*'
    path clinvar_vcf
    path clinvar_vcf_tbi
    val gnomad_cache_dir

    output:
    tuple val(meta), path("annotation_${meta.partition_id}"), emit: partition_dirs

    script:
    def resultDir = "annotation_${meta.partition_id}"
    """
    echo "INFO: Annotation resources partition=${meta.partition_id} support_rows=${meta.annotation_event_ortholog_support_count} attempt=${task.attempt} memory=${task.memory}" >&2
    GAPH_GNOMAD_CACHE_DIR="${gnomad_cache_dir}" \\
    python3 "${annotate_script}" \\
        --events-tsv "${alignment_partition}/alignment_events.tsv.gz" \\
        --event-ortholog-support-tsv "${alignment_partition}/event_ortholog_support.tsv.gz" \\
        --snv-site-depth-tsv "${alignment_partition}/snv_site_depth.tsv.gz" \\
        --snv-taxonomic-depth-tsv "${alignment_partition}/snv_taxonomic_depth.tsv.gz" \\
        --snv-alt-taxonomic-support-tsv "${alignment_partition}/snv_alt_taxonomic_support.tsv.gz" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir targets \\
        --target-features-dir target_features \\
        --outdir "${resultDir}" \\
        --partition-id "${meta.partition_id}" \\
        --clinvar-vcf "${clinvar_vcf}"
    """
}
