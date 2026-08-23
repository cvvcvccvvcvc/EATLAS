process ANNOTATE_EVENTS_PARTITION {
    tag { "${meta.partition_id} events=${meta.annotation_event_count} base_memory=${meta.annotation_memory_gb}GB" }

    input:
    tuple val(meta), path(alignment_partition), path(genes_tsv), path(target_fastas, stageAs: 'targets/*')
    path annotate_script, stageAs: 'bin/annotate_events.py'
    path bin_package_init, stageAs: 'bin/__init__.py'
    path genomics_sources, stageAs: 'genomics/*'
    path clinvar_vcf
    path clinvar_vcf_tbi
    val gnomad_cache_dir

    output:
    tuple val(meta), path("annotation_${meta.partition_id}"), emit: partition_dirs

    script:
    def resultDir = "annotation_${meta.partition_id}"
    """
    echo "INFO: Annotation resources partition=${meta.partition_id} events=${meta.annotation_event_count} attempt=${task.attempt} memory=${task.memory}" >&2
    GAPH_GNOMAD_CACHE_DIR="${gnomad_cache_dir}" \\
    python3 -m bin.annotate_events \\
        --alignment-manifest "${alignment_partition}/manifest.json" \\
        --events-tsv "${alignment_partition}/alignment_events.tsv.gz" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir targets \\
        --outdir "${resultDir}" \\
        --partition-id "${meta.partition_id}" \\
        --clinvar-vcf "${clinvar_vcf}"
    """
}
