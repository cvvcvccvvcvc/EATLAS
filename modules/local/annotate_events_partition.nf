process ANNOTATE_EVENTS_PARTITION {
    tag { "${meta.partition_id} support_rows=${meta.annotation_event_ortholog_support_count} base_memory=${meta.annotation_memory_gb}GB" }

    input:
    tuple val(meta), path(alignment_partition), path(genes_tsv), path(target_fastas, stageAs: 'targets/*')
    path annotate_script
    path genomics_sources, stageAs: 'genomics/*'
    path clinvar_vcf
    path clinvar_vcf_tbi
    val gnomad_cache_dir

    output:
    tuple val(meta), path("annotation_${meta.partition_id}"), emit: partition_dirs

    script:
    def resultDir = "annotation_${meta.partition_id}"
    def duckdbMemoryGb = Math.max(2, (task.memory.toGiga() * 0.25) as int)
    """
    echo "INFO: Annotation resources partition=${meta.partition_id} support_rows=${meta.annotation_event_ortholog_support_count} attempt=${task.attempt} memory=${task.memory}" >&2
    GAPH_GNOMAD_CACHE_DIR="${gnomad_cache_dir}" \\
    GAPH_ANNOTATION_DUCKDB_MEMORY_LIMIT="${duckdbMemoryGb}GB" \\
    GAPH_ANNOTATION_DUCKDB_THREADS="${task.cpus}" \\
    python3 "${annotate_script}" \\
        --alignment-manifest "${alignment_partition}/manifest.json" \\
        --events-tsv "${alignment_partition}/alignment_events.tsv.gz" \\
        --event-ortholog-support-tsv "${alignment_partition}/event_ortholog_support.tsv.gz" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir targets \\
        --outdir "${resultDir}" \\
        --partition-id "${meta.partition_id}" \\
        --clinvar-vcf "${clinvar_vcf}"
    """
}
