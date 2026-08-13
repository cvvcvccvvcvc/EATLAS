process ANNOTATE_EVENTS {
    tag "annotate"

    input:
    path alignment_manifest, stageAs: 'alignment_manifest.json'
    path events_tsv
    path event_ortholog_support_tsv, stageAs: 'support/*'
    path snv_site_depth_tsv
    path snv_taxonomic_depth_tsv
    path snv_alt_taxonomic_support_tsv
    path genes_tsv
    path target_features
    path target_sequences_dir
    path annotate_script
    path annotation_helpers
    path genomics_sources, stageAs: 'genomics/*'
    path clinvar_vcf
    path clinvar_vcf_tbi
    val gnomad_cache_dir

    output:
    path "variant_annotations.tsv.gz", emit: variant_annotations
    path "variant_strategy_support.tsv.gz", emit: variant_strategy_support
    path "variant_ortholog_support", emit: variant_ortholog_support
    path "ortholog_evidence_summary.tsv.gz", emit: ortholog_evidence_summary
    path "manifest.json", emit: manifest
    path "failures.tsv.gz", emit: failures

    script:
    def duckdbMemoryGb = Math.max(2, (task.memory.toGiga() * 0.25) as int)
    """
    GAPH_GNOMAD_CACHE_DIR="${gnomad_cache_dir}" \\
    GAPH_ANNOTATION_DUCKDB_MEMORY_LIMIT="${duckdbMemoryGb}GB" \\
    GAPH_ANNOTATION_DUCKDB_THREADS="${task.cpus}" \\
    python3 "${annotate_script}" \\
        --alignment-manifest alignment_manifest.json \\
        --events-tsv "${events_tsv}" \\
        --event-ortholog-support-tsv "${event_ortholog_support_tsv}" \\
        --snv-site-depth-tsv "${snv_site_depth_tsv}" \\
        --snv-taxonomic-depth-tsv "${snv_taxonomic_depth_tsv}" \\
        --snv-alt-taxonomic-support-tsv "${snv_alt_taxonomic_support_tsv}" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir "${target_sequences_dir}" \\
        --target-features "${target_features}" \\
        --outdir . \\
        --clinvar-vcf "${clinvar_vcf}"
    """
}
