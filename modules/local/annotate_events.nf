process ANNOTATE_EVENTS {
    tag "annotate"

    input:
    tuple val(has_event_support), path(annotation_inputs, stageAs: 'annotation_inputs/*')
    path segments_tsv
    path genes_tsv
    path sequences_dir
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
    path "manifest.json", emit: manifest
    path "failures.tsv.gz", emit: failures

    script:
    def inputs = annotation_inputs instanceof List ? annotation_inputs : [annotation_inputs]
    if (inputs.size() != (has_event_support ? 2 : 1)) {
        error "Annotation input bundle has ${inputs.size()} file(s), expected ${has_event_support ? 2 : 1}"
    }
    def eventsTsv = inputs[0]
    def eventSupportArg = has_event_support ? "--event-ortholog-support-tsv \"${inputs[1]}\"" : ""
    def duckdbMemoryGb = Math.max(2, (task.memory.toGiga() * 0.25) as int)
    """
    GAPH_GNOMAD_CACHE_DIR="${gnomad_cache_dir}" \\
    GAPH_ANNOTATION_DUCKDB_MEMORY_LIMIT="${duckdbMemoryGb}GB" \\
    GAPH_ANNOTATION_DUCKDB_THREADS="${task.cpus}" \\
    python3 "${annotate_script}" ${eventSupportArg} \\
        --events-tsv "${eventsTsv}" \\
        --segments-tsv "${segments_tsv}" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir "${sequences_dir}/targets" \\
        --outdir . \\
        --clinvar-vcf "${clinvar_vcf}"
    """
}
