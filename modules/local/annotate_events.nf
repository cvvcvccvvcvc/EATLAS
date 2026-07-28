process ANNOTATE_EVENTS {
    tag "annotate"

    input:
    path events_tsv
    path segments_tsv
    path genes_tsv
    path sequences_dir
    path annotate_script
    path clinvar_vcf
    path clinvar_vcf_tbi
    val gnomad_cache_dir

    output:
    path "variant_annotations.tsv.gz", emit: variant_annotations
    path "variant_strategy_support.tsv.gz", emit: variant_strategy_support
    path "manifest.json", emit: manifest
    path "failures.tsv.gz", emit: failures

    script:
    """
    GAPH_GNOMAD_CACHE_DIR="${gnomad_cache_dir}" \\
    python3 "${annotate_script}" \\
        --events-tsv "${events_tsv}" \\
        --segments-tsv "${segments_tsv}" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir "${sequences_dir}/targets" \\
        --outdir . \\
        --clinvar-vcf "${clinvar_vcf}"
    """
}
