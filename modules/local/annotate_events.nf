process ANNOTATE_EVENTS {
    tag "annotate"

    input:
    path events_tsv
    path genes_tsv
    path sequences_dir
    path annotate_script
    path clinvar_vcf
    path clinvar_vcf_tbi

    output:
    path "variant_annotations.tsv.gz", emit: variant_annotations
    path "manifest.json", emit: manifest
    path "failures.tsv.gz", emit: failures

    script:
    """
    python3 "${annotate_script}" \\
        --events-tsv "${events_tsv}" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir "${sequences_dir}/targets" \\
        --outdir . \\
        --clinvar-vcf "${clinvar_vcf}"
    """
}
