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
    path "alignment_events_annotated.tsv.gz", emit: annotated_events

    script:
    def clinvarArg = clinvar_vcf.name != 'no_clinvar.vcf' ? "--clinvar-vcf \"${clinvar_vcf}\"" : ""
    """
    python3 "${annotate_script}" \\
        --events-tsv "${events_tsv}" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir "${sequences_dir}/targets" \\
        --outdir . \\
        ${clinvarArg}
    """
}
