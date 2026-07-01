process ANNOTATE_EVENTS {
    tag "annotate"

    input:
    path events_tsv
    path genes_tsv
    path sequences_dir
    path annotate_script

    output:
    path "alignment_events_annotated.tsv.gz", emit: annotated_events

    script:
    """
    python3 "${annotate_script}" \\
        --events-tsv "${events_tsv}" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir "${sequences_dir}/targets" \\
        --outdir .
    """
}

process ANNOTATE_EVENTS_WITH_CLINVAR {
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
    """
    python3 "${annotate_script}" \\
        --events-tsv "${events_tsv}" \\
        --genes-tsv "${genes_tsv}" \\
        --target-sequences-dir "${sequences_dir}/targets" \\
        --outdir . \\
        --clinvar-vcf "${clinvar_vcf}"
    """
}
