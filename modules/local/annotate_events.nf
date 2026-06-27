process ANNOTATE_EVENTS {
    tag "annotate"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path events_tsv
    path annotate_script
    path clinvar_vcf
    path clinvar_vcf_tbi

    output:
    path "alignment_events_annotated.tsv.gz", emit: annotated_events

    script:
    def clinvarArg = clinvar_vcf.name != 'NO_CLINVAR' ? "--clinvar-vcf \"${clinvar_vcf}\"" : ""
    """
    python3 "${annotate_script}" \\
        --events-tsv "${events_tsv}" \\
        --outdir . \\
        ${clinvarArg}
    """
}
