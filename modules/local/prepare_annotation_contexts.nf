process PREPARE_ANNOTATION_CONTEXTS {
    tag "annotation_contexts"

    input:
    path alignment_evidence
    path genes_tsv
    path target_sequences_dir
    path prepare_script

    output:
    path "contexts/*", emit: context_dirs

    script:
    """
    python3 "${prepare_script}" \
        --alignment-evidence-dir "${alignment_evidence}" \
        --genes-tsv "${genes_tsv}" \
        --target-sequences-dir "${target_sequences_dir}" \
        --outdir contexts
    """
}
