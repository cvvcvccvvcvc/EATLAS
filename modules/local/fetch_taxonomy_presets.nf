process FETCH_TAXONOMY_PRESETS {
    tag "taxonomy"

    input:
    path orthologs_tsv
    path taxonomy_script
    path taxonomy_classes

    output:
    path "taxonomy_presets.tsv.gz", emit: taxonomy_presets
    path "taxonomy_failures.tsv.gz", emit: taxonomy_failures
    path "taxonomy_summary.tsv.gz", emit: taxonomy_summary

    script:
    """
    if [[ -f "${projectDir}/.env" ]]; then
        set -a
        source "${projectDir}/.env"
        set +a
    fi

    python3 "${taxonomy_script}" \\
        --orthologs-tsv "${orthologs_tsv}" \\
        --outdir . \\
        --taxonomy-classes "${taxonomy_classes}" \\
        --datasets-bin "${params.datasets_bin}"
    """
}
