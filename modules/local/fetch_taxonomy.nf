process FETCH_TAXONOMY {
    tag "taxonomy"

    input:
    path orthologs_tsv
    path taxonomy_script

    output:
    path "taxonomy.tsv.gz", emit: taxonomy
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
        --outdir .
    """
}
