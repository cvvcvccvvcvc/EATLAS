process FETCH_TAXONOMY {
    tag "taxonomy"

    input:
    path orthologs_tsv
    path taxonomy_script, stageAs: 'bin/fetch_taxonomy.py'
    path bin_package_init, stageAs: 'bin/__init__.py'
    path taxonomy_sources, stageAs: 'genomics/*'

    output:
    path "taxonomy.tsv.gz", emit: taxonomy
    path "taxonomy_failures.tsv.gz", emit: taxonomy_failures

    script:
    """
    if [[ -f "${projectDir}/.env" ]]; then
        set -a
        source "${projectDir}/.env"
        set +a
    fi

    python3 -m bin.fetch_taxonomy \\
        --orthologs-tsv "${orthologs_tsv}" \\
        --outdir .
    """
}
