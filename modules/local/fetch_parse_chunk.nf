process FETCH_PARSE_CHUNK {
    tag { chunk_file.baseName }

    input:
    path chunk_file
    path fetch_script

    output:
    path "fetch_*", emit: chunk_dirs

    script:
    """
    chunk_name=\$(basename "${chunk_file}" .ids.txt)

    python3 "${fetch_script}" \\
        --ids-file "${chunk_file}" \\
        --outdir "fetch_\${chunk_name}" \\
        --datasets-bin "${params.datasets_bin}" \\
        --target-assembly-accession "${params.target_assembly_accession}" \\
        --target-assembly-name "${params.target_assembly_name}" \\
        --target-tax-id "${params.target_tax_id}"
    """
}
