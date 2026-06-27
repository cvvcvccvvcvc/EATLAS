process VALIDATE_IDS {
    tag "ids"

    input:
    path ids_file
    path normalize_script

    output:
    path "input.ids.tsv", emit: ids_tsv
    path "chunks.tsv", emit: chunks_tsv
    path "chunks/*.ids.txt", emit: chunk_files

    script:
    """
    python3 "${normalize_script}" \\
        --ids-file "${ids_file}" \\
        --chunk-size "${params.chunk_size}" \\
        --outdir .
    """
}
