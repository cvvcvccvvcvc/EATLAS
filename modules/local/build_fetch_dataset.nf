process BUILD_FETCH_DATASET {
    tag "fetch_dataset"

    input:
    path ids_tsv
    path chunks_tsv
    path chunk_dirs
    path target_annotation_gff3
    path build_fetch_dataset_script

    output:
    path "manifest.json", emit: manifest
    path "input.ids.tsv", emit: input_ids
    path "chunks.tsv", emit: chunks
    path "genes.tsv.gz", emit: genes
    path "target_features.tsv.gz", emit: target_features
    path "orthologs.selected.tsv.gz", emit: orthologs_selected
    path "orthologs.candidates.tsv.gz", emit: orthologs_candidates
    path "failures.tsv.gz", emit: failures
    path "sequences", optional: true, emit: sequences

    script:
    def chunkArgs = chunk_dirs.collect { "--chunk-dir \"${it}\"" }.join(' ')
    """
    python3 "${build_fetch_dataset_script}" \\
        --ids-tsv "${ids_tsv}" \\
        --chunks-tsv "${chunks_tsv}" \\
        --outdir . \\
        --target-assembly-accession "${params.target_assembly_accession}" \\
        --target-assembly-name "${params.target_assembly_name}" \\
        --target-tax-id "${params.target_tax_id}" \\
        --target-annotation-gff3 "${target_annotation_gff3}" \\
        ${chunkArgs}
    """
}
