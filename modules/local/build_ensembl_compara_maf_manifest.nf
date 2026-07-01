process BUILD_ENSEMBL_COMPARA_MAF_MANIFEST {
    tag "ensembl_compara_maf_manifest"

    input:
    path genes
    path manifest_script

    output:
    path "ensembl_compara_maf_manifest.tsv.gz", emit: maf_manifest
    path "ensembl_compara_maf_manifest.json", emit: manifest
    path "ensembl_compara_maf_manifest_failures.tsv.gz", emit: failures

    script:
    """
    python3 "${manifest_script}" \\
        --genes-tsv "${genes}" \\
        --outdir . \\
        --release "${params.ensembl_compara_maf_release}" \\
        --species-set "${params.ensembl_compara_maf_species_set}" \\
        --base-url "${params.ensembl_compara_maf_base_url}"
    """
}
