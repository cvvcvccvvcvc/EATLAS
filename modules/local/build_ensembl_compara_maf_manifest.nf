process BUILD_ENSEMBL_COMPARA_MAF_MANIFEST {
    tag "ensembl_compara_maf_manifest"

    input:
    path genes
    path manifest_script, stageAs: 'bin/build_ensembl_compara_maf_manifest.py'
    path bin_sources, stageAs: 'bin/*'

    output:
    path "ensembl_compara_maf_manifest.tsv.gz", emit: maf_manifest
    path "ensembl_compara_maf_manifest.json", emit: manifest
    path "ensembl_compara_maf_manifest_failures.tsv.gz", emit: failures

    script:
    """
    python3 -m bin.build_ensembl_compara_maf_manifest \\
        --genes-tsv "${genes}" \\
        --outdir .
    """
}
