process CHECK_RUNTIME {
    tag "runtime"

    input:
    path check_script
    val stage
    val alignment_strategies

    output:
    path "runtime_check.json", emit: runtime_check

    script:
    """
    PYTHONPATH="${projectDir}/bin:\${PYTHONPATH:-}" python3 "${check_script}" \\
        --stage "${stage}" \\
        --alignment-strategies "${alignment_strategies}" \\
        --datasets-bin "${params.datasets_bin}" \\
        --minimap2-bin "${params.minimap2_bin}" \\
        --nucmer-bin "${params.nucmer_bin}" \\
        --show-coords-bin "${params.show_coords_bin}" \\
        --show-snps-bin "${params.show_snps_bin}" \\
        --bwa-bin "${params.bwa_bin}" \\
        --samtools-bin "${params.samtools_bin}" \\
        --varscan-jar "${params.varscan_jar}" \\
        --out-json runtime_check.json
    """
}
