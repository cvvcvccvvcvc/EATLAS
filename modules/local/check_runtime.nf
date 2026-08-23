process CHECK_RUNTIME {
    tag "runtime"

    input:
    path check_script
    val alignment_strategies

    output:
    path "runtime_check.json", emit: runtime_check

    script:
    """
    PYTHONPATH="${projectDir}/bin:\${PYTHONPATH:-}" python3 "${check_script}" \\
        --alignment-strategies "${alignment_strategies}" \\
        --out-json runtime_check.json
    """
}
