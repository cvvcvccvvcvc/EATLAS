process CHECK_RUNTIME {
    tag "runtime"

    input:
    path check_script, stageAs: 'bin/check_runtime.py'
    path bin_package_init, stageAs: 'bin/__init__.py'
    val alignment_strategies

    output:
    path "runtime_check.json", emit: runtime_check

    script:
    """
    python3 -m bin.check_runtime \\
        --alignment-strategies "${alignment_strategies}" \\
        --out-json runtime_check.json
    """
}
