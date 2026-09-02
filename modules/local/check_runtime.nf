process CHECK_RUNTIME {
    tag "runtime"

    input:
    path check_script, stageAs: 'bin/check_runtime.py'
    path bin_package_init, stageAs: 'bin/__init__.py'
    path genomics_package_init, stageAs: 'genomics/__init__.py'
    path vep_package_init, stageAs: 'genomics/vep/__init__.py'
    path vep_runtime, stageAs: 'genomics/vep/local_runtime.py'
    val alignment_strategies
    val vep_backend
    val vep_release
    val vep_executable
    val vep_cache_dir

    output:
    path "runtime_check.json", emit: runtime_check

    script:
    """
    python3 -m bin.check_runtime \\
        --alignment-strategies "${alignment_strategies}" \\
        --vep-backend "${vep_backend}" \\
        --vep-release "${vep_release}" \\
        --vep-executable "${vep_executable}" \\
        --vep-cache-dir "${vep_cache_dir}" \\
        --out-json runtime_check.json
    """
}
