process ANNOTATE_VEP_PARTITION {
    tag { "${meta.partition_id}/${meta.shard_id} rows=${meta.vep_row_count}" }

    input:
    tuple val(meta), path(input_tsv)
    path annotate_script, stageAs: 'bin/annotate_vep_partition.py'
    path bin_package_init, stageAs: 'bin/__init__.py'
    path genomics_package_init, stageAs: 'genomics/__init__.py'
    path variants_source, stageAs: 'genomics/variants.py'
    path vep_sources, stageAs: 'genomics/vep/*'
    val vep_backend
    val vep_release
    val vep_executable
    val vep_cache_dir
    val vep_result_cache_dir
    val vep_result_cache_tile_size_bp
    val vep_forks

    output:
    tuple val(meta), path("vep_${meta.partition_id}_${meta.shard_id}"), emit: shard_dirs

    script:
    def resultDir = "vep_${meta.partition_id}_${meta.shard_id}"
    def releaseArg = vep_release ? "--vep-release \"${vep_release}\"" : ""
    def cacheArg = vep_cache_dir ? "--vep-cache-dir \"${vep_cache_dir}\"" : ""
    def resultCacheArg = vep_result_cache_dir ? "--vep-result-cache-dir \"${vep_result_cache_dir}\"" : ""
    """
    python3 -m bin.annotate_vep_partition \\
        --input-tsv "${input_tsv}" \\
        --outdir "${resultDir}" \\
        --partition-id "${meta.partition_id}" \\
        --shard-id "${meta.shard_id}" \\
        --expected-row-count "${meta.vep_row_count}" \\
        --vep-backend "${vep_backend}" \\
        ${releaseArg} \\
        --vep-executable "${vep_executable}" \\
        ${cacheArg} \\
        --vep-forks "${vep_forks}" \\
        --rest-workers "${task.cpus}" \\
        ${resultCacheArg} \\
        --vep-result-cache-tile-size-bp "${vep_result_cache_tile_size_bp}"
    """
}
