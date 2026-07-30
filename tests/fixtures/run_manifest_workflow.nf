nextflow.enable.dsl = 2

import RunManifest

params.outdir = null
params.source_dir = null
params.schema_path = null
params.api_token = null
params.endpoint = null
params.fail = false
params.hold_seconds = 0
params.stage = 'test'

process TERMINAL_TASK {
    input:
    val fail
    val hold_seconds

    script:
    """
    sleep ${hold_seconds}
    if [ "${fail}" = "true" ]; then
        exit 1
    fi
    """
}

workflow {
    run_manifest_path = file("${params.outdir}/run_manifest.json")
    RunManifest.start(
        run_manifest_path,
        file(params.source_dir),
        workflow,
        params,
        file(params.schema_path)
    )
    TERMINAL_TASK(params.fail, params.hold_seconds)
}

workflow.onComplete {
    RunManifest.finish(file("${params.outdir}/run_manifest.json"), workflow)
}
