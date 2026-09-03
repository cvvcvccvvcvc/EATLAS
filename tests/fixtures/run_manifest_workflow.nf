nextflow.enable.dsl = 2

import RunManifest

include { BUILD_EVIDENCE_INVENTORY } from '../../modules/local/build_evidence_inventory.nf'

params.outdir = null
params.source_dir = null
params.schema_path = null
params.api_token = null
params.endpoint = null
params.fail = false
params.hold_seconds = 0
process TERMINAL_TASK {
    publishDir "${params.outdir}", mode: 'move', pattern: '{fetch,alignment,annotation}'

    input:
    val fail
    val hold_seconds

    output:
    path 'fetch'
    path 'alignment'
    path 'annotation'
    path 'fetch.inventory.json', emit: fetch_inventory
    path 'alignment.inventory.json', emit: alignment_inventory
    path 'annotation.inventory.json', emit: annotation_inventory

    script:
    """
    sleep ${hold_seconds}
    if [ "${fail}" = "true" ]; then
        exit 1
    fi
    mkdir -p fetch alignment annotation
    printf '%s\n' '{"stage":"fetch"}' > fetch/manifest.json
    printf '%s\n' '{"stage":"alignment"}' > alignment/manifest.json
    printf '%s\n' '{"stage":"annotation"}' > annotation/manifest.json
    task_dir=\$PWD
    cd "${params.source_dir}"
    python3 -m provenance.evidence_inventory create \
        --scope fetch --output "\${task_dir}/fetch.inventory.json" --input "fetch=\${task_dir}/fetch"
    python3 -m provenance.evidence_inventory create \
        --scope alignment --output "\${task_dir}/alignment.inventory.json" --input "alignment=\${task_dir}/alignment"
    python3 -m provenance.evidence_inventory create \
        --scope annotation --output "\${task_dir}/annotation.inventory.json" --input "annotation=\${task_dir}/annotation"
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
    provenance_sources = [
        file("${params.source_dir}/provenance/__init__.py"),
        file("${params.source_dir}/provenance/evidence_inventory.py"),
    ]
    BUILD_EVIDENCE_INVENTORY(
        TERMINAL_TASK.out.fetch_inventory,
        TERMINAL_TASK.out.alignment_inventory,
        TERMINAL_TASK.out.annotation_inventory,
        provenance_sources
    )
}

workflow.onComplete {
    RunManifest.finish(file("${params.outdir}/run_manifest.json"), workflow)
}
