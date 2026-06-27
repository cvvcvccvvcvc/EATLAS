import sys

with open("main.nf", "r") as f:
    content = f.read()

ANNOTATE_PROCESS = """
process ANNOTATE_EVENTS {
    tag "annotate"

    input:
    path events_tsv
    path annotate_script
    path clinvar_vcf

    output:
    path "alignment_events_annotated.tsv.gz", emit: annotated_events

    script:
    def clinvarArg = clinvar_vcf.name != 'NO_CLINVAR' ? "--clinvar-vcf \"${clinvar_vcf}\"" : ""
    """
    python3 "${annotate_script}" \\
        --events-tsv "${events_tsv}" \\
        --outdir . \\
        ${clinvarArg}
    """
}

workflow ANNOTATION_STAGE {
    take:
    events_tsv
    clinvar_vcf

    main:
    annotate_script = file("${projectDir}/bin/annotate_events.py")
    ANNOTATE_EVENTS(events_tsv, annotate_script, clinvar_vcf)

    emit:
    annotated_events = ANNOTATE_EVENTS.out.annotated_events
}
"""

# Insert before ALIGNMENT_STAGE_FROM_DIR
idx = content.find("workflow ALIGNMENT_STAGE_FROM_DIR")
content = content[:idx] + ANNOTATE_PROCESS + "\n" + content[idx:]

# Modify main workflow
old_main = """workflow {
    if (params.stage == 'all') {
        FETCH_STAGE(file(params.ids_file))
        ALIGNMENT_STAGE(
            FETCH_STAGE.out.genes,
            FETCH_STAGE.out.orthologs_selected,
            FETCH_STAGE.out.sequences
        )
    } else if (params.stage == 'fetch') {
        FETCH_STAGE(file(params.ids_file))
    } else if (params.stage == 'align') {
        ALIGNMENT_STAGE_FROM_DIR()
    }
}"""

new_main = """workflow {
    clinvar_vcf = params.clinvar_vcf ? file(params.clinvar_vcf) : file('NO_CLINVAR')

    if (params.stage == 'all') {
        FETCH_STAGE(file(params.ids_file))
        ALIGNMENT_STAGE(
            FETCH_STAGE.out.genes,
            FETCH_STAGE.out.orthologs_selected,
            FETCH_STAGE.out.sequences
        )
        ANNOTATION_STAGE(ALIGNMENT_STAGE.out.events, clinvar_vcf)
    } else if (params.stage == 'fetch') {
        FETCH_STAGE(file(params.ids_file))
    } else if (params.stage == 'align') {
        ALIGNMENT_STAGE_FROM_DIR()
        ANNOTATION_STAGE(ALIGNMENT_STAGE_FROM_DIR.out.events, clinvar_vcf)
    } else if (params.stage == 'annotate') {
        ANNOTATION_STAGE(file(params.events_tsv), clinvar_vcf)
    }
}"""

content = content.replace(old_main, new_main)

# We also need ALIGNMENT_STAGE_FROM_DIR to emit events
old_align_dir = """workflow ALIGNMENT_STAGE_FROM_DIR {
    main:
    fetch_dir = file(params.fetch_dir)
    ALIGNMENT_STAGE(
        Channel.value(file("${fetch_dir}/genes.tsv.gz")),
        Channel.value(file("${fetch_dir}/orthologs.selected.tsv.gz")),
        Channel.value(file("${fetch_dir}/sequences"))
    )
}"""
new_align_dir = """workflow ALIGNMENT_STAGE_FROM_DIR {
    main:
    fetch_dir = file(params.fetch_dir)
    ALIGNMENT_STAGE(
        Channel.value(file("${fetch_dir}/genes.tsv.gz")),
        Channel.value(file("${fetch_dir}/orthologs.selected.tsv.gz")),
        Channel.value(file("${fetch_dir}/sequences"))
    )
    emit:
    events = ALIGNMENT_STAGE.out.events
}"""
content = content.replace(old_align_dir, new_align_dir)

with open("main.nf", "w") as f:
    f.write(content)
