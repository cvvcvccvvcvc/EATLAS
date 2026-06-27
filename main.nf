nextflow.enable.dsl = 2

allowed_stages = ['all', 'fetch', 'align']
if (!allowed_stages.contains(params.stage)) {
    error "Invalid --stage '${params.stage}'. Expected one of: ${allowed_stages.join(', ')}"
}

if (['all', 'fetch'].contains(params.stage) && !params.ids_file) {
    error "Missing required parameter: --ids_file"
}

if (params.stage == 'align' && !params.fetch_dir) {
    error "Missing required parameter for --stage align: --fetch_dir"
}

process VALIDATE_IDS {
    tag "ids"

    input:
    path ids_file
    path normalize_script

    output:
    path "input.ids.tsv", emit: ids_tsv
    path "chunks.tsv", emit: chunks_tsv
    path "chunks/*.ids.txt", emit: chunk_files

    script:
    """
    python3 "${normalize_script}" \\
        --ids-file "${ids_file}" \\
        --chunk-size "${params.chunk_size}" \\
        --outdir .
    """
}

process FETCH_PARSE_CHUNK {
    tag { chunk_file.baseName }

    input:
    path chunk_file
    path fetch_script

    output:
    path "fetch_*", emit: chunk_dirs

    script:
    """
    chunk_name=\$(basename "${chunk_file}" .ids.txt)

    python3 "${fetch_script}" \\
        --ids-file "${chunk_file}" \\
        --outdir "fetch_\${chunk_name}" \\
        --datasets-bin "${params.datasets_bin}" \\
        --target-assembly-accession "${params.target_assembly_accession}" \\
        --target-assembly-name "${params.target_assembly_name}" \\
        --target-tax-id "${params.target_tax_id}"
    """
}

process MERGE_FETCH_RESULTS {
    tag "merge"

    input:
    path ids_tsv
    path chunks_tsv
    path chunk_dirs
    path merge_script

    output:
    path "manifest.json", emit: manifest
    path "input.ids.tsv", emit: input_ids
    path "chunks.tsv", emit: chunks
    path "genes.tsv.gz", emit: genes
    path "orthologs.selected.tsv.gz", emit: orthologs_selected
    path "orthologs.candidates.tsv.gz", emit: orthologs_candidates
    path "failures.tsv.gz", emit: failures
    path "sequences", optional: true, emit: sequences

    script:
    def chunkArgs = chunk_dirs.collect { "--chunk-dir \"${it}\"" }.join(' ')
    """
    python3 "${merge_script}" \\
        --ids-tsv "${ids_tsv}" \\
        --chunks-tsv "${chunks_tsv}" \\
        --outdir . \\
        --target-assembly-accession "${params.target_assembly_accession}" \\
        --target-assembly-name "${params.target_assembly_name}" \\
        --target-tax-id "${params.target_tax_id}" \\
        ${chunkArgs}
    """
}

process FETCH_TAXONOMY_PRESETS {
    tag "taxonomy"

    input:
    path orthologs_tsv
    path taxonomy_script
    path taxonomy_classes

    output:
    path "taxonomy_presets.tsv.gz", emit: taxonomy_presets
    path "taxonomy_failures.tsv.gz", emit: taxonomy_failures

    script:
    """
    python3 "${taxonomy_script}" \\
        --orthologs-tsv "${orthologs_tsv}" \\
        --outdir . \\
        --taxonomy-classes "${taxonomy_classes}"
    """
}

process BUILD_ALIGNMENT_TASKS {
    tag "alignment_tasks"

    input:
    path genes_tsv
    path orthologs_tsv
    path "sequences/*", stageAs: 'sequences/*'
    path taxonomy_presets
    path prepare_script

    output:
    path "alignment_tasks.tsv.gz", emit: alignment_tasks
    path "tasks/task_*", emit: task_dirs

    script:
    """
    targetArgs=\$(find sequences/*/targets -name "*.fa.gz" | sed 's/^/--target-fasta /' | tr '\\n' ' ')
    orthologArgs=\$(find sequences/*/orthologs -name "*.fa.gz" | sed 's/^/--ortholog-fasta /' | tr '\\n' ' ')
    python3 "${prepare_script}" \\
        --genes-tsv "${genes_tsv}" \\
        --orthologs-tsv "${orthologs_tsv}" \\
        --taxonomy-presets "${taxonomy_presets}" \\
        --outdir . \\
        \$targetArgs \\
        \$orthologArgs
    """
}

process ALIGN_MINIMAP2_ASM10 {
    tag { task_dir.baseName }

    input:
    path task_dir
    path minimap2_script

    output:
    path "align_minimap2_asm10", emit: asm10_result_dirs

    script:
    """
    python3 "${minimap2_script}" \\
        --task-dir "${task_dir}" \\
        --outdir "align_minimap2_asm10" \\
        --strategy minimap2_asm10 \\
        --mode fixed \\
        --fixed-preset asm10 \\
        --minimap2-bin "${params.minimap2_bin}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}

process ALIGN_MINIMAP2_ASM20 {
    tag { task_dir.baseName }

    input:
    path task_dir
    path minimap2_script

    output:
    path "align_minimap2_asm20", emit: asm20_result_dirs

    script:
    """
    python3 "${minimap2_script}" \\
        --task-dir "${task_dir}" \\
        --outdir "align_minimap2_asm20" \\
        --strategy minimap2_asm20 \\
        --mode fixed \\
        --fixed-preset asm20 \\
        --minimap2-bin "${params.minimap2_bin}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}

process ALIGN_MINIMAP2_ADAPTIVE {
    tag { task_dir.baseName }

    input:
    path task_dir
    path minimap2_script

    output:
    path "align_minimap2_taxonomy_adaptive", emit: adaptive_result_dirs

    script:
    """
    python3 "${minimap2_script}" \\
        --task-dir "${task_dir}" \\
        --outdir "align_minimap2_taxonomy_adaptive" \\
        --strategy minimap2_taxonomy_adaptive \\
        --mode adaptive \\
        --minimap2-bin "${params.minimap2_bin}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}

process MERGE_ALIGNMENT {
    tag "merge_alignment"

    input:
    path alignment_tasks
    path taxonomy_presets
    path taxonomy_failures
    path(minimap2_asm10_dirs, stageAs: 'minimap2_asm10_dirs/*')
    path(minimap2_asm20_dirs, stageAs: 'minimap2_asm20_dirs/*')
    path(minimap2_adaptive_dirs, stageAs: 'minimap2_adaptive_dirs/*')
    path(nucmer_dirs, stageAs: 'nucmer_dirs/*')
    path(bwa_dirs, stageAs: 'bwa_dirs/*')
    path merge_script

    output:
    path "manifest.json", emit: manifest
    path "alignment_tasks.tsv.gz", emit: alignment_tasks
    path "taxonomy_presets.tsv.gz", emit: taxonomy_presets
    path "taxonomy_failures.tsv.gz", emit: taxonomy_failures
    path "ortholog_alignment_summary.tsv.gz", emit: summaries
    path "alignment_segments.tsv.gz", emit: segments
    path "alignment_events.tsv.gz", emit: events
    path "failures.tsv.gz", emit: failures
    path "native", optional: true, emit: native_outputs

    script:
    def asm10Args = minimap2_asm10_dirs.collect { "--result-dir \"${it}\"" }.join(' ')
    def asm20Args = minimap2_asm20_dirs.collect { "--result-dir \"${it}\"" }.join(' ')
    def adaptiveArgs = minimap2_adaptive_dirs.collect { "--result-dir \"${it}\"" }.join(' ')
    def nucmerArgs = nucmer_dirs.collect { "--result-dir \"${it}\"" }.join(' ')
    def bwaArgs = bwa_dirs.collect { "--result-dir \"${it}\"" }.join(' ')
    """
    python3 "${merge_script}" \\
        --alignment-tasks "${alignment_tasks}" \\
        --taxonomy-presets "${taxonomy_presets}" \\
        --taxonomy-failures "${taxonomy_failures}" \\
        --outdir . \\
        ${asm10Args} \\
        ${asm20Args} \\
        ${adaptiveArgs} \\
        ${nucmerArgs} \\
        ${bwaArgs}
    """
}

process ALIGN_NUCMER_COMPARATOR {
    tag { task_dir.baseName }

    input:
    path task_dir
    path nucmer_script

    output:
    path "align_nucmer", emit: nucmer_result_dirs

    script:
    """
    python3 "${nucmer_script}" \\
        --task-dir "${task_dir}" \\
        --outdir "align_nucmer" \\
        --nucmer-bin "${params.nucmer_bin}" \\
        --show-coords-bin "${params.show_coords_bin}" \\
        --show-snps-bin "${params.show_snps_bin}" \\
        --keep-native "${params.keep_native_alignments}"
    """
}

process ALIGN_BWA_PSEUDOREADS {
    tag { task_dir.baseName }

    input:
    path task_dir
    path bwa_script

    output:
    path "align_bwa", emit: bwa_result_dirs

    script:
    """
    cp -r "${task_dir}" align_bwa
    python3 "${bwa_script}" \\
        --task-dir "align_bwa" \\
        --bwa-bin "${params.bwa_bin}" \\
        --samtools-bin "${params.samtools_bin}" \\
        --varscan-jar "${params.varscan_jar}"
    """
}

workflow FETCH_STAGE {
    take:
    ids

    main:
    normalize_script = file("${projectDir}/bin/normalize_ids.py")
    fetch_script = file("${projectDir}/bin/fetch_parse_chunk.py")
    merge_script = file("${projectDir}/bin/merge_fetch_results.py")

    VALIDATE_IDS(ids, normalize_script)

    chunk_files = VALIDATE_IDS.out.chunk_files.flatten()
    FETCH_PARSE_CHUNK(chunk_files, fetch_script)

    MERGE_FETCH_RESULTS(
        VALIDATE_IDS.out.ids_tsv,
        VALIDATE_IDS.out.chunks_tsv,
        FETCH_PARSE_CHUNK.out.chunk_dirs.collect(),
        merge_script
    )

    emit:
    manifest = MERGE_FETCH_RESULTS.out.manifest
    input_ids = MERGE_FETCH_RESULTS.out.input_ids
    chunks = MERGE_FETCH_RESULTS.out.chunks
    genes = MERGE_FETCH_RESULTS.out.genes
    orthologs_selected = MERGE_FETCH_RESULTS.out.orthologs_selected
    orthologs_candidates = MERGE_FETCH_RESULTS.out.orthologs_candidates
    failures = MERGE_FETCH_RESULTS.out.failures
    sequences = MERGE_FETCH_RESULTS.out.sequences
}

workflow ALIGNMENT_STAGE {
    take:
    genes
    orthologs_selected
    sequences

    main:
    taxonomy_script = file("${projectDir}/bin/fetch_taxonomy_presets.py")
    taxonomy_classes = file("${projectDir}/assets/taxonomy_classes.json.gz")
    prepare_script = file("${projectDir}/bin/prepare_alignment_tasks.py")
    minimap2_script = file("${projectDir}/bin/run_minimap2_alignment.py")
    nucmer_script = file("${projectDir}/bin/run_nucmer_alignment.py")
    bwa_script = file("${projectDir}/bin/run_bwa_pseudoreads.py")
    merge_script = file("${projectDir}/bin/merge_alignment_results.py")

    target_fastas = sequences.map { seq_dir -> file("${seq_dir}/targets/*.fa.gz") }.flatten()
    ortholog_fastas = sequences.map { seq_dir -> file("${seq_dir}/orthologs/*.fa.gz") }.flatten().unique { it.name }

    FETCH_TAXONOMY_PRESETS(orthologs_selected, taxonomy_script, taxonomy_classes)
    BUILD_ALIGNMENT_TASKS(
        genes,
        orthologs_selected,
        sequences.collect(),
        FETCH_TAXONOMY_PRESETS.out.taxonomy_presets,
        prepare_script
    )

    task_dirs = BUILD_ALIGNMENT_TASKS.out.task_dirs.flatten()
    ALIGN_MINIMAP2_ASM10(task_dirs, minimap2_script)
    ALIGN_MINIMAP2_ASM20(task_dirs, minimap2_script)
    ALIGN_MINIMAP2_ADAPTIVE(task_dirs, minimap2_script)
    ALIGN_NUCMER_COMPARATOR(task_dirs, nucmer_script)
    ALIGN_BWA_PSEUDOREADS(task_dirs, bwa_script)

    MERGE_ALIGNMENT(
        BUILD_ALIGNMENT_TASKS.out.alignment_tasks,
        FETCH_TAXONOMY_PRESETS.out.taxonomy_presets,
        FETCH_TAXONOMY_PRESETS.out.taxonomy_failures,
        ALIGN_MINIMAP2_ASM10.out.asm10_result_dirs.collect(),
        ALIGN_MINIMAP2_ASM20.out.asm20_result_dirs.collect(),
        ALIGN_MINIMAP2_ADAPTIVE.out.adaptive_result_dirs.collect(),
        ALIGN_NUCMER_COMPARATOR.out.nucmer_result_dirs.collect(),
        ALIGN_BWA_PSEUDOREADS.out.bwa_result_dirs.collect(),
        merge_script
    )

    emit:
    manifest = MERGE_ALIGNMENT.out.manifest
    tasks = MERGE_ALIGNMENT.out.alignment_tasks
    taxonomy_presets = MERGE_ALIGNMENT.out.taxonomy_presets
    taxonomy_failures = MERGE_ALIGNMENT.out.taxonomy_failures
    summaries = MERGE_ALIGNMENT.out.summaries
    segments = MERGE_ALIGNMENT.out.segments
    events = MERGE_ALIGNMENT.out.events
    failures = MERGE_ALIGNMENT.out.failures
}

workflow ALIGNMENT_STAGE_FROM_DIR {
    main:
    fetch_dir = file(params.fetch_dir)
    ALIGNMENT_STAGE(
        Channel.value(file("${fetch_dir}/genes.tsv.gz")),
        Channel.value(file("${fetch_dir}/orthologs.selected.tsv.gz")),
        Channel.value(file("${fetch_dir}/sequences"))
    )
    emit:
    events = ALIGNMENT_STAGE.out.events
}

process ANNOTATE_EVENTS {
    tag "annotate"
    publishDir "${params.outdir}", mode: 'copy'

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

workflow {
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
}
