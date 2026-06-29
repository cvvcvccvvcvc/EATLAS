nextflow.enable.dsl = 2

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

AVAILABLE_ALIGNMENT_STRATEGIES = [
    'minimap2_asm10',
    'minimap2_asm20',
    'minimap2_taxonomy_adaptive',
    'nucmer',
    'bwa_pseudoreads',
]

def parseAlignmentStrategies(rawValue) {
    def raw = rawValue == null ? 'all' : rawValue.toString().trim()
    def selected = []
    if (!raw || raw == 'all') {
        selected = AVAILABLE_ALIGNMENT_STRATEGIES
    } else {
        selected = raw.split(',')
            .collect { it.trim() }
            .findAll { it }
            .unique()
    }
    def unknown = selected.findAll { !AVAILABLE_ALIGNMENT_STRATEGIES.contains(it) }
    if (unknown) {
        throw new IllegalArgumentException(
            "Unknown alignment_strategies value(s): ${unknown.join(', ')}. Available: ${AVAILABLE_ALIGNMENT_STRATEGIES.join(', ')}"
        )
    }
    if (!selected) {
        throw new IllegalArgumentException("--alignment_strategies must select at least one strategy")
    }
    return selected
}

// Print help message
if (params.help) {
    log.info paramsHelp("gaph_v2")
    exit 0
}

// Validate parameters against nextflow_schema.json
validateParameters()

if (['all', 'fetch'].contains(params.stage) && !params.ids_file) {
    error "Missing required parameter: --ids_file"
}

if (params.stage == 'align' && !params.fetch_dir) {
    error "Missing required parameter for --stage align: --fetch_dir"
}

if (params.stage == 'annotate' && !params.events_tsv) {
    error "Missing required parameter for --stage annotate: --events_tsv"
}

if (params.stage == 'annotate' && !params.fetch_dir) {
    error "Missing required parameter for --stage annotate: --fetch_dir"
}

SELECTED_ALIGNMENT_STRATEGIES = parseAlignmentStrategies(params.alignment_strategies)

include { VALIDATE_IDS } from './modules/local/validate_ids.nf'
include { FETCH_PARSE_CHUNK } from './modules/local/fetch_parse_chunk.nf'
include { MERGE_FETCH_RESULTS } from './modules/local/merge_fetch_results.nf'
include { FETCH_TAXONOMY_PRESETS } from './modules/local/fetch_taxonomy_presets.nf'
include { BUILD_ALIGNMENT_TASKS } from './modules/local/build_alignment_tasks.nf'
include { ALIGN_MINIMAP2_ASM10 } from './modules/local/align_minimap2_asm10.nf'
include { ALIGN_MINIMAP2_ASM20 } from './modules/local/align_minimap2_asm20.nf'
include { ALIGN_MINIMAP2_ADAPTIVE } from './modules/local/align_minimap2_adaptive.nf'
include { MERGE_ALIGNMENT } from './modules/local/merge_alignment.nf'
include { ALIGN_NUCMER_COMPARATOR } from './modules/local/align_nucmer_comparator.nf'
include { ALIGN_BWA_PSEUDOREADS } from './modules/local/align_bwa_pseudoreads.nf'
include { ANNOTATE_EVENTS } from './modules/local/annotate_events.nf'

workflow FETCH_STAGE {
    take:
    ids

    main:
    normalize_script = file("${projectDir}/bin/normalize_ids.py")
    fetch_script = file("${projectDir}/bin/fetch_parse_chunk.py")
    merge_script = file("${projectDir}/bin/merge_fetch_results.py")

    VALIDATE_IDS(ids, normalize_script)

    chunk_files = VALIDATE_IDS.out.chunk_files.flatten().map { file -> tuple([id: file.baseName], file) }
    FETCH_PARSE_CHUNK(chunk_files, fetch_script)

    MERGE_FETCH_RESULTS(
        VALIDATE_IDS.out.ids_tsv,
        VALIDATE_IDS.out.chunks_tsv,
        FETCH_PARSE_CHUNK.out.chunk_dirs.map { meta, dir -> dir }.collect(),
        merge_script
    )

    emit:
    manifest = MERGE_FETCH_RESULTS.out.manifest
    input_ids = MERGE_FETCH_RESULTS.out.input_ids
    chunks = MERGE_FETCH_RESULTS.out.chunks
    genes = MERGE_FETCH_RESULTS.out.genes
    target_features = MERGE_FETCH_RESULTS.out.target_features
    orthologs_selected = MERGE_FETCH_RESULTS.out.orthologs_selected
    orthologs_candidates = MERGE_FETCH_RESULTS.out.orthologs_candidates
    failures = MERGE_FETCH_RESULTS.out.failures
    sequences = MERGE_FETCH_RESULTS.out.sequences
}

workflow ALIGNMENT_STAGE {
    take:
    genes
    target_features
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

    FETCH_TAXONOMY_PRESETS(orthologs_selected, taxonomy_script, taxonomy_classes)
    BUILD_ALIGNMENT_TASKS(
        genes,
        orthologs_selected,
        sequences.collect(),
        FETCH_TAXONOMY_PRESETS.out.taxonomy_presets,
        prepare_script
    )

    task_dirs = BUILD_ALIGNMENT_TASKS.out.task_dirs.flatten().map { dir -> tuple([id: dir.baseName], dir) }
    alignment_result_dirs = Channel.empty()

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('minimap2_asm10')) {
        ALIGN_MINIMAP2_ASM10(task_dirs, minimap2_script)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_MINIMAP2_ASM10.out.asm10_result_dirs.map { meta, dir -> dir })
    }

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('minimap2_asm20')) {
        ALIGN_MINIMAP2_ASM20(task_dirs, minimap2_script)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_MINIMAP2_ASM20.out.asm20_result_dirs.map { meta, dir -> dir })
    }

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('minimap2_taxonomy_adaptive')) {
        ALIGN_MINIMAP2_ADAPTIVE(task_dirs, minimap2_script)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_MINIMAP2_ADAPTIVE.out.adaptive_result_dirs.map { meta, dir -> dir })
    }

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('nucmer')) {
        ALIGN_NUCMER_COMPARATOR(task_dirs, nucmer_script)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_NUCMER_COMPARATOR.out.nucmer_result_dirs.map { meta, dir -> dir })
    }

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('bwa_pseudoreads')) {
        ALIGN_BWA_PSEUDOREADS(task_dirs, bwa_script)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_BWA_PSEUDOREADS.out.bwa_result_dirs.map { meta, dir -> dir })
    }

    MERGE_ALIGNMENT(
        BUILD_ALIGNMENT_TASKS.out.alignment_tasks,
        FETCH_TAXONOMY_PRESETS.out.taxonomy_presets,
        FETCH_TAXONOMY_PRESETS.out.taxonomy_failures,
        target_features,
        alignment_result_dirs.collect(),
        merge_script
    )

    emit:
    manifest = MERGE_ALIGNMENT.out.manifest
    tasks = MERGE_ALIGNMENT.out.alignment_tasks
    taxonomy_presets = MERGE_ALIGNMENT.out.taxonomy_presets
    taxonomy_failures = MERGE_ALIGNMENT.out.taxonomy_failures
    summaries = MERGE_ALIGNMENT.out.summaries
    segments = MERGE_ALIGNMENT.out.segments
    feature_coverage = MERGE_ALIGNMENT.out.feature_coverage
    events = MERGE_ALIGNMENT.out.events
    failures = MERGE_ALIGNMENT.out.failures
}

workflow ALIGNMENT_STAGE_FROM_DIR {
    main:
    fetch_dir = file(params.fetch_dir)
    genes = Channel.value(file("${fetch_dir}/genes.tsv.gz"))
    target_features = Channel.value(file("${fetch_dir}/target_features.tsv.gz"))
    orthologs_selected = Channel.value(file("${fetch_dir}/orthologs.selected.tsv.gz"))
    sequences = Channel.value(file("${fetch_dir}/sequences"))
    ALIGNMENT_STAGE(
        genes,
        target_features,
        orthologs_selected,
        sequences
    )
    emit:
    events = ALIGNMENT_STAGE.out.events
    genes = genes
    sequences = sequences
}

workflow ANNOTATION_STAGE {
    take:
    events_tsv
    genes_tsv
    sequences_dir
    clinvar_vcf
    clinvar_vcf_tbi

    main:
    annotate_script = file("${projectDir}/bin/annotate_events.py")
    ANNOTATE_EVENTS(events_tsv, genes_tsv, sequences_dir, annotate_script, clinvar_vcf, clinvar_vcf_tbi)

    emit:
    annotated_events = ANNOTATE_EVENTS.out.annotated_events
}

workflow {
    clinvar_vcf = params.clinvar_vcf ? file(params.clinvar_vcf) : file('NO_CLINVAR')
    clinvar_vcf_tbi = params.clinvar_vcf ? file("${params.clinvar_vcf}.tbi") : file('NO_CLINVAR_TBI')

    if (params.stage == 'all') {
        FETCH_STAGE(file(params.ids_file))
        ALIGNMENT_STAGE(
            FETCH_STAGE.out.genes,
            FETCH_STAGE.out.target_features,
            FETCH_STAGE.out.orthologs_selected,
            FETCH_STAGE.out.sequences
        )
        ANNOTATION_STAGE(ALIGNMENT_STAGE.out.events, FETCH_STAGE.out.genes, FETCH_STAGE.out.sequences, clinvar_vcf, clinvar_vcf_tbi)
    } else if (params.stage == 'fetch') {
        FETCH_STAGE(file(params.ids_file))
    } else if (params.stage == 'align') {
        ALIGNMENT_STAGE_FROM_DIR()
        ANNOTATION_STAGE(
            ALIGNMENT_STAGE_FROM_DIR.out.events,
            ALIGNMENT_STAGE_FROM_DIR.out.genes,
            ALIGNMENT_STAGE_FROM_DIR.out.sequences,
            clinvar_vcf,
            clinvar_vcf_tbi
        )
    } else if (params.stage == 'annotate') {
        fetch_dir = file(params.fetch_dir)
        ANNOTATION_STAGE(
            file(params.events_tsv),
            file("${fetch_dir}/genes.tsv.gz"),
            file("${fetch_dir}/sequences"),
            clinvar_vcf,
            clinvar_vcf_tbi
        )
    }
}
