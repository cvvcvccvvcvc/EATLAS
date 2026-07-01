nextflow.enable.dsl = 2

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

AVAILABLE_ALIGNMENT_STRATEGIES = [
    'minimap2_asm10',
    'minimap2_asm20',
    'minimap2_taxonomy_adaptive',
    'nucmer',
    'bwa_pseudoreads',
    'precomputed_ensembl_92_mammals_epo_extended',
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

def resolveClinvarInputs() {
    def selectedVcf = params.clinvar_vcf
    if (!selectedVcf) {
        def assetClinvar = "${projectDir}/assets/reference/clinvar/clinvar.vcf.gz"
        selectedVcf = file(assetClinvar).exists() ? assetClinvar : null
    }
    if (!selectedVcf) {
        return [enabled: false, vcf: null, tbi: null, path: null]
    }
    def selectedTbi = "${selectedVcf}.tbi"
    if (!file(selectedVcf).exists()) {
        error "ClinVar VCF not found: ${selectedVcf}. Pass --clinvar_vcf, set CLINVAR_VCF, or place the file at assets/reference/clinvar/clinvar.vcf.gz"
    }
    if (!file(selectedTbi).exists()) {
        error "ClinVar VCF index not found: ${selectedTbi}. Place the .tbi next to the VCF."
    }
    return [enabled: true, vcf: file(selectedVcf), tbi: file(selectedTbi), path: selectedVcf]
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

if (['all', 'fetch'].contains(params.stage) && !file(params.target_annotation_gff3).exists()) {
    error "Target annotation GFF3 not found for --stage ${params.stage}: ${params.target_annotation_gff3}. Pass --target_annotation_gff3, set GAPH_TARGET_ANNOTATION_GFF3, or place the file at assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz"
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
include { CHECK_RUNTIME } from './modules/local/check_runtime.nf'
include { FETCH_PARSE_CHUNK } from './modules/local/fetch_parse_chunk.nf'
include { BUILD_FETCH_DATASET } from './modules/local/build_fetch_dataset.nf'
include { FETCH_TAXONOMY_PRESETS } from './modules/local/fetch_taxonomy_presets.nf'
include { BUILD_ALIGNMENT_TASKS } from './modules/local/build_alignment_tasks.nf'
include { ALIGN_MINIMAP2_ASM10 } from './modules/local/align_minimap2_asm10.nf'
include { ALIGN_MINIMAP2_ASM20 } from './modules/local/align_minimap2_asm20.nf'
include { ALIGN_MINIMAP2_ADAPTIVE } from './modules/local/align_minimap2_adaptive.nf'
include { MERGE_ALIGNMENT } from './modules/local/merge_alignment.nf'
include { ALIGN_NUCMER_COMPARATOR } from './modules/local/align_nucmer_comparator.nf'
include { ALIGN_BWA_PSEUDOREADS } from './modules/local/align_bwa_pseudoreads.nf'
include { BUILD_ENSEMBL_COMPARA_MAF_MANIFEST } from './modules/local/build_ensembl_compara_maf_manifest.nf'
include { ALIGN_ENSEMBL_COMPARA_MAF } from './modules/local/align_ensembl_compara_maf.nf'
include { ANNOTATE_EVENTS; ANNOTATE_EVENTS_WITH_CLINVAR } from './modules/local/annotate_events.nf'

workflow FETCH_STAGE {
    take:
    ids

    main:
    normalize_script = file("${projectDir}/bin/normalize_ids.py")
    fetch_script = file("${projectDir}/bin/fetch_parse_chunk.py")
    build_fetch_dataset_script = file("${projectDir}/bin/build_fetch_dataset.py")
    target_annotation_gff3 = file(params.target_annotation_gff3)

    VALIDATE_IDS(ids, normalize_script)

    chunk_files = VALIDATE_IDS.out.chunk_files.flatten().map { file -> tuple([id: file.baseName], file) }
    FETCH_PARSE_CHUNK(chunk_files, fetch_script)

    BUILD_FETCH_DATASET(
        VALIDATE_IDS.out.ids_tsv,
        VALIDATE_IDS.out.chunks_tsv,
        FETCH_PARSE_CHUNK.out.chunk_dirs.map { meta, dir -> dir }.collect(),
        target_annotation_gff3,
        build_fetch_dataset_script
    )

    emit:
    manifest = BUILD_FETCH_DATASET.out.manifest
    input_ids = BUILD_FETCH_DATASET.out.input_ids
    chunks = BUILD_FETCH_DATASET.out.chunks
    genes = BUILD_FETCH_DATASET.out.genes
    target_features = BUILD_FETCH_DATASET.out.target_features
    orthologs_selected = BUILD_FETCH_DATASET.out.orthologs_selected
    orthologs_candidates = BUILD_FETCH_DATASET.out.orthologs_candidates
    failures = BUILD_FETCH_DATASET.out.failures
    sequences = BUILD_FETCH_DATASET.out.sequences
}

workflow ALIGNMENT_STAGE {
    take:
    genes
    target_features
    orthologs_selected
    sequences

    main:
    taxonomy_script = file("${projectDir}/bin/fetch_taxonomy_presets.py")
    taxonomy_classes = file("${projectDir}/assets/reference/ncbi/taxonomy/taxonomy_classes.json.gz")
    prepare_script = file("${projectDir}/bin/prepare_alignment_tasks.py")
    minimap2_script = file("${projectDir}/bin/run_minimap2_alignment.py")
    nucmer_script = file("${projectDir}/bin/run_nucmer_alignment.py")
    bwa_script = file("${projectDir}/bin/run_bwa_pseudoreads.py")
    bam_filtering_script = file("${projectDir}/bin/bam_filtering_v1.py")
    ensembl_compara_maf_manifest_script = file("${projectDir}/bin/build_ensembl_compara_maf_manifest.py")
    ensembl_compara_maf_script = file("${projectDir}/bin/run_ensembl_compara_maf_alignment.py")
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
        ALIGN_BWA_PSEUDOREADS(task_dirs, bwa_script, bam_filtering_script)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_BWA_PSEUDOREADS.out.bwa_result_dirs.map { meta, dir -> dir })
    }

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('precomputed_ensembl_92_mammals_epo_extended')) {
        default_maf_manifest = file("${projectDir}/assets/reference/ensembl/compara/release-${params.ensembl_compara_maf_release}/${params.ensembl_compara_maf_species_set}/ensembl_compara_maf_manifest.tsv.gz")
        configured_maf_manifest = params.ensembl_compara_maf_manifest ? file(params.ensembl_compara_maf_manifest) : null
        if (configured_maf_manifest && !configured_maf_manifest.exists()) {
            error "Configured Ensembl Compara MAF manifest not found: ${params.ensembl_compara_maf_manifest}"
        }
        if (configured_maf_manifest) {
            maf_manifest = Channel.value(configured_maf_manifest)
        } else if (default_maf_manifest.exists()) {
            maf_manifest = Channel.value(default_maf_manifest)
        } else {
            BUILD_ENSEMBL_COMPARA_MAF_MANIFEST(genes, ensembl_compara_maf_manifest_script)
            maf_manifest = BUILD_ENSEMBL_COMPARA_MAF_MANIFEST.out.maf_manifest
        }
        ALIGN_ENSEMBL_COMPARA_MAF(
            task_dirs,
            maf_manifest,
            ensembl_compara_maf_script
        )
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_ENSEMBL_COMPARA_MAF.out.ensembl_compara_maf_result_dirs.map { meta, dir -> dir })
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

    main:
    annotate_script = file("${projectDir}/bin/annotate_events.py")
    ANNOTATE_EVENTS(events_tsv, genes_tsv, sequences_dir, annotate_script)

    emit:
    annotated_events = ANNOTATE_EVENTS.out.annotated_events
}

workflow ANNOTATION_STAGE_WITH_CLINVAR {
    take:
    events_tsv
    genes_tsv
    sequences_dir
    clinvar_vcf
    clinvar_vcf_tbi

    main:
    annotate_script = file("${projectDir}/bin/annotate_events.py")
    ANNOTATE_EVENTS_WITH_CLINVAR(events_tsv, genes_tsv, sequences_dir, annotate_script, clinvar_vcf, clinvar_vcf_tbi)

    emit:
    annotated_events = ANNOTATE_EVENTS_WITH_CLINVAR.out.annotated_events
}

workflow {
    runtime_check_script = file("${projectDir}/bin/check_runtime.py")
    CHECK_RUNTIME(runtime_check_script, params.stage, SELECTED_ALIGNMENT_STRATEGIES.join(','))

    if (params.stage == 'all') {
        clinvar_inputs = resolveClinvarInputs()
        if (clinvar_inputs.enabled) {
            log.info "Using ClinVar VCF: ${clinvar_inputs.path}"
        } else {
            log.info "ClinVar VCF is not configured; ClinVar annotation will be skipped"
        }
        FETCH_STAGE(file(params.ids_file))
        ALIGNMENT_STAGE(
            FETCH_STAGE.out.genes,
            FETCH_STAGE.out.target_features,
            FETCH_STAGE.out.orthologs_selected,
            FETCH_STAGE.out.sequences
        )
        if (clinvar_inputs.enabled) {
            ANNOTATION_STAGE_WITH_CLINVAR(
                ALIGNMENT_STAGE.out.events,
                FETCH_STAGE.out.genes,
                FETCH_STAGE.out.sequences,
                clinvar_inputs.vcf,
                clinvar_inputs.tbi
            )
        } else {
            ANNOTATION_STAGE(
                ALIGNMENT_STAGE.out.events,
                FETCH_STAGE.out.genes,
                FETCH_STAGE.out.sequences
            )
        }
    } else if (params.stage == 'fetch') {
        FETCH_STAGE(file(params.ids_file))
    } else if (params.stage == 'align') {
        ALIGNMENT_STAGE_FROM_DIR()
    } else if (params.stage == 'annotate') {
        clinvar_inputs = resolveClinvarInputs()
        if (clinvar_inputs.enabled) {
            log.info "Using ClinVar VCF: ${clinvar_inputs.path}"
        } else {
            log.info "ClinVar VCF is not configured; ClinVar annotation will be skipped"
        }
        fetch_dir = file(params.fetch_dir)
        if (clinvar_inputs.enabled) {
            ANNOTATION_STAGE_WITH_CLINVAR(
                file(params.events_tsv),
                file("${fetch_dir}/genes.tsv.gz"),
                file("${fetch_dir}/sequences"),
                clinvar_inputs.vcf,
                clinvar_inputs.tbi
            )
        } else {
            ANNOTATION_STAGE(
                file(params.events_tsv),
                file("${fetch_dir}/genes.tsv.gz"),
                file("${fetch_dir}/sequences")
            )
        }
    }
}
