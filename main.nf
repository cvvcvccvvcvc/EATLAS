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


def alignmentResultProcessCount(selectedStrategies) {
    return selectedStrategies.size()
}

def geneIdFromFastaPath(value) {
    def name = value instanceof java.nio.file.Path ? value.getFileName().toString() : new File(value.toString()).name
    return name
        .replaceFirst(/\.fasta\.gz$/, '')
        .replaceFirst(/\.fa\.gz$/, '')
        .replaceFirst(/\.fasta$/, '')
        .replaceFirst(/\.fa$/, '')
}

def fastaFilesByGene(seqDir, subdir) {
    def matches = file("${seqDir}/${subdir}/*.fa.gz")
    def paths = matches instanceof List ? matches : [matches]
    return paths
        .sort { left, right -> left.toString() <=> right.toString() }
        .collect { path -> tuple(geneIdFromFastaPath(path), path) }
}

def resolveClinvarInputs() {
    def selectedVcf = params.clinvar_vcf
    if (!selectedVcf) {
        def assetClinvar = "${projectDir}/assets/reference/clinvar/clinvar.vcf.gz"
        selectedVcf = file(assetClinvar).exists() ? assetClinvar : null
    }
    if (!selectedVcf) {
        error "ClinVar VCF is required for annotation. Pass --clinvar_vcf, set CLINVAR_VCF, or place the file at assets/reference/clinvar/clinvar.vcf.gz"
    }
    def selectedTbi = "${selectedVcf}.tbi"
    if (!file(selectedVcf).exists()) {
        error "ClinVar VCF not found: ${selectedVcf}. Pass --clinvar_vcf, set CLINVAR_VCF, or place the file at assets/reference/clinvar/clinvar.vcf.gz"
    }
    if (!file(selectedTbi).exists()) {
        error "ClinVar VCF index not found: ${selectedTbi}. Place the .tbi next to the VCF."
    }
    return [vcf: file(selectedVcf), tbi: file(selectedTbi), path: selectedVcf]
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
ALIGNMENT_RESULT_PROCESS_COUNT = alignmentResultProcessCount(SELECTED_ALIGNMENT_STRATEGIES)

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
include { MERGE_ALIGNMENT_PARTITION } from './modules/local/merge_alignment_partition.nf'
include { ALIGN_NUCMER_COMPARATOR } from './modules/local/align_nucmer_comparator.nf'
include { ALIGN_BWA_PSEUDOREADS } from './modules/local/align_bwa_pseudoreads.nf'
include { BUILD_ENSEMBL_COMPARA_MAF_MANIFEST } from './modules/local/build_ensembl_compara_maf_manifest.nf'
include { ALIGN_ENSEMBL_COMPARA_MAF } from './modules/local/align_ensembl_compara_maf.nf'
include { BUILD_ENSEMBL_COMPARA_MAF_CHUNK_TASKS } from './modules/local/build_ensembl_compara_maf_chunk_tasks.nf'
include { ALIGN_ENSEMBL_COMPARA_MAF_CHUNK } from './modules/local/align_ensembl_compara_maf_chunk.nf'
include { MERGE_ENSEMBL_COMPARA_MAF_GENE } from './modules/local/merge_ensembl_compara_maf_gene.nf'
include { ANNOTATE_EVENTS } from './modules/local/annotate_events.nf'
include { ANNOTATE_EVENTS_PARTITION } from './modules/local/annotate_events_partition.nf'
include { FINALIZE_ANNOTATION } from './modules/local/finalize_annotation.nf'

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
    chunk_metrics = BUILD_FETCH_DATASET.out.chunk_metrics
    genes = BUILD_FETCH_DATASET.out.genes
    target_features = BUILD_FETCH_DATASET.out.target_features
    orthologs_selected = BUILD_FETCH_DATASET.out.orthologs_selected
    orthologs_candidates = BUILD_FETCH_DATASET.out.orthologs_candidates
    failures = BUILD_FETCH_DATASET.out.failures
    sequences = BUILD_FETCH_DATASET.out.sequences
}

workflow ALIGNMENT_STAGE {
    take:
    fetch_manifest
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
    ensembl_compara_maf_chunk_tasks_script = file("${projectDir}/bin/prepare_ensembl_compara_maf_chunk_tasks.py")
    ensembl_compara_maf_chunk_script = file("${projectDir}/bin/run_ensembl_compara_maf_chunk_alignment.py")
    ensembl_compara_maf_gene_merge_script = file("${projectDir}/bin/merge_ensembl_compara_maf_gene.py")
    merge_script = file("${projectDir}/bin/merge_alignment_results.py")

    FETCH_TAXONOMY_PRESETS(orthologs_selected, taxonomy_script, taxonomy_classes)
    BUILD_ALIGNMENT_TASKS(
        genes,
        orthologs_selected,
        fetch_manifest,
        target_features,
        sequences,
        FETCH_TAXONOMY_PRESETS.out.taxonomy_presets,
        prepare_script
    )

    task_dirs_by_gene_unpartitioned = BUILD_ALIGNMENT_TASKS.out.task_dirs.flatten().map { dir ->
        gene_id = dir.baseName.replaceFirst(/^task_/, '')
        tuple(gene_id, dir)
    }
    task_partitions = BUILD_ALIGNMENT_TASKS.out.alignment_tasks
        .splitCsv(header: true, sep: '\t', decompress: true)
        .filter { row -> row.status == 'ready' }
        .map { row -> tuple(row.gene_id as String, row.partition_id as String) }
    ready_genes_by_partition = task_partitions
        .map { gene_id, partition_id -> tuple(partition_id, gene_id) }
        .groupTuple()
        .map { partition_id, gene_ids ->
            tuple(partition_id, gene_ids.unique().sort())
        }
    task_dirs_by_gene = task_dirs_by_gene_unpartitioned
        .join(task_partitions)
        .map { gene_id, dir, partition_id -> tuple(gene_id, partition_id, dir) }
    target_fastas_by_gene = sequences.flatMap { seq_dir -> fastaFilesByGene(seq_dir, 'targets') }
    ortholog_fastas_by_gene = sequences.flatMap { seq_dir -> fastaFilesByGene(seq_dir, 'orthologs') }
    partition_genes = BUILD_ALIGNMENT_TASKS.out.partition_genes.flatten().map { path ->
        tuple(path.baseName.replaceFirst(/\.tsv$/, ''), path)
    }
    target_fastas_by_partition = task_dirs_by_gene
        .map { gene_id, partition_id, dir -> tuple(gene_id, partition_id) }
        .join(target_fastas_by_gene)
        .map { gene_id, partition_id, fasta -> tuple(partition_id, fasta) }
        .groupTuple()
    alignment_inputs = task_dirs_by_gene
        .join(target_fastas_by_gene)
        .join(ortholog_fastas_by_gene)
        .map { gene_id, partition_id, task_dir, source_target_fasta, source_ortholog_fasta ->
            tuple(
                [id: "task_${gene_id}", gene_id: gene_id, partition_id: partition_id],
                task_dir,
                source_target_fasta,
                source_ortholog_fasta
            )
        }
    task_dirs = task_dirs_by_gene.map { gene_id, partition_id, dir ->
        tuple([id: "task_${gene_id}", gene_id: gene_id, partition_id: partition_id], dir)
    }
    alignment_result_dirs = Channel.empty()

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('minimap2_asm10')) {
        ALIGN_MINIMAP2_ASM10(alignment_inputs, minimap2_script)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_MINIMAP2_ASM10.out.asm10_result_dirs)
    }

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('minimap2_asm20')) {
        ALIGN_MINIMAP2_ASM20(alignment_inputs, minimap2_script)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_MINIMAP2_ASM20.out.asm20_result_dirs)
    }

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('minimap2_taxonomy_adaptive')) {
        ALIGN_MINIMAP2_ADAPTIVE(alignment_inputs, minimap2_script)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_MINIMAP2_ADAPTIVE.out.adaptive_result_dirs)
    }

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('nucmer')) {
        ALIGN_NUCMER_COMPARATOR(alignment_inputs, nucmer_script)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_NUCMER_COMPARATOR.out.nucmer_result_dirs)
    }

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('bwa_pseudoreads')) {
        ALIGN_BWA_PSEUDOREADS(alignment_inputs, bwa_script, bam_filtering_script)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_BWA_PSEUDOREADS.out.bwa_result_dirs)
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
        BUILD_ENSEMBL_COMPARA_MAF_CHUNK_TASKS(
            genes,
            maf_manifest,
            ensembl_compara_maf_chunk_tasks_script
        )
        maf_chunk_task_dirs = BUILD_ENSEMBL_COMPARA_MAF_CHUNK_TASKS.out.chunk_task_dirs.flatten().map { dir ->
            tuple([id: dir.baseName], dir)
        }
        ALIGN_ENSEMBL_COMPARA_MAF_CHUNK(
            maf_chunk_task_dirs,
            ensembl_compara_maf_chunk_script
        )
        maf_fragments_by_gene = ALIGN_ENSEMBL_COMPARA_MAF_CHUNK.out.ensembl_compara_maf_gene_fragments
            .flatMap { meta, dirs ->
                def fragmentDirs = dirs instanceof List ? dirs : [dirs]
                fragmentDirs.collect { dir ->
                    def geneId = dir.baseName.replaceFirst(/^gene_/, '').replaceFirst(/__chunk_.*$/, '')
                    tuple(geneId, dir)
                }
            }
            .groupTuple()
        maf_gene_merge_inputs = task_dirs_by_gene
            .join(maf_fragments_by_gene)
            .map { gene_id, partition_id, task_dir, fragment_dirs ->
                tuple(
                    [id: "task_${gene_id}", gene_id: gene_id, partition_id: partition_id],
                    task_dir,
                    fragment_dirs
                )
            }
        MERGE_ENSEMBL_COMPARA_MAF_GENE(maf_gene_merge_inputs, ensembl_compara_maf_gene_merge_script)
        alignment_result_dirs = alignment_result_dirs.mix(
            MERGE_ENSEMBL_COMPARA_MAF_GENE.out.gene_result_dirs
        )
    }

    gene_result_dirs = alignment_result_dirs
        .ifEmpty { error "No alignment result directories were produced" }
        .map { meta, dir ->
            def key = tuple(meta.partition_id as String, meta.gene_id as String)
            tuple(groupKey(key, ALIGNMENT_RESULT_PROCESS_COUNT), dir)
        }
        .groupTuple(remainder: true)
        .map { key, dirs ->
            def target = key.getGroupTarget()
            if (dirs.size() != ALIGNMENT_RESULT_PROCESS_COUNT) {
                error(
                    "Incomplete alignment outputs for gene ${target[1]} in partition ${target[0]}: " +
                    "expected ${ALIGNMENT_RESULT_PROCESS_COUNT} result directories, observed ${dirs.size()}"
                )
            }
            tuple(target[0], target[1], dirs)
        }
    partition_merge_inputs = gene_result_dirs
        .combine(ready_genes_by_partition, by: 0)
        .map { partition_id, gene_id, dirs, expected_gene_ids ->
            def key = tuple(partition_id, expected_gene_ids)
            tuple(groupKey(key, expected_gene_ids.size()), gene_id, dirs)
        }
        .groupTuple(remainder: true)
        .map { key, observed_gene_ids, dirs_by_gene ->
            def target = key.getGroupTarget()
            def partition_id = target[0]
            def expected_gene_ids = target[1]
            def missing_gene_ids = expected_gene_ids.findAll { !observed_gene_ids.contains(it) }
            def unexpected_gene_ids = observed_gene_ids.findAll { !expected_gene_ids.contains(it) }
            if (
                observed_gene_ids.size() != expected_gene_ids.size() ||
                missing_gene_ids ||
                unexpected_gene_ids
            ) {
                error(
                    "Incomplete alignment partition ${partition_id}: " +
                    "expected ${expected_gene_ids.size()} genes, observed ${observed_gene_ids.size()}; " +
                    "missing=${missing_gene_ids}; unexpected=${unexpected_gene_ids}"
                )
            }
            tuple(
                [id: partition_id, partition_id: partition_id, gene_ids: expected_gene_ids],
                dirs_by_gene.flatten()
            )
        }
    MERGE_ALIGNMENT_PARTITION(
        partition_merge_inputs,
        SELECTED_ALIGNMENT_STRATEGIES.join(','),
        merge_script
    )

    MERGE_ALIGNMENT(
        BUILD_ALIGNMENT_TASKS.out.alignment_tasks,
        FETCH_TAXONOMY_PRESETS.out.taxonomy_presets,
        FETCH_TAXONOMY_PRESETS.out.taxonomy_failures,
        target_features,
        MERGE_ALIGNMENT_PARTITION.out.partition_dirs.map { meta, dir -> dir }.collect(),
        SELECTED_ALIGNMENT_STRATEGIES.join(','),
        merge_script
    )

    emit:
    manifest = MERGE_ALIGNMENT.out.manifest
    tasks = MERGE_ALIGNMENT.out.alignment_tasks
    taxonomy_presets = MERGE_ALIGNMENT.out.taxonomy_presets
    taxonomy_failures = MERGE_ALIGNMENT.out.taxonomy_failures
    summaries = MERGE_ALIGNMENT.out.summaries
    strategy_summary = MERGE_ALIGNMENT.out.strategy_summary
    segments = MERGE_ALIGNMENT.out.segments
    feature_coverage = MERGE_ALIGNMENT.out.feature_coverage
    events = MERGE_ALIGNMENT.out.events
    failures = MERGE_ALIGNMENT.out.failures
    partitions = MERGE_ALIGNMENT_PARTITION.out.partition_dirs
    partition_genes = partition_genes
    partition_target_fastas = target_fastas_by_partition
}

workflow ALIGNMENT_STAGE_FROM_DIR {
    main:
    fetch_dir = file(params.fetch_dir)
    genes = Channel.value(file("${fetch_dir}/genes.tsv.gz"))
    target_features = Channel.value(file("${fetch_dir}/target_features.tsv.gz"))
    orthologs_selected = Channel.value(file("${fetch_dir}/orthologs.selected.tsv.gz"))
    sequences = Channel.value(file("${fetch_dir}/sequences"))
    ALIGNMENT_STAGE(
        Channel.value(file("${fetch_dir}/manifest.json")),
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
    variant_annotations = ANNOTATE_EVENTS.out.variant_annotations
    manifest = ANNOTATE_EVENTS.out.manifest
    failures = ANNOTATE_EVENTS.out.failures
}

workflow PARTITIONED_ANNOTATION_STAGE {
    take:
    alignment_partitions
    partition_genes
    partition_target_fastas
    clinvar_vcf
    clinvar_vcf_tbi

    main:
    annotate_script = file("${projectDir}/bin/annotate_events.py")
    finalize_script = file("${projectDir}/bin/finalize_annotation_partitions.py")
    annotation_inputs = alignment_partitions
        .map { meta, dir -> tuple(meta.partition_id as String, meta, dir) }
        .join(partition_genes)
        .join(partition_target_fastas)
        .map { partition_id, meta, alignment_partition, genes_tsv, target_fastas ->
            tuple(meta, alignment_partition, genes_tsv, target_fastas)
        }
    ANNOTATE_EVENTS_PARTITION(
        annotation_inputs,
        annotate_script,
        clinvar_vcf,
        clinvar_vcf_tbi
    )
    FINALIZE_ANNOTATION(
        ANNOTATE_EVENTS_PARTITION.out.partition_dirs.map { meta, dir -> dir }.collect(),
        finalize_script
    )

    emit:
    variant_annotations = FINALIZE_ANNOTATION.out.variant_annotations
    manifest = FINALIZE_ANNOTATION.out.manifest
    failures = FINALIZE_ANNOTATION.out.failures
}

workflow {
    runtime_check_script = file("${projectDir}/bin/check_runtime.py")
    CHECK_RUNTIME(runtime_check_script, params.stage, SELECTED_ALIGNMENT_STRATEGIES.join(','))

    if (params.stage == 'all') {
        clinvar_inputs = resolveClinvarInputs()
        log.info "Using ClinVar VCF: ${clinvar_inputs.path}"
        FETCH_STAGE(file(params.ids_file))
        ALIGNMENT_STAGE(
            FETCH_STAGE.out.manifest,
            FETCH_STAGE.out.genes,
            FETCH_STAGE.out.target_features,
            FETCH_STAGE.out.orthologs_selected,
            FETCH_STAGE.out.sequences
        )
        PARTITIONED_ANNOTATION_STAGE(
            ALIGNMENT_STAGE.out.partitions,
            ALIGNMENT_STAGE.out.partition_genes,
            ALIGNMENT_STAGE.out.partition_target_fastas,
            clinvar_inputs.vcf,
            clinvar_inputs.tbi
        )
    } else if (params.stage == 'fetch') {
        FETCH_STAGE(file(params.ids_file))
    } else if (params.stage == 'align') {
        ALIGNMENT_STAGE_FROM_DIR()
    } else if (params.stage == 'annotate') {
        clinvar_inputs = resolveClinvarInputs()
        log.info "Using ClinVar VCF: ${clinvar_inputs.path}"
        fetch_dir = file(params.fetch_dir)
        ANNOTATION_STAGE(
            file(params.events_tsv),
            file("${fetch_dir}/genes.tsv.gz"),
            file("${fetch_dir}/sequences"),
            clinvar_inputs.vcf,
            clinvar_inputs.tbi
        )
    }
}
