nextflow.enable.dsl = 2

import groovy.json.JsonSlurper
import RunManifest

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

ENSEMBL_COMPARA_STRATEGY = 'precomputed_ensembl_92_mammals_epo_extended'
ALIGNMENT_STRATEGY_REGISTRY = [
    [name: 'minimap2_asm10', default_enabled: true, minimap2_preset: 'asm10'],
    [name: 'minimap2_asm20', default_enabled: true, minimap2_preset: 'asm20'],
    [
        name: 'minimap2_map_ont_pseudoreads_30000_15000',
        default_enabled: false,
        minimap2_preset: 'map-ont',
        minimap2_pseudoread_len: 30000,
        minimap2_pseudoread_step: 15000,
    ],
    [name: 'nucmer', default_enabled: true],
    [
        name: 'bwa_pseudoreads_150_75',
        default_enabled: true,
        bwa_pseudoread_len: 150,
        bwa_pseudoread_step: 75,
        bwa_pseudoread_phred: 30,
    ],
    [name: ENSEMBL_COMPARA_STRATEGY, default_enabled: false],
]
AVAILABLE_ALIGNMENT_STRATEGIES = ALIGNMENT_STRATEGY_REGISTRY.collect { it.name }
DEFAULT_ALIGNMENT_STRATEGIES = ALIGNMENT_STRATEGY_REGISTRY
    .findAll { it.default_enabled }
    .collect { it.name }

def parseAlignmentStrategies(rawValue) {
    def raw = rawValue == null ? 'all' : rawValue.toString().trim()
    def selected = []
    if (!raw || raw == 'all') {
        selected = DEFAULT_ALIGNMENT_STRATEGIES
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


def annotationBaseMemoryGbForSupportRows(rawSupportRowCount) {
    def supportRowCount = rawSupportRowCount as Long
    if (supportRowCount < 0) {
        throw new IllegalArgumentException(
            "event_ortholog_support_count must be non-negative: ${supportRowCount}"
        )
    }
    if (supportRowCount <= 1_000_000L) {
        return 8
    }
    if (supportRowCount <= 5_000_000L) {
        return 16
    }
    if (supportRowCount <= 15_000_000L) {
        return 32
    }
    if (supportRowCount <= 30_000_000L) {
        return 48
    }
    if (supportRowCount <= 40_000_000L) {
        return 64
    }
    return 96
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

def removedExecutionParameters = ['stage', 'fetch_dir', 'alignment_dir'].findAll {
    params.containsKey(it)
}
if (removedExecutionParameters) {
    error(
        "Removed parameter(s): " +
        removedExecutionParameters.collect { "--${it}" }.join(', ') +
        ". The pipeline has one end-to-end execution path; use -resume to continue an interrupted run."
    )
}

// Validate parameters against nextflow_schema.json
validateParameters()

if (!params.ids_file) {
    error "Missing required parameter: --ids_file"
}

if (!file(params.target_annotation_gff3).exists()) {
    error "Target annotation GFF3 not found: ${params.target_annotation_gff3}. Pass --target_annotation_gff3, set GAPH_TARGET_ANNOTATION_GFF3, or place the file at assets/reference/ncbi/refseq/GCF_000001405.40_GRCh38.p14/genomic.gff.gz"
}

SELECTED_ALIGNMENT_STRATEGIES = parseAlignmentStrategies(params.alignment_strategies)
SELECTED_ORTHOLOG_ALIGNMENT_STRATEGIES = SELECTED_ALIGNMENT_STRATEGIES.findAll {
    it != ENSEMBL_COMPARA_STRATEGY
}
SELECTED_MINIMAP2_STRATEGIES = ALIGNMENT_STRATEGY_REGISTRY.findAll {
    it.minimap2_preset && SELECTED_ALIGNMENT_STRATEGIES.contains(it.name)
}
SELECTED_BWA_STRATEGIES = ALIGNMENT_STRATEGY_REGISTRY.findAll {
    it.bwa_pseudoread_len && SELECTED_ALIGNMENT_STRATEGIES.contains(it.name)
}

include { VALIDATE_IDS } from './modules/local/validate_ids.nf'
include { CHECK_RUNTIME } from './modules/local/check_runtime.nf'
include { FETCH_PARSE_CHUNK } from './modules/local/fetch_parse_chunk.nf'
include { BUILD_FETCH_DATASET } from './modules/local/build_fetch_dataset.nf'
include { FINALIZE_FETCH_OUTPUT } from './modules/local/finalize_fetch_output.nf'
include { FETCH_TAXONOMY } from './modules/local/fetch_taxonomy.nf'
include { BUILD_ALIGNMENT_TASKS } from './modules/local/build_alignment_tasks.nf'
include { ALIGN_MINIMAP2 } from './modules/local/align_minimap2.nf'
include { MERGE_ALIGNMENT } from './modules/local/merge_alignment.nf'
include { MERGE_ALIGNMENT_PARTITION } from './modules/local/merge_alignment_partition.nf'
include { ALIGN_NUCMER_COMPARATOR } from './modules/local/align_nucmer_comparator.nf'
include { ALIGN_BWA_PSEUDOREADS } from './modules/local/align_bwa_pseudoreads.nf'
include { BUILD_ENSEMBL_COMPARA_MAF_MANIFEST } from './modules/local/build_ensembl_compara_maf_manifest.nf'
include { BUILD_ENSEMBL_COMPARA_MAF_CHUNK_TASKS } from './modules/local/build_ensembl_compara_maf_chunk_tasks.nf'
include { ALIGN_ENSEMBL_COMPARA_MAF_CHUNK } from './modules/local/align_ensembl_compara_maf_chunk.nf'
include { MERGE_ENSEMBL_COMPARA_MAF_GENE } from './modules/local/merge_ensembl_compara_maf_gene.nf'
include { ANNOTATE_EVENTS_PARTITION } from './modules/local/annotate_events_partition.nf'
include { FINALIZE_ANNOTATION } from './modules/local/finalize_annotation.nf'
include { PREPARE_ANNOTATION_CONTEXTS } from './modules/local/prepare_annotation_contexts.nf'

workflow FETCH_STAGE {
    take:
    ids

    main:
    normalize_script = file("${projectDir}/bin/normalize_ids.py")
    fetch_script = file("${projectDir}/bin/fetch_parse_chunk.py")
    build_fetch_dataset_script = file("${projectDir}/bin/build_fetch_dataset.py")
    taxonomy_script = file("${projectDir}/bin/fetch_taxonomy.py")
    target_annotation_gff3 = file(params.target_annotation_gff3)
    request_throttle_dir = "${workflow.workDir}/.gaph/ncbi_fetch_throttle"

    VALIDATE_IDS(ids, normalize_script)

    chunk_files = VALIDATE_IDS.out.chunk_files.flatten().map { file -> tuple([id: file.baseName], file) }
    FETCH_PARSE_CHUNK(chunk_files, fetch_script, request_throttle_dir)

    BUILD_FETCH_DATASET(
        VALIDATE_IDS.out.ids_tsv,
        VALIDATE_IDS.out.chunks_tsv,
        FETCH_PARSE_CHUNK.out.chunk_dirs.map { meta, dir -> dir }.collect(),
        target_annotation_gff3,
        build_fetch_dataset_script
    )
    FETCH_TAXONOMY(BUILD_FETCH_DATASET.out.orthologs_selected, taxonomy_script)
    FINALIZE_FETCH_OUTPUT(
        BUILD_FETCH_DATASET.out.manifest,
        BUILD_FETCH_DATASET.out.input_ids,
        BUILD_FETCH_DATASET.out.genes,
        BUILD_FETCH_DATASET.out.target_features,
        BUILD_FETCH_DATASET.out.orthologs_selected,
        BUILD_FETCH_DATASET.out.failures,
        BUILD_FETCH_DATASET.out.sequences,
        FETCH_TAXONOMY.out.taxonomy,
        FETCH_TAXONOMY.out.taxonomy_failures
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
    taxonomy = FETCH_TAXONOMY.out.taxonomy
    taxonomy_failures = FETCH_TAXONOMY.out.taxonomy_failures
}

workflow ALIGNMENT_STAGE {
    take:
    genes
    target_features
    orthologs_selected
    sequences

    main:
    prepare_script = file("${projectDir}/bin/prepare_alignment_tasks.py")
    minimap2_script = file("${projectDir}/bin/run_minimap2_alignment.py")
    nucmer_script = file("${projectDir}/bin/run_nucmer_alignment.py")
    bwa_script = file("${projectDir}/bin/run_bwa_pseudoreads.py")
    bam_filtering_script = file("${projectDir}/bin/bam_filtering_v1.py")
    ensembl_compara_maf_manifest_script = file("${projectDir}/bin/build_ensembl_compara_maf_manifest.py")
    ensembl_compara_maf_chunk_tasks_script = file("${projectDir}/bin/prepare_ensembl_compara_maf_chunk_tasks.py")
    ensembl_compara_maf_chunk_script = file("${projectDir}/bin/run_ensembl_compara_maf_chunk_alignment.py")
    ensembl_compara_maf_gene_merge_script = file("${projectDir}/bin/merge_ensembl_compara_maf_gene.py")
    merge_script = file("${projectDir}/bin/merge_alignment_results.py")
    alignment_table_schema = file("${projectDir}/bin/alignment_table_schema.py")

    BUILD_ALIGNMENT_TASKS(
        genes,
        orthologs_selected,
        sequences,
        prepare_script
    )

    task_dirs_by_gene_unpartitioned = BUILD_ALIGNMENT_TASKS.out.task_dirs.flatten().map { dir ->
        gene_id = dir.baseName.replaceFirst(/^task_/, '')
        tuple(gene_id, dir)
    }
    task_capabilities = BUILD_ALIGNMENT_TASKS.out.alignment_tasks
        .splitCsv(header: true, sep: '\t', decompress: true)
        .map { row ->
            tuple(
                row.gene_id as String,
                row.partition_id as String,
                (row.target_ready as String) == 'true',
                (row.ortholog_ready as String) == 'true'
            )
        }
    target_task_partitions = task_capabilities
        .filter { gene_id, partition_id, target_ready, ortholog_ready -> target_ready }
        .map { gene_id, partition_id, target_ready, ortholog_ready -> tuple(gene_id, partition_id) }
    ortholog_task_partitions = task_capabilities
        .filter { gene_id, partition_id, target_ready, ortholog_ready -> ortholog_ready }
        .map { gene_id, partition_id, target_ready, ortholog_ready -> tuple(gene_id, partition_id) }
    eligible_task_partitions = task_capabilities
        .filter { gene_id, partition_id, target_ready, ortholog_ready ->
            (SELECTED_ALIGNMENT_STRATEGIES.contains(ENSEMBL_COMPARA_STRATEGY) && target_ready) ||
            (!SELECTED_ORTHOLOG_ALIGNMENT_STRATEGIES.isEmpty() && ortholog_ready)
        }
        .map { gene_id, partition_id, target_ready, ortholog_ready -> tuple(gene_id, partition_id) }
    eligible_genes_by_partition = eligible_task_partitions
        .map { gene_id, partition_id -> tuple(partition_id, gene_id) }
        .groupTuple()
        .map { partition_id, gene_ids ->
            tuple(partition_id, gene_ids.unique().sort())
        }
    target_task_dirs_by_gene = task_dirs_by_gene_unpartitioned
        .join(target_task_partitions)
        .map { gene_id, dir, partition_id -> tuple(gene_id, partition_id, dir) }
    ortholog_task_dirs_by_gene = task_dirs_by_gene_unpartitioned
        .join(ortholog_task_partitions)
        .map { gene_id, dir, partition_id -> tuple(gene_id, partition_id, dir) }
    eligible_task_dirs_by_gene = task_dirs_by_gene_unpartitioned
        .join(eligible_task_partitions)
        .map { gene_id, dir, partition_id -> tuple(gene_id, partition_id, dir) }
    target_fastas_by_gene = sequences.flatMap { seq_dir -> fastaFilesByGene(seq_dir, 'targets') }
    ortholog_fastas_by_gene = sequences.flatMap { seq_dir -> fastaFilesByGene(seq_dir, 'orthologs') }
    alignment_inputs = ortholog_task_dirs_by_gene
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
    minimap2_inputs = alignment_inputs.flatMap {
        meta, task_dir, source_target_fasta, source_ortholog_fasta ->
            SELECTED_MINIMAP2_STRATEGIES.collect { strategy ->
                tuple(
                    meta + [
                        strategy: strategy.name,
                        preset: strategy.minimap2_preset,
                        pseudoread_len: strategy.minimap2_pseudoread_len ?: 0,
                        pseudoread_step: strategy.minimap2_pseudoread_step ?: 0,
                    ],
                    task_dir,
                    source_target_fasta,
                    source_ortholog_fasta
                )
            }
    }
    bwa_inputs = alignment_inputs.flatMap {
        meta, task_dir, source_target_fasta, source_ortholog_fasta ->
            SELECTED_BWA_STRATEGIES.collect { strategy ->
                tuple(
                    meta + [
                        strategy: strategy.name,
                        pseudoread_len: strategy.bwa_pseudoread_len,
                        pseudoread_step: strategy.bwa_pseudoread_step,
                        pseudoread_phred: strategy.bwa_pseudoread_phred,
                    ],
                    task_dir,
                    source_target_fasta,
                    source_ortholog_fasta
                )
            }
    }
    alignment_result_dirs = Channel.empty()

    if (!SELECTED_MINIMAP2_STRATEGIES.isEmpty()) {
        ALIGN_MINIMAP2(minimap2_inputs, minimap2_script, alignment_table_schema)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_MINIMAP2.out.result_dirs)
    }

    if (SELECTED_ALIGNMENT_STRATEGIES.contains('nucmer')) {
        ALIGN_NUCMER_COMPARATOR(alignment_inputs, nucmer_script, alignment_table_schema)
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_NUCMER_COMPARATOR.out.nucmer_result_dirs)
    }

    if (!SELECTED_BWA_STRATEGIES.isEmpty()) {
        ALIGN_BWA_PSEUDOREADS(
            bwa_inputs,
            bwa_script,
            bam_filtering_script,
            alignment_table_schema
        )
        alignment_result_dirs = alignment_result_dirs.mix(ALIGN_BWA_PSEUDOREADS.out.bwa_result_dirs)
    }

    if (SELECTED_ALIGNMENT_STRATEGIES.contains(ENSEMBL_COMPARA_STRATEGY)) {
        default_maf_manifest = file("${projectDir}/assets/reference/ensembl/compara/release-116/92_mammals.epo_extended/ensembl_compara_maf_manifest.tsv.gz")
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
            ensembl_compara_maf_chunk_script,
            alignment_table_schema
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
        maf_gene_merge_inputs = target_task_dirs_by_gene
            .join(maf_fragments_by_gene)
            .map { gene_id, partition_id, task_dir, fragment_dirs ->
                tuple(
                    [id: "task_${gene_id}", gene_id: gene_id, partition_id: partition_id],
                    task_dir,
                    fragment_dirs
                )
            }
        MERGE_ENSEMBL_COMPARA_MAF_GENE(
            maf_gene_merge_inputs,
            ensembl_compara_maf_gene_merge_script,
            alignment_table_schema
        )
        alignment_result_dirs = alignment_result_dirs.mix(
            MERGE_ENSEMBL_COMPARA_MAF_GENE.out.gene_result_dirs
        )
    }

    gene_result_dirs = alignment_result_dirs
        .ifEmpty { error "No alignment result directories were produced" }
        .map { meta, dir ->
            tuple(meta.partition_id as String, meta.gene_id as String, dir)
        }
        .groupTuple(by: [0, 1])
    partition_merge_inputs = gene_result_dirs
        .combine(eligible_genes_by_partition, by: 0)
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
        BUILD_ALIGNMENT_TASKS.out.alignment_tasks,
        SELECTED_ALIGNMENT_STRATEGIES.join(','),
        merge_script,
        alignment_table_schema
    )

    MERGE_ALIGNMENT(
        BUILD_ALIGNMENT_TASKS.out.alignment_tasks,
        genes,
        target_features,
        MERGE_ALIGNMENT_PARTITION.out.partition_dirs.map { meta, dir -> dir }.collect(),
        SELECTED_ALIGNMENT_STRATEGIES.join(','),
        merge_script,
        alignment_table_schema
    )

    emit:
    manifest = MERGE_ALIGNMENT.out.manifest
    evidence = MERGE_ALIGNMENT.out.evidence
    failures = MERGE_ALIGNMENT.out.failures
}

workflow PARTITIONED_ANNOTATION_STAGE {
    take:
    alignment_evidence
    genes_tsv
    target_sequences_dir
    clinvar_vcf
    clinvar_vcf_tbi
    gnomad_cache_dir

    main:
    annotate_script = file("${projectDir}/bin/annotate_events.py")
    prepare_contexts_script = file("${projectDir}/bin/prepare_annotation_contexts.py")
    genomics_sources = [
        file("${projectDir}/genomics/__init__.py"),
        file("${projectDir}/genomics/clinvar.py"),
        file("${projectDir}/genomics/gnomad.py"),
        file("${projectDir}/genomics/gnomad_cache.py"),
        file("${projectDir}/genomics/variants.py"),
    ]
    finalize_script = file("${projectDir}/bin/finalize_annotation_partitions.py")
    PREPARE_ANNOTATION_CONTEXTS(
        alignment_evidence,
        genes_tsv,
        target_sequences_dir,
        prepare_contexts_script
    )
    alignment_partitions = alignment_evidence.flatMap { evidence_dir ->
        def partitionsRoot = evidence_dir.resolve('partitions').toFile()
        def directories = partitionsRoot.listFiles()
            .findAll { it.isDirectory() }
            .sort { left, right -> left.name <=> right.name }
        directories.collect { directory ->
            def partition_id = directory.name
            tuple(partition_id, [id: partition_id, partition_id: partition_id], file(directory))
        }
    }
    partition_contexts = PREPARE_ANNOTATION_CONTEXTS.out.context_dirs.flatten().map { context_dir ->
        tuple(context_dir.baseName, context_dir)
    }
    annotation_inputs = alignment_partitions
        .join(partition_contexts)
        .map { partition_id, meta, alignment_partition, context_dir ->
            def manifestPath = alignment_partition.resolve('manifest.json')
            if (!manifestPath.exists()) {
                error "Alignment partition ${partition_id} is missing manifest.json: ${alignment_partition}"
            }
            def partitionManifest = new JsonSlurper().parse(manifestPath.toFile())
            if (partitionManifest.schema != 'normalized_alignment_evidence_partition_v1') {
                error "Alignment partition ${partition_id} has unsupported schema=${partitionManifest.schema}"
            }
            if (partitionManifest.event_ortholog_support_count == null) {
                error "Alignment partition ${partition_id} manifest is missing event_ortholog_support_count"
            }
            def supportRowCount = partitionManifest.event_ortholog_support_count as Long
            def annotationMeta = meta + [
                annotation_event_ortholog_support_count: supportRowCount,
                annotation_memory_gb: annotationBaseMemoryGbForSupportRows(supportRowCount),
            ]
            def contextGenes = context_dir.resolve('genes.tsv.gz')
            def targetFastas = context_dir.resolve('targets').toFile().listFiles()
                .findAll { it.isFile() && it.name.endsWith('.fa.gz') }
                .sort { left, right -> left.name <=> right.name }
                .collect { file(it) }
            if (!contextGenes.exists() || !targetFastas) {
                error "Annotation context is incomplete for partition ${partition_id}: ${context_dir}"
            }
            tuple(annotationMeta, alignment_partition, contextGenes, targetFastas)
        }
    ANNOTATE_EVENTS_PARTITION(
        annotation_inputs,
        annotate_script,
        genomics_sources,
        clinvar_vcf,
        clinvar_vcf_tbi,
        gnomad_cache_dir
    )
    FINALIZE_ANNOTATION(
        ANNOTATE_EVENTS_PARTITION.out.partition_dirs.map { meta, dir -> dir }.collect(),
        finalize_script
    )

    emit:
    variant_annotations = FINALIZE_ANNOTATION.out.variant_annotations
    event_variant_map = FINALIZE_ANNOTATION.out.event_variant_map
    manifest = FINALIZE_ANNOTATION.out.manifest
    failures = FINALIZE_ANNOTATION.out.failures
}

workflow {
    run_manifest_path = file("${params.outdir}/run_manifest.json")
    RunManifest.start(
        run_manifest_path,
        projectDir,
        workflow,
        params,
        file("${projectDir}/nextflow_schema.json")
    )

    if (params.gnomad_cache_dir) {
        log.info "Using shared gnomAD cache: ${params.gnomad_cache_dir}"
    } else {
        log.warn "Shared gnomAD cache is disabled; set --gnomad_cache_dir or GAPH_GNOMAD_CACHE_DIR to reuse regional responses"
    }

    runtime_check_script = file("${projectDir}/bin/check_runtime.py")
    CHECK_RUNTIME(runtime_check_script, SELECTED_ALIGNMENT_STRATEGIES.join(','))

    clinvar_inputs = resolveClinvarInputs()
    log.info "Using ClinVar VCF: ${clinvar_inputs.path}"
    FETCH_STAGE(file(params.ids_file))
    ALIGNMENT_STAGE(
        FETCH_STAGE.out.genes,
        FETCH_STAGE.out.target_features,
        FETCH_STAGE.out.orthologs_selected,
        FETCH_STAGE.out.sequences
    )
    PARTITIONED_ANNOTATION_STAGE(
        ALIGNMENT_STAGE.out.evidence,
        FETCH_STAGE.out.genes,
        FETCH_STAGE.out.sequences.map { sequences_dir -> file("${sequences_dir}/targets") },
        clinvar_inputs.vcf,
        clinvar_inputs.tbi,
        params.gnomad_cache_dir ?: ''
    )
}

workflow.onComplete {
    RunManifest.finish(file("${params.outdir}/run_manifest.json"), workflow)
}
