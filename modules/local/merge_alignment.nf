process MERGE_ALIGNMENT {
    tag "merge_alignment"

    input:
    path alignment_tasks
    path taxonomy_presets
    path taxonomy_failures
    path minimap2_asm10_dirs
    path minimap2_asm20_dirs
    path minimap2_adaptive_dirs
    path nucmer_dirs
    path bwa_dirs
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
