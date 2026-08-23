process FINALIZE_FETCH_OUTPUT {
    tag "fetch_output"

    input:
    path manifest, stageAs: "source/manifest.json"
    path input_ids, stageAs: "source/input.ids.tsv"
    path genes, stageAs: "source/genes.tsv.gz"
    path target_features, stageAs: "source/target_features.tsv.gz"
    path orthologs_selected, stageAs: "source/orthologs.selected.tsv.gz"
    path failures, stageAs: "source/failures.tsv.gz"
    path source_sequences, stageAs: "source/sequences"
    path taxonomy, stageAs: "source/taxonomy.tsv.gz"
    path taxonomy_failures, stageAs: "source/taxonomy_failures.tsv.gz"

    output:
    path "manifest.json"
    path "input.ids.tsv"
    path "genes.tsv.gz"
    path "target_features.tsv.gz"
    path "orthologs.selected.tsv.gz"
    path "failures.tsv.gz"
    path "sequences"
    path "taxonomy.tsv.gz"
    path "taxonomy_failures.tsv.gz"

    script:
    """
    mkdir -p sequences
    cp "${manifest}" manifest.json
    cp "${input_ids}" input.ids.tsv
    cp "${genes}" genes.tsv.gz
    cp "${target_features}" target_features.tsv.gz
    cp "${orthologs_selected}" orthologs.selected.tsv.gz
    cp "${failures}" failures.tsv.gz
    cp -R "${source_sequences}/targets" sequences/targets
    cp "${taxonomy}" taxonomy.tsv.gz
    cp "${taxonomy_failures}" taxonomy_failures.tsv.gz
    """
}
