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
    path provenance_sources, stageAs: "provenance/*"

    output:
    path "manifest.json", emit: manifest
    path "input.ids.tsv", emit: input_ids
    path "genes.tsv.gz", emit: genes
    path "target_features.tsv.gz", emit: target_features
    path "orthologs.selected.tsv.gz", emit: orthologs_selected
    path "failures.tsv.gz", emit: failures
    path "sequences", emit: sequences
    path "taxonomy.tsv.gz", emit: taxonomy
    path "taxonomy_failures.tsv.gz", emit: taxonomy_failures
    path "fetch.inventory.json", emit: inventory

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
    python3 -m provenance.evidence_inventory create --scope fetch \\
        --output fetch.inventory.json \\
        --input fetch/manifest.json=manifest.json \\
        --input fetch/input.ids.tsv=input.ids.tsv \\
        --input fetch/genes.tsv.gz=genes.tsv.gz \\
        --input fetch/target_features.tsv.gz=target_features.tsv.gz \\
        --input fetch/orthologs.selected.tsv.gz=orthologs.selected.tsv.gz \\
        --input fetch/failures.tsv.gz=failures.tsv.gz \\
        --input fetch/sequences=sequences \\
        --input fetch/taxonomy.tsv.gz=taxonomy.tsv.gz \\
        --input fetch/taxonomy_failures.tsv.gz=taxonomy_failures.tsv.gz
    """
}
