process BUILD_EVIDENCE_INVENTORY {
    tag "evidence_inventory"
    input:
    path fetch_inventory
    path alignment_inventory
    path annotation_inventory
    path provenance_sources, stageAs: "provenance/*"

    output:
    path "evidence_inventory.json", emit: inventory

    script:
    """
    python3 -m provenance.evidence_inventory combine \\
        --output evidence_inventory.json \\
        --fetch "${fetch_inventory}" \\
        --alignment "${alignment_inventory}" \\
        --annotation "${annotation_inventory}"
    """
}
