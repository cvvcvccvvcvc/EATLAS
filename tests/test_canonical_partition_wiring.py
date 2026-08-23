from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]


def test_all_align_and_annotate_share_one_partition_evidence_contract() -> None:
    main = (PROJECT_DIR / "main.nf").read_text()
    partition_merge = (PROJECT_DIR / "modules/local/merge_alignment_partition.nf").read_text()
    final_merge = (PROJECT_DIR / "modules/local/merge_alignment.nf").read_text()

    assert "workflow ANNOTATION_STAGE" not in main
    assert main.count("PARTITIONED_ANNOTATION_STAGE(") == 2
    assert "params.stage" not in partition_merge
    assert "params.stage" not in final_merge
    assert "--output-profile" not in partition_merge
    assert "--output-profile" not in final_merge


def test_pipeline_modules_do_not_publish_derived_alignment_or_annotation_tables() -> None:
    module_paths = [
        PROJECT_DIR / "modules/local/align_minimap2.nf",
        PROJECT_DIR / "modules/local/align_nucmer_comparator.nf",
        PROJECT_DIR / "modules/local/align_bwa_pseudoreads.nf",
        PROJECT_DIR / "modules/local/merge_alignment_partition.nf",
        PROJECT_DIR / "modules/local/merge_alignment.nf",
        PROJECT_DIR / "modules/local/annotate_events_partition.nf",
        PROJECT_DIR / "modules/local/finalize_annotation.nf",
    ]
    module_text = "\n".join(path.read_text() for path in module_paths)
    forbidden = {
        "snv_site_depth.tsv.gz",
        "snv_taxonomic_depth.tsv.gz",
        "snv_alt_taxonomic_support.tsv.gz",
        "feature_coverage.tsv.gz",
        "variant_strategy_support.tsv.gz",
        "variant_ortholog_support",
        "ortholog_evidence_summary.tsv.gz",
    }
    for filename in forbidden:
        assert filename not in module_text


def test_small_annotation_partitions_do_not_request_large_run_memory() -> None:
    main = (PROJECT_DIR / "main.nf").read_text()

    assert "supportRowCount <= 1_000_000L" in main
    assert "supportRowCount <= 5_000_000L" in main
    assert "return 8" in main
    assert "return 16" in main

    annotation_runtime = "\n".join(
        (PROJECT_DIR / path).read_text()
        for path in (
            "bin/annotate_events.py",
            "bin/finalize_annotation_partitions.py",
        )
    )
    assert "variant_ortholog_support" not in annotation_runtime

    aligner_scripts = [
        PROJECT_DIR / "bin/run_minimap2_alignment.py",
        PROJECT_DIR / "bin/run_nucmer_alignment.py",
        PROJECT_DIR / "bin/run_bwa_pseudoreads.py",
        PROJECT_DIR / "bin/merge_ensembl_compara_maf_gene.py",
    ]
    aligner_text = "\n".join(path.read_text() for path in aligner_scripts)
    assert "feature_coverage.tsv.gz" not in aligner_text
    assert "feature_coverage_count" not in aligner_text
