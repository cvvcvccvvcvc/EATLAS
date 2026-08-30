from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]


def test_pipeline_has_one_end_to_end_partition_evidence_path() -> None:
    main = (PROJECT_DIR / "main.nf").read_text()
    partition_merge = (PROJECT_DIR / "modules/local/merge_alignment_partition.nf").read_text()
    final_merge = (PROJECT_DIR / "modules/local/merge_alignment.nf").read_text()

    assert "workflow ANNOTATION_STAGE" not in main
    assert main.count("PARTITIONED_ANNOTATION_STAGE(") == 1
    assert "ALIGNMENT_STAGE_FROM_DIR" not in main
    for parameter in ("params.stage", "params.fetch_dir", "params.alignment_dir"):
        assert parameter not in main
    assert "params.stage" not in partition_merge
    assert "params.stage" not in final_merge
    assert "--output-profile" not in partition_merge
    assert "--output-profile" not in final_merge


def test_removed_stage_modes_are_not_public_parameters() -> None:
    main = (PROJECT_DIR / "main.nf").read_text()
    config = (PROJECT_DIR / "nextflow.config").read_text()
    schema = (PROJECT_DIR / "nextflow_schema.json").read_text()

    for parameter in ("stage", "fetch_dir", "alignment_dir"):
        assert f"{parameter} =" not in config
        assert f'"{parameter}"' not in schema
    assert "removedExecutionParameters = ['stage', 'fetch_dir', 'alignment_dir']" in main
    assert "The pipeline has one end-to-end execution path" in main


def test_strategy_default_and_native_output_contract_are_unambiguous() -> None:
    main = (PROJECT_DIR / "main.nf").read_text()
    config = (PROJECT_DIR / "nextflow.config").read_text()
    schema = (PROJECT_DIR / "nextflow_schema.json").read_text()
    runtime = "\n".join(
        path.read_text()
        for path in (
            PROJECT_DIR / "bin/run_minimap2_alignment.py",
            PROJECT_DIR / "bin/run_nucmer_alignment.py",
            PROJECT_DIR / "bin/run_bwa_pseudoreads.py",
        )
    )

    assert "alignment_strategies = 'default'" in config
    assert "raw == 'default'" in main
    assert "raw == 'all'" not in main
    assert '"default": "default"' in schema
    assert "keep_native_alignments" not in config + schema
    assert "--keep-native" not in runtime


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


def test_analytics_derivations_do_not_import_pipeline_entrypoints() -> None:
    support_io = (PROJECT_DIR / "analytics/io/annotation_support.py").read_text()
    aggregate_io = (PROJECT_DIR / "analytics/io/alignment_aggregates.py").read_text()

    for source in (support_io, aggregate_io):
        assert "sys.path" not in source
        assert "from bin" not in source
    assert "from analytics.derivations" in support_io
    assert "from analytics.derivations" in aggregate_io


def test_pipeline_python_modules_are_staged_packages_without_path_bridges() -> None:
    module_text = "\n".join(
        path.read_text() for path in (PROJECT_DIR / "modules/local").glob("*.nf")
    )

    assert "PYTHONPATH" not in module_text
    for module_name in (
        "bin.check_runtime",
        "bin.run_minimap2_alignment",
        "bin.run_nucmer_alignment",
        "bin.run_bwa_pseudoreads",
        "bin.merge_alignment_results",
        "bin.finalize_annotation_partitions",
    ):
        assert f"python3 -m {module_name}" in module_text


def test_small_annotation_partitions_do_not_request_large_run_memory() -> None:
    main = (PROJECT_DIR / "main.nf").read_text()

    assert "eventRowCount <= 1_000_000L" in main
    assert "eventRowCount <= 5_000_000L" in main
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
    ]
    aligner_text = "\n".join(path.read_text() for path in aligner_scripts)
    assert "feature_coverage.tsv.gz" not in aligner_text
    assert "feature_coverage_count" not in aligner_text


def test_alignment_memory_is_sized_from_ortholog_sequence_volume() -> None:
    main = (PROJECT_DIR / "main.nf").read_text()
    config = (PROJECT_DIR / "nextflow.config").read_text()

    assert "row.ortholog_sequence_bp as Long" in main
    assert "orthologSequenceBp < 150_000_000L" in main
    assert "orthologSequenceBp < 600_000_000L" in main
    assert "longPseudoreadMinimapBaseMemoryGbForOrthologBp" in main
    assert "bwaBaseMemoryGbForOrthologBp" in main
    assert "meta.alignment_memory_gb" in config
    assert "retryStepGb" in config
    assert "attempt,cpus,memory,time" in config
