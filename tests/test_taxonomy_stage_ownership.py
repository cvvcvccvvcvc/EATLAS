from __future__ import annotations

from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]


def source_between(path: Path, start: str, end: str) -> str:
    source = path.read_text()
    return source.split(start, 1)[1].split(end, 1)[0]


def test_taxonomy_fetch_is_owned_by_stage_one() -> None:
    workflow = PROJECT_DIR / "main.nf"
    fetch_stage = source_between(
        workflow,
        "workflow FETCH_STAGE {",
        "workflow ALIGNMENT_STAGE {",
    )
    alignment_stage = source_between(
        workflow,
        "workflow ALIGNMENT_STAGE {",
        "workflow ALIGNMENT_STAGE_FROM_DIR {",
    )

    assert fetch_stage.count("FETCH_TAXONOMY(") == 1
    assert "FETCH_TAXONOMY(" not in alignment_stage
    for output in ("taxonomy", "taxonomy_failures", "taxonomy_summary"):
        assert f"{output} = FETCH_TAXONOMY.out.{output}" in fetch_stage
        assert output in alignment_stage.split("main:", 1)[0]


def test_standalone_alignment_requires_published_taxonomy() -> None:
    workflow = source_between(
        PROJECT_DIR / "main.nf",
        "workflow ALIGNMENT_STAGE_FROM_DIR {",
        "workflow ANNOTATION_STAGE {",
    )

    for filename in (
        "taxonomy.tsv.gz",
        "taxonomy_failures.tsv.gz",
        "taxonomy_summary.tsv.gz",
    ):
        assert filename in workflow
    assert "alignment does not fetch taxonomy metadata" in workflow


def test_fetch_publication_contains_taxonomy_handoff() -> None:
    finalizer = (PROJECT_DIR / "modules" / "local" / "finalize_fetch_output.nf").read_text()
    for filename in (
        "taxonomy.tsv.gz",
        "taxonomy_failures.tsv.gz",
        "taxonomy_summary.tsv.gz",
    ):
        assert finalizer.count(f'path "{filename}"') == 1
        assert filename in finalizer.split("script:", 1)[1]

    config = source_between(
        PROJECT_DIR / "nextflow.config",
        "withName: FETCH_TAXONOMY {",
        "withName: BUILD_ALIGNMENT_TASKS {",
    )
    assert 'conda = "${projectDir}/envs/fetch.yml"' in config
    assert "enabled: params.stage == 'fetch'" in config

    build_config = source_between(
        PROJECT_DIR / "nextflow.config",
        "withName: BUILD_FETCH_DATASET {",
        "withName: FINALIZE_FETCH_OUTPUT {",
    )
    assert "enabled: params.stage == 'fetch'" in build_config
