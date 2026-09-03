from pathlib import Path

from analytics.io.calculation_identity import (
    build_calculation_identity,
    calculation_cache_versions,
    repository_provenance,
)


PROJECT_DIR = Path(__file__).resolve().parents[1]


def test_calculation_identity_records_loaded_computation_runtime() -> None:
    identity = build_calculation_identity(
        firth_runtime={"R": "4.5.2", "logistf": "1.26.1"},
    )

    assert identity["cache_versions"] == calculation_cache_versions()
    assert identity["runtime"]["python"]["version"]
    assert set(identity["runtime"]["packages"]) == {
        "duckdb",
        "numpy",
        "pandas",
        "pyBigWig",
        "scipy",
        "statsmodels",
    }
    assert identity["runtime"]["R"] == {"R": "4.5.2", "logistf": "1.26.1"}


def test_repository_provenance_records_revision_and_dirty_state() -> None:
    provenance = repository_provenance(PROJECT_DIR)

    assert len(provenance["git_commit"]) == 40
    assert isinstance(provenance["git_dirty"], bool)
