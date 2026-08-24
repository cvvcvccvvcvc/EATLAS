from __future__ import annotations

import json
from pathlib import Path

import pytest

from analytics.io.performance import PerformanceProfile


def test_performance_profile_persists_nested_stages_and_resources(tmp_path: Path) -> None:
    analytics_dir = tmp_path / "analytics"
    analytics_dir.mkdir()
    report_path = tmp_path / "reports" / "report.html"
    profile_path = analytics_dir / "performance" / "report.json"
    profile = PerformanceProfile(
        profile_path,
        analysis_dir=analytics_dir,
        analysis_id="analysis-123",
        report_path=report_path,
        tracked_directory=analytics_dir,
        command=["strategy_report", "--run-dir", str(tmp_path)],
    )

    with profile.stage("Parent") as parent:
        parent["details"] = "cache miss"
        with profile.stage("Child"):
            (analytics_dir / "artifact.txt").write_text("payload")
            profile.add_metric("temporary_sqlite_bytes", 123)
    profile.disabled_stage("Optional", "disabled for test")
    report_path.parent.mkdir()
    report_path.write_text("<html></html>")
    profile.finish(artifacts=[report_path])

    payload = json.loads(profile_path.read_text())
    parent_record, child_record, disabled_record = payload["stages"]
    assert payload["status"] == "completed"
    assert payload["schema_version"] == 2
    assert payload["analysis_id"] == "analysis-123"
    assert payload["command"][0] == "strategy_report"
    assert parent_record["status"] == "completed"
    assert child_record["parent_id"] == parent_record["id"]
    assert child_record["metrics"]["temporary_sqlite_bytes"] == 123
    assert child_record["wall_seconds"] >= 0
    assert child_record["cpu_seconds"] >= 0
    assert child_record["process_peak_rss_bytes"] > 0
    assert disabled_record["status"] == "disabled"
    assert payload["summary"]["tracked_directory_bytes_after"] >= (
        payload["summary"]["tracked_directory_bytes_before"]
    )
    assert payload["artifacts"][0]["size_bytes"] == report_path.stat().st_size
    assert [row["Stage"] for row in profile.table_rows()] == ["Parent", "Optional"]


def test_performance_profile_flushes_failed_stage(tmp_path: Path) -> None:
    profile_path = tmp_path / "performance.json"
    profile = PerformanceProfile(
        profile_path,
        analysis_dir=tmp_path,
        analysis_id="analysis-failed",
        report_path=tmp_path / "report.html",
    )

    with pytest.raises(ValueError, match="broken stage"):
        with profile.stage("Failure"):
            raise ValueError("broken stage")

    payload = json.loads(profile_path.read_text())
    assert payload["status"] == "failed"
    assert payload["finished_at_utc"]
    assert payload["stages"][0]["status"] == "failed"
    assert payload["stages"][0]["error_type"] == "ValueError"


def test_performance_profile_records_source_run_provenance(tmp_path: Path) -> None:
    source_runs = [tmp_path / "run-a", tmp_path / "run-b"]
    profile_path = tmp_path / "analytics" / "performance.json"
    profile = PerformanceProfile(
        profile_path,
        analysis_dir=tmp_path / "analytics" / "analyses" / "analysis-123",
        analysis_id="analysis-123",
        report_path=tmp_path / "analytics" / "reports" / "report.html",
        source_run_dirs=source_runs,
    )
    profile.finish()

    payload = json.loads(profile_path.read_text())
    assert payload["analysis_id"] == "analysis-123"
    assert payload["source_run_dirs"] == [
        str(path.resolve()) for path in source_runs
    ]
