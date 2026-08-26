"""Low-overhead runtime profiling for long analytics report jobs."""

from __future__ import annotations

import os
import platform
import resource
import sys
import time
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from .artifacts import path_metadata, write_json_atomic


PROFILE_SCHEMA_VERSION = 2


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _rss_bytes(value: int | float) -> int:
    # macOS reports bytes; Linux and the ITMO Slurm nodes report KiB.
    return int(value) if sys.platform == "darwin" else int(value) * 1024


def _usage_snapshot() -> dict[str, float | int]:
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "user_cpu_seconds": float(own.ru_utime + children.ru_utime),
        "system_cpu_seconds": float(own.ru_stime + children.ru_stime),
        "process_peak_rss_bytes": _rss_bytes(own.ru_maxrss),
        "children_peak_rss_bytes": _rss_bytes(children.ru_maxrss),
        "input_blocks": int(own.ru_inblock + children.ru_inblock),
        "output_blocks": int(own.ru_oublock + children.ru_oublock),
    }


def _directory_size(path: Path | None) -> int | None:
    if path is None:
        return None
    if not path.exists():
        return 0
    total = 0
    pending = [path]
    while pending:
        directory = pending.pop()
        try:
            with os.scandir(directory) as entries:
                for entry in entries:
                    try:
                        if entry.is_dir(follow_symlinks=False):
                            pending.append(Path(entry.path))
                        elif entry.is_file(follow_symlinks=False):
                            total += entry.stat(follow_symlinks=False).st_size
                    except FileNotFoundError:
                        continue
        except FileNotFoundError:
            continue
    return total


class PerformanceProfile:
    """Collect nested stage metrics and atomically persist progress after each stage."""

    def __init__(
        self,
        path: Path,
        *,
        analysis_dir: Path,
        analysis_id: str,
        report_path: Path,
        tracked_directory: Path | None = None,
        command: Sequence[str] | None = None,
        source_run_dirs: Sequence[Path] = (),
    ) -> None:
        self.path = path
        self.analysis_dir = analysis_dir
        self.analysis_id = analysis_id
        self.report_path = report_path
        self.tracked_directory = tracked_directory
        self.command = list(command or sys.argv)
        self.source_run_dirs = tuple(source_run_dirs)
        self.status = "running"
        self.started_at_utc = _utc_now()
        self.finished_at_utc: str | None = None
        self.stages: list[dict[str, object]] = []
        self._stack: list[str] = []
        self._counter = 0
        self._started = time.perf_counter()
        self._start_usage = _usage_snapshot()
        self._tracked_bytes_before = _directory_size(tracked_directory)
        self._tracked_bytes_after: int | None = None
        self._artifacts: list[dict[str, object]] = []
        self._flush()

    @contextmanager
    def stage(self, name: str) -> Iterator[dict[str, object]]:
        self._counter += 1
        stage_id = f"stage-{self._counter:03d}"
        record: dict[str, object] = {
            "id": stage_id,
            "name": name,
            "parent_id": self._stack[-1] if self._stack else None,
            "depth": len(self._stack),
            "status": "running",
            "details": "",
            "started_at_utc": _utc_now(),
            "metrics": {},
        }
        started = time.perf_counter()
        start_usage = _usage_snapshot()
        self.stages.append(record)
        self._stack.append(stage_id)
        self._flush()
        try:
            yield record
        except BaseException as error:
            record["status"] = "failed"
            record["error_type"] = type(error).__name__
            if not record["details"]:
                record["details"] = str(error)
            self.status = "failed"
            self.finished_at_utc = _utc_now()
            raise
        else:
            record["status"] = "completed"
        finally:
            end_usage = _usage_snapshot()
            wall_seconds = time.perf_counter() - started
            user_cpu = float(end_usage["user_cpu_seconds"]) - float(
                start_usage["user_cpu_seconds"]
            )
            system_cpu = float(end_usage["system_cpu_seconds"]) - float(
                start_usage["system_cpu_seconds"]
            )
            cpu_seconds = user_cpu + system_cpu
            record.update(
                {
                    "finished_at_utc": _utc_now(),
                    "wall_seconds": round(wall_seconds, 6),
                    "user_cpu_seconds": round(user_cpu, 6),
                    "system_cpu_seconds": round(system_cpu, 6),
                    "cpu_seconds": round(cpu_seconds, 6),
                    "cpu_to_wall_ratio": round(cpu_seconds / wall_seconds, 4)
                    if wall_seconds > 0
                    else None,
                    "process_peak_rss_bytes": end_usage["process_peak_rss_bytes"],
                    "children_peak_rss_bytes": end_usage["children_peak_rss_bytes"],
                    "input_blocks": int(end_usage["input_blocks"])
                    - int(start_usage["input_blocks"]),
                    "output_blocks": int(end_usage["output_blocks"])
                    - int(start_usage["output_blocks"]),
                }
            )
            if self._stack and self._stack[-1] == stage_id:
                self._stack.pop()
            else:
                self._stack.remove(stage_id)
            self._flush()
            indent = "  " * int(record["depth"])
            print(
                f"{indent}{name}: {record['status']} in "
                f"{float(record['wall_seconds']):.3f} s"
            )

    def disabled_stage(self, name: str, details: str) -> None:
        usage = _usage_snapshot()
        self._counter += 1
        self.stages.append(
            {
                "id": f"stage-{self._counter:03d}",
                "name": name,
                "parent_id": self._stack[-1] if self._stack else None,
                "depth": len(self._stack),
                "status": "disabled",
                "details": details,
                "started_at_utc": _utc_now(),
                "finished_at_utc": _utc_now(),
                "wall_seconds": 0.0,
                "user_cpu_seconds": 0.0,
                "system_cpu_seconds": 0.0,
                "cpu_seconds": 0.0,
                "cpu_to_wall_ratio": None,
                "process_peak_rss_bytes": usage["process_peak_rss_bytes"],
                "children_peak_rss_bytes": usage["children_peak_rss_bytes"],
                "input_blocks": 0,
                "output_blocks": 0,
                "metrics": {},
            }
        )
        self._flush()

    def add_metric(self, name: str, value: object) -> None:
        if not self._stack:
            raise RuntimeError("Performance metrics require an active stage")
        stage_id = self._stack[-1]
        record = next(record for record in reversed(self.stages) if record["id"] == stage_id)
        metrics = record["metrics"]
        assert isinstance(metrics, dict)
        metrics[name] = value

    def checkpoint(self, metrics: dict[str, object] | None = None) -> None:
        """Persist current-stage progress without closing the stage."""

        if not self._stack:
            raise RuntimeError("Performance checkpoints require an active stage")
        for name, value in (metrics or {}).items():
            self.add_metric(name, value)
        self._flush()

    def table_rows(self) -> list[dict[str, object]]:
        rows = []
        for stage in self.stages:
            if stage["depth"] != 0 or stage["status"] == "running":
                continue
            rows.append(
                {
                    "Stage": stage["name"],
                    "Status": stage["status"],
                    "Wall s": stage["wall_seconds"],
                    "CPU s": stage["cpu_seconds"],
                    "CPU / wall": stage["cpu_to_wall_ratio"],
                    "Process peak RSS MiB": round(
                        int(stage["process_peak_rss_bytes"]) / (1024 * 1024), 1
                    ),
                    "Details": stage["details"],
                }
            )
        return rows

    def finish(self, *, artifacts: Sequence[Path] = ()) -> None:
        if self.status != "failed":
            self.status = "completed"
            self.finished_at_utc = _utc_now()
        self._tracked_bytes_after = _directory_size(self.tracked_directory)
        self._artifacts = [
            path_metadata(path)
            for path in artifacts
            if path.exists() and path.is_file()
        ]
        self._flush()

    @property
    def total_wall_seconds(self) -> float:
        return time.perf_counter() - self._started

    def _summary(self) -> dict[str, object]:
        usage = _usage_snapshot()
        user_cpu = float(usage["user_cpu_seconds"]) - float(
            self._start_usage["user_cpu_seconds"]
        )
        system_cpu = float(usage["system_cpu_seconds"]) - float(
            self._start_usage["system_cpu_seconds"]
        )
        wall_seconds = self.total_wall_seconds
        cpu_seconds = user_cpu + system_cpu
        tracked_delta = None
        if self._tracked_bytes_before is not None and self._tracked_bytes_after is not None:
            tracked_delta = self._tracked_bytes_after - self._tracked_bytes_before
        return {
            "wall_seconds": round(wall_seconds, 6),
            "user_cpu_seconds": round(user_cpu, 6),
            "system_cpu_seconds": round(system_cpu, 6),
            "cpu_seconds": round(cpu_seconds, 6),
            "cpu_to_wall_ratio": round(cpu_seconds / wall_seconds, 4)
            if wall_seconds > 0
            else None,
            "process_peak_rss_bytes": usage["process_peak_rss_bytes"],
            "children_peak_rss_bytes": usage["children_peak_rss_bytes"],
            "input_blocks": int(usage["input_blocks"])
            - int(self._start_usage["input_blocks"]),
            "output_blocks": int(usage["output_blocks"])
            - int(self._start_usage["output_blocks"]),
            "tracked_directory_bytes_before": self._tracked_bytes_before,
            "tracked_directory_bytes_after": self._tracked_bytes_after,
            "tracked_directory_bytes_delta": tracked_delta,
        }

    def _payload(self) -> dict[str, object]:
        payload = {
            "schema_version": PROFILE_SCHEMA_VERSION,
            "status": self.status,
            "started_at_utc": self.started_at_utc,
            "finished_at_utc": self.finished_at_utc,
            "analysis_id": self.analysis_id,
            "analysis_dir": str(self.analysis_dir.resolve()),
            "source_run_dirs": [
                str(path.resolve()) for path in self.source_run_dirs
            ],
            "report_path": str(self.report_path.resolve()),
            "profile_path": str(self.path.resolve()),
            "command": self.command,
            "host": platform.node(),
            "platform": platform.platform(),
            "python_version": platform.python_version(),
            "summary": self._summary(),
            "stages": self.stages,
            "artifacts": self._artifacts,
        }
        return payload

    def _flush(self) -> None:
        write_json_atomic(self.path, self._payload())


@contextmanager
def profile_stage(
    profile: PerformanceProfile | None,
    name: str,
) -> Iterator[dict[str, object]]:
    """Open a persisted stage when profiling is enabled."""
    if profile is None:
        yield {}
        return
    with profile.stage(name) as record:
        yield record
