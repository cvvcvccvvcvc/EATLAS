"""Shared CPU and memory policy for analytics DuckDB calculations."""

from __future__ import annotations

import os
import re

from analytics.io.variant_source import sql_string


DUCKDB_MEMORY_LIMIT_ENV = "GAPH_DUCKDB_MEMORY_LIMIT"
DUCKDB_MEMORY_FRACTION = 0.5

_MEMORY_SETTING = re.compile(r"^([0-9]+(?:\.[0-9]+)?)\s*([KMGT]i?B)$", re.IGNORECASE)
_MEMORY_UNITS = {
    "KB": 1000,
    "MB": 1000**2,
    "GB": 1000**3,
    "TB": 1000**4,
    "KIB": 1024,
    "MIB": 1024**2,
    "GIB": 1024**3,
    "TIB": 1024**4,
}


def available_cpu_count() -> int:
    """Return CPUs available to this process, respecting Slurm allocation."""

    allocated = os.environ.get("SLURM_CPUS_PER_TASK")
    if allocated and allocated.isdigit() and int(allocated) > 0:
        return int(allocated)
    return os.cpu_count() or 1


def configure_duckdb_memory(connection, thread_count: int) -> dict[str, object]:
    override = os.environ.get(DUCKDB_MEMORY_LIMIT_ENV, "").strip()
    if override:
        requested = override
        source = DUCKDB_MEMORY_LIMIT_ENV
    else:
        slurm_bytes, source = _slurm_memory_bytes(thread_count)
        if slurm_bytes is not None:
            requested = _memory_limit_setting(slurm_bytes * DUCKDB_MEMORY_FRACTION)
        else:
            current = str(
                connection.execute("SELECT current_setting('memory_limit')").fetchone()[0]
            )
            requested = _memory_limit_setting(
                _parse_memory_setting(current) * DUCKDB_MEMORY_FRACTION
            )
            source = "duckdb_default_fraction"

    connection.execute(f"SET memory_limit={sql_string(requested)}")
    return {
        "memory_limit": str(
            connection.execute("SELECT current_setting('memory_limit')").fetchone()[0]
        ),
        "memory_limit_source": source,
    }


def _slurm_memory_bytes(thread_count: int) -> tuple[int | None, str]:
    per_node_mb = _positive_int_environment("SLURM_MEM_PER_NODE")
    if per_node_mb is not None:
        return per_node_mb * 1024**2, "SLURM_MEM_PER_NODE"

    per_cpu_mb = _positive_int_environment("SLURM_MEM_PER_CPU")
    if per_cpu_mb is None:
        return None, "duckdb_default_fraction"
    allocated_cpus = (
        _positive_int_environment("SLURM_CPUS_PER_TASK") or thread_count
    )
    return per_cpu_mb * allocated_cpus * 1024**2, "SLURM_MEM_PER_CPU"


def _positive_int_environment(name: str) -> int | None:
    value = os.environ.get(name, "").strip()
    if not value.isdigit() or int(value) < 1:
        return None
    return int(value)


def _parse_memory_setting(value: str) -> int:
    match = _MEMORY_SETTING.fullmatch(value.strip())
    if match is None:
        raise ValueError(f"Could not parse DuckDB memory limit: {value!r}")
    amount, unit = match.groups()
    return int(float(amount) * _MEMORY_UNITS[unit.upper()])


def _memory_limit_setting(value: float) -> str:
    mebibytes = max(128, int(value) // (1024**2))
    return f"{mebibytes}MiB"
