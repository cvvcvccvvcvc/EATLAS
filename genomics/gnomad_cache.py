#!/usr/bin/env python3
"""Shared, resumable cache for gnomAD regional GraphQL responses."""

from __future__ import annotations

import gzip
import json
import logging
import os
import tempfile
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Callable
from urllib.error import HTTPError, URLError

from .gnomad import (
    GNOMAD_MAX_ATTEMPTS,
    fetch_region_variants_recursive,
    is_retryable_network_error,
)


logger = logging.getLogger(__name__)

CACHE_SCHEMA_VERSION = 1
DEFAULT_TILE_SIZE_BP = 25_000
DEFAULT_GROUP_ATTEMPTS = 2
GNOMAD_DATASET = "gnomad_r4"
GNOMAD_REFERENCE_GENOME = "GRCh38"

FetchRegion = Callable[..., list[dict]]


@dataclass(frozen=True, slots=True)
class Tile:
    chrom: str
    start: int
    end: int


def normalize_chrom(chrom: object) -> str:
    value = str(chrom or "").strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    if value == "M":
        value = "MT"
    if not value or any(character in value for character in "/\\\0"):
        raise ValueError(f"Invalid chromosome: {chrom!r}")
    return value


def tiles_for_region(
    chrom: object,
    start: int,
    end: int,
    tile_size_bp: int = DEFAULT_TILE_SIZE_BP,
) -> list[Tile]:
    if tile_size_bp < 1:
        raise ValueError("tile_size_bp must be >= 1")
    normalized_chrom = normalize_chrom(chrom)
    start = int(start)
    end = int(end)
    if start < 1 or end < start:
        raise ValueError(f"Invalid region: {normalized_chrom}:{start}-{end}")

    first_start = ((start - 1) // tile_size_bp) * tile_size_bp + 1
    last_start = ((end - 1) // tile_size_bp) * tile_size_bp + 1
    return [
        Tile(normalized_chrom, tile_start, tile_start + tile_size_bp - 1)
        for tile_start in range(first_start, last_start + 1, tile_size_bp)
    ]


class GnomadRegionCache:
    """Fetch gnomAD regions while persisting complete fixed-size tiles."""

    def __init__(
        self,
        cache_dir: Path | str | None,
        *,
        fetcher: FetchRegion = fetch_region_variants_recursive,
        tile_size_bp: int = DEFAULT_TILE_SIZE_BP,
        group_attempts: int = DEFAULT_GROUP_ATTEMPTS,
        max_attempts: int = GNOMAD_MAX_ATTEMPTS,
    ) -> None:
        if tile_size_bp < 1:
            raise ValueError("tile_size_bp must be >= 1")
        if group_attempts < 1:
            raise ValueError("group_attempts must be >= 1")
        if max_attempts < group_attempts:
            raise ValueError("max_attempts must be >= group_attempts")

        self.cache_dir = Path(cache_dir).expanduser().resolve() if cache_dir else None
        self.namespace_dir = (
            self.cache_dir
            / GNOMAD_DATASET
            / GNOMAD_REFERENCE_GENOME
            / f"schema_v{CACHE_SCHEMA_VERSION}"
            / f"tiles_{tile_size_bp}bp"
            if self.cache_dir is not None
            else None
        )
        self.fetcher = fetcher
        self.tile_size_bp = tile_size_bp
        self.group_attempts = group_attempts
        self.max_attempts = max_attempts
        self._stats_lock = threading.Lock()
        self._stats = {
            "tile_hit_count": 0,
            "tile_miss_count": 0,
            "tile_write_count": 0,
            "corrupt_tile_count": 0,
            "fetch_batch_count": 0,
            "split_count": 0,
        }

    @property
    def enabled(self) -> bool:
        return self.namespace_dir is not None

    def fetch_region(self, chrom: object, start: int, end: int) -> list[dict]:
        normalized_chrom = normalize_chrom(chrom)
        start = int(start)
        end = int(end)
        if start < 1 or end < start:
            raise ValueError(f"Invalid region: {normalized_chrom}:{start}-{end}")

        if not self.enabled:
            return self._call_fetcher(normalized_chrom, start, end, self.max_attempts)

        tiles = tiles_for_region(normalized_chrom, start, end, self.tile_size_bp)
        records_by_tile = self.fetch_tiles(tiles)

        records = [record for tile in tiles for record in records_by_tile[tile]]
        records = [
            record
            for record in records
            if normalize_chrom(record.get("chrom")) == normalized_chrom
            and start <= int(record.get("pos", 0)) <= end
        ]
        records.sort(key=_record_sort_key)
        return records

    def fetch_tiles(self, tiles: list[Tile]) -> dict[Tile, list[dict]]:
        """Return complete fixed-size tiles, fetching only missing groups."""

        if not tiles:
            return {}
        if not self.enabled:
            raise RuntimeError("Fixed-tile access requires an enabled gnomAD cache")
        for tile in tiles:
            if tile.end - tile.start + 1 != self.tile_size_bp:
                raise ValueError(f"Unexpected gnomAD tile size: {tile}")

        records_by_tile: dict[Tile, list[dict]] = {}
        missing_tiles: list[Tile] = []
        for tile in tiles:
            records = self.read_cached_tile(tile)
            if records is None:
                missing_tiles.append(tile)
            else:
                records_by_tile[tile] = records

        for group in _consecutive_groups(missing_tiles):
            records_by_tile.update(self._fetch_and_store(group))
        return records_by_tile

    def read_cached_tile(self, tile: Tile) -> list[dict] | None:
        """Read and validate one tile without using the network."""

        if not self.enabled:
            return None
        records = self._read_tile(tile)
        if records is None:
            self._increment("tile_miss_count")
        else:
            self._increment("tile_hit_count")
        return records

    def tile_path(self, tile: Tile) -> Path:
        """Return the durable JSON path for a fixed tile."""

        if not self.enabled:
            raise RuntimeError("Tile paths require an enabled gnomAD cache")
        return self._tile_path(tile)

    def snapshot(self) -> dict[str, object]:
        with self._stats_lock:
            stats = dict(self._stats)
        return {
            "enabled": self.enabled,
            "directory": str(self.namespace_dir) if self.namespace_dir else "",
            "schema_version": CACHE_SCHEMA_VERSION,
            "dataset": GNOMAD_DATASET,
            "reference_genome": GNOMAD_REFERENCE_GENOME,
            "tile_size_bp": self.tile_size_bp,
            **stats,
        }

    def _fetch_and_store(self, tiles: list[Tile]) -> dict[Tile, list[dict]]:
        if not tiles:
            return {}
        first = tiles[0]
        last = tiles[-1]
        attempts = self.max_attempts if len(tiles) == 1 else self.group_attempts
        try:
            records = self._call_fetcher(first.chrom, first.start, last.end, attempts)
        except Exception as exc:
            if len(tiles) > 1 and (
                _is_split_worthy(exc) or _is_retryable_fetch_error(exc)
            ):
                self._increment("split_count")
                midpoint = len(tiles) // 2
                return {
                    **self._fetch_and_store(tiles[:midpoint]),
                    **self._fetch_and_store(tiles[midpoint:]),
                }
            raise

        records_by_tile = _partition_records(records, tiles)
        for tile in tiles:
            self._write_tile(tile, records_by_tile[tile])
        return records_by_tile

    def _call_fetcher(
        self,
        chrom: str,
        start: int,
        end: int,
        max_attempts: int,
    ) -> list[dict]:
        self._increment("fetch_batch_count")
        return self.fetcher(chrom, start, end, max_attempts=max_attempts)

    def _tile_path(self, tile: Tile) -> Path:
        assert self.namespace_dir is not None
        return self.namespace_dir / tile.chrom / f"{tile.start:012d}-{tile.end:012d}.json.gz"

    def _read_tile(self, tile: Tile) -> list[dict] | None:
        path = self._tile_path(tile)
        if not path.exists():
            return None
        try:
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                payload = json.load(handle)
            self._validate_payload(payload, tile)
            return payload["variants"]
        except (EOFError, OSError, UnicodeDecodeError, json.JSONDecodeError, TypeError, ValueError) as exc:
            self._increment("corrupt_tile_count")
            logger.warning("Ignoring invalid gnomAD cache tile %s: %s", path, exc)
            return None

    def _validate_payload(self, payload: object, tile: Tile) -> None:
        if not isinstance(payload, dict):
            raise ValueError("cache payload is not an object")
        expected = {
            "schema_version": CACHE_SCHEMA_VERSION,
            "dataset": GNOMAD_DATASET,
            "reference_genome": GNOMAD_REFERENCE_GENOME,
            "tile_size_bp": self.tile_size_bp,
            "chrom": tile.chrom,
            "start": tile.start,
            "end": tile.end,
        }
        for key, value in expected.items():
            if payload.get(key) != value:
                raise ValueError(f"unexpected {key}: {payload.get(key)!r}")
        variants = payload.get("variants")
        if not isinstance(variants, list):
            raise ValueError("variants is not a list")
        for record in variants:
            if not isinstance(record, dict):
                raise ValueError("variant record is not an object")
            if normalize_chrom(record.get("chrom")) != tile.chrom:
                raise ValueError("variant chromosome is outside the tile")
            position = int(record.get("pos", 0))
            if not tile.start <= position <= tile.end:
                raise ValueError("variant position is outside the tile")

    def _write_tile(self, tile: Tile, records: list[dict]) -> None:
        path = self._tile_path(tile)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema_version": CACHE_SCHEMA_VERSION,
            "dataset": GNOMAD_DATASET,
            "reference_genome": GNOMAD_REFERENCE_GENOME,
            "tile_size_bp": self.tile_size_bp,
            "chrom": tile.chrom,
            "start": tile.start,
            "end": tile.end,
            "variants": records,
        }
        temporary: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(
                dir=path.parent,
                prefix=f".{path.name}.",
                suffix=".tmp",
                delete=False,
            ) as handle:
                temporary = Path(handle.name)
            with gzip.open(temporary, "wt", encoding="utf-8", compresslevel=6) as handle:
                json.dump(payload, handle, separators=(",", ":"), sort_keys=True)
            os.replace(temporary, path)
            self._increment("tile_write_count")
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    def _increment(self, key: str) -> None:
        with self._stats_lock:
            self._stats[key] += 1


def _consecutive_groups(tiles: list[Tile]) -> list[list[Tile]]:
    groups: list[list[Tile]] = []
    for tile in tiles:
        if (
            not groups
            or groups[-1][-1].chrom != tile.chrom
            or groups[-1][-1].end + 1 != tile.start
        ):
            groups.append([tile])
        else:
            groups[-1].append(tile)
    return groups


def _partition_records(records: list[dict], tiles: list[Tile]) -> dict[Tile, list[dict]]:
    records_by_tile = {tile: {} for tile in tiles}
    tile_by_start = {tile.start: tile for tile in tiles}
    tile_size_bp = tiles[0].end - tiles[0].start + 1
    first_start = tiles[0].start
    chrom = tiles[0].chrom
    for record in records:
        record_chrom = normalize_chrom(record.get("chrom"))
        position = int(record.get("pos", 0))
        if record_chrom != chrom or position < first_start:
            continue
        tile_start = first_start + ((position - first_start) // tile_size_bp) * tile_size_bp
        tile = tile_by_start.get(tile_start)
        if tile is None or position > tile.end:
            continue
        key = (
            record_chrom,
            position,
            str(record.get("ref", "")),
            str(record.get("alt", "")),
        )
        records_by_tile[tile][key] = record
    return {
        tile: sorted(tile_records.values(), key=_record_sort_key)
        for tile, tile_records in records_by_tile.items()
    }


def _record_sort_key(record: dict) -> tuple[int, str, str, str]:
    return (
        int(record.get("pos", 0)),
        str(record.get("ref", "")),
        str(record.get("alt", "")),
        str(record.get("variant_id", "")),
    )


def _is_split_worthy(exc: Exception) -> bool:
    if isinstance(exc, TimeoutError):
        return True
    if isinstance(exc, HTTPError):
        return exc.code in {408, 504}
    if isinstance(exc, URLError):
        reason = exc.reason
        if isinstance(reason, TimeoutError):
            return True
    message = str(exc).lower()
    return any(
        marker in message
        for marker in (
            "timed out",
            "timeout",
            "unexpected eof",
            "stream error",
            "response ended prematurely",
        )
    )


def _is_retryable_fetch_error(exc: Exception) -> bool:
    return is_retryable_network_error(exc) or "rate limit" in str(exc).lower()
