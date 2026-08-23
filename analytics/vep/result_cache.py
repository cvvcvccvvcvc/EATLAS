"""Shared immutable cache for completed VEP variant/gene annotations."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import re
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import duckdb
import pandas as pd

from genomics.variants import normalize_chrom, parse_variant_key


CACHE_SCHEMA_VERSION = 1
DEFAULT_TILE_SIZE_BP = 1_000_000
KEY_COLUMNS = ["chrom", "pos", "ref", "alt", "gene_id"]
RESULT_COLUMNS = [
    "status",
    "primary_consequence",
    "consequence_terms",
    "transcript_id",
    "mane_select",
    "canonical",
    "impact",
    "variant_class",
]
ANNOTATION_COLUMNS = ["variant_key", "gene_id", *RESULT_COLUMNS]
STORAGE_COLUMNS = [*KEY_COLUMNS, *RESULT_COLUMNS]
CACHEABLE_STATUSES = {"ok", "no_target_gene", "no_consequence"}

_SAFE_COMPONENT = re.compile(r"^[A-Za-z0-9_.-]+$")


@dataclass(frozen=True, slots=True)
class CacheTile:
    chrom: str
    start: int
    end: int


class VepResultCache:
    """Read and append immutable regional Parquet fragments."""

    def __init__(
        self,
        cache_dir: Path | str,
        *,
        config: dict[str, Any],
        tile_size_bp: int = DEFAULT_TILE_SIZE_BP,
    ) -> None:
        if tile_size_bp < 1:
            raise ValueError("VEP result cache tile size must be >= 1")
        self.cache_dir = Path(cache_dir).expanduser().resolve()
        self.config = _normalized_config(config)
        self.config_hash = _json_hash(self.config)
        self.tile_size_bp = int(tile_size_bp)
        release = _safe_component(
            str(self.config.get("release") or "unknown"),
            "release",
        )
        backend = _safe_component(
            str(self.config.get("backend") or "rest"),
            "backend",
        )
        self.namespace_dir = (
            self.cache_dir
            / f"schema_v{CACHE_SCHEMA_VERSION}"
            / f"release={release}"
            / f"backend={backend}"
            / f"config={self.config_hash}"
            / f"tiles_{self.tile_size_bp}bp"
        )
        self.metadata_path = self.namespace_dir / "metadata.json"

    def lookup(self, requests: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
        """Return cached annotations for unique, validated requests."""

        wanted = _normalize_requests(requests)
        if wanted.empty:
            return _empty_annotations(), self._lookup_summary(0, 0, 0, 0)

        self._validate_metadata_for_read()
        files_by_tile = {
            tile: sorted(self._tile_dir(tile).glob("part-*.parquet"))
            for tile in _request_tiles(wanted, self.tile_size_bp)
        }
        parquet_files = [
            str(path)
            for paths in files_by_tile.values()
            for path in paths
        ]
        if not parquet_files:
            return (
                _empty_annotations(),
                self._lookup_summary(len(wanted), 0, len(files_by_tile), 0),
            )

        with duckdb.connect() as connection:
            connection.register("vep_requests", wanted)
            connection.read_parquet(parquet_files).create_view("vep_cache_rows")
            selected_results = ", ".join(f"c.{column}" for column in RESULT_COLUMNS)
            hits = connection.execute(
                f"""
                SELECT
                    r.variant_key,
                    r.gene_id,
                    {selected_results}
                FROM vep_requests AS r
                JOIN vep_cache_rows AS c
                  ON c.chrom = r.chrom
                 AND c.pos = r.pos
                 AND c.ref = r.ref
                 AND c.alt = r.alt
                 AND c.gene_id = r.gene_id
                """
            ).df()

        hits = _deduplicate_annotations(hits, context="shared VEP cache")
        return (
            hits,
            self._lookup_summary(
                len(wanted),
                len(hits),
                len(files_by_tile),
                len(parquet_files),
            ),
        )

    def publish(self, annotations: pd.DataFrame) -> dict[str, object]:
        """Append cache misses as validated, content-addressed fragments."""

        rows, skipped_count = _normalize_annotations(
            annotations,
            tile_size_bp=self.tile_size_bp,
        )
        if rows.empty:
            return self._publish_summary(0, skipped_count, 0, 0, 0)

        requests = rows[
            ["variant_key", "gene_id", "chrom", "pos", "ref", "alt"]
        ]
        existing, _ = self.lookup(requests)
        missing_records, existing_count = _partition_cached_annotations(
            rows,
            existing,
        )

        if not missing_records:
            return self._publish_summary(
                len(rows),
                skipped_count,
                existing_count,
                0,
                0,
            )

        missing = pd.DataFrame.from_records(missing_records)
        self._ensure_metadata()
        fragment_count = 0
        published_count = 0
        for (chrom, tile_start), group in missing.groupby(
            ["chrom", "tile_start"],
            sort=True,
        ):
            tile = CacheTile(
                chrom=str(chrom),
                start=int(tile_start),
                end=int(tile_start) + self.tile_size_bp - 1,
            )
            with self._tile_publish_lock(tile):
                current, _ = self.lookup(
                    group[["variant_key", "gene_id", "chrom", "pos", "ref", "alt"]]
                )
                pending_records, concurrent_existing = _partition_cached_annotations(
                    group,
                    current,
                )
                existing_count += concurrent_existing
                if not pending_records:
                    continue
                pending = pd.DataFrame.from_records(pending_records)
                storage = pending[STORAGE_COLUMNS].sort_values(
                    KEY_COLUMNS,
                    kind="mergesort",
                )
                created = self._write_fragment(tile, storage)
                fragment_count += int(created)
                if created:
                    published_count += len(storage)

        return self._publish_summary(
            len(rows),
            skipped_count,
            existing_count,
            published_count,
            fragment_count,
        )

    def _tile_dir(self, tile: CacheTile) -> Path:
        chrom = _safe_component(tile.chrom, "chromosome")
        return (
            self.namespace_dir
            / f"chrom={chrom}"
            / f"tile={tile.start:012d}-{tile.end:012d}"
        )

    @contextmanager
    def _tile_publish_lock(self, tile: CacheTile) -> Iterator[None]:
        directory = self._tile_dir(tile)
        directory.mkdir(parents=True, exist_ok=True)
        lock_path = directory / ".publish.lock"
        with lock_path.open("a+b") as handle:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    def _write_fragment(self, tile: CacheTile, rows: pd.DataFrame) -> bool:
        directory = self._tile_dir(tile)
        directory.mkdir(parents=True, exist_ok=True)
        temporary: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(
                dir=directory,
                prefix=".vep-cache.",
                suffix=".parquet.tmp",
                delete=False,
            ) as handle:
                temporary = Path(handle.name)

            with duckdb.connect() as connection:
                connection.register("cache_fragment", rows)
                connection.table("cache_fragment").write_parquet(
                    str(temporary),
                    compression="zstd",
                )
                relation = connection.read_parquet(str(temporary))
                if relation.columns != STORAGE_COLUMNS:
                    raise ValueError(
                        "VEP cache fragment columns changed: "
                        f"{', '.join(relation.columns)}"
                    )
                if relation.count("*").fetchone()[0] != len(rows):
                    raise ValueError("VEP cache fragment row count changed while writing")

            digest = _file_sha256(temporary)
            destination = directory / f"part-{digest}.parquet"
            if destination.exists():
                if _file_sha256(destination) != digest:
                    raise ValueError(f"VEP cache fragment hash mismatch: {destination}")
                temporary.unlink()
                temporary = None
                return False
            os.replace(temporary, destination)
            temporary = None
            return True
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    def _metadata(self) -> dict[str, object]:
        return {
            "schema_version": CACHE_SCHEMA_VERSION,
            "config": self.config,
            "config_hash": self.config_hash,
            "tile_size_bp": self.tile_size_bp,
            "key_columns": KEY_COLUMNS,
            "result_columns": RESULT_COLUMNS,
            "format": "parquet",
            "compression": "zstd",
            "layout": "immutable_fragments",
        }

    def _ensure_metadata(self) -> None:
        expected = self._metadata()
        if self.metadata_path.exists():
            observed = json.loads(self.metadata_path.read_text())
            if observed != expected:
                raise ValueError(f"VEP result cache metadata changed: {self.metadata_path}")
            return

        self.namespace_dir.mkdir(parents=True, exist_ok=True)
        temporary: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(
                "w",
                dir=self.namespace_dir,
                prefix=".metadata.",
                suffix=".json.tmp",
                delete=False,
            ) as handle:
                temporary = Path(handle.name)
                json.dump(expected, handle, indent=2, sort_keys=True)
                handle.write("\n")
            os.replace(temporary, self.metadata_path)
            temporary = None
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    def _validate_metadata_for_read(self) -> None:
        if not self.namespace_dir.exists():
            return
        if not self.metadata_path.exists():
            if any(self.namespace_dir.rglob("part-*.parquet")):
                raise ValueError(
                    f"VEP result cache contains data without metadata: {self.namespace_dir}"
                )
            return
        observed = json.loads(self.metadata_path.read_text())
        if observed != self._metadata():
            raise ValueError(f"VEP result cache metadata changed: {self.metadata_path}")

    def _lookup_summary(
        self,
        requested_count: int,
        hit_count: int,
        tile_count: int,
        fragment_count: int,
    ) -> dict[str, object]:
        return {
            "enabled": True,
            "directory": str(self.namespace_dir),
            "config_hash": self.config_hash,
            "tile_size_bp": self.tile_size_bp,
            "requested_count": requested_count,
            "hit_count": hit_count,
            "miss_count": requested_count - hit_count,
            "tile_count": tile_count,
            "fragment_count": fragment_count,
        }

    def _publish_summary(
        self,
        accepted_count: int,
        skipped_count: int,
        existing_count: int,
        published_count: int,
        fragment_count: int,
    ) -> dict[str, object]:
        return {
            "directory": str(self.namespace_dir),
            "config_hash": self.config_hash,
            "tile_size_bp": self.tile_size_bp,
            "accepted_count": accepted_count,
            "skipped_count": skipped_count,
            "existing_count": existing_count,
            "published_count": published_count,
            "fragment_count": fragment_count,
        }


def _partition_cached_annotations(
    rows: pd.DataFrame,
    cached: pd.DataFrame,
) -> tuple[list[dict[str, object]], int]:
    cached_by_key = {
        (str(row.variant_key), str(row.gene_id)): tuple(
            getattr(row, column) for column in RESULT_COLUMNS
        )
        for row in cached.itertuples(index=False)
    }
    missing_records = []
    existing_count = 0
    for row in rows.itertuples(index=False):
        key = (str(row.variant_key), str(row.gene_id))
        values = tuple(getattr(row, column) for column in RESULT_COLUMNS)
        cached_values = cached_by_key.get(key)
        if cached_values is None:
            missing_records.append(row._asdict())
        elif cached_values != values:
            raise ValueError(
                "Conflicting VEP result for cached key "
                f"{row.variant_key} / gene {row.gene_id}"
            )
        else:
            existing_count += 1
    return missing_records, existing_count


def _normalize_requests(frame: pd.DataFrame) -> pd.DataFrame:
    required = {"variant_key", "gene_id", "chrom", "pos", "ref", "alt"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(
            f"VEP cache requests missing columns: {', '.join(sorted(missing))}"
        )

    records = []
    coordinates_by_key: dict[tuple[str, str], tuple[str, int, str, str]] = {}
    for row in frame[
        ["variant_key", "gene_id", "chrom", "pos", "ref", "alt"]
    ].itertuples(index=False):
        parsed = parse_variant_key(row.variant_key)
        supplied = (
            normalize_chrom(str(row.chrom)),
            int(row.pos),
            str(row.ref).upper(),
            str(row.alt).upper(),
        )
        if parsed is None or parsed != supplied:
            raise ValueError(
                f"VEP cache request coordinates do not match {row.variant_key!r}"
            )
        key = (str(row.variant_key), str(row.gene_id))
        previous = coordinates_by_key.setdefault(key, supplied)
        if previous != supplied:
            raise ValueError(
                f"VEP cache request has conflicting coordinates: {row.variant_key}"
            )
        records.append(
            {
                "variant_key": key[0],
                "gene_id": key[1],
                "chrom": supplied[0],
                "pos": supplied[1],
                "ref": supplied[2],
                "alt": supplied[3],
            }
        )

    if not records:
        return pd.DataFrame(
            columns=["variant_key", "gene_id", "chrom", "pos", "ref", "alt"]
        )
    return (
        pd.DataFrame.from_records(records)
        .drop_duplicates(["variant_key", "gene_id"])
        .sort_values(["chrom", "pos", "variant_key", "gene_id"], kind="mergesort")
        .reset_index(drop=True)
    )


def _normalize_annotations(
    frame: pd.DataFrame,
    *,
    tile_size_bp: int,
) -> tuple[pd.DataFrame, int]:
    missing = set(ANNOTATION_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(
            f"VEP cache annotations missing columns: {', '.join(sorted(missing))}"
        )

    records = []
    skipped = 0
    for row in frame[ANNOTATION_COLUMNS].itertuples(index=False):
        status = _as_text(row.status)
        if status not in CACHEABLE_STATUSES:
            skipped += 1
            continue
        parsed = parse_variant_key(row.variant_key)
        if parsed is None:
            raise ValueError(f"Invalid VEP cache variant key: {row.variant_key!r}")
        chrom, pos, ref, alt = parsed
        records.append(
            {
                "variant_key": str(row.variant_key),
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene_id": str(row.gene_id),
                "status": status,
                "primary_consequence": _as_text(row.primary_consequence),
                "consequence_terms": _as_text(row.consequence_terms),
                "transcript_id": _as_text(row.transcript_id),
                "mane_select": _as_text(row.mane_select),
                "canonical": _as_bool(row.canonical),
                "impact": _as_text(row.impact),
                "variant_class": _as_text(row.variant_class),
                "tile_start": ((pos - 1) // tile_size_bp) * tile_size_bp + 1,
            }
        )
    if not records:
        return (
            pd.DataFrame(columns=["variant_key", *STORAGE_COLUMNS, "tile_start"]),
            skipped,
        )

    normalized = pd.DataFrame.from_records(records)
    normalized = _deduplicate_storage(normalized)
    return normalized, skipped


def _request_tiles(frame: pd.DataFrame, tile_size_bp: int) -> list[CacheTile]:
    tiles = {
        CacheTile(
            chrom=str(row.chrom),
            start=((int(row.pos) - 1) // tile_size_bp) * tile_size_bp + 1,
            end=((int(row.pos) - 1) // tile_size_bp) * tile_size_bp + tile_size_bp,
        )
        for row in frame[["chrom", "pos"]].itertuples(index=False)
    }
    return sorted(tiles, key=lambda tile: (tile.chrom, tile.start))


def _deduplicate_storage(frame: pd.DataFrame) -> pd.DataFrame:
    unique = frame.drop_duplicates(["variant_key", *STORAGE_COLUMNS, "tile_start"])
    conflicts = unique.groupby(KEY_COLUMNS, sort=False, dropna=False).size()
    if not conflicts.empty and int(conflicts.max()) > 1:
        raise ValueError("Conflicting VEP annotations for the same cache key")
    return (
        unique.drop_duplicates(KEY_COLUMNS)
        .sort_values(KEY_COLUMNS, kind="mergesort")
        .reset_index(drop=True)
    )


def _deduplicate_annotations(frame: pd.DataFrame, *, context: str) -> pd.DataFrame:
    if frame.empty:
        return _empty_annotations()
    frame = frame[ANNOTATION_COLUMNS].copy()
    frame["canonical"] = frame["canonical"].map(_as_bool)
    unique = frame.drop_duplicates(ANNOTATION_COLUMNS)
    conflicts = unique.groupby(["variant_key", "gene_id"], sort=False).size()
    if not conflicts.empty and int(conflicts.max()) > 1:
        raise ValueError(f"Conflicting annotations found in {context}")
    return (
        unique.drop_duplicates(["variant_key", "gene_id"])
        .sort_values(["variant_key", "gene_id"], kind="mergesort")
        .reset_index(drop=True)
    )


def _empty_annotations() -> pd.DataFrame:
    return pd.DataFrame(columns=ANNOTATION_COLUMNS)


def _normalized_config(config: dict[str, Any]) -> dict[str, Any]:
    try:
        return json.loads(json.dumps(config, sort_keys=True))
    except (TypeError, ValueError) as exc:
        raise ValueError("VEP cache configuration must be JSON serializable") from exc


def _safe_component(value: str, label: str) -> str:
    if not value or value in {".", ".."} or not _SAFE_COMPONENT.fullmatch(value):
        raise ValueError(f"Unsafe VEP cache {label}: {value!r}")
    return value


def _as_bool(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    if value in {0, 1}:
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "yes", "1"}:
        return True
    if text in {"false", "no", "0", ""}:
        return False
    raise ValueError(f"Invalid VEP canonical value: {value!r}")


def _as_text(value: object) -> str:
    return "" if pd.isna(value) else str(value)


def _json_hash(payload: dict[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
