"""Derived Parquet index for exact-allele gnomAD lookups."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import tempfile
from collections import defaultdict
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import duckdb
import pandas as pd

from .gnomad import GNOMAD_DATASET, select_af_metrics
from .gnomad_cache import (
    CACHE_SCHEMA_VERSION as REGION_CACHE_SCHEMA_VERSION,
    GNOMAD_REFERENCE_GENOME,
    GnomadRegionCache,
    Tile,
)
from .variants import normalize_chrom, parse_variant_key


INDEX_SCHEMA_VERSION = 1
AF_POLICY_VERSION = 1
FRAGMENT_TILE_COUNT = 64
INDEX_COLUMNS = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "tile_start",
    "gnomad_af",
    "gnomad_af_source",
    "af_exome",
    "af_genome",
    "af_joint",
    "an_joint",
    "ac_joint",
    "consequence",
]
EVIDENCE_COLUMNS = [
    "variant_key",
    "gnomad_status",
    "gnomad_found",
    "gnomad_af",
]
REQUEST_COLUMNS = ["variant_key", "chrom", "pos", "ref", "alt"]


class GnomadAlleleIndex:
    """Build and query immutable allele-index fragments from regional JSON."""

    def __init__(
        self,
        cache_dir: Path | str,
        *,
        region_cache: GnomadRegionCache,
    ) -> None:
        if not region_cache.enabled:
            raise ValueError("gnomAD allele indexing requires an enabled regional cache")
        self.cache_dir = Path(cache_dir).expanduser().resolve()
        self.region_cache = region_cache
        self.tile_size_bp = region_cache.tile_size_bp
        self.namespace_dir = (
            self.cache_dir
            / GNOMAD_DATASET
            / GNOMAD_REFERENCE_GENOME
            / f"allele_index_v{INDEX_SCHEMA_VERSION}"
            / f"region_schema_v{REGION_CACHE_SCHEMA_VERSION}"
            / f"tiles_{self.tile_size_bp}bp"
        )
        self.metadata_path = self.namespace_dir / "metadata.json"

    def lookup(
        self,
        requests: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
        """Resolve requests covered by complete index tiles.

        Requests whose regional JSON tile is not locally available remain
        unresolved so the caller can use the existing network/retry path.
        """

        wanted = _normalize_requests(requests, self.tile_size_bp)
        if wanted.empty:
            return (
                _empty_evidence(),
                wanted[REQUEST_COLUMNS],
                self._summary(0, 0, 0, 0, 0, 0, 0),
            )

        tiles = _request_tiles(wanted, self.tile_size_bp)
        existing_coverage = self._coverage(tiles)
        preparation = self.prepare_tiles(tiles)
        coverage = self._coverage(tiles)
        covered_tiles = set(coverage)
        covered_mask = [
            Tile(str(row.chrom), int(row.tile_start), int(row.tile_end))
            in covered_tiles
            for row in wanted.itertuples(index=False)
        ]
        covered = wanted.loc[covered_mask].copy()
        unresolved = wanted.loc[
            [not value for value in covered_mask], REQUEST_COLUMNS
        ].copy()
        if covered.empty:
            return (
                _empty_evidence(),
                unresolved.reset_index(drop=True),
                self._summary(
                    len(wanted),
                    0,
                    0,
                    len(tiles),
                    len(existing_coverage),
                    int(preparation["tile_build_count"]),
                    0,
                    preparation=preparation,
                ),
            )

        parquet_files = sorted({str(coverage[tile]) for tile in covered_tiles})
        try:
            with duckdb.connect() as connection:
                connection.register("gnomad_requests", covered[REQUEST_COLUMNS])
                connection.read_parquet(parquet_files).create_view(
                    "gnomad_index_rows"
                )
                evidence = connection.execute(
                    """
                    SELECT
                        r.variant_key,
                        'ok' AS gnomad_status,
                        c.pos IS NOT NULL AS gnomad_found,
                        c.gnomad_af
                    FROM gnomad_requests AS r
                    LEFT JOIN gnomad_index_rows AS c
                      ON c.chrom = r.chrom
                     AND c.pos = r.pos
                     AND c.ref = r.ref
                     AND c.alt = r.alt
                    ORDER BY r.chrom, r.pos, r.variant_key
                    """
                ).df()
        except Exception as exc:
            raise ValueError(
                f"Unable to read gnomAD allele-index under {self.namespace_dir}"
            ) from exc

        if len(evidence) != len(covered) or evidence["variant_key"].duplicated().any():
            raise ValueError("gnomAD allele index returned duplicate or missing request rows")
        evidence = evidence[EVIDENCE_COLUMNS]
        evidence["gnomad_found"] = evidence["gnomad_found"].astype(bool)
        evidence["gnomad_af"] = pd.to_numeric(
            evidence["gnomad_af"],
            errors="coerce",
        )
        return (
            evidence.reset_index(drop=True),
            unresolved.reset_index(drop=True),
            self._summary(
                len(wanted),
                len(evidence),
                int(evidence["gnomad_found"].sum()),
                len(tiles),
                len(existing_coverage),
                int(preparation["tile_build_count"]),
                len(parquet_files),
                preparation=preparation,
            ),
        )

    def prepare(self, requests: pd.DataFrame) -> dict[str, object]:
        """Build any missing index tiles whose regional JSON is available."""

        wanted = _normalize_requests(requests, self.tile_size_bp)
        return self.prepare_tiles(_request_tiles(wanted, self.tile_size_bp))

    def prepare_tiles(self, tiles: list[Tile]) -> dict[str, object]:
        unique_tiles = sorted(set(tiles), key=_tile_sort_key)
        covered = self._coverage(unique_tiles)
        candidates = [tile for tile in unique_tiles if tile not in covered]
        if not candidates:
            return {
                "requested_tile_count": len(unique_tiles),
                "tile_build_count": 0,
                "raw_tile_missing_count": 0,
                "indexed_variant_count": 0,
                "fragment_build_count": 0,
            }

        self._ensure_metadata()
        candidates_by_chrom: dict[str, list[Tile]] = defaultdict(list)
        for tile in candidates:
            candidates_by_chrom[tile.chrom].append(tile)

        built_tiles = 0
        raw_missing = 0
        indexed_variants = 0
        built_fragments = 0
        with duckdb.connect() as connection:
            for chrom in sorted(candidates_by_chrom):
                chrom_candidates = candidates_by_chrom[chrom]
                with self._chrom_publish_lock(chrom):
                    current = self._coverage(chrom_candidates)
                    pending = [tile for tile in chrom_candidates if tile not in current]
                    for tile_batch in _batches(pending, FRAGMENT_TILE_COUNT):
                        frames = []
                        manifest_tiles = []
                        for tile in tile_batch:
                            raw_path = self.region_cache.tile_path(tile)
                            if not raw_path.exists():
                                raw_missing += 1
                                continue
                            records = self.region_cache.read_cached_tile(tile)
                            if records is None:
                                raw_missing += 1
                                continue
                            frames.append(_index_frame(records, tile))
                            manifest_tiles.append(
                                {
                                    "start": tile.start,
                                    "end": tile.end,
                                    "row_count": len(records),
                                }
                            )
                            built_tiles += 1
                            indexed_variants += len(records)
                        if not manifest_tiles:
                            continue
                        fragment = pd.concat(frames, ignore_index=True)
                        self._write_fragment(
                            connection,
                            chrom,
                            fragment,
                            manifest_tiles,
                        )
                        built_fragments += 1
        return {
            "requested_tile_count": len(unique_tiles),
            "tile_build_count": built_tiles,
            "raw_tile_missing_count": raw_missing,
            "indexed_variant_count": indexed_variants,
            "fragment_build_count": built_fragments,
        }

    def fragment_paths(self) -> list[Path]:
        return sorted(self.namespace_dir.glob("chrom=*/fragment-*.parquet"))

    def _chrom_dir(self, chrom: str) -> Path:
        if not chrom or chrom in {".", ".."} or any(char in chrom for char in "/\\\0"):
            raise ValueError(f"Unsafe gnomAD chromosome: {chrom!r}")
        return self.namespace_dir / f"chrom={chrom}"

    @contextmanager
    def _chrom_publish_lock(self, chrom: str) -> Iterator[None]:
        directory = self._chrom_dir(chrom)
        directory.mkdir(parents=True, exist_ok=True)
        with (directory / ".publish.lock").open("a+b") as handle:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    def _write_fragment(
        self,
        connection,
        chrom: str,
        frame: pd.DataFrame,
        tiles: list[dict[str, int]],
    ) -> None:
        directory = self._chrom_dir(chrom)
        parquet_temporary: Path | None = None
        manifest_temporary: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(
                dir=directory,
                prefix=".fragment.",
                suffix=".parquet.tmp",
                delete=False,
            ) as handle:
                parquet_temporary = Path(handle.name)
            connection.register("gnomad_fragment_rows", frame)
            try:
                connection.table("gnomad_fragment_rows").write_parquet(
                    str(parquet_temporary),
                    compression="zstd",
                )
            finally:
                connection.unregister("gnomad_fragment_rows")
            relation = connection.read_parquet(str(parquet_temporary))
            if relation.columns != INDEX_COLUMNS:
                raise ValueError(
                    "gnomAD allele-index columns changed: "
                    + ", ".join(relation.columns)
                )
            if relation.count("*").fetchone()[0] != len(frame):
                raise ValueError("gnomAD allele-index row count changed while writing")

            parquet_sha256 = _file_sha256(parquet_temporary)
            identity = {
                "chrom": chrom,
                "tiles": tiles,
                "parquet_sha256": parquet_sha256,
            }
            digest = _json_sha256(identity)
            stem = f"fragment-{digest}"
            parquet_path = directory / f"{stem}.parquet"
            manifest_path = directory / f"{stem}.json"
            os.replace(parquet_temporary, parquet_path)
            parquet_temporary = None

            manifest = {
                "schema_version": INDEX_SCHEMA_VERSION,
                "chrom": chrom,
                "parquet": parquet_path.name,
                "parquet_sha256": parquet_sha256,
                "row_count": len(frame),
                "tiles": tiles,
            }
            with tempfile.NamedTemporaryFile(
                "w",
                dir=directory,
                prefix=f".{stem}.",
                suffix=".json.tmp",
                delete=False,
            ) as handle:
                manifest_temporary = Path(handle.name)
                json.dump(manifest, handle, indent=2, sort_keys=True)
                handle.write("\n")
            os.replace(manifest_temporary, manifest_path)
            manifest_temporary = None
        finally:
            if parquet_temporary is not None:
                parquet_temporary.unlink(missing_ok=True)
            if manifest_temporary is not None:
                manifest_temporary.unlink(missing_ok=True)

    def _coverage(self, requested_tiles: list[Tile]) -> dict[Tile, Path]:
        if not requested_tiles or not self.namespace_dir.exists():
            return {}
        self._validate_metadata_for_read()
        requested = set(requested_tiles)
        coverage: dict[Tile, Path] = {}
        for chrom in sorted({tile.chrom for tile in requested}):
            directory = self._chrom_dir(chrom)
            for manifest_path in sorted(directory.glob("fragment-*.json")):
                try:
                    manifest = json.loads(manifest_path.read_text())
                except (OSError, json.JSONDecodeError) as exc:
                    raise ValueError(
                        f"Invalid gnomAD allele-index manifest: {manifest_path}"
                    ) from exc
                if (
                    manifest.get("schema_version") != INDEX_SCHEMA_VERSION
                    or manifest.get("chrom") != chrom
                    or not isinstance(manifest.get("tiles"), list)
                ):
                    raise ValueError(
                        f"Invalid gnomAD allele-index manifest: {manifest_path}"
                    )
                parquet_name = str(manifest.get("parquet") or "")
                parquet_sha256 = str(manifest.get("parquet_sha256") or "")
                expected_digest = _json_sha256(
                    {
                        "chrom": chrom,
                        "tiles": manifest["tiles"],
                        "parquet_sha256": parquet_sha256,
                    }
                )
                if (
                    manifest_path.name != f"fragment-{expected_digest}.json"
                    or parquet_name != f"fragment-{expected_digest}.parquet"
                ):
                    raise ValueError(
                        f"gnomAD allele-index manifest identity changed: {manifest_path}"
                    )
                parquet_path = directory / parquet_name
                if not parquet_path.is_file():
                    raise ValueError(
                        f"gnomAD allele-index fragment is missing: {parquet_path}"
                    )
                manifest_row_count = 0
                for item in manifest["tiles"]:
                    try:
                        tile = Tile(chrom, int(item["start"]), int(item["end"]))
                        tile_row_count = int(item["row_count"])
                    except (KeyError, TypeError, ValueError) as exc:
                        raise ValueError(
                            f"Invalid gnomAD allele-index tile: {manifest_path}"
                        ) from exc
                    if tile.end - tile.start + 1 != self.tile_size_bp:
                        raise ValueError(
                            f"Invalid gnomAD allele-index tile size: {manifest_path}"
                        )
                    if tile_row_count < 0:
                        raise ValueError(
                            f"Invalid gnomAD allele-index row count: {manifest_path}"
                        )
                    manifest_row_count += tile_row_count
                    if tile not in requested:
                        continue
                    previous = coverage.setdefault(tile, parquet_path)
                    if previous != parquet_path:
                        raise ValueError(
                            f"Duplicate gnomAD allele-index coverage for {tile}"
                        )
                if manifest_row_count != int(manifest.get("row_count", -1)):
                    raise ValueError(
                        f"Invalid gnomAD allele-index row total: {manifest_path}"
                    )
        return coverage

    def _metadata(self) -> dict[str, object]:
        return {
            "schema_version": INDEX_SCHEMA_VERSION,
            "region_cache_schema_version": REGION_CACHE_SCHEMA_VERSION,
            "dataset": GNOMAD_DATASET,
            "reference_genome": GNOMAD_REFERENCE_GENOME,
            "tile_size_bp": self.tile_size_bp,
            "af_policy_version": AF_POLICY_VERSION,
            "columns": INDEX_COLUMNS,
            "format": "parquet",
            "compression": "zstd",
            "fragment_tile_count": FRAGMENT_TILE_COUNT,
            "layout": "immutable_tile_fragments",
        }

    def _ensure_metadata(self) -> None:
        expected = self._metadata()
        if self.metadata_path.exists():
            observed = json.loads(self.metadata_path.read_text())
            if observed != expected:
                raise ValueError(
                    f"gnomAD allele-index metadata changed: {self.metadata_path}"
                )
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
        if not self.metadata_path.exists():
            raise ValueError(
                f"gnomAD allele index contains fragments without metadata: {self.namespace_dir}"
            )
        observed = json.loads(self.metadata_path.read_text())
        if observed != self._metadata():
            raise ValueError(f"gnomAD allele-index metadata changed: {self.metadata_path}")

    def _summary(
        self,
        requested_count: int,
        resolved_count: int,
        found_count: int,
        requested_tile_count: int,
        tile_hit_count: int,
        tile_build_count: int,
        fragment_count: int,
        *,
        preparation: dict[str, object] | None = None,
    ) -> dict[str, object]:
        preparation = preparation or {}
        return {
            "enabled": True,
            "directory": str(self.namespace_dir),
            "schema_version": INDEX_SCHEMA_VERSION,
            "requested_count": requested_count,
            "resolved_count": resolved_count,
            "unresolved_count": requested_count - resolved_count,
            "found_count": found_count,
            "requested_tile_count": requested_tile_count,
            "tile_hit_count": tile_hit_count,
            "tile_build_count": tile_build_count,
            "fragment_count": fragment_count,
            "fragment_build_count": int(
                preparation.get("fragment_build_count", 0)
            ),
            "raw_tile_missing_count": int(
                preparation.get("raw_tile_missing_count", 0)
            ),
            "indexed_variant_count": int(
                preparation.get("indexed_variant_count", 0)
            ),
        }


def _normalize_requests(frame: pd.DataFrame, tile_size_bp: int) -> pd.DataFrame:
    missing = set(REQUEST_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(
            f"gnomAD allele-index requests missing columns: {', '.join(sorted(missing))}"
        )
    records = []
    coordinates_by_key: dict[str, tuple[str, int, str, str]] = {}
    for row in frame[REQUEST_COLUMNS].itertuples(index=False):
        supplied = (
            normalize_chrom(str(row.chrom)) or "",
            int(row.pos),
            str(row.ref).upper(),
            str(row.alt).upper(),
        )
        if parse_variant_key(row.variant_key) != supplied:
            raise ValueError(
                f"gnomAD allele-index coordinates do not match {row.variant_key!r}"
            )
        key = str(row.variant_key)
        previous = coordinates_by_key.setdefault(key, supplied)
        if previous != supplied:
            raise ValueError(f"Conflicting coordinates for gnomAD key {key}")
        tile_start = ((supplied[1] - 1) // tile_size_bp) * tile_size_bp + 1
        records.append(
            {
                "variant_key": key,
                "chrom": supplied[0],
                "pos": supplied[1],
                "ref": supplied[2],
                "alt": supplied[3],
                "tile_start": tile_start,
                "tile_end": tile_start + tile_size_bp - 1,
            }
        )
    if not records:
        return pd.DataFrame(columns=[*REQUEST_COLUMNS, "tile_start", "tile_end"])
    return (
        pd.DataFrame.from_records(records)
        .drop_duplicates("variant_key")
        .sort_values(["chrom", "pos", "variant_key"], kind="mergesort")
        .reset_index(drop=True)
    )


def _request_tiles(frame: pd.DataFrame, tile_size_bp: int) -> list[Tile]:
    return sorted(
        {
            Tile(
                str(row.chrom),
                ((int(row.pos) - 1) // tile_size_bp) * tile_size_bp + 1,
                ((int(row.pos) - 1) // tile_size_bp) * tile_size_bp + tile_size_bp,
            )
            for row in frame[["chrom", "pos"]].itertuples(index=False)
        },
        key=_tile_sort_key,
    )


def _index_frame(records: list[dict], tile: Tile) -> pd.DataFrame:
    rows = []
    for record in records:
        (
            af,
            source,
            af_exome,
            af_genome,
            af_joint,
            an_joint,
            ac_joint,
        ) = select_af_metrics(record)
        rows.append(
            {
                "chrom": normalize_chrom(record.get("chrom")) or "",
                "pos": int(record.get("pos", 0)),
                "ref": str(record.get("ref", "")).upper(),
                "alt": str(record.get("alt", "")).upper(),
                "tile_start": tile.start,
                "gnomad_af": af,
                "gnomad_af_source": source or "",
                "af_exome": af_exome,
                "af_genome": af_genome,
                "af_joint": af_joint,
                "an_joint": an_joint,
                "ac_joint": ac_joint,
                "consequence": str(record.get("consequence") or ""),
            }
        )
    frame = pd.DataFrame.from_records(rows, columns=INDEX_COLUMNS)
    text_columns = ["chrom", "ref", "alt", "gnomad_af_source", "consequence"]
    float_columns = ["gnomad_af", "af_exome", "af_genome", "af_joint"]
    integer_columns = ["pos", "tile_start", "an_joint", "ac_joint"]
    for column in text_columns:
        frame[column] = frame[column].astype("string")
    for column in float_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce").astype("Float64")
    for column in integer_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce").astype("Int64")
    return frame


def _empty_evidence() -> pd.DataFrame:
    return pd.DataFrame(columns=EVIDENCE_COLUMNS)


def _tile_sort_key(tile: Tile) -> tuple[str, int]:
    return tile.chrom, tile.start


def _batches(values: list[Tile], size: int) -> Iterator[list[Tile]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(payload: dict[str, object]) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
