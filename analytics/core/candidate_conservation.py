"""Candidate-wide phyloP distributions for gnomAD-stratified reporting."""

from __future__ import annotations

import json
import sqlite3
import tempfile
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd

from .clinvar_validation import split_strategies
from .conservation import (
    DEFAULT_TRACK_NAMES,
    PositionScores,
    format_chrom,
    parse_tracks,
    path_metadata,
    read_position_scores,
    score_positions,
)


CACHE_VERSION = 1
QUANTILES = np.linspace(0.0, 1.0, 101)
REQUIRED_COLUMNS = {
    "variant_key",
    "lookup_chrom",
    "lookup_pos",
    "lookup_ref",
    "lookup_alt",
    "lookup_status",
    "strategies",
    "gnomad_af",
}


@dataclass(frozen=True)
class CandidateConservation:
    distributions_path: Path
    manifest_path: Path
    distributions: pd.DataFrame
    manifest: dict
    position_scores: PositionScores | None = None

    def without_position_scores(self) -> "CandidateConservation":
        return replace(self, position_scores=None)


def build_candidate_conservation(
    *,
    variant_annotations_tsv: Path,
    analytics_dir: Path,
    additional_rows: list[dict[str, str]] | None = None,
    track_names: str = DEFAULT_TRACK_NAMES,
    max_block_bp: int = 250_000,
    max_gap_bp: int = 50_000,
    remote_retries: int = 3,
    retry_sleep_seconds: float = 5.0,
    precision: int = 6,
    chunk_size: int = 100_000,
) -> CandidateConservation:
    """Compute compact exact percentile curves without persisting allele-level scores."""
    tracks = parse_tracks(track_names)
    if len(tracks) != 1:
        raise ValueError("Candidate-wide conservation currently requires exactly one track.")
    track = tracks[0]
    analytics_dir.mkdir(parents=True, exist_ok=True)
    distributions_path = analytics_dir / "candidate_variants.phyloP100way.distributions.tsv.gz"
    manifest_path = analytics_dir / "candidate_variants.phyloP100way.manifest.json"
    expected_inputs = {
        "cache_version": CACHE_VERSION,
        "variant_annotations": path_metadata(variant_annotations_tsv),
        "track": asdict(track),
        "max_block_bp": max_block_bp,
        "max_gap_bp": max_gap_bp,
        "remote_retries": remote_retries,
        "retry_sleep_seconds": retry_sleep_seconds,
        "precision": precision,
    }
    cached = _load_cache(distributions_path, manifest_path, expected_inputs)
    if cached is not None:
        return cached

    positions_by_chrom, scan_summary = _candidate_positions(
        variant_annotations_tsv,
        track.chrom_style,
        chunk_size,
    )
    _add_positions(positions_by_chrom, additional_rows or [], track.chrom_style)
    position_scores = read_position_scores(
        positions_by_chrom=positions_by_chrom,
        track=track,
        max_block_bp=max_block_bp,
        max_gap_bp=max_gap_bp,
        remote_retries=remote_retries,
        retry_sleep_seconds=retry_sleep_seconds,
        precision=precision,
    )
    distributions, groups, membership_summary = _aggregate_distributions(
        variant_annotations_tsv,
        position_scores,
        analytics_dir,
        chunk_size,
    )
    _write_frame(distributions_path, distributions)
    manifest = {
        "inputs": expected_inputs,
        "complete": position_scores.summary.get("status") == "complete",
        "candidate_scan": scan_summary,
        "position_read": position_scores.summary,
        "memberships": membership_summary,
        "groups": groups,
        "quantile_count": len(QUANTILES),
        "distributions_tsv": str(distributions_path),
    }
    _write_json(manifest_path, manifest)
    return CandidateConservation(
        distributions_path,
        manifest_path,
        distributions,
        manifest,
        position_scores,
    )


def _load_cache(
    distributions_path: Path,
    manifest_path: Path,
    expected_inputs: dict,
) -> CandidateConservation | None:
    if not distributions_path.exists() or not manifest_path.exists():
        return None
    try:
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("inputs") != expected_inputs or manifest.get("complete") is not True:
            return None
        distributions = pd.read_csv(distributions_path, sep="\t", compression="gzip")
        return CandidateConservation(distributions_path, manifest_path, distributions, manifest)
    except (OSError, json.JSONDecodeError, ValueError):
        return None


def _candidate_positions(
    path: Path,
    chrom_style: str,
    chunk_size: int,
) -> tuple[dict[str, set[int]], dict[str, int]]:
    header = pd.read_csv(path, sep="\t", compression="gzip", nrows=0).columns.tolist()
    missing = REQUIRED_COLUMNS - set(header)
    if missing:
        raise ValueError(f"Variant annotations missing conservation columns: {', '.join(sorted(missing))}")
    positions_by_chrom: dict[str, set[int]] = {}
    row_count = 0
    usable_allele_count = 0
    unsupported_allele_count = 0
    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        keep_default_na=False,
        usecols=["lookup_chrom", "lookup_pos", "lookup_ref", "lookup_alt", "lookup_status"],
        chunksize=chunk_size,
    ):
        row_count += len(chunk)
        chunk = chunk[chunk["lookup_status"].astype(str).eq("ok")]
        for chrom, pos, ref, alt in chunk[
            ["lookup_chrom", "lookup_pos", "lookup_ref", "lookup_alt"]
        ].itertuples(index=False, name=None):
            positions, _basis = score_positions(int(pos), str(ref), str(alt))
            positions = [position for position in positions if position >= 0]
            if not positions:
                unsupported_allele_count += 1
                continue
            usable_allele_count += 1
            formatted_chrom = format_chrom(chrom, chrom_style)
            positions_by_chrom.setdefault(formatted_chrom, set()).update(positions)
    return positions_by_chrom, {
        "variant_context_row_count": row_count,
        "usable_allele_context_count": usable_allele_count,
        "unsupported_allele_context_count": unsupported_allele_count,
        "candidate_unique_position_count": sum(len(values) for values in positions_by_chrom.values()),
    }


def _add_positions(
    positions_by_chrom: dict[str, set[int]],
    rows: list[dict[str, str]],
    chrom_style: str,
) -> None:
    for row in rows:
        chrom = format_chrom(row.get("chrom", ""), chrom_style)
        positions, _basis = score_positions(int(row["pos"]), row.get("ref", ""), row.get("alt", ""))
        positions_by_chrom.setdefault(chrom, set()).update(position for position in positions if position >= 0)


def _aggregate_distributions(
    path: Path,
    position_scores: PositionScores,
    analytics_dir: Path,
    chunk_size: int,
) -> tuple[pd.DataFrame, list[dict[str, object]], dict[str, int]]:
    with tempfile.NamedTemporaryFile(
        prefix=".candidate_phylop.", suffix=".sqlite3", dir=analytics_dir, delete=False
    ) as handle:
        database_path = Path(handle.name)
    connection = sqlite3.connect(database_path)
    connection.executescript(
        """
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous = OFF;
        PRAGMA temp_store = FILE;
        CREATE TABLE scores (
            variant_id TEXT NOT NULL,
            strategy TEXT NOT NULL,
            gnomad_status TEXT NOT NULL,
            score REAL,
            PRIMARY KEY (variant_id, strategy)
        ) WITHOUT ROWID;
        """
    )
    try:
        for chunk in pd.read_csv(
            path,
            sep="\t",
            compression="gzip",
            keep_default_na=False,
            usecols=list(REQUIRED_COLUMNS),
            chunksize=chunk_size,
        ):
            records = []
            gnomad_af = pd.to_numeric(chunk["gnomad_af"], errors="coerce")
            for row, af in zip(chunk.itertuples(index=False), gnomad_af):
                if str(row.lookup_status) != "ok":
                    continue
                required = _required_positions(
                    row.lookup_chrom,
                    row.lookup_pos,
                    row.lookup_ref,
                    row.lookup_alt,
                    position_scores.track.chrom_style,
                )
                if not required:
                    continue
                values = [position_scores.values.get(position) for position in required]
                score = float(np.mean(values)) if all(value is not None for value in values) else None
                status = "found" if not pd.isna(af) else "not_found"
                for strategy in split_strategies(str(row.strategies)):
                    records.append((str(row.variant_key), strategy, status, score))
            connection.executemany(
                "INSERT OR IGNORE INTO scores (variant_id, strategy, gnomad_status, score) VALUES (?, ?, ?, ?)",
                records,
            )
            connection.commit()

        groups_frame = pd.read_sql_query(
            """
            SELECT strategy, gnomad_status, COUNT(*) AS variant_count,
                   SUM(score IS NOT NULL) AS scored_count
            FROM scores
            GROUP BY strategy, gnomad_status
            ORDER BY strategy, gnomad_status
            """,
            connection,
        )
        rows = []
        for group in groups_frame.itertuples(index=False):
            values = np.asarray(
                [row[0] for row in connection.execute(
                    "SELECT score FROM scores WHERE strategy = ? AND gnomad_status = ? AND score IS NOT NULL",
                    (group.strategy, group.gnomad_status),
                )],
                dtype=float,
            )
            if values.size == 0:
                continue
            quantile_values = np.quantile(values, QUANTILES)
            rows.extend(
                {
                    "strategy": str(group.strategy),
                    "gnomad_status": str(group.gnomad_status),
                    "quantile": float(quantile),
                    "phyloP100way": float(score),
                    "variant_count": int(group.variant_count),
                    "scored_count": int(group.scored_count),
                }
                for quantile, score in zip(QUANTILES, quantile_values)
            )
        unique_alleles = int(connection.execute("SELECT COUNT(DISTINCT variant_id) FROM scores").fetchone()[0])
        membership_count = int(connection.execute("SELECT COUNT(*) FROM scores").fetchone()[0])
    finally:
        connection.close()
        database_path.unlink(missing_ok=True)

    groups = groups_frame.to_dict(orient="records")
    for group in groups:
        variant_count = int(group["variant_count"])
        scored_count = int(group["scored_count"])
        group["variant_count"] = variant_count
        group["scored_count"] = scored_count
        group["score_coverage"] = scored_count / variant_count if variant_count else 0.0
    return pd.DataFrame(rows), groups, {
        "unique_usable_allele_count": unique_alleles,
        "strategy_variant_membership_count": membership_count,
    }


def _required_positions(
    chrom: object,
    pos: object,
    ref: object,
    alt: object,
    chrom_style: str,
) -> list[tuple[str, int]]:
    positions, _basis = score_positions(int(pos), str(ref), str(alt))
    formatted_chrom = format_chrom(chrom, chrom_style)
    return [(formatted_chrom, position) for position in positions if position >= 0]


def _write_frame(path: Path, frame: pd.DataFrame) -> None:
    with tempfile.NamedTemporaryFile(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        frame.to_csv(temporary, sep="\t", index=False, compression="gzip", lineterminator="\n")
        temporary.chmod(0o644)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _write_json(path: Path, payload: dict) -> None:
    with tempfile.NamedTemporaryFile(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, mode="w", delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    try:
        temporary.chmod(0o644)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
