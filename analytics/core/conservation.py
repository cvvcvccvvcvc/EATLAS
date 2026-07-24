"""Conservation score cache for ClinVar validation alleles."""

from __future__ import annotations

import csv
import gzip
import json
import math
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd

try:
    import pyBigWig
except ImportError:  # pragma: no cover - depends on the local analytics env
    pyBigWig = None


@dataclass(frozen=True)
class Track:
    name: str
    url: str
    chrom_style: str


@dataclass(frozen=True)
class ConservationAnnotations:
    annotations_path: Path
    manifest_path: Path
    annotations: pd.DataFrame
    manifest: dict
    score_columns: list[str]


DEFAULT_TRACKS = {
    "phyloP100way": Track(
        name="phyloP100way",
        url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw",
        chrom_style="ucsc",
    ),
    "phastCons100way": Track(
        name="phastCons100way",
        url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw",
        chrom_style="ucsc",
    ),
    "GERP_RS_92mammals": Track(
        name="GERP_RS_92mammals",
        url=(
            "https://ftp.ensembl.org/pub/current/compara/conservation_scores/"
            "92_mammals.gerp_conservation_score/"
            "gerp_conservation_scores.homo_sapiens.GRCh38.bw"
        ),
        chrom_style="ensembl",
    ),
}
DEFAULT_TRACK_NAMES = "phyloP100way"
CACHE_VERSION = 2
CONSERVATION_FIELDS = ["variant_key", "chrom", "pos"]


def build_conservation_annotations(
    *,
    universe: pd.DataFrame,
    universe_path: Path,
    analytics_dir: Path,
    track_names: str = DEFAULT_TRACK_NAMES,
    max_block_bp: int = 250_000,
    max_gap_bp: int = 50_000,
    remote_retries: int = 3,
    retry_sleep_seconds: float = 5.0,
    precision: int = 6,
) -> ConservationAnnotations:
    tracks = parse_tracks(track_names)
    analytics_dir.mkdir(parents=True, exist_ok=True)
    annotations_path = analytics_dir / "clinvar_universe.snv.conservation.tsv.gz"
    manifest_path = analytics_dir / "clinvar_universe.snv.conservation.manifest.json"
    score_columns = [track.name for track in tracks]
    expected_inputs = {
        "cache_version": CACHE_VERSION,
        "universe": path_metadata(universe_path),
        "tracks": [asdict(track) for track in tracks],
        "max_block_bp": max_block_bp,
        "max_gap_bp": max_gap_bp,
        "remote_retries": remote_retries,
        "retry_sleep_seconds": retry_sleep_seconds,
        "precision": precision,
    }
    previous_manifest = None
    if annotations_path.exists() and manifest_path.exists():
        candidate = json.loads(manifest_path.read_text())
        if candidate.get("inputs") == expected_inputs:
            if candidate.get("complete") is True:
                annotations = read_annotations(annotations_path, score_columns)
                return ConservationAnnotations(annotations_path, manifest_path, annotations, candidate, score_columns)
            previous_manifest = candidate

    rows = (
        read_annotation_rows(annotations_path, score_columns)
        if previous_manifest is not None
        else snv_universe_rows(universe)
    )
    previous_summaries = {
        str(summary.get("track")): summary
        for summary in (previous_manifest or {}).get("tracks", [])
    }
    summaries = []
    for track in tracks:
        previous_summary = previous_summaries.get(track.name)
        if previous_summary and previous_summary.get("status") == "complete":
            summaries.append(previous_summary)
            continue
        print(f"Annotating ClinVar SNVs with {track.name}: {track.url}", file=sys.stderr)
        try:
            summary = annotate_track(
                rows=rows,
                track=track,
                max_block_bp=max_block_bp,
                max_gap_bp=max_gap_bp,
                remote_retries=remote_retries,
                retry_sleep_seconds=retry_sleep_seconds,
                precision=precision,
            )
        except RuntimeError as exc:
            summary = failed_track_summary(rows, track, str(exc), remote_retries)
            print(f"warning: {track.name}: {exc}", file=sys.stderr)
        summaries.append(summary)
        print(
            f"{track.name}: {summary['annotated_positions']}/{summary['unique_positions']} "
            f"positions annotated in {summary['block_count']} blocks "
            f"(open {summary['open_seconds']}s, read {summary['read_seconds']}s)",
            file=sys.stderr,
        )

    write_annotations(annotations_path, rows, score_columns)
    manifest = {
        "inputs": expected_inputs,
        "complete": all(summary.get("status") == "complete" for summary in summaries),
        "row_count": len(rows),
        "score_columns": score_columns,
        "tracks": summaries,
        "annotation_tsv": str(annotations_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    annotations = read_annotations(annotations_path, score_columns)
    return ConservationAnnotations(annotations_path, manifest_path, annotations, manifest, score_columns)


def parse_tracks(raw_names: str) -> list[Track]:
    tracks = []
    for raw_name in raw_names.split(","):
        name = raw_name.strip()
        if not name:
            continue
        if name not in DEFAULT_TRACKS:
            raise ValueError(f"Unknown conservation track {name!r}; available: {', '.join(DEFAULT_TRACKS)}")
        tracks.append(DEFAULT_TRACKS[name])
    if not tracks:
        raise ValueError("At least one conservation track is required.")
    return tracks


def snv_universe_rows(universe: pd.DataFrame) -> list[dict[str, str]]:
    required = {"variant_key", "variant_type", "chrom", "pos"}
    missing = required - set(universe.columns)
    if missing:
        raise ValueError(f"ClinVar universe missing required columns: {', '.join(sorted(missing))}")
    snv = universe[universe["variant_type"].astype(str) == "snv"][["variant_key", "chrom", "pos"]].copy()
    snv = snv.drop_duplicates("variant_key").sort_values(["chrom", "pos", "variant_key"], kind="mergesort")
    return [
        {"variant_key": str(row.variant_key), "chrom": str(row.chrom), "pos": str(int(row.pos))}
        for row in snv.itertuples(index=False)
    ]


def open_text(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, newline="") if str(path).endswith(".gz") else path.open(mode, newline="")


def read_annotations(path: Path, score_columns: list[str]) -> pd.DataFrame:
    annotations = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False, low_memory=False)
    for column in score_columns:
        if column not in annotations.columns:
            annotations[column] = ""
        annotations[column] = pd.to_numeric(annotations[column], errors="coerce")
    return annotations


def read_annotation_rows(path: Path, score_columns: list[str]) -> list[dict[str, str]]:
    fields = [*CONSERVATION_FIELDS, *score_columns]
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(fields) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Conservation cache missing columns: {', '.join(sorted(missing))}")
        return [{field: row.get(field, "") for field in fields} for row in reader]


def write_annotations(path: Path, rows: list[dict[str, str]], score_columns: list[str]) -> None:
    fields = [*CONSERVATION_FIELDS, *score_columns]
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def base_chrom(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if text.startswith("chr"):
        text = text[3:]
    if text == "MT":
        return "M"
    return str(int(text)) if text.isdigit() else text


def format_chrom(value: object, style: str) -> str:
    chrom = base_chrom(value)
    if not chrom:
        return ""
    if style == "ucsc":
        return "chrM" if chrom == "M" else f"chr{chrom}"
    if style == "ensembl":
        return "MT" if chrom == "M" else chrom
    raise ValueError(f"Unsupported chromosome style: {style}")


def make_blocks(positions: set[int], max_block_bp: int, max_gap_bp: int) -> list[tuple[int, int, list[int]]]:
    if not positions:
        return []
    sorted_positions = sorted(positions)
    blocks = []
    block_positions = []
    block_start = sorted_positions[0]
    previous = sorted_positions[0]
    for pos in sorted_positions:
        span = (pos + 1) - block_start
        gap = pos - previous
        if block_positions and (span > max_block_bp or gap > max_gap_bp):
            blocks.append((block_start, previous + 1, block_positions))
            block_start = pos
            block_positions = []
        block_positions.append(pos)
        previous = pos
    if block_positions:
        blocks.append((block_start, previous + 1, block_positions))
    return blocks


def annotate_track(
    *,
    rows: list[dict[str, str]],
    track: Track,
    max_block_bp: int,
    max_gap_bp: int,
    remote_retries: int,
    retry_sleep_seconds: float,
    precision: int,
) -> dict[str, object]:
    if pyBigWig is None:
        raise RuntimeError("pyBigWig is required for conservation annotation.")

    positions_by_chrom, row_indices_by_position = build_position_index(rows, track)
    for row in rows:
        row[track.name] = ""

    bw, open_attempts, open_seconds = open_bigwig_with_retries(track, remote_retries, retry_sleep_seconds)
    chrom_sizes = bw.chroms()
    block_count = 0
    annotated_positions = 0
    missing_positions = 0
    read_seconds = 0.0
    read_retry_count = 0
    failed_block_count = 0
    first_error = ""

    for chrom in sorted(positions_by_chrom):
        if chrom not in chrom_sizes:
            missing_positions += len(positions_by_chrom[chrom])
            print(f"warning: {track.name}: chromosome {chrom} not found in bigWig", file=sys.stderr)
            continue
        chrom_size = int(chrom_sizes[chrom])
        for start0, end0, block_positions in make_blocks(positions_by_chrom[chrom], max_block_bp, max_gap_bp):
            start0 = max(0, start0)
            end0 = min(chrom_size, end0)
            if end0 <= start0:
                missing_positions += len(block_positions)
                continue
            read_start = time.perf_counter()
            try:
                values, read_attempts = read_values_with_retries(
                    bw, track.name, chrom, start0, end0, remote_retries, retry_sleep_seconds
                )
            except RuntimeError as exc:
                print(f"warning: {track.name}: failed {chrom}:{start0}-{end0}: {exc}", file=sys.stderr)
                missing_positions += len(block_positions)
                failed_block_count += 1
                if not first_error:
                    first_error = str(exc)
                continue
            read_seconds += time.perf_counter() - read_start
            read_retry_count += read_attempts - 1
            block_count += 1
            for pos0 in block_positions:
                if pos0 < start0 or pos0 >= end0:
                    missing_positions += 1
                    continue
                text = value_to_text(values[pos0 - start0], precision)
                if text == "":
                    missing_positions += 1
                    continue
                annotated_positions += 1
                for row_index in row_indices_by_position.get((chrom, pos0), []):
                    rows[row_index][track.name] = text
    bw.close()
    return {
        "track": track.name,
        "status": "complete" if failed_block_count == 0 else "partial",
        "error": first_error,
        "url": track.url,
        "chrom_style": track.chrom_style,
        "open_seconds": round(open_seconds, 3),
        "open_attempts": open_attempts,
        "read_seconds": round(read_seconds, 3),
        "read_retry_count": read_retry_count,
        "block_count": block_count,
        "failed_block_count": failed_block_count,
        "unique_positions": sum(len(items) for items in positions_by_chrom.values()),
        "annotated_positions": annotated_positions,
        "missing_positions": missing_positions,
    }


def failed_track_summary(
    rows: list[dict[str, str]],
    track: Track,
    error: str,
    open_attempts: int,
) -> dict[str, object]:
    positions_by_chrom, _row_indices = build_position_index(rows, track)
    unique_positions = sum(len(items) for items in positions_by_chrom.values())
    for row in rows:
        row[track.name] = ""
    return {
        "track": track.name,
        "status": "failed",
        "error": error,
        "url": track.url,
        "chrom_style": track.chrom_style,
        "open_seconds": 0.0,
        "open_attempts": open_attempts,
        "read_seconds": 0.0,
        "read_retry_count": 0,
        "block_count": 0,
        "failed_block_count": 0,
        "unique_positions": unique_positions,
        "annotated_positions": 0,
        "missing_positions": unique_positions,
    }


def build_position_index(
    rows: list[dict[str, str]], track: Track
) -> tuple[dict[str, set[int]], dict[tuple[str, int], list[int]]]:
    positions_by_chrom: dict[str, set[int]] = {}
    row_indices_by_position: dict[tuple[str, int], list[int]] = {}
    for index, row in enumerate(rows):
        chrom = format_chrom(row.get("chrom", ""), track.chrom_style)
        pos1 = int(row["pos"])
        pos0 = pos1 - 1
        positions_by_chrom.setdefault(chrom, set()).add(pos0)
        row_indices_by_position.setdefault((chrom, pos0), []).append(index)
    return positions_by_chrom, row_indices_by_position


def value_to_text(value: float | None, precision: int) -> str:
    if value is None:
        return ""
    try:
        if math.isnan(value):
            return ""
    except TypeError:
        return ""
    return f"{float(value):.{precision}g}"


def open_bigwig_with_retries(track: Track, remote_retries: int, retry_sleep_seconds: float) -> tuple[object, int, float]:
    start_time = time.perf_counter()
    last_error: Exception | None = None
    for attempt in range(1, remote_retries + 1):
        try:
            bw = pyBigWig.open(track.url)
            if bw is None:
                raise RuntimeError("pyBigWig.open returned None")
            return bw, attempt, time.perf_counter() - start_time
        except (RuntimeError, OSError) as exc:
            last_error = exc
            if attempt >= remote_retries:
                break
            sleep_before_retry(retry_sleep_seconds, attempt)
    raise RuntimeError(f"{track.name}: failed to open remote bigWig after {remote_retries} attempts: {last_error}")


def read_values_with_retries(
    bw: object,
    track_name: str,
    chrom: str,
    start0: int,
    end0: int,
    remote_retries: int,
    retry_sleep_seconds: float,
) -> tuple[list[float | None], int]:
    last_error: Exception | None = None
    for attempt in range(1, remote_retries + 1):
        try:
            values = bw.values(chrom, start0, end0)
            if values is None:
                raise RuntimeError("pyBigWig returned no values")
            return values, attempt
        except (RuntimeError, OSError) as exc:
            last_error = exc
            if attempt >= remote_retries:
                break
            sleep_before_retry(retry_sleep_seconds, attempt)
    raise RuntimeError(f"{track_name}: failed {chrom}:{start0}-{end0} after {remote_retries} attempts: {last_error}")


def sleep_before_retry(base_seconds: float, attempt: int) -> None:
    if base_seconds > 0:
        time.sleep(base_seconds * (2 ** (attempt - 1)))


def path_metadata(path: Path) -> dict[str, object]:
    if not path.exists():
        raise FileNotFoundError(path)
    stat = path.stat()
    return {"path": str(path.resolve()), "size_bytes": stat.st_size, "mtime": int(stat.st_mtime)}
