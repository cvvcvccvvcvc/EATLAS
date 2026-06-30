#!/usr/bin/env python3
"""Annotate variant rows with remote bigWig conservation scores.

This script is standalone and optimized for low-disk validation runs. It does
not download whole bigWig files. Instead, it opens remote bigWig URLs with
pyBigWig, reads genomic intervals in blocks, extracts values for variant
positions, and writes only the final annotated TSV.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

try:
    import pyBigWig
except ImportError:  # pragma: no cover - exercised only in missing dependency envs
    pyBigWig = None


@dataclass(frozen=True)
class Track:
    name: str
    url: str
    chrom_style: str


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants-tsv", required=True, type=Path, help="Input variant TSV/TSV.GZ.")
    parser.add_argument("--out-tsv", required=True, type=Path, help="Output annotated TSV/TSV.GZ.")
    parser.add_argument("--summary-json", type=Path, help="Optional JSON run summary.")
    parser.add_argument("--chrom-column", default="genomic_accession")
    parser.add_argument("--pos-column", default="genomic_start1", help="1-based variant position column.")
    parser.add_argument("--gene-id", action="append", help="Optional gene_id filter. Can be supplied multiple times.")
    parser.add_argument("--strategy", help="Optional strategy filter.")
    parser.add_argument(
        "--tracks",
        default="phyloP100way,phastCons100way,GERP_RS_92mammals",
        help=f"Comma-separated built-in tracks. Available: {', '.join(DEFAULT_TRACKS)}",
    )
    parser.add_argument(
        "--track",
        action="append",
        default=[],
        metavar="NAME=URL=STYLE",
        help="Custom track spec. STYLE must be ucsc or ensembl. Can be supplied multiple times.",
    )
    parser.add_argument("--max-block-bp", type=int, default=250_000)
    parser.add_argument("--max-gap-bp", type=int, default=50_000)
    parser.add_argument(
        "--remote-retries",
        type=int,
        default=3,
        help="Total attempts for each remote bigWig open/read operation.",
    )
    parser.add_argument(
        "--retry-sleep-seconds",
        type=float,
        default=5.0,
        help="Initial sleep before retrying a failed remote bigWig operation. Uses exponential backoff.",
    )
    parser.add_argument("--precision", type=int, default=6)
    return parser.parse_args()


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
        return rows, list(reader.fieldnames or [])


def write_rows(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
            count += 1
    return count


def parse_tracks(raw_names: str, custom_specs: list[str]) -> list[Track]:
    tracks: list[Track] = []
    for raw_name in raw_names.split(","):
        name = raw_name.strip()
        if not name:
            continue
        if name not in DEFAULT_TRACKS:
            raise ValueError(f"Unknown built-in track {name!r}; available: {', '.join(DEFAULT_TRACKS)}")
        tracks.append(DEFAULT_TRACKS[name])

    for spec in custom_specs:
        parts = spec.split("=", 2)
        if len(parts) != 3:
            raise ValueError("--track must use NAME=URL=STYLE")
        name, url, style = [part.strip() for part in parts]
        if not name or not url or style not in {"ucsc", "ensembl"}:
            raise ValueError("--track must use non-empty NAME/URL and STYLE of ucsc or ensembl")
        tracks.append(Track(name=name, url=url, chrom_style=style))

    if not tracks:
        raise ValueError("At least one conservation track must be selected")
    seen = set()
    for track in tracks:
        if track.name in seen:
            raise ValueError(f"Duplicate output track name: {track.name}")
        seen.add(track.name)
    return tracks


def base_chrom(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if text.startswith("chr"):
        text = text[3:]
    if text == "MT":
        return "M"
    if text in {"X", "Y", "M"}:
        return text
    if text.isdigit():
        return str(int(text))
    if text.startswith("NC_"):
        base = text.split(".", 1)[0]
        try:
            number = int(base.split("_", 1)[1])
        except (IndexError, ValueError):
            return text
        if 1 <= number <= 22:
            return str(number)
        if number == 23:
            return "X"
        if number == 24:
            return "Y"
        if number == 12920:
            return "M"
    return text


def format_chrom(value: object, style: str) -> str:
    chrom = base_chrom(value)
    if not chrom:
        return ""
    if style == "ucsc":
        return "chrM" if chrom == "M" else f"chr{chrom}"
    if style == "ensembl":
        return "MT" if chrom == "M" else chrom
    raise ValueError(f"Unsupported chromosome style: {style}")


def to_int(value: object) -> int | None:
    try:
        if value is None or value == "":
            return None
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def row_included(row: dict[str, str], gene_ids: set[str] | None, strategy: str | None) -> bool:
    if gene_ids is not None and row.get("gene_id", "") not in gene_ids:
        return False
    if strategy is not None and row.get("strategy", "") != strategy:
        return False
    return True


def build_position_index(
    rows: list[dict[str, str]],
    chrom_column: str,
    pos_column: str,
    track: Track,
) -> tuple[dict[str, set[int]], dict[tuple[str, int], list[int]], int]:
    positions_by_chrom: dict[str, set[int]] = {}
    row_indices_by_position: dict[tuple[str, int], list[int]] = {}
    skipped = 0
    for index, row in enumerate(rows):
        chrom = format_chrom(row.get(chrom_column, ""), track.chrom_style)
        pos1 = to_int(row.get(pos_column, ""))
        if not chrom or pos1 is None or pos1 <= 0:
            skipped += 1
            continue
        pos0 = pos1 - 1
        positions_by_chrom.setdefault(chrom, set()).add(pos0)
        row_indices_by_position.setdefault((chrom, pos0), []).append(index)
    return positions_by_chrom, row_indices_by_position, skipped


def make_blocks(positions: set[int], max_block_bp: int, max_gap_bp: int) -> list[tuple[int, int, list[int]]]:
    if not positions:
        return []
    sorted_positions = sorted(positions)
    blocks: list[tuple[int, int, list[int]]] = []
    block_positions: list[int] = []
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


def value_to_text(value: float | None, precision: int) -> str:
    if value is None:
        return ""
    try:
        if math.isnan(value):
            return ""
    except TypeError:
        return ""
    return f"{float(value):.{precision}g}"


def sleep_before_retry(base_seconds: float, attempt: int) -> None:
    if base_seconds <= 0:
        return
    time.sleep(base_seconds * (2 ** (attempt - 1)))


def open_bigwig_with_retries(
    track: Track,
    remote_retries: int,
    retry_sleep_seconds: float,
) -> tuple[object, int, float]:
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
            print(
                f"warning: {track.name}: open failed on attempt {attempt}/{remote_retries}: {exc}; retrying",
                file=sys.stderr,
            )
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
            print(
                f"warning: {track_name}: failed {chrom}:{start0}-{end0} "
                f"on attempt {attempt}/{remote_retries}: {exc}; retrying",
                file=sys.stderr,
            )
            sleep_before_retry(retry_sleep_seconds, attempt)
    raise RuntimeError(f"failed {chrom}:{start0}-{end0} after {remote_retries} attempts: {last_error}")


def annotate_track(
    rows: list[dict[str, str]],
    track: Track,
    chrom_column: str,
    pos_column: str,
    max_block_bp: int,
    max_gap_bp: int,
    remote_retries: int,
    retry_sleep_seconds: float,
    precision: int,
) -> dict[str, object]:
    if pyBigWig is None:
        raise RuntimeError("pyBigWig is required. Install it with: python -m pip install pyBigWig")

    positions_by_chrom, row_indices_by_position, skipped = build_position_index(
        rows=rows,
        chrom_column=chrom_column,
        pos_column=pos_column,
        track=track,
    )
    for row in rows:
        row[track.name] = ""

    bw, open_attempts, open_seconds = open_bigwig_with_retries(
        track=track,
        remote_retries=remote_retries,
        retry_sleep_seconds=retry_sleep_seconds,
    )
    chrom_sizes = bw.chroms()

    block_count = 0
    annotated_positions = 0
    missing_positions = 0
    read_seconds = 0.0
    read_retry_count = 0

    for chrom in sorted(positions_by_chrom):
        if chrom not in chrom_sizes:
            missing = len(positions_by_chrom[chrom])
            missing_positions += missing
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
                    bw=bw,
                    track_name=track.name,
                    chrom=chrom,
                    start0=start0,
                    end0=end0,
                    remote_retries=remote_retries,
                    retry_sleep_seconds=retry_sleep_seconds,
                )
            except RuntimeError as exc:
                print(f"warning: {track.name}: failed {chrom}:{start0}-{end0}: {exc}", file=sys.stderr)
                missing_positions += len(block_positions)
                continue
            read_seconds += time.perf_counter() - read_start
            read_retry_count += read_attempts - 1
            block_count += 1
            for pos0 in block_positions:
                if pos0 < start0 or pos0 >= end0:
                    missing_positions += 1
                    continue
                value = values[pos0 - start0]
                text = value_to_text(value, precision)
                if text == "":
                    missing_positions += 1
                    continue
                annotated_positions += 1
                for row_index in row_indices_by_position.get((chrom, pos0), []):
                    rows[row_index][track.name] = text
    bw.close()

    return {
        "track": track.name,
        "url": track.url,
        "chrom_style": track.chrom_style,
        "open_seconds": round(open_seconds, 3),
        "open_attempts": open_attempts,
        "read_seconds": round(read_seconds, 3),
        "read_retry_count": read_retry_count,
        "block_count": block_count,
        "unique_positions": sum(len(items) for items in positions_by_chrom.values()),
        "annotated_positions": annotated_positions,
        "missing_positions": missing_positions,
        "skipped_bad_coordinate_rows": skipped,
    }


def main() -> None:
    args = parse_args()
    if args.remote_retries < 1:
        raise ValueError("--remote-retries must be >= 1")
    if args.retry_sleep_seconds < 0:
        raise ValueError("--retry-sleep-seconds must be >= 0")
    tracks = parse_tracks(args.tracks, args.track)
    rows, fields = read_rows(args.variants_tsv)
    required = {args.chrom_column, args.pos_column}
    missing = required - set(fields)
    if missing:
        raise ValueError(f"Missing required input column(s): {', '.join(sorted(missing))}")

    gene_ids = set(args.gene_id) if args.gene_id else None
    rows = [row for row in rows if row_included(row, gene_ids, args.strategy)]
    if not rows:
        raise ValueError("No input rows remain after filters")

    print(f"Annotating {len(rows)} rows with {len(tracks)} conservation tracks", file=sys.stderr)
    summaries = []
    for track in tracks:
        print(f"Opening {track.name}: {track.url}", file=sys.stderr)
        summary = annotate_track(
            rows=rows,
            track=track,
            chrom_column=args.chrom_column,
            pos_column=args.pos_column,
            max_block_bp=args.max_block_bp,
            max_gap_bp=args.max_gap_bp,
            remote_retries=args.remote_retries,
            retry_sleep_seconds=args.retry_sleep_seconds,
            precision=args.precision,
        )
        summaries.append(summary)
        print(
            f"{track.name}: {summary['annotated_positions']}/{summary['unique_positions']} "
            f"positions annotated in {summary['block_count']} blocks "
            f"(open {summary['open_seconds']}s, read {summary['read_seconds']}s)",
            file=sys.stderr,
        )

    output_fields = list(fields)
    for track in tracks:
        if track.name not in output_fields:
            output_fields.append(track.name)
    row_count = write_rows(args.out_tsv, rows, output_fields)
    run_summary = {
        "input": str(args.variants_tsv),
        "output": str(args.out_tsv),
        "row_count": row_count,
        "tracks": summaries,
        "filters": {
            "gene_id": sorted(gene_ids) if gene_ids else [],
            "strategy": args.strategy or "",
        },
        "max_block_bp": args.max_block_bp,
        "max_gap_bp": args.max_gap_bp,
        "remote_retries": args.remote_retries,
        "retry_sleep_seconds": args.retry_sleep_seconds,
    }
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(run_summary, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
