#!/usr/bin/env python3
"""Run nucmer for one gene task and normalize comparator alignment evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable

from alignment_task_io import load_task_context, materialize_task_fastas
from feature_coverage import summarize_feature_coverage


TSV_NULL = ""

SEGMENT_FIELDS = [
    "gene_id",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "strategy",
    "tool",
    "preset",
    "sequence_id",
    "target_id",
    "query_id",
    "target_start0",
    "target_end0",
    "query_start0",
    "query_end0",
    "strand",
    "matches",
    "block_length",
    "identity",
    "mapq",
    "is_primary",
    "divergence",
    "gap_compressed_divergence",
    "native_record_id",
    "qc_flags",
]

EVENT_FIELDS = [
    "gene_id",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "strategy",
    "tool",
    "preset",
    "event_id",
    "event_type",
    "target_start0",
    "target_end0",
    "genomic_accession",
    "genomic_start1",
    "genomic_end1",
    "ref",
    "alt",
    "query_id",
    "strand",
    "native_record_id",
    "qc_flags",
]

SUMMARY_FIELDS = [
    "gene_id",
    "ortholog_gene_id",
    "tax_id",
    "taxname",
    "strategy",
    "tool",
    "preset",
    "status",
    "target_length",
    "query_length",
    "segment_count",
    "primary_segment_count",
    "secondary_segment_count",
    "aligned_target_bp",
    "aligned_query_bp",
    "target_coverage",
    "query_coverage",
    "best_identity",
    "mean_identity",
    "event_count",
    "qc_flags",
]

FAILURE_FIELDS = ["gene_id", "ortholog_gene_id", "strategy", "tool", "failure_type", "message"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-dir", required=True, type=Path)
    parser.add_argument("--source-target-fasta", required=True, type=Path)
    parser.add_argument("--source-ortholog-fasta", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--nucmer-bin", required=True)
    parser.add_argument("--show-coords-bin", required=True)
    parser.add_argument("--show-snps-bin", required=True)
    parser.add_argument("--target-features", type=Path)
    parser.add_argument("--keep-native", default="false")
    return parser.parse_args()


def truthy(value: str) -> bool:
    return str(value).lower() in {"1", "true", "yes", "y"}


def write_tsv_gz(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, TSV_NULL) for field in fields})
            count += 1
    return count


def interval_union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def genomic_coords(target_meta: dict[str, str], start0: int, end0: int) -> tuple[str, str]:
    begin_text = target_meta.get("genomic_begin") or ""
    if not begin_text:
        return "", ""
    begin = int(begin_text)
    genomic_start = begin + start0
    if end0 > start0:
        genomic_end = begin + end0 - 1
    else:
        genomic_end = genomic_start
    return str(genomic_start), str(genomic_end)


def run_command(cmd: list[str], stdout_path: Path | None = None) -> str:
    if stdout_path:
        with stdout_path.open("w") as handle:
            result = subprocess.run(cmd, text=True, stdout=handle, stderr=subprocess.PIPE)
    else:
        result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {(result.stderr or result.stdout or '').strip()}")
    return " ".join(cmd)


def empty_summary(gene_id: str, meta: dict[str, str], target_length: int) -> dict[str, object]:
    return {
        "gene_id": gene_id,
        "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
        "tax_id": meta.get("tax_id", ""),
        "taxname": meta.get("taxname", ""),
        "strategy": "nucmer",
        "tool": "nucmer",
        "preset": "default",
        "status": "not_run",
        "target_length": target_length,
        "query_length": int(meta.get("sequence_length") or 0),
        "segment_count": 0,
        "primary_segment_count": 0,
        "secondary_segment_count": 0,
        "target_intervals": [],
        "query_intervals": [],
        "best_identity": 0.0,
        "identities": [],
        "event_count": 0,
        "qc_flags": {"unfiltered_nucmer"},
    }


def parse_coords(
    coords_path: Path,
    gene_id: str,
    target_meta: dict[str, str],
    meta_by_sequence: dict[str, dict[str, str]],
    summaries: dict[str, dict[str, object]],
) -> tuple[list[dict[str, object]], dict[str, list[dict[str, object]]]]:
    segments: list[dict[str, object]] = []
    segments_by_query: dict[str, list[dict[str, object]]] = {sequence_id: [] for sequence_id in meta_by_sequence}
    target_id = f"target_{gene_id}"

    with coords_path.open() as handle:
        for line_index, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 13:
                continue
            ref_id = fields[-2]
            query_id = fields[-1]
            meta = meta_by_sequence.get(query_id)
            if meta is None:
                continue

            s1, e1, s2, e2 = (int(fields[index]) for index in range(4))
            len1 = int(fields[4])
            len2 = int(fields[5])
            identity_pct = float(fields[6])
            len_r = int(float(fields[7]))
            len_q = int(float(fields[8]))
            target_start0 = min(s1, e1) - 1
            target_end0 = max(s1, e1)
            query_start0 = min(s2, e2) - 1
            query_end0 = max(s2, e2)
            strand = "-" if s2 > e2 else "+"
            identity = identity_pct / 100.0
            block_length = max(len1, len2)
            matches = round(identity * block_length)

            segment = {
                "gene_id": gene_id,
                "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
                "tax_id": meta.get("tax_id", ""),
                "taxname": meta.get("taxname", ""),
                "strategy": "nucmer",
                "tool": "nucmer",
                "preset": "default",
                "sequence_id": query_id,
                "target_id": ref_id or target_id,
                "query_id": query_id,
                "target_start0": target_start0,
                "target_end0": target_end0,
                "query_start0": query_start0,
                "query_end0": query_end0,
                "strand": strand,
                "matches": matches,
                "block_length": block_length,
                "identity": f"{identity:.6f}",
                "mapq": "",
                "is_primary": "true",
                "divergence": "",
                "gap_compressed_divergence": "",
                "native_record_id": line_index,
                "qc_flags": "unfiltered_nucmer",
            }
            segments.append(segment)
            segments_by_query[query_id].append(segment)

            summary = summaries[query_id]
            summary["status"] = "aligned"
            summary["target_length"] = len_r or summary["target_length"]
            summary["query_length"] = len_q or summary["query_length"]
            summary["segment_count"] += 1
            summary["primary_segment_count"] += 1
            summary["target_intervals"].append((target_start0, target_end0))
            summary["query_intervals"].append((query_start0, query_end0))
            summary["identities"].append(identity)
            summary["best_identity"] = max(summary["best_identity"], identity)

    return segments, segments_by_query


def segment_for_event(segments: list[dict[str, object]], target_pos0: int) -> dict[str, object] | None:
    for segment in segments:
        if int(segment["target_start0"]) <= target_pos0 < int(segment["target_end0"]):
            return segment
    return segments[0] if segments else None


def parse_snps(
    snps_path: Path,
    gene_id: str,
    target_meta: dict[str, str],
    meta_by_sequence: dict[str, dict[str, str]],
    segments_by_query: dict[str, list[dict[str, object]]],
    summaries: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    events: list[dict[str, object]] = []
    with snps_path.open() as handle:
        for event_index, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            query_id = fields[-1]
            meta = meta_by_sequence.get(query_id)
            if meta is None:
                continue
            try:
                p1 = int(fields[0])
            except ValueError:
                continue
            ref = fields[1].upper()
            alt = fields[2].upper()
            if ref == ".":
                event_type = "ins"
                target_start0 = max(p1 - 1, 0)
                target_end0 = target_start0
                ref = ""
            elif alt == ".":
                event_type = "del"
                target_start0 = p1 - 1
                target_end0 = target_start0 + 1
                alt = ""
            else:
                event_type = "snv"
                target_start0 = p1 - 1
                target_end0 = target_start0 + 1
            segment = segment_for_event(segments_by_query.get(query_id, []), target_start0)
            strand = segment.get("strand", "") if segment else ""
            genomic_start, genomic_end = genomic_coords(target_meta, target_start0, target_end0)
            events.append(
                {
                    "gene_id": gene_id,
                    "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
                    "tax_id": meta.get("tax_id", ""),
                    "taxname": meta.get("taxname", ""),
                    "strategy": "nucmer",
                    "tool": "nucmer",
                    "preset": "default",
                    "event_id": event_index,
                    "event_type": event_type,
                    "target_start0": target_start0,
                    "target_end0": target_end0,
                    "genomic_accession": target_meta.get("genomic_accession", ""),
                    "genomic_start1": genomic_start,
                    "genomic_end1": genomic_end,
                    "ref": ref,
                    "alt": alt,
                    "query_id": query_id,
                    "strand": strand,
                    "native_record_id": event_index,
                    "qc_flags": "unfiltered_nucmer",
                }
            )
            summaries[query_id]["event_count"] += 1
    return events


def finalize_summary(row: dict[str, object]) -> dict[str, object]:
    target_length = int(row.get("target_length") or 0)
    query_length = int(row.get("query_length") or 0)
    aligned_target = interval_union_length(row.pop("target_intervals"))
    aligned_query = interval_union_length(row.pop("query_intervals"))
    identities = row.pop("identities")
    flags = row.pop("qc_flags")
    if row["status"] == "not_run":
        row["status"] = "no_alignment"
        flags.add("no_alignment")
    mean_identity = sum(identities) / len(identities) if identities else 0.0
    row.update(
        {
            "aligned_target_bp": aligned_target,
            "aligned_query_bp": aligned_query,
            "target_coverage": f"{(aligned_target / target_length) if target_length else 0.0:.6f}",
            "query_coverage": f"{(aligned_query / query_length) if query_length else 0.0:.6f}",
            "best_identity": f"{float(row['best_identity']):.6f}",
            "mean_identity": f"{mean_identity:.6f}",
            "qc_flags": ",".join(sorted(flags)),
        }
    )
    return row


def gzip_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    with src.open("rb") as inp, gzip.open(dst, "wb") as out:
        shutil.copyfileobj(inp, out)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    keep_native = truthy(args.keep_native)

    task, target_meta, ortholog_meta = load_task_context(args.task_dir)
    gene_id = task["gene_id"]
    meta_by_sequence = {row["sequence_id"]: row for row in ortholog_meta}
    target_length = int(target_meta.get("sequence_length") or task.get("target_length") or 0)
    summaries = {
        sequence_id: empty_summary(gene_id, meta, target_length)
        for sequence_id, meta in meta_by_sequence.items()
    }

    failures: list[dict[str, object]] = []
    commands: list[str] = []

    with tempfile.TemporaryDirectory(prefix="nucmer_", dir=args.outdir) as tmp_name:
        work_dir = Path(tmp_name)
        target_fasta, orthologs_fasta = materialize_task_fastas(
            args.source_target_fasta,
            args.source_ortholog_fasta,
            task,
            ortholog_meta,
            work_dir,
        )
        prefix = str(work_dir / "nucmer")
        delta_path = work_dir / "nucmer.delta"
        coords_path = work_dir / "nucmer.coords"
        snps_path = work_dir / "nucmer.snps"

        try:
            commands.append(
                run_command(
                    [
                        args.nucmer_bin,
                        "--prefix",
                        prefix,
                        str(target_fasta),
                        str(orthologs_fasta),
                    ]
                )
            )
            commands.append(
                run_command(
                    [args.show_coords_bin, "-THrcl", str(delta_path)],
                    stdout_path=coords_path,
                )
            )
            commands.append(
                run_command(
                    [args.show_snps_bin, "-THrl", str(delta_path)],
                    stdout_path=snps_path,
                )
            )
            segments, segments_by_query = parse_coords(coords_path, gene_id, target_meta, meta_by_sequence, summaries)
            events = parse_snps(snps_path, gene_id, target_meta, meta_by_sequence, segments_by_query, summaries)
            if keep_native:
                gzip_copy(delta_path, args.outdir / "native" / f"{gene_id}.delta.gz")
                gzip_copy(coords_path, args.outdir / "native" / f"{gene_id}.coords.gz")
                gzip_copy(snps_path, args.outdir / "native" / f"{gene_id}.snps.gz")
        except Exception as exc:
            failures.append(
                {
                    "gene_id": gene_id,
                    "ortholog_gene_id": "",
                    "strategy": "nucmer",
                    "tool": "nucmer",
                    "failure_type": "nucmer_failed",
                    "message": str(exc),
                }
            )
            raise

    summary_rows = [finalize_summary(row) for row in summaries.values()]
    write_tsv_gz(args.outdir / "alignment_segments.tsv.gz", SEGMENT_FIELDS, segments)
    write_tsv_gz(args.outdir / "alignment_events.tsv.gz", EVENT_FIELDS, events)
    write_tsv_gz(args.outdir / "ortholog_alignment_summary.tsv.gz", SUMMARY_FIELDS, summary_rows)
    write_tsv_gz(args.outdir / "failures.tsv.gz", FAILURE_FIELDS, failures)
    feature_coverage_count = None
    if args.target_features:
        feature_coverage_count = summarize_feature_coverage(
            args.target_features,
            args.outdir / "ortholog_alignment_summary.tsv.gz",
            args.outdir / "alignment_segments.tsv.gz",
            args.outdir / "feature_coverage.tsv.gz",
        )
    manifest = {
        "gene_id": gene_id,
        "strategy": "nucmer",
        "tool": "nucmer",
        "commands": commands,
        "segment_count": len(segments),
        "event_count": len(events),
        "feature_coverage_count": feature_coverage_count,
        "ortholog_count": len(ortholog_meta),
        "keep_native": keep_native,
        "filtering": "no global delta-filter; downstream parser evaluates records per ortholog",
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
