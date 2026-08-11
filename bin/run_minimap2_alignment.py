#!/usr/bin/env python3
"""Run minimap2 for one gene task and normalize alignment evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from alignment_task_io import load_task_context, materialize_task_fastas
from feature_coverage import summarize_feature_coverage_rows


CS_OP_RE = re.compile(r"(:\d+|=[A-Za-z]+|\*[A-Za-z][A-Za-z]|[+\-][A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2})")
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
    parser.add_argument("--strategy", required=True)
    parser.add_argument("--preset", choices=["asm10", "asm20"], required=True)
    parser.add_argument("--minimap2-bin", default="minimap2")
    parser.add_argument("--threads", default=1, type=int)
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


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def is_primary(tags: dict[str, str]) -> bool:
    return tags.get("tp", "P") == "P"


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


def paf_record_digest(line: str) -> str:
    return hashlib.sha256(line.encode("utf-8")).hexdigest()[:32]


def stable_paf_record_id(record_digest: str, occurrence: int) -> str:
    return f"paf:{record_digest}:{occurrence}"


def cs_events(
    cs: str,
    target_start0: int,
    record: dict[str, object],
    target_meta: dict[str, str],
    event_id_prefix: str,
) -> list[dict[str, object]]:
    events: list[dict[str, object]] = []
    target_pos = target_start0
    event_index = 1

    for match in CS_OP_RE.finditer(cs):
        op = match.group(0)
        if op.startswith(":"):
            target_pos += int(op[1:])
        elif op.startswith("="):
            target_pos += len(op) - 1
        elif op.startswith("*"):
            ref = op[1].upper()
            alt = op[2].upper()
            genomic_start, genomic_end = genomic_coords(target_meta, target_pos, target_pos + 1)
            events.append(
                {
                    **record,
                    "event_id": f"{event_id_prefix}:{event_index}",
                    "event_type": "snv",
                    "target_start0": target_pos,
                    "target_end0": target_pos + 1,
                    "genomic_accession": target_meta.get("genomic_accession", ""),
                    "genomic_start1": genomic_start,
                    "genomic_end1": genomic_end,
                    "ref": ref,
                    "alt": alt,
                }
            )
            event_index += 1
            target_pos += 1
        elif op.startswith("+"):
            alt = op[1:].upper()
            genomic_start, genomic_end = genomic_coords(target_meta, target_pos, target_pos)
            events.append(
                {
                    **record,
                    "event_id": f"{event_id_prefix}:{event_index}",
                    "event_type": "ins",
                    "target_start0": target_pos,
                    "target_end0": target_pos,
                    "genomic_accession": target_meta.get("genomic_accession", ""),
                    "genomic_start1": genomic_start,
                    "genomic_end1": genomic_end,
                    "ref": "",
                    "alt": alt,
                }
            )
            event_index += 1
        elif op.startswith("-"):
            ref = op[1:].upper()
            genomic_start, genomic_end = genomic_coords(target_meta, target_pos, target_pos + len(ref))
            events.append(
                {
                    **record,
                    "event_id": f"{event_id_prefix}:{event_index}",
                    "event_type": "del",
                    "target_start0": target_pos,
                    "target_end0": target_pos + len(ref),
                    "genomic_accession": target_meta.get("genomic_accession", ""),
                    "genomic_start1": genomic_start,
                    "genomic_end1": genomic_end,
                    "ref": ref,
                    "alt": "",
                }
            )
            event_index += 1
            target_pos += len(ref)
        elif op.startswith("~"):
            intron_len = int(re.search(r"\d+", op).group(0))
            target_pos += intron_len

    return events


def run_minimap2(
    minimap2_bin: str,
    preset: str,
    target_fa: Path,
    query_fa: Path,
    paf_path: Path,
    threads: int,
) -> str:
    cmd = [
        minimap2_bin,
        "-t",
        str(threads),
        "-x",
        preset,
        "-c",
        "--cs=long",
        "--paf-no-hit",
        str(target_fa),
        str(query_fa),
    ]
    with paf_path.open("w") as handle:
        result = subprocess.run(cmd, text=True, stdout=handle, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(f"minimap2 failed for preset={preset}: {result.stderr.strip()}")
    return " ".join(cmd)


def parse_paf(
    paf_path: Path,
    gene_id: str,
    strategy: str,
    preset: str,
    target_meta: dict[str, str],
    meta_by_sequence: dict[str, dict[str, str]],
    summaries: dict[str, dict[str, object]],
    event_start_index: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], int]:
    segments: list[dict[str, object]] = []
    event_by_key: dict[tuple[object, ...], dict[str, object]] = {}
    record_occurrences: dict[str, int] = defaultdict(int)

    with paf_path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            qname = fields[0]
            meta = meta_by_sequence.get(qname)
            if meta is None:
                continue
            summary = summaries[qname]

            if fields[5] == "*" or fields[4] == "*":
                summary["status"] = "no_alignment"
                summary["qc_flags"].add("no_alignment")
                continue

            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            strand = fields[4]
            target_id = fields[5]
            target_length = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])
            matches = int(fields[9])
            block_length = int(fields[10])
            mapq = int(fields[11])
            tags = parse_tags(fields)
            primary = is_primary(tags)
            record_digest = paf_record_digest(line)
            record_occurrences[record_digest] += 1
            native_record_id = stable_paf_record_id(record_digest, record_occurrences[record_digest])
            identity = matches / block_length if block_length else 0.0
            flags = []
            if not primary:
                flags.append("non_primary")
            if mapq < 10:
                flags.append("low_mapq")

            segment = {
                "gene_id": gene_id,
                "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
                "tax_id": meta.get("tax_id", ""),
                "taxname": meta.get("taxname", ""),
                "strategy": strategy,
                "tool": "minimap2",
                "preset": preset,
                "sequence_id": qname,
                "target_id": target_id,
                "query_id": qname,
                "target_start0": target_start,
                "target_end0": target_end,
                "query_start0": qstart,
                "query_end0": qend,
                "strand": strand,
                "matches": matches,
                "block_length": block_length,
                "identity": f"{identity:.6f}",
                "mapq": mapq,
                "is_primary": str(primary).lower(),
                "divergence": tags.get("dv", ""),
                "gap_compressed_divergence": tags.get("de", ""),
                "native_record_id": native_record_id,
                "qc_flags": ",".join(flags),
            }
            segments.append(segment)

            summary["status"] = "aligned"
            summary["segment_count"] += 1
            summary["target_length"] = target_length
            summary["query_length"] = qlen
            summary["identities"].append(identity)
            summary["best_identity"] = max(summary["best_identity"], identity)
            if primary:
                summary["primary_segment_count"] += 1
                summary["target_intervals"].append((target_start, target_end))
                summary["query_intervals"].append((qstart, qend))
            else:
                summary["secondary_segment_count"] += 1
                summary["qc_flags"].add("has_secondary")

            cs = tags.get("cs")
            if cs:
                event_record = {
                    "gene_id": gene_id,
                    "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
                    "tax_id": meta.get("tax_id", ""),
                    "taxname": meta.get("taxname", ""),
                    "strategy": strategy,
                    "tool": "minimap2",
                    "preset": preset,
                    "query_id": qname,
                    "strand": strand,
                    "native_record_id": native_record_id,
                    "qc_flags": ",".join(flags),
                }
                new_events = cs_events(
                    cs,
                    target_start,
                    event_record,
                    target_meta,
                    f"{strategy}:{native_record_id}",
                )
                for event in new_events:
                    event["_is_primary"] = primary
                    key = (
                        event["query_id"],
                        event["event_type"],
                        event["target_start0"],
                        event["target_end0"],
                        event["ref"],
                        event["alt"],
                    )
                    current = event_by_key.get(key)
                    candidate_rank = (not primary, str(event["native_record_id"]))
                    current_rank = (
                        not bool(current["_is_primary"]),
                        str(current["native_record_id"]),
                    ) if current is not None else None
                    if current_rank is None or candidate_rank < current_rank:
                        event_by_key[key] = event

    events = sorted(
        event_by_key.values(),
        key=lambda row: (
            str(row["query_id"]),
            int(row["target_start0"]),
            str(row["event_type"]),
            str(row["ref"]),
            str(row["alt"]),
        ),
    )
    for event in events:
        event.pop("_is_primary", None)
        summaries[str(event["query_id"])]["event_count"] += 1

    return segments, events, event_start_index + len(events)


def empty_summary(gene_id: str, strategy: str, preset: str, meta: dict[str, str], target_length: int) -> dict[str, object]:
    return {
        "gene_id": gene_id,
        "ortholog_gene_id": meta.get("ortholog_gene_id", ""),
        "tax_id": meta.get("tax_id", ""),
        "taxname": meta.get("taxname", ""),
        "strategy": strategy,
        "tool": "minimap2",
        "preset": preset,
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
        "qc_flags": set(),
    }


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
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")
    args.outdir.mkdir(parents=True, exist_ok=True)
    keep_native = truthy(args.keep_native)

    task, target_meta, ortholog_meta = load_task_context(args.task_dir)
    gene_id = task["gene_id"]
    meta_by_sequence = {row["sequence_id"]: row for row in ortholog_meta}
    target_length = int(target_meta["sequence_length"])

    summaries = {
        sequence_id: empty_summary(
            gene_id,
            args.strategy,
            args.preset,
            meta,
            target_length,
        )
        for sequence_id, meta in meta_by_sequence.items()
    }

    all_segments: list[dict[str, object]] = []
    all_events: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    commands: list[str] = []
    event_index = 1

    with tempfile.TemporaryDirectory(prefix=f"{args.strategy}_", dir=args.outdir) as tmp_name:
        work_dir = Path(tmp_name)
        target_fasta, orthologs_fasta = materialize_task_fastas(
            args.source_target_fasta,
            args.source_ortholog_fasta,
            task,
            ortholog_meta,
            work_dir,
        )
        paf_path = Path(f"{args.strategy}.{args.preset}.paf")
        try:
            command = run_minimap2(
                args.minimap2_bin,
                args.preset,
                target_fasta,
                orthologs_fasta,
                paf_path,
                args.threads,
            )
            commands.append(command)
            segments, events, event_index = parse_paf(
                paf_path,
                gene_id,
                args.strategy,
                args.preset,
                target_meta,
                meta_by_sequence,
                summaries,
                event_index,
            )
            all_segments.extend(segments)
            all_events.extend(events)
            if keep_native:
                gzip_copy(
                    paf_path,
                    args.outdir / "native" / f"{gene_id}.{args.preset}.paf.gz",
                )
        except Exception as exc:
            failures.append(
                {
                    "gene_id": gene_id,
                    "ortholog_gene_id": "",
                    "strategy": args.strategy,
                    "tool": "minimap2",
                    "failure_type": "minimap2_failed",
                    "message": str(exc),
                }
            )
            raise
        finally:
            if paf_path.exists() and not keep_native:
                paf_path.unlink()

    summary_rows = [finalize_summary(row) for row in summaries.values()]
    write_tsv_gz(args.outdir / "alignment_segments.tsv.gz", SEGMENT_FIELDS, all_segments)
    write_tsv_gz(args.outdir / "alignment_events.tsv.gz", EVENT_FIELDS, all_events)
    write_tsv_gz(args.outdir / "ortholog_alignment_summary.tsv.gz", SUMMARY_FIELDS, summary_rows)
    write_tsv_gz(args.outdir / "failures.tsv.gz", FAILURE_FIELDS, failures)
    feature_coverage_count = None
    if args.target_features:
        feature_coverage_count = summarize_feature_coverage_rows(
            args.target_features,
            summary_rows,
            all_segments,
            args.outdir / "feature_coverage.tsv.gz",
        )
    manifest = {
        "gene_id": gene_id,
        "strategy": args.strategy,
        "tool": "minimap2",
        "preset": args.preset,
        "commands": commands,
        "segment_count": len(all_segments),
        "event_count": len(all_events),
        "feature_coverage_count": feature_coverage_count,
        "ortholog_count": len(ortholog_meta),
        "keep_native": keep_native,
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
