#!/usr/bin/env python3
"""Convert Ensembl Compara REST MSA JSON into GAPH-like Stage 2 tables."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


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

DNA_BASES = {"A", "C", "G", "T"}
DEFAULT_TOOL_NAME = "ensembl_compara_rest"


@dataclass(frozen=True)
class Cursor:
    start1: int
    end1: int
    strand: int

    def first(self) -> int:
        return self.start1 if self.strand >= 0 else self.end1

    def advance(self, pos: int) -> int:
        return pos + 1 if self.strand >= 0 else pos - 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--strategy", default="ensembl_compara_rest_epo")
    parser.add_argument("--tool-name", default=DEFAULT_TOOL_NAME)
    parser.add_argument("--method", default="EPO")
    parser.add_argument("--species-set-group", default="mammals")
    parser.add_argument("--human-species", default="homo_sapiens")
    parser.add_argument(
        "--target-origin1",
        type=int,
        help="One-based genomic origin for target_start0. Defaults to the first human base in the payload.",
    )
    parser.add_argument(
        "--genomic-accession",
        help="Override human seq_region in event genomic_accession, useful for RefSeq accession normalization.",
    )
    parser.add_argument("--include-ancestral", action="store_true")
    return parser.parse_args()


def write_tsv_gz(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
            count += 1
    return count


def alignment_blocks(payload: object) -> list[dict[str, object]]:
    if isinstance(payload, dict) and "alignments" in payload:
        return [payload]
    if isinstance(payload, list):
        blocks = [item for item in payload if isinstance(item, dict) and "alignments" in item]
        if blocks:
            return blocks
    raise ValueError("Input JSON does not look like an Ensembl alignment/region response")


def row_species(row: dict[str, object]) -> str:
    return str(row.get("species") or "")


def is_ancestral(row: dict[str, object]) -> bool:
    species = row_species(row)
    return (
        "[" in species
        or "]" in species
        or "-" in species
        or species.startswith("ancestral")
        or species == "ancestral_sequences"
    )


def row_cursor(row: dict[str, object]) -> Cursor:
    start = int(row["start"])
    end = int(row["end"])
    strand = int(row.get("strand") or 1)
    return Cursor(start, end, strand)


def query_id(row: dict[str, object]) -> str:
    species = row_species(row)
    seq_region = row.get("seq_region") or ""
    strand = row.get("strand") or ""
    return f"{species}:{seq_region}:{row.get('start')}:{row.get('end')}:{strand}"


def genomic_to_target0(genomic_pos1: int, origin1: int) -> int:
    return genomic_pos1 - origin1


def target0_to_genomic(target0: int, origin1: int) -> int:
    return origin1 + target0


def interval_union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def empty_summary(args: argparse.Namespace, row: dict[str, object], target_length: int) -> dict[str, object]:
    query_sequence = str(row.get("seq") or "")
    return {
        "gene_id": args.gene_id,
        "ortholog_gene_id": row_species(row),
        "tax_id": "",
        "taxname": row_species(row),
        "strategy": args.strategy,
        "tool": args.tool_name,
        "preset": f"{args.method}:{args.species_set_group}",
        "status": "not_run",
        "target_length": target_length,
        "query_length": sum(1 for char in query_sequence if char != "-"),
        "segment_count": 0,
        "primary_segment_count": 0,
        "secondary_segment_count": 0,
        "target_intervals": [],
        "query_intervals": [],
        "identities": [],
        "best_identity": 0.0,
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
            "target_coverage": f"{aligned_target / target_length if target_length else 0.0:.6f}",
            "query_coverage": f"{aligned_query / query_length if query_length else 0.0:.6f}",
            "best_identity": f"{float(row['best_identity']):.6f}",
            "mean_identity": f"{mean_identity:.6f}",
            "qc_flags": ",".join(sorted(flags)),
        }
    )
    return row


def append_segment(
    segments: list[dict[str, object]],
    summary: dict[str, object],
    args: argparse.Namespace,
    human_row: dict[str, object],
    query_row: dict[str, object],
    native_record_id: str,
    target_start0: int,
    target_end0: int,
    query_positions: list[int],
    matches: int,
    block_length: int,
    qc_flags: set[str],
) -> None:
    if target_end0 <= target_start0 or block_length <= 0:
        return
    query_start0 = min(query_positions) - 1 if query_positions else ""
    query_end0 = max(query_positions) if query_positions else ""
    identity = matches / block_length if block_length else 0.0
    sequence_id = query_id(query_row)
    row = {
        "gene_id": args.gene_id,
        "ortholog_gene_id": row_species(query_row),
        "tax_id": "",
        "taxname": row_species(query_row),
        "strategy": args.strategy,
        "tool": args.tool_name,
        "preset": f"{args.method}:{args.species_set_group}",
        "sequence_id": sequence_id,
        "target_id": str(human_row.get("seq_region") or ""),
        "query_id": sequence_id,
        "target_start0": target_start0,
        "target_end0": target_end0,
        "query_start0": query_start0,
        "query_end0": query_end0,
        "strand": str(query_row.get("strand") or ""),
        "matches": matches,
        "block_length": block_length,
        "identity": f"{identity:.6f}",
        "mapq": "",
        "is_primary": "true",
        "divergence": "",
        "gap_compressed_divergence": "",
        "native_record_id": native_record_id,
        "qc_flags": ",".join(sorted(qc_flags)),
    }
    segments.append(row)
    summary["status"] = "aligned"
    summary["segment_count"] += 1
    summary["primary_segment_count"] += 1
    summary["target_intervals"].append((target_start0, target_end0))
    if query_positions:
        summary["query_intervals"].append((min(query_positions) - 1, max(query_positions)))
    summary["identities"].append(identity)
    summary["best_identity"] = max(float(summary["best_identity"]), identity)


def append_event(
    events: list[dict[str, object]],
    event_id: int,
    args: argparse.Namespace,
    human_row: dict[str, object],
    query_row: dict[str, object],
    event_type: str,
    target_start0: int,
    target_end0: int,
    ref: str,
    alt: str,
    native_record_id: str,
    qc_flags: set[str],
) -> None:
    genomic_accession = args.genomic_accession or str(human_row.get("seq_region") or "")
    if target_end0 > target_start0:
        genomic_start1 = target0_to_genomic(target_start0, args.target_origin1)
        genomic_end1 = target0_to_genomic(target_end0, args.target_origin1) - 1
    else:
        genomic_start1 = target0_to_genomic(target_start0, args.target_origin1)
        genomic_end1 = genomic_start1
    events.append(
        {
            "gene_id": args.gene_id,
            "ortholog_gene_id": row_species(query_row),
            "tax_id": "",
            "taxname": row_species(query_row),
            "strategy": args.strategy,
            "tool": args.tool_name,
            "preset": f"{args.method}:{args.species_set_group}",
            "event_id": event_id,
            "event_type": event_type,
            "target_start0": target_start0,
            "target_end0": target_end0,
            "genomic_accession": genomic_accession,
            "genomic_start1": genomic_start1,
            "genomic_end1": genomic_end1,
            "ref": ref,
            "alt": alt,
            "query_id": query_id(query_row),
            "strand": str(query_row.get("strand") or ""),
            "native_record_id": native_record_id,
            "qc_flags": ",".join(sorted(qc_flags)),
        }
    )


def emit_pending_indel(
    pending: dict[str, object] | None,
    events: list[dict[str, object]],
    event_id: int,
    args: argparse.Namespace,
    human_row: dict[str, object],
    query_row: dict[str, object],
    native_record_id: str,
    summary: dict[str, object],
) -> int:
    if not pending:
        return event_id
    event_type = str(pending["event_type"])
    ref = str(pending.get("ref") or "")
    alt = str(pending.get("alt") or "")
    append_event(
        events,
        event_id,
        args,
        human_row,
        query_row,
        event_type,
        int(pending["target_start0"]),
        int(pending["target_end0"]),
        ref,
        alt,
        native_record_id,
        set(pending.get("qc_flags") or set()),
    )
    summary["event_count"] += 1
    return event_id + 1


def convert_pair(
    args: argparse.Namespace,
    human_row: dict[str, object],
    query_row: dict[str, object],
    block_index: int,
    query_index: int,
    summary: dict[str, object],
    event_id: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], int]:
    human_seq = str(human_row.get("seq") or "").upper()
    query_seq = str(query_row.get("seq") or "").upper()
    if len(human_seq) != len(query_seq):
        raise ValueError(f"MSA row length mismatch for {row_species(query_row)}")

    human_cursor = row_cursor(human_row)
    query_cursor = row_cursor(query_row)
    human_pos = human_cursor.first()
    query_pos = query_cursor.first()
    next_target0 = genomic_to_target0(human_pos, args.target_origin1)

    segments: list[dict[str, object]] = []
    events: list[dict[str, object]] = []
    native_record_id = f"block{block_index}:row{query_index}"

    active_segment: dict[str, object] | None = None
    pending_indel: dict[str, object] | None = None

    def close_segment() -> None:
        nonlocal active_segment
        if not active_segment:
            return
        append_segment(
            segments,
            summary,
            args,
            human_row,
            query_row,
            native_record_id,
            int(active_segment["target_start0"]),
            int(active_segment["target_end0"]),
            list(active_segment["query_positions"]),
            int(active_segment["matches"]),
            int(active_segment["block_length"]),
            set(active_segment["qc_flags"]),
        )
        active_segment = None

    def close_indel() -> None:
        nonlocal pending_indel, event_id
        event_id = emit_pending_indel(
            pending_indel,
            events,
            event_id,
            args,
            human_row,
            query_row,
            native_record_id,
            summary,
        )
        pending_indel = None

    for human_base, query_base in zip(human_seq, query_seq):
        human_has_base = human_base != "-"
        query_has_base = query_base != "-"
        target0 = genomic_to_target0(human_pos, args.target_origin1) if human_has_base else next_target0
        current_query_pos = query_pos if query_has_base else None

        if human_has_base and query_has_base:
            close_indel()
            qc_flags: set[str] = set()
            if human_base not in DNA_BASES or query_base not in DNA_BASES:
                qc_flags.add("ambiguous_base")
            if active_segment is None:
                active_segment = {
                    "target_start0": target0,
                    "target_end0": target0,
                    "query_positions": [],
                    "matches": 0,
                    "block_length": 0,
                    "qc_flags": set(),
                }
            active_segment["target_end0"] = target0 + 1
            active_segment["query_positions"].append(current_query_pos)
            active_segment["block_length"] += 1
            active_segment["qc_flags"].update(qc_flags)
            if human_base == query_base:
                active_segment["matches"] += 1
            elif not qc_flags:
                append_event(
                    events,
                    event_id,
                    args,
                    human_row,
                    query_row,
                    "snv",
                    target0,
                    target0 + 1,
                    human_base,
                    query_base,
                    native_record_id,
                    set(),
                )
                event_id += 1
                summary["event_count"] += 1

        elif human_has_base and not query_has_base:
            close_segment()
            if pending_indel and pending_indel["event_type"] == "del" and pending_indel["target_end0"] == target0:
                pending_indel["target_end0"] = target0 + 1
                pending_indel["ref"] += human_base
            else:
                close_indel()
                pending_indel = {
                    "event_type": "del",
                    "target_start0": target0,
                    "target_end0": target0 + 1,
                    "ref": human_base,
                    "alt": "",
                    "qc_flags": set(),
                }

        elif not human_has_base and query_has_base:
            close_segment()
            if pending_indel and pending_indel["event_type"] == "ins" and pending_indel["target_start0"] == target0:
                pending_indel["alt"] += query_base
            else:
                close_indel()
                pending_indel = {
                    "event_type": "ins",
                    "target_start0": target0,
                    "target_end0": target0,
                    "ref": "",
                    "alt": query_base,
                    "qc_flags": set(),
                }
        else:
            close_segment()
            close_indel()

        if human_has_base:
            human_pos = human_cursor.advance(human_pos)
            next_target0 = genomic_to_target0(human_pos, args.target_origin1)
        if query_has_base:
            query_pos = query_cursor.advance(query_pos)

    close_segment()
    close_indel()
    return segments, events, event_id


def infer_target_origin(blocks: list[dict[str, object]], human_species: str) -> int:
    starts = []
    for block in blocks:
        for row in block.get("alignments", []):
            if isinstance(row, dict) and row_species(row) == human_species:
                starts.append(int(row["start"]))
    if not starts:
        raise ValueError(f"No human row found for species={human_species}")
    return min(starts)


def human_span_length(blocks: list[dict[str, object]], human_species: str, origin1: int) -> int:
    ends = []
    for block in blocks:
        for row in block.get("alignments", []):
            if isinstance(row, dict) and row_species(row) == human_species:
                ends.append(int(row["end"]))
    if not ends:
        return 0
    return max(ends) - origin1 + 1


def main() -> None:
    args = parse_args()
    payload = json.loads(args.input.read_text())
    blocks = alignment_blocks(payload)
    if args.target_origin1 is None:
        args.target_origin1 = infer_target_origin(blocks, args.human_species)
    target_length = human_span_length(blocks, args.human_species, args.target_origin1)

    all_segments: list[dict[str, object]] = []
    all_events: list[dict[str, object]] = []
    summaries: dict[str, dict[str, object]] = {}
    event_id = 1
    skipped_rows = 0

    for block_index, block in enumerate(blocks, start=1):
        rows = [row for row in block.get("alignments", []) if isinstance(row, dict)]
        human_rows = [row for row in rows if row_species(row) == args.human_species]
        if not human_rows:
            raise ValueError(f"Block {block_index} has no {args.human_species} row")
        human_row = human_rows[0]
        for query_index, row in enumerate(rows, start=1):
            if row is human_row:
                continue
            if is_ancestral(row) and not args.include_ancestral:
                skipped_rows += 1
                continue
            key = row_species(row)
            summaries.setdefault(key, empty_summary(args, row, target_length))
            segments, events, event_id = convert_pair(
                args,
                human_row,
                row,
                block_index,
                query_index,
                summaries[key],
                event_id,
            )
            all_segments.extend(segments)
            all_events.extend(events)

    summary_rows = [finalize_summary(row) for row in summaries.values()]
    failures: list[dict[str, object]] = []
    args.outdir.mkdir(parents=True, exist_ok=True)
    segment_count = write_tsv_gz(args.outdir / "alignment_segments.tsv.gz", SEGMENT_FIELDS, all_segments)
    event_count = write_tsv_gz(args.outdir / "alignment_events.tsv.gz", EVENT_FIELDS, all_events)
    summary_count = write_tsv_gz(args.outdir / "ortholog_alignment_summary.tsv.gz", SUMMARY_FIELDS, summary_rows)
    write_tsv_gz(args.outdir / "failures.tsv.gz", FAILURE_FIELDS, failures)
    manifest = {
        "source": str(args.input),
        "strategy": args.strategy,
        "tool": args.tool_name,
        "method": args.method,
        "species_set_group": args.species_set_group,
        "target_origin1": args.target_origin1,
        "target_length": target_length,
        "segment_count": segment_count,
        "event_count": event_count,
        "summary_count": summary_count,
        "skipped_ancestral_rows": skipped_rows,
        "include_ancestral": args.include_ancestral,
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.outdir}")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
