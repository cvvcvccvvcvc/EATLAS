#!/usr/bin/env python3
"""Normalize one Ensembl Compara MAF source chunk for all overlapping genes."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import sys
import time
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path

from run_ensembl_compara_maf_alignment import (
    EVENT_FIELDS,
    FAILURE_FIELDS,
    OUTPUT_GZIP_COMPRESSLEVEL,
    SEGMENT_FIELDS,
    SUMMARY_FIELDS,
    TOOL_NAME,
    convert_pair,
    empty_summary,
    finalize_summary,
    is_ancestral,
    iter_maf_blocks,
    maf_source_name,
    open_maf_text,
    overlaps,
    retry_sleep_seconds,
    retryable_maf_error,
    resolve_maf_dots,
    source_read_failure,
    to_alignment_row,
    write_tsv_gz,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-task-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--strategy", default="precomputed_ensembl_92_mammals_epo_extended")
    parser.add_argument("--release", default="116")
    parser.add_argument("--species-set", default="92_mammals.epo_extended")
    parser.add_argument("--method", default="EPO_EXTENDED")
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--retry-base-seconds", type=float, default=2.0)
    parser.add_argument("--retry-max-seconds", type=float, default=30.0)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


@dataclass(frozen=True)
class GeneIntervalIndex:
    genes_by_human_src: dict[str, list[dict[str, str]]]
    starts_by_human_src: dict[str, list[int]]

    @classmethod
    def build(cls, genes: list[dict[str, str]]) -> "GeneIntervalIndex":
        genes_by_human_src: dict[str, list[dict[str, str]]] = {}
        starts_by_human_src: dict[str, list[int]] = {}
        for gene in genes:
            genes_by_human_src.setdefault(gene["human_src"], []).append(gene)
        for human_src, rows in genes_by_human_src.items():
            rows.sort(key=lambda row: int(row["target_origin1"]))
            starts_by_human_src[human_src] = [int(row["target_origin1"]) for row in rows]
        return cls(genes_by_human_src, starts_by_human_src)

    def overlapping(self, human_src: str, start1: int, end1: int) -> list[dict[str, str]]:
        genes = self.genes_by_human_src.get(human_src, [])
        if not genes:
            return []
        starts = self.starts_by_human_src[human_src]
        limit = bisect.bisect_right(starts, end1)
        return [gene for gene in genes[:limit] if int(gene["target_end1"]) >= start1]


class PartitionedTsvGzWriter:
    """Route rows to bounded-open gzip writers keyed by gene_id."""

    def __init__(
        self,
        root: Path,
        filename: str,
        fields: list[str],
        chunk_id: str,
        max_open: int = 32,
    ) -> None:
        self.root = root
        self.filename = filename
        self.fields = fields
        self.chunk_id = chunk_id
        self.max_open = max_open
        self.handles: OrderedDict[str, tuple[object, csv.DictWriter]] = OrderedDict()
        self.created: set[str] = set()
        self.counts: dict[str, int] = defaultdict(int)
        self.count = 0

    def path_for(self, gene_id: str) -> Path:
        return self.root / fragment_dir_name(gene_id, self.chunk_id) / self.filename

    def _writer(self, gene_id: str) -> csv.DictWriter:
        if gene_id in self.handles:
            handle, writer = self.handles.pop(gene_id)
            self.handles[gene_id] = (handle, writer)
            return writer
        if len(self.handles) >= self.max_open:
            _, (handle, _) = self.handles.popitem(last=False)
            handle.close()
        path = self.path_for(gene_id)
        path.parent.mkdir(parents=True, exist_ok=True)
        mode = "at" if gene_id in self.created else "wt"
        handle = gzip.open(path, mode, newline="", compresslevel=OUTPUT_GZIP_COMPRESSLEVEL)
        writer = csv.DictWriter(handle, fieldnames=self.fields, delimiter="\t", extrasaction="ignore")
        if gene_id not in self.created:
            writer.writeheader()
            self.created.add(gene_id)
        self.handles[gene_id] = (handle, writer)
        return writer

    def write(self, row: dict[str, object]) -> None:
        gene_id = str(row.get("gene_id") or "")
        if not gene_id:
            raise ValueError(f"Cannot partition {self.filename} row without gene_id")
        writer = self._writer(gene_id)
        writer.writerow({field: row.get(field, "") for field in self.fields})
        self.counts[gene_id] += 1
        self.count += 1

    def ensure_gene(self, gene_id: str) -> None:
        self._writer(gene_id)

    def close(self) -> None:
        for handle, _ in self.handles.values():
            handle.close()
        self.handles.clear()


def write_gene_fragment_outputs(
    args: argparse.Namespace,
    chunk_manifest: dict,
    gene_ids: list[str],
    summary_rows: list[dict[str, object]],
    failures: list[dict[str, object]],
    segment_writer: PartitionedTsvGzWriter,
    event_writer: PartitionedTsvGzWriter,
) -> None:
    summaries_by_gene: dict[str, list[dict[str, object]]] = defaultdict(list)
    failures_by_gene: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in summary_rows:
        summaries_by_gene[str(row.get("gene_id") or "")].append(row)
    for row in failures:
        failures_by_gene[str(row.get("gene_id") or "")].append(row)

    result_root = args.outdir / "gene_results"
    for gene_id in gene_ids:
        segment_writer.ensure_gene(gene_id)
        event_writer.ensure_gene(gene_id)
    segment_writer.close()
    event_writer.close()

    for gene_id in gene_ids:
        gene_dir = result_root / fragment_dir_name(gene_id, str(chunk_manifest["chunk_id"]))
        gene_summaries = summaries_by_gene.get(gene_id, [])
        gene_failures = failures_by_gene.get(gene_id, [])
        write_tsv_gz(gene_dir / "ortholog_alignment_summary.tsv.gz", SUMMARY_FIELDS, gene_summaries)
        write_tsv_gz(gene_dir / "failures.tsv.gz", FAILURE_FIELDS, gene_failures)
        manifest = {
            "task_type": "maf_gene_fragment",
            "chunk_id": chunk_manifest["chunk_id"],
            "gene_id": gene_id,
            "gene_ids": [gene_id],
            "strategy": args.strategy,
            "strategies": [args.strategy],
            "tool": TOOL_NAME,
            "release": args.release,
            "species_set": args.species_set,
            "method": args.method,
            "source": chunk_manifest.get("source", ""),
            "seq_region": chunk_manifest.get("seq_region", ""),
            "chunk_order": chunk_manifest.get("chunk_order", ""),
            "output_gzip_compresslevel": OUTPUT_GZIP_COMPRESSLEVEL,
            "summary_count": len(gene_summaries),
            "segment_count": segment_writer.counts.get(gene_id, 0),
            "event_count": event_writer.counts.get(gene_id, 0),
            "failure_count": len(gene_failures),
        }
        (gene_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def fragment_dir_name(gene_id: str, chunk_id: str) -> str:
    return f"gene_{gene_id}__{chunk_id}"


def scan_chunk_source(
    args: argparse.Namespace,
    source: str,
    gene_index: GeneIntervalIndex,
    gene_ids: list[str],
    summaries: dict[tuple[str, str], dict[str, object]],
    segment_writer: TsvGzWriter,
    event_writer: TsvGzWriter,
) -> tuple[int, int, int, list[dict[str, object]]]:
    source_name = maf_source_name(source)
    completed_block_count = 0
    used_block_count = 0
    alignment_row_count = 0
    event_id = 1
    failures: list[dict[str, object]] = []
    attempts = max(int(args.retries), 1)
    last_error: Exception | None = None

    for attempt in range(1, attempts + 1):
        current_block_count = 0
        attempt_error: Exception | None = None
        try:
            handle = open_maf_text(source, args.timeout)
        except Exception as exc:
            if not retryable_maf_error(exc):
                raise
            attempt_error = exc
        else:
            with handle:
                block_iter = iter_maf_blocks(handle)
                while True:
                    try:
                        block = next(block_iter)
                    except StopIteration:
                        return event_id, used_block_count, alignment_row_count, failures
                    except Exception as exc:
                        if not retryable_maf_error(exc):
                            raise
                        attempt_error = exc
                        break

                    current_block_count += 1
                    if current_block_count <= completed_block_count:
                        continue

                    block_used = False
                    human_rows = [row for row in block if row.src in gene_index.genes_by_human_src]
                    for human_maf in human_rows:
                        human_start0, human_end0 = human_maf.forward_interval0()
                        genes = gene_index.overlapping(human_maf.src, human_start0 + 1, human_end0)
                        if not genes:
                            continue
                        flip_orientation = human_maf.strand == "-"
                        human_row = to_alignment_row(human_maf, flip_orientation)
                        for gene in genes:
                            target_origin1 = int(gene["target_origin1"])
                            target_end1 = int(gene["target_end1"])
                            if not overlaps(human_start0 + 1, human_end0, target_origin1, target_end1):
                                continue
                            gene_id = str(gene["gene_id"])
                            for query_index, maf_row in enumerate(block, start=1):
                                if maf_row.src == human_maf.src:
                                    continue
                                query_row = to_alignment_row(maf_row, flip_orientation)
                                if is_ancestral(query_row):
                                    continue
                                query_row = resolve_maf_dots(human_row, query_row)
                                summary_key = (gene_id, query_row.species)
                                summary = summaries.setdefault(
                                    summary_key,
                                    empty_summary(args, query_row, int(gene["target_length"])),
                                )
                                summary["gene_id"] = gene_id
                                native_record_id = f"{source_name}:block{current_block_count}:row{query_index}"
                                event_id = convert_pair(
                                    args,
                                    gene_id,
                                    str(gene["genomic_accession"]),
                                    target_origin1,
                                    target_end1,
                                    human_row,
                                    query_row,
                                    native_record_id,
                                    summary,
                                    event_id,
                                    segment_writer,
                                    event_writer,
                                )
                                block_used = True
                                alignment_row_count += 1
                    if block_used:
                        used_block_count += 1
                    completed_block_count = current_block_count

        if attempt_error is None:
            continue
        last_error = attempt_error
        print(
            f"MAF source read attempt {attempt}/{attempts} failed for {source_name} "
            f"after {completed_block_count} committed blocks: "
            f"{type(attempt_error).__name__}: {attempt_error}",
            file=sys.stderr,
        )
        if attempt < attempts:
            time.sleep(retry_sleep_seconds(args, attempt))
            continue
        for gene_id in gene_ids:
            failures.append(
                source_read_failure(
                    args,
                    gene_id,
                    source,
                    attempts,
                    completed_block_count,
                    used_block_count,
                    attempt_error,
                )
            )
        return event_id, used_block_count, alignment_row_count, failures

    for gene_id in gene_ids:
        failures.append(
            source_read_failure(args, gene_id, source, attempts, completed_block_count, used_block_count, last_error)
        )
    return event_id, used_block_count, alignment_row_count, failures


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    chunk_manifest = json.loads((args.chunk_task_dir / "chunk.json").read_text())
    task_type = chunk_manifest.get("task_type", "maf_chunk")
    genes_path = args.chunk_task_dir / str(chunk_manifest["genes_tsv"])
    result_root = args.outdir / "gene_results"
    chunk_id = str(chunk_manifest["chunk_id"])
    segment_writer = PartitionedTsvGzWriter(
        result_root,
        "alignment_segments.tsv.gz",
        SEGMENT_FIELDS,
        chunk_id,
    )
    event_writer = PartitionedTsvGzWriter(
        result_root,
        "alignment_events.tsv.gz",
        EVENT_FIELDS,
        chunk_id,
    )

    if task_type == "failures":
        failures = [
            {
                "gene_id": row.get("gene_id", ""),
                "ortholog_gene_id": "",
                "strategy": args.strategy,
                "tool": TOOL_NAME,
                "failure_type": row.get("failure_type", ""),
                "message": row.get("message", ""),
            }
            for row in read_tsv(genes_path)
        ]
        gene_ids = sorted(
            {str(row["gene_id"]) for row in failures if row.get("gene_id")},
            key=lambda value: int(value) if value.isdigit() else value,
        )
        write_gene_fragment_outputs(
            args,
            chunk_manifest,
            gene_ids,
            [],
            failures,
            segment_writer,
            event_writer,
        )
        manifest = {
            "chunk_id": chunk_manifest["chunk_id"],
            "strategy": args.strategy,
            "tool": TOOL_NAME,
            "task_type": task_type,
            "output_gzip_compresslevel": OUTPUT_GZIP_COMPRESSLEVEL,
            "gene_count": len(gene_ids),
            "gene_ids": gene_ids,
            "summary_count": 0,
            "segment_count": 0,
            "event_count": 0,
            "failure_count": len(failures),
        }
        (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        return

    genes = read_tsv(genes_path)
    gene_index = GeneIntervalIndex.build(genes)
    gene_ids = sorted({row["gene_id"] for row in genes}, key=lambda value: int(value) if value.isdigit() else value)
    summaries: dict[tuple[str, str], dict[str, object]] = {}
    failures: list[dict[str, object]] = []
    used_block_count = 0
    alignment_row_count = 0

    try:
        _, used_block_count, alignment_row_count, failures = scan_chunk_source(
            args,
            str(chunk_manifest["source"]),
            gene_index,
            gene_ids,
            summaries,
            segment_writer,
            event_writer,
        )
    except Exception as exc:
        failures.extend(
            {
                "gene_id": gene_id,
                "ortholog_gene_id": "",
                "strategy": args.strategy,
                "tool": TOOL_NAME,
                "failure_type": "ensembl_compara_maf_failed",
                "message": str(exc),
            }
            for gene_id in gene_ids
        )
        raise
    finally:
        segment_writer.close()
        event_writer.close()

    if failures:
        for summary in summaries.values():
            summary["qc_flags"].add("maf_source_read_failed")

    summary_rows = [
        finalize_summary(row)
        for _, row in sorted(
            summaries.items(),
            key=lambda item: (
                int(item[0][0]) if item[0][0].isdigit() else item[0][0],
                item[0][1],
            ),
        )
    ]
    write_gene_fragment_outputs(
        args,
        chunk_manifest,
        gene_ids,
        summary_rows,
        failures,
        segment_writer,
        event_writer,
    )

    manifest = {
        "chunk_id": chunk_manifest["chunk_id"],
        "strategy": args.strategy,
        "tool": TOOL_NAME,
        "release": args.release,
        "species_set": args.species_set,
        "method": args.method,
        "task_type": task_type,
        "output_gzip_compresslevel": OUTPUT_GZIP_COMPRESSLEVEL,
        "source": chunk_manifest["source"],
        "seq_region": chunk_manifest.get("seq_region", ""),
        "chunk_order": chunk_manifest.get("chunk_order", ""),
        "gene_count": len(gene_ids),
        "gene_ids": gene_ids,
        "used_block_count": used_block_count,
        "alignment_row_count": alignment_row_count,
        "summary_count": len(summary_rows),
        "segment_count": segment_writer.count,
        "event_count": event_writer.count,
        "failure_count": len(failures),
    }
    (args.outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    if failures:
        print(f"{TOOL_NAME} completed with failures for chunk_id={chunk_manifest['chunk_id']}", file=sys.stderr)


if __name__ == "__main__":
    main()
