#!/usr/bin/env python3
"""Add bounded Ensembl VEP evidence to one variant-annotation shard."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd

from genomics.variants import parse_variant_key
from genomics.vep.annotator import annotate_vep_consequences, vep_result_cache_config


csv.field_size_limit(sys.maxsize)


SCHEMA = "normalized_vep_annotation_shard_v1"
REQUIRED_INPUT_FIELDS = {"variant_key", "gene_id"}
VEP_FIELDS = [
    "vep_status",
    "vep_primary_consequence",
    "vep_consequence_terms",
    "vep_transcript_id",
    "vep_mane_select",
    "vep_canonical",
    "vep_impact",
    "vep_variant_class",
]
VEP_RENAME = {
    "status": "vep_status",
    "primary_consequence": "vep_primary_consequence",
    "consequence_terms": "vep_consequence_terms",
    "transcript_id": "vep_transcript_id",
    "mane_select": "vep_mane_select",
    "canonical": "vep_canonical",
    "impact": "vep_impact",
    "variant_class": "vep_variant_class",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-tsv", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--partition-id", required=True)
    parser.add_argument("--shard-id", required=True)
    parser.add_argument("--expected-row-count", required=True, type=int)
    parser.add_argument("--vep-backend", choices=("rest", "local"), required=True)
    parser.add_argument("--vep-release")
    parser.add_argument("--vep-executable", default="vep")
    parser.add_argument("--vep-cache-dir", type=Path)
    parser.add_argument("--vep-forks", type=int, default=1)
    parser.add_argument("--rest-workers", type=int, default=2)
    parser.add_argument("--vep-result-cache-dir", type=Path)
    parser.add_argument("--vep-result-cache-tile-size-bp", type=int, default=1_000_000)
    return parser.parse_args()


def vep_requests(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    rows = []
    invalid_count = 0
    for row in frame[["variant_key", "gene_id"]].itertuples(index=False):
        parsed = parse_variant_key(row.variant_key)
        if parsed is None:
            invalid_count += 1
            continue
        chrom, pos, ref, alt = parsed
        rows.append(
            {
                "variant_key": str(row.variant_key),
                "gene_id": str(row.gene_id),
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
            }
        )
    fields = ["variant_key", "gene_id", "chrom", "pos", "ref", "alt"]
    return pd.DataFrame(rows, columns=fields), invalid_count


def merge_vep_annotations(
    frame: pd.DataFrame,
    annotations: pd.DataFrame,
    invalid_key_count: int,
) -> pd.DataFrame:
    renamed = annotations.rename(columns=VEP_RENAME)
    keep = ["variant_key", "gene_id", *VEP_FIELDS]
    enriched = frame.assign(_source_order=range(len(frame))).merge(
        renamed[keep],
        on=["variant_key", "gene_id"],
        how="left",
        sort=False,
        validate="many_to_one",
    )
    enriched = enriched.sort_values("_source_order", kind="stable").drop(
        columns="_source_order"
    )
    missing = enriched["vep_status"].isna()
    if int(missing.sum()) != invalid_key_count:
        raise ValueError("VEP output did not cover every valid variant/gene pair")
    enriched.loc[missing, "vep_status"] = "invalid_variant_key"
    for field in VEP_FIELDS:
        if field == "vep_canonical":
            enriched[field] = enriched[field].map(
                lambda value: bool(value) if pd.notna(value) else False
            )
        else:
            enriched[field] = enriched[field].fillna("").astype(str)
    return enriched


def semantic_config(summary: dict[str, object]) -> dict[str, object]:
    backend = str(summary.get("backend") or "")
    release = str(summary.get("release") or "")
    if backend not in {"rest", "local"} or not release:
        raise ValueError("VEP did not report a complete backend/release contract")
    return {
        **vep_result_cache_config(backend=backend, release=release),
        "backend": backend,
    }


def write_tsv_atomic(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
        frame.to_csv(
            temporary,
            sep="\t",
            index=False,
            lineterminator="\n",
            compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
        )
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main() -> None:
    args = parse_args()
    if args.expected_row_count < 0:
        raise ValueError("Expected VEP shard row count must be >= 0")
    if args.vep_forks < 1 or args.rest_workers < 1:
        raise ValueError("VEP worker counts must be >= 1")
    if args.vep_result_cache_tile_size_bp < 1:
        raise ValueError("VEP result-cache tile size must be >= 1")
    if not args.input_tsv.is_file():
        raise FileNotFoundError(args.input_tsv)

    frame = pd.read_csv(
        args.input_tsv,
        sep="\t",
        compression="gzip",
        dtype=str,
        keep_default_na=False,
    )
    missing = REQUIRED_INPUT_FIELDS - set(frame.columns)
    if missing:
        raise ValueError(
            "Variant annotation shard missing fields: " + ", ".join(sorted(missing))
        )
    conflicts = set(VEP_FIELDS) & set(frame.columns)
    if conflicts:
        raise ValueError(
            "Variant annotation shard already contains VEP fields: "
            + ", ".join(sorted(conflicts))
        )
    if len(frame) != args.expected_row_count:
        raise ValueError(
            "Variant annotation shard row count changed: "
            f"observed={len(frame)}, expected={args.expected_row_count}"
        )

    requests, invalid_key_count = vep_requests(frame)
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="gaph_vep_") as temporary:
        annotations, summary = annotate_vep_consequences(
            requests,
            Path(temporary) / "partition.sqlite",
            backend=args.vep_backend,
            release=args.vep_release,
            max_workers=args.rest_workers,
            vep_executable=args.vep_executable,
            vep_cache_dir=args.vep_cache_dir,
            vep_forks=args.vep_forks,
            vep_result_cache_dir=args.vep_result_cache_dir,
            vep_result_cache_tile_size_bp=args.vep_result_cache_tile_size_bp,
        )
    enriched = merge_vep_annotations(frame, annotations, invalid_key_count)

    args.outdir.mkdir(parents=True, exist_ok=True)
    output = args.outdir / "variant_annotations.tsv.gz"
    write_tsv_atomic(enriched, output)
    status_counts = {
        str(status): int(count)
        for status, count in enriched["vep_status"].value_counts().sort_index().items()
    }
    execution_summary = {
        key: value
        for key, value in summary.items()
        if key not in {"base_url", "cache_path", "vep_cache_dir", "vep_executable"}
    }
    manifest = {
        "stage": "annotation",
        "schema": SCHEMA,
        "partition_id": args.partition_id,
        "shard_id": args.shard_id,
        "row_count": len(enriched),
        "invalid_variant_key_count": invalid_key_count,
        "status_counts": status_counts,
        "input": {
            "name": args.input_tsv.name,
            "size_bytes": args.input_tsv.stat().st_size,
            "fields": list(frame.columns),
        },
        "output": {
            "name": output.name,
            "size_bytes": output.stat().st_size,
            "fields": list(enriched.columns),
        },
        "config": semantic_config(summary),
        "vep": execution_summary,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }
    (args.outdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
