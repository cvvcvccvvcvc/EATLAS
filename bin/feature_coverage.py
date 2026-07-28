"""Feature coverage summarization backed by bedtools."""

from __future__ import annotations

import csv
import gzip
import os
import shutil
import subprocess
import tempfile
from collections.abc import Iterable, Iterator
from pathlib import Path


FEATURE_COVERAGE_FIELDS = [
    "gene_id",
    "strategy",
    "feature_type",
    "feature_id",
    "genomic_accession",
    "genomic_start1",
    "genomic_end1",
    "target_start0",
    "target_end0",
    "length_bp",
    "ortholog_count",
    "orthologs_covered",
    "covered_bases",
    "coverage_breadth",
    "depth_bases",
    "mean_depth",
]

KEY_SEPARATOR = "|"
SORT_MEMORY = "128M"
FEATURE_BED_FIELD_COUNT = 12


def _iter_tsv_gz(path: Path) -> Iterator[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _format_fraction(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return "0.000000"
    return f"{numerator / denominator:.6f}"


def _key_part(value: object, field: str) -> str:
    text = str(value or "")
    if not text:
        raise ValueError(f"Feature coverage row has an empty {field}")
    if KEY_SEPARATOR in text or any(character.isspace() for character in text):
        raise ValueError(f"Feature coverage {field} contains unsupported whitespace or '{KEY_SEPARATOR}': {text!r}")
    return text


def _tsv_value(value: object, field: str) -> str:
    text = str(value or "")
    if "\t" in text or "\n" in text or "\r" in text:
        raise ValueError(f"Feature coverage {field} contains a tab or newline")
    return text


def _group_key(gene_id: str, strategy: str) -> str:
    return f"{gene_id}{KEY_SEPARATOR}{strategy}"


def _ortholog_key(gene_id: str, strategy: str, ortholog_gene_id: str) -> str:
    return f"{gene_id}{KEY_SEPARATOR}{strategy}{KEY_SEPARATOR}{ortholog_gene_id}"


def _run_to_file(command: list[str], output: Path, *, env: dict[str, str]) -> None:
    with output.open("w") as handle:
        result = subprocess.run(
            command,
            text=True,
            stdout=handle,
            stderr=subprocess.PIPE,
            env=env,
        )
    if result.returncode != 0:
        detail = (result.stderr or "").strip()
        raise RuntimeError(f"{command[0]} failed: {detail}")


def _sort_file(
    source: Path,
    output: Path,
    temp_dir: Path,
    keys: list[str],
    *,
    unique: bool = False,
    env: dict[str, str],
) -> None:
    command = ["sort", "-S", SORT_MEMORY, "-T", str(temp_dir)]
    if unique:
        command.append("-u")
    command.extend(keys)
    command.extend(["-o", str(output), str(source)])
    result = subprocess.run(command, text=True, capture_output=True, env=env)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"sort failed: {detail}")


def _write_inputs(
    summary_rows: Iterable[dict[str, object]],
    segment_rows: Iterable[dict[str, object]],
    summaries_raw: Path,
    segments_raw: Path,
) -> dict[str, set[str]]:
    strategies_by_gene: dict[str, set[str]] = {}
    with summaries_raw.open("w") as handle:
        for row in summary_rows:
            gene_id = str(row.get("gene_id") or "")
            strategy = str(row.get("strategy") or "")
            ortholog_gene_id = str(row.get("ortholog_gene_id") or "")
            if not gene_id or not strategy or not ortholog_gene_id:
                continue
            gene_id = _key_part(gene_id, "gene_id")
            strategy = _key_part(strategy, "strategy")
            ortholog_gene_id = _key_part(ortholog_gene_id, "ortholog_gene_id")
            strategies_by_gene.setdefault(gene_id, set()).add(strategy)
            handle.write(f"{gene_id}\t{strategy}\t{ortholog_gene_id}\n")

    with segments_raw.open("w") as handle:
        for row in segment_rows:
            gene_id = str(row.get("gene_id") or "")
            strategy = str(row.get("strategy") or "")
            ortholog_gene_id = str(row.get("ortholog_gene_id") or "")
            if not gene_id or not strategy or not ortholog_gene_id:
                continue
            gene_id = _key_part(gene_id, "gene_id")
            strategy = _key_part(strategy, "strategy")
            ortholog_gene_id = _key_part(ortholog_gene_id, "ortholog_gene_id")
            start0 = int(row.get("target_start0") or 0)
            end0 = int(row.get("target_end0") or 0)
            if end0 <= start0:
                continue
            strategies_by_gene.setdefault(gene_id, set()).add(strategy)
            handle.write(f"{_ortholog_key(gene_id, strategy, ortholog_gene_id)}\t{start0}\t{end0}\n")
    return strategies_by_gene


def _count_orthologs(
    summaries_raw: Path,
    summaries_unique: Path,
    temp_dir: Path,
    env: dict[str, str],
) -> dict[tuple[str, str], int]:
    _sort_file(
        summaries_raw,
        summaries_unique,
        temp_dir,
        ["-k1,1", "-k2,2", "-k3,3"],
        unique=True,
        env=env,
    )
    counts: dict[tuple[str, str], int] = {}
    with summaries_unique.open() as handle:
        for line in handle:
            gene_id, strategy, _ortholog_gene_id = line.rstrip("\n").split("\t")
            key = (gene_id, strategy)
            counts[key] = counts.get(key, 0) + 1
    return counts


def _merge_ortholog_intervals(
    segments_raw: Path,
    segments_sorted: Path,
    merged_raw: Path,
    merged_sorted: Path,
    temp_dir: Path,
    env: dict[str, str],
) -> None:
    _sort_file(
        segments_raw,
        segments_sorted,
        temp_dir,
        ["-k1,1", "-k2,2n", "-k3,3n"],
        env=env,
    )
    merged_by_ortholog = temp_dir / "segments.merged_by_ortholog.bed"
    _run_to_file(
        ["bedtools", "merge", "-i", str(segments_sorted)],
        merged_by_ortholog,
        env=env,
    )
    with merged_by_ortholog.open() as source, merged_raw.open("w") as output:
        for line in source:
            ortholog_group, start0, end0 = line.rstrip("\n").split("\t")
            gene_id, strategy, ortholog_gene_id = ortholog_group.split(KEY_SEPARATOR, 2)
            output.write(
                f"{_group_key(gene_id, strategy)}\t{start0}\t{end0}\t{ortholog_gene_id}\n"
            )
    _sort_file(
        merged_raw,
        merged_sorted,
        temp_dir,
        ["-k1,1", "-k2,2n", "-k3,3n", "-k4,4"],
        env=env,
    )


def _write_sorted_features(
    target_features: Path,
    strategies_by_gene: dict[str, set[str]],
    features_raw: Path,
    features_sorted_rows: Path,
    expanded_features_raw: Path,
    expanded_features_sorted: Path,
    temp_dir: Path,
    env: dict[str, str],
) -> int:
    with features_raw.open("w") as handle:
        for input_order, row in enumerate(_iter_tsv_gz(target_features), start=1):
            gene_id = _key_part(row.get("gene_id"), "gene_id")
            if not gene_id.isdigit():
                raise ValueError(f"Feature coverage requires numeric Entrez gene IDs, got {gene_id!r}")
            start0 = int(row.get("target_start0") or 0)
            end0 = int(row.get("target_end0") or 0)
            if end0 <= start0:
                raise ValueError(
                    f"Target feature {row.get('feature_id', '')!r} has invalid interval {start0}-{end0}"
                )
            fields = [
                gene_id,
                str(start0),
                str(end0),
                str(input_order),
                _tsv_value(row.get("feature_type"), "feature_type"),
                _tsv_value(row.get("feature_id"), "feature_id"),
                _tsv_value(row.get("genomic_accession"), "genomic_accession"),
                _tsv_value(row.get("genomic_start1"), "genomic_start1"),
                _tsv_value(row.get("genomic_end1"), "genomic_end1"),
                str(int(row.get("length_bp") or 0)),
            ]
            handle.write("\t".join(fields) + "\n")

    _sort_file(
        features_raw,
        features_sorted_rows,
        temp_dir,
        ["-k1,1n", "-k2,2n", "-k4,4n"],
        env=env,
    )

    serial = 0
    current_gene = ""
    gene_features: list[list[str]] = []

    def write_gene(output, gene_id: str, features: list[list[str]]) -> None:
        nonlocal serial
        for strategy in sorted(strategies_by_gene.get(gene_id, set())):
            for feature in features:
                serial += 1
                (
                    _gene_id,
                    start0,
                    end0,
                    _input_order,
                    feature_type,
                    feature_id,
                    genomic_accession,
                    genomic_start1,
                    genomic_end1,
                    length_bp,
                ) = feature
                output.write(
                    "\t".join(
                        [
                            _group_key(gene_id, strategy),
                            start0,
                            end0,
                            str(serial),
                            gene_id,
                            strategy,
                            feature_type,
                            feature_id,
                            genomic_accession,
                            genomic_start1,
                            genomic_end1,
                            length_bp,
                        ]
                    )
                    + "\n"
                )

    with features_sorted_rows.open() as source, expanded_features_raw.open("w") as output:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            gene_id = fields[0]
            if current_gene and gene_id != current_gene:
                write_gene(output, current_gene, gene_features)
                gene_features = []
            current_gene = gene_id
            gene_features.append(fields)
        if current_gene:
            write_gene(output, current_gene, gene_features)

    _sort_file(
        expanded_features_raw,
        expanded_features_sorted,
        temp_dir,
        ["-k1,1", "-k2,2n", "-k3,3n", "-k4,4n"],
        env=env,
    )
    return serial


def _write_coverage_metrics(
    features_bed: Path,
    merged_segments_bed: Path,
    output: Path,
    env: dict[str, str],
) -> None:
    process = subprocess.Popen(
        [
            "bedtools",
            "coverage",
            "-a",
            str(features_bed),
            "-b",
            str(merged_segments_bed),
            "-hist",
            "-sorted",
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
    )
    assert process.stdout is not None
    current_serial = ""
    metadata: list[str] = []
    covered_bases = 0
    depth_bases = 0
    histogram_bases = 0

    def flush(handle) -> None:
        if not current_serial:
            return
        length_bp = int(metadata[7] or 0)
        interval_length = int(metadata[6]) - int(metadata[5])
        if histogram_bases != interval_length:
            raise ValueError(
                f"bedtools histogram length mismatch for feature {metadata[3]!r}: "
                f"{histogram_bases} != {interval_length}"
            )
        handle.write(
            "\t".join(
                [
                    current_serial,
                    *metadata[:5],
                    metadata[5],
                    metadata[6],
                    metadata[7],
                    str(covered_bases),
                    _format_fraction(covered_bases, length_bp),
                    str(depth_bases),
                    _format_fraction(depth_bases, length_bp),
                    metadata[8],
                    metadata[9],
                ]
            )
            + "\n"
        )

    try:
        with output.open("w") as handle:
            for line in process.stdout:
                fields = line.rstrip("\n").split("\t")
                if fields[0] == "all":
                    continue
                if len(fields) != FEATURE_BED_FIELD_COUNT + 4:
                    raise ValueError(f"Unexpected bedtools coverage row with {len(fields)} fields")
                serial = fields[3]
                if serial != current_serial:
                    flush(handle)
                    current_serial = serial
                    metadata = [
                        fields[4],
                        fields[5],
                        fields[6],
                        fields[7],
                        fields[8],
                        fields[1],
                        fields[2],
                        fields[11],
                        fields[9],
                        fields[10],
                    ]
                    covered_bases = 0
                    depth_bases = 0
                    histogram_bases = 0
                depth = int(fields[-4])
                bases = int(fields[-3])
                histogram_bases += bases
                if depth > 0:
                    covered_bases += bases
                    depth_bases += depth * bases
            flush(handle)
    except Exception:
        process.kill()
        process.wait()
        raise

    stderr = process.stderr.read() if process.stderr is not None else ""
    if process.wait() != 0:
        raise RuntimeError(f"bedtools coverage failed: {stderr.strip()}")


def _write_feature_ortholog_pairs(
    features_bed: Path,
    merged_segments_bed: Path,
    output: Path,
    env: dict[str, str],
) -> None:
    process = subprocess.Popen(
        [
            "bedtools",
            "intersect",
            "-a",
            str(features_bed),
            "-b",
            str(merged_segments_bed),
            "-wa",
            "-wb",
            "-sorted",
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
    )
    assert process.stdout is not None
    try:
        with output.open("w") as handle:
            for line in process.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != FEATURE_BED_FIELD_COUNT + 4:
                    raise ValueError(f"Unexpected bedtools intersect row with {len(fields)} fields")
                handle.write(f"{fields[3]}\t{fields[-1]}\n")
    except Exception:
        process.kill()
        process.wait()
        raise

    stderr = process.stderr.read() if process.stderr is not None else ""
    if process.wait() != 0:
        raise RuntimeError(f"bedtools intersect failed: {stderr.strip()}")


def _count_feature_orthologs(unique_pairs: Path, output: Path) -> None:
    current_serial = ""
    count = 0
    with unique_pairs.open() as source, output.open("w") as destination:
        for line in source:
            serial, _ortholog_gene_id = line.rstrip("\n").split("\t", 1)
            if current_serial and serial != current_serial:
                destination.write(f"{current_serial}\t{count}\n")
                count = 0
            current_serial = serial
            count += 1
        if current_serial:
            destination.write(f"{current_serial}\t{count}\n")


def _write_final_output(
    coverage_metrics: Path,
    ortholog_counts: Path,
    orthologs_by_group: dict[tuple[str, str], int],
    output: Path,
) -> int:
    count_handle = ortholog_counts.open()
    count_line = count_handle.readline()
    count_serial = int(count_line.split("\t", 1)[0]) if count_line else None
    row_count = 0
    try:
        with coverage_metrics.open() as source, gzip.open(output, "wt", newline="") as destination:
            writer = csv.DictWriter(destination, fieldnames=FEATURE_COVERAGE_FIELDS, delimiter="\t")
            writer.writeheader()
            for line in source:
                fields = line.rstrip("\n").split("\t")
                serial = int(fields[0])
                while count_serial is not None and count_serial < serial:
                    raise ValueError(f"Ortholog coverage refers to unknown feature serial {count_serial}")
                orthologs_covered = 0
                if count_serial == serial:
                    orthologs_covered = int(count_line.rstrip("\n").split("\t")[1])
                    count_line = count_handle.readline()
                    count_serial = int(count_line.split("\t", 1)[0]) if count_line else None
                (
                    _serial,
                    gene_id,
                    strategy,
                    feature_type,
                    feature_id,
                    genomic_accession,
                    target_start0,
                    target_end0,
                    length_bp,
                    covered_bases,
                    coverage_breadth,
                    depth_bases,
                    mean_depth,
                    genomic_start1,
                    genomic_end1,
                ) = fields
                writer.writerow(
                    {
                        "gene_id": gene_id,
                        "strategy": strategy,
                        "feature_type": feature_type,
                        "feature_id": feature_id,
                        "genomic_accession": genomic_accession,
                        "genomic_start1": genomic_start1,
                        "genomic_end1": genomic_end1,
                        "target_start0": target_start0,
                        "target_end0": target_end0,
                        "length_bp": length_bp,
                        "ortholog_count": orthologs_by_group.get((gene_id, strategy), 0),
                        "orthologs_covered": orthologs_covered,
                        "covered_bases": covered_bases,
                        "coverage_breadth": coverage_breadth,
                        "depth_bases": depth_bases,
                        "mean_depth": mean_depth,
                    }
                )
                row_count += 1
        if count_serial is not None:
            raise ValueError(f"Ortholog coverage refers to unknown feature serial {count_serial}")
    finally:
        count_handle.close()
    return row_count


def summarize_feature_coverage_rows(
    target_features: Path,
    summary_rows: Iterable[dict[str, object]],
    segment_rows: Iterable[dict[str, object]],
    output: Path,
) -> int:
    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools is required for feature coverage summarization")
    output.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["LC_ALL"] = "C"

    with tempfile.TemporaryDirectory(prefix=".feature_coverage_", dir=output.parent) as temp_name:
        temp_dir = Path(temp_name)
        summaries_raw = temp_dir / "summaries.raw.tsv"
        summaries_unique = temp_dir / "summaries.unique.tsv"
        segments_raw = temp_dir / "segments.raw.bed"
        segments_sorted = temp_dir / "segments.sorted.bed"
        merged_raw = temp_dir / "segments.merged.raw.bed"
        merged_sorted = temp_dir / "segments.merged.sorted.bed"
        features_raw = temp_dir / "features.raw.tsv"
        features_sorted_rows = temp_dir / "features.sorted.tsv"
        expanded_features_raw = temp_dir / "features.expanded.raw.bed"
        expanded_features_sorted = temp_dir / "features.expanded.sorted.bed"
        coverage_raw = temp_dir / "coverage.raw.tsv"
        coverage_sorted = temp_dir / "coverage.sorted.tsv"
        pairs_raw = temp_dir / "feature_ortholog_pairs.raw.tsv"
        pairs_unique = temp_dir / "feature_ortholog_pairs.unique.tsv"
        ortholog_counts = temp_dir / "feature_ortholog_counts.tsv"
        temporary_output = temp_dir / "feature_coverage.tsv.gz"

        strategies_by_gene = _write_inputs(
            summary_rows,
            segment_rows,
            summaries_raw,
            segments_raw,
        )
        orthologs_by_group = _count_orthologs(
            summaries_raw,
            summaries_unique,
            temp_dir,
            env,
        )
        _merge_ortholog_intervals(
            segments_raw,
            segments_sorted,
            merged_raw,
            merged_sorted,
            temp_dir,
            env,
        )
        feature_count = _write_sorted_features(
            target_features,
            strategies_by_gene,
            features_raw,
            features_sorted_rows,
            expanded_features_raw,
            expanded_features_sorted,
            temp_dir,
            env,
        )
        if feature_count == 0:
            with gzip.open(temporary_output, "wt", newline="") as handle:
                csv.DictWriter(handle, fieldnames=FEATURE_COVERAGE_FIELDS, delimiter="\t").writeheader()
            temporary_output.replace(output)
            return 0

        _write_coverage_metrics(expanded_features_sorted, merged_sorted, coverage_raw, env)
        _sort_file(
            coverage_raw,
            coverage_sorted,
            temp_dir,
            ["-k1,1n"],
            env=env,
        )
        _write_feature_ortholog_pairs(
            expanded_features_sorted,
            merged_sorted,
            pairs_raw,
            env,
        )
        _sort_file(
            pairs_raw,
            pairs_unique,
            temp_dir,
            ["-k1,1n", "-k2,2"],
            unique=True,
            env=env,
        )
        _count_feature_orthologs(pairs_unique, ortholog_counts)
        row_count = _write_final_output(
            coverage_sorted,
            ortholog_counts,
            orthologs_by_group,
            temporary_output,
        )
        if row_count != feature_count:
            raise ValueError(f"Feature coverage row count mismatch: {row_count} != {feature_count}")
        temporary_output.replace(output)
        return row_count


def summarize_feature_coverage(
    target_features: Path,
    summaries_path: Path,
    segments_path: Path,
    output: Path,
) -> int:
    return summarize_feature_coverage_rows(
        target_features,
        _iter_tsv_gz(summaries_path),
        _iter_tsv_gz(segments_path),
        output,
    )


def site_aligned_ortholog_counts(
    segments_path: Path,
    site_rows: Iterable[dict[str, object]],
    temp_parent: Path,
) -> dict[tuple[str, str, str], int]:
    """Count distinct ortholog intervals covering each variant-strategy SNV site."""
    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools is required for site ortholog-depth calculation")
    env = os.environ.copy()
    env["LC_ALL"] = "C"

    with tempfile.TemporaryDirectory(prefix=".site_ortholog_depth_", dir=temp_parent) as temp_name:
        temp_dir = Path(temp_name)
        segments_raw = temp_dir / "segments.raw.bed"
        segments_sorted = temp_dir / "segments.sorted.bed"
        merged_raw = temp_dir / "segments.merged.raw.bed"
        merged_sorted = temp_dir / "segments.merged.sorted.bed"
        sites_raw = temp_dir / "sites.raw.bed"
        sites_sorted = temp_dir / "sites.sorted.bed"
        coverage = temp_dir / "sites.coverage.tsv"

        with segments_raw.open("w") as handle:
            for row in _iter_tsv_gz(segments_path):
                gene_id = str(row.get("gene_id") or "")
                strategy = str(row.get("strategy") or "")
                ortholog_gene_id = str(row.get("ortholog_gene_id") or "")
                if not gene_id or not strategy or not ortholog_gene_id:
                    continue
                if str(row.get("is_primary") or "").lower() == "false":
                    continue
                start0 = int(row.get("target_start0") or 0)
                end0 = int(row.get("target_end0") or 0)
                if end0 <= start0:
                    continue
                ortholog_key = _ortholog_key(
                    _key_part(gene_id, "gene_id"),
                    _key_part(strategy, "strategy"),
                    _key_part(ortholog_gene_id, "ortholog_gene_id"),
                )
                handle.write(
                    f"{ortholog_key}\t{start0}\t{end0}\n"
                )

        site_count = 0
        with sites_raw.open("w") as handle:
            for row in site_rows:
                gene_id = _key_part(row.get("gene_id"), "gene_id")
                strategy = _key_part(row.get("strategy"), "strategy")
                variant_key = _key_part(row.get("variant_key"), "variant_key")
                start0 = int(row.get("target_start0") or 0)
                handle.write(
                    f"{_group_key(gene_id, strategy)}\t{start0}\t{start0 + 1}\t{variant_key}\n"
                )
                site_count += 1
        if site_count == 0:
            return {}

        _merge_ortholog_intervals(
            segments_raw,
            segments_sorted,
            merged_raw,
            merged_sorted,
            temp_dir,
            env,
        )
        _sort_file(
            sites_raw,
            sites_sorted,
            temp_dir,
            ["-k1,1", "-k2,2n", "-k3,3n", "-k4,4"],
            env=env,
        )
        _run_to_file(
            [
                "bedtools",
                "coverage",
                "-a",
                str(sites_sorted),
                "-b",
                str(merged_sorted),
                "-counts",
                "-sorted",
            ],
            coverage,
            env=env,
        )

        counts: dict[tuple[str, str, str], int] = {}
        with coverage.open() as handle:
            for line in handle:
                group, _start0, _end0, variant_key, depth = line.rstrip("\n").split("\t")
                gene_id, strategy = group.split(KEY_SEPARATOR, 1)
                key = (gene_id, strategy, variant_key)
                if key in counts:
                    raise ValueError(f"Duplicate variant-strategy site: {key}")
                counts[key] = int(depth)
        if len(counts) != site_count:
            raise ValueError(f"Site ortholog-depth row count mismatch: {len(counts)} != {site_count}")
        return counts
