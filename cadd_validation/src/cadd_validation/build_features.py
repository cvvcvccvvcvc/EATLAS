#!/usr/bin/env python3
"""Build variant-level GAPH ortholog features from alignment evidence."""

from __future__ import annotations

import argparse
import bisect
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from .io import iter_tsv, read_tsv, write_tsv


GROUPS = ["all", "primates", "other_mammals", "non_mammal_vertebrates", "other_or_unknown"]
COUNT_CLASSES = ["ref", "alt", "other", "indel"]
FEATURE_TYPE_PRIORITY = {"cds": 0, "utr": 1, "intron": 2, "exon": 3, "gene": 4}
PASSTHROUGH_COLUMNS = [
    "variant_id",
    "label",
    "gene_id",
    "genomic_accession",
    "genomic_start1",
    "ref",
    "alt",
    "target_start0",
    "target_end0",
    "event_type",
    "target_feature_types",
]


@dataclass(frozen=True)
class Variant:
    index: int
    variant_id: str
    label: str
    gene_id: str
    genomic_accession: str
    genomic_start1: str
    ref: str
    alt: str
    target_start0: int
    target_end0: int
    event_type: str
    target_feature_types: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants-tsv", required=True, type=Path)
    parser.add_argument("--segments-tsv", required=True, type=Path)
    parser.add_argument("--events-tsv", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--summaries-tsv", type=Path)
    parser.add_argument("--taxonomy-presets-tsv", type=Path)
    parser.add_argument("--target-features-tsv", type=Path)
    parser.add_argument("--feature-coverage-tsv", type=Path)
    parser.add_argument("--strategies", help="Comma-separated strategy allow-list. Default: all observed strategies.")
    return parser.parse_args()


def clean_allele(value: str) -> str:
    return (value or "").strip().upper()


def to_int(value: str | int | None) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def event_type_for(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) > len(alt):
        return "del"
    if len(alt) > len(ref):
        return "ins"
    return "complex"


def entropy(counts: list[int]) -> float:
    total = sum(counts)
    if total <= 0:
        return 0.0
    value = 0.0
    for count in counts:
        if count:
            p = count / total
            value -= p * math.log2(p)
    return value


def parse_strategy_filter(raw: str | None) -> set[str] | None:
    if not raw or raw.strip().lower() == "all":
        return None
    values = {part.strip() for part in raw.split(",") if part.strip()}
    if not values:
        raise ValueError("--strategies must contain at least one non-empty strategy")
    return values


def load_gene_coordinate_map(target_features_tsv: Path | None) -> dict[str, dict[str, str]]:
    if target_features_tsv is None:
        return {}
    mapping: dict[str, dict[str, str]] = {}
    for row in iter_tsv(target_features_tsv):
        if row.get("feature_type") != "gene":
            continue
        gene_id = row.get("gene_id", "")
        if gene_id:
            mapping[gene_id] = row
    return mapping


def load_target_feature_intervals(target_features_tsv: Path | None) -> dict[str, list[dict[str, str]]]:
    if target_features_tsv is None:
        return {}
    intervals: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in iter_tsv(target_features_tsv):
        gene_id = row.get("gene_id", "")
        start0 = to_int(row.get("target_start0"))
        end0 = to_int(row.get("target_end0"))
        if not gene_id or start0 is None or end0 is None or end0 <= start0:
            continue
        intervals[gene_id].append(row)
    for rows in intervals.values():
        rows.sort(
            key=lambda row: (
                int(row.get("target_start0") or 0),
                int(row.get("target_end0") or 0),
                row.get("feature_type", ""),
            )
        )
    return intervals


def feature_types_for_variant(
    row: dict[str, str],
    feature_intervals: dict[str, list[dict[str, str]]],
    target_start0: int,
    target_end0: int,
) -> str:
    gene_id = row.get("gene_id", "")
    if not gene_id:
        return ""
    types = set()
    for feature in feature_intervals.get(gene_id, []):
        start0 = int(feature.get("target_start0") or 0)
        end0 = int(feature.get("target_end0") or 0)
        if end0 <= target_start0:
            continue
        if start0 >= target_end0:
            break
        if end0 > target_start0 and start0 < target_end0:
            feature_type = feature.get("feature_type", "")
            if feature_type:
                types.add(feature_type)
    return ",".join(sorted(types, key=lambda value: FEATURE_TYPE_PRIORITY.get(value, 99)))


def derive_target_start0(row: dict[str, str], gene_features: dict[str, dict[str, str]]) -> int | None:
    direct = to_int(row.get("target_start0"))
    if direct is not None:
        return direct
    gene_id = row.get("gene_id", "")
    pos1 = to_int(row.get("genomic_start1") or row.get("pos"))
    if not gene_id or pos1 is None:
        return None
    feature = gene_features.get(gene_id)
    if not feature:
        return None
    feature_start1 = to_int(feature.get("genomic_start1"))
    feature_end1 = to_int(feature.get("genomic_end1"))
    if feature_start1 is None or feature_end1 is None:
        return None
    low = min(feature_start1, feature_end1)
    high = max(feature_start1, feature_end1)
    if pos1 < low or pos1 > high:
        return None
    genomic_accession = row.get("genomic_accession", "")
    if genomic_accession and feature.get("genomic_accession") and genomic_accession != feature["genomic_accession"]:
        return None
    return pos1 - low


def make_variant_id(row: dict[str, str], gene_id: str, target_start0: int, ref: str, alt: str) -> str:
    existing = row.get("variant_id") or row.get("id")
    if existing:
        return existing
    accession = row.get("genomic_accession") or row.get("chrom") or gene_id
    pos = row.get("genomic_start1") or row.get("pos") or str(target_start0 + 1)
    return f"{accession}:{pos}:{ref}>{alt}"


def load_variants(path: Path, target_features_tsv: Path | None) -> tuple[list[Variant], dict[str, int]]:
    gene_features = load_gene_coordinate_map(target_features_tsv)
    feature_intervals = load_target_feature_intervals(target_features_tsv)
    variants: list[Variant] = []
    skipped = defaultdict(int)
    for index, row in enumerate(read_tsv(path), start=1):
        ref = clean_allele(row.get("ref", ""))
        alt = clean_allele(row.get("alt", ""))
        if not ref or not alt:
            skipped["missing_ref_alt"] += 1
            continue
        gene_id = str(row.get("gene_id", "")).strip()
        if not gene_id:
            skipped["missing_gene_id"] += 1
            continue
        target_start0 = derive_target_start0(row, gene_features)
        if target_start0 is None:
            skipped["unresolved_target_coordinate"] += 1
            continue
        target_end0 = to_int(row.get("target_end0"))
        if target_end0 is None:
            target_end0 = target_start0 + max(1, len(ref))
        feature_types = feature_types_for_variant(row, feature_intervals, target_start0, target_end0)
        variant_id = make_variant_id(row, gene_id, target_start0, ref, alt)
        variants.append(
            Variant(
                index=len(variants),
                variant_id=variant_id,
                label=row.get("label", ""),
                gene_id=gene_id,
                genomic_accession=row.get("genomic_accession", "") or row.get("chrom", ""),
                genomic_start1=str(row.get("genomic_start1", "") or row.get("pos", "")),
                ref=ref,
                alt=alt,
                target_start0=target_start0,
                target_end0=target_end0,
                event_type=event_type_for(ref, alt),
                target_feature_types=feature_types,
            )
        )
    return variants, dict(skipped)


def load_feature_coverage(path: Path | None) -> dict[tuple[str, str], list[dict[str, str]]]:
    if path is None:
        return {}
    by_gene_strategy: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in iter_tsv(path):
        gene_id = row.get("gene_id", "")
        strategy = row.get("strategy", "")
        start0 = to_int(row.get("target_start0"))
        end0 = to_int(row.get("target_end0"))
        if not gene_id or not strategy or start0 is None or end0 is None or end0 <= start0:
            continue
        by_gene_strategy[(gene_id, strategy)].append(row)
    for rows in by_gene_strategy.values():
        rows.sort(
            key=lambda row: (
                FEATURE_TYPE_PRIORITY.get(row.get("feature_type", ""), 99),
                int(row.get("target_end0") or 0) - int(row.get("target_start0") or 0),
            )
        )
    return by_gene_strategy


def feature_coverage_for_variant(
    feature_coverage: dict[tuple[str, str], list[dict[str, str]]],
    variant: Variant,
    strategy: str,
) -> dict[str, object]:
    for row in feature_coverage.get((variant.gene_id, strategy), []):
        start0 = int(row.get("target_start0") or 0)
        end0 = int(row.get("target_end0") or 0)
        if start0 <= variant.target_start0 and end0 >= variant.target_end0:
            return {
                "gaph_feature_context_type": row.get("feature_type", ""),
                "gaph_feature_context_id": row.get("feature_id", ""),
                "gaph_feature_context_coverage_breadth": row.get("coverage_breadth", ""),
                "gaph_feature_context_mean_depth": row.get("mean_depth", ""),
                "gaph_feature_context_orthologs_covered": row.get("orthologs_covered", ""),
            }
    return {
        "gaph_feature_context_type": "",
        "gaph_feature_context_id": "",
        "gaph_feature_context_coverage_breadth": "",
        "gaph_feature_context_mean_depth": "",
        "gaph_feature_context_orthologs_covered": "",
    }


def load_taxonomy_groups(path: Path | None) -> dict[str, str]:
    groups: dict[str, str] = {}
    if path is None:
        return groups
    for row in iter_tsv(path):
        tax_id = row.get("tax_id", "")
        if not tax_id:
            continue
        if truthy(row.get("is_primate", "")):
            group = "primates"
        elif truthy(row.get("is_mammal", "")):
            group = "other_mammals"
        elif truthy(row.get("is_vertebrate", "")):
            group = "non_mammal_vertebrates"
        else:
            group = "other_or_unknown"
        groups[tax_id] = group
    return groups


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def variant_index_by_gene(variants: list[Variant]) -> dict[str, tuple[list[int], list[Variant]]]:
    grouped: dict[str, list[Variant]] = defaultdict(list)
    for variant in variants:
        grouped[variant.gene_id].append(variant)
    indexed = {}
    for gene_id, rows in grouped.items():
        rows.sort(key=lambda item: (item.target_start0, item.target_end0, item.variant_id))
        starts = [item.target_start0 for item in rows]
        indexed[gene_id] = (starts, rows)
    return indexed


def overlapping_variants(
    gene_index: dict[str, tuple[list[int], list[Variant]]],
    gene_id: str,
    start0: int,
    end0: int,
) -> list[Variant]:
    if gene_id not in gene_index:
        return []
    starts, rows = gene_index[gene_id]
    stop = bisect.bisect_left(starts, end0)
    out = []
    for variant in rows[:stop]:
        if variant.target_end0 > start0:
            out.append(variant)
    return out


def group_for_tax(tax_id: str, taxonomy_groups: dict[str, str]) -> str:
    return taxonomy_groups.get(tax_id, "other_or_unknown")


def add_grouped(
    store: dict[tuple[int, str, str], set[str]],
    variant: Variant,
    strategy: str,
    group: str,
    ortholog_gene_id: str,
) -> None:
    store[(variant.index, strategy, "all")].add(ortholog_gene_id)
    store[(variant.index, strategy, group)].add(ortholog_gene_id)


def load_ortholog_universe(
    summaries_tsv: Path | None,
    taxonomy_groups: dict[str, str],
    strategy_filter: set[str] | None,
) -> dict[tuple[str, str, str], set[str]]:
    universe: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    if summaries_tsv is None:
        return universe
    for row in iter_tsv(summaries_tsv):
        strategy = row.get("strategy", "")
        if strategy_filter is not None and strategy not in strategy_filter:
            continue
        gene_id = row.get("gene_id", "")
        ortholog_gene_id = row.get("ortholog_gene_id", "")
        if not gene_id or not strategy or not ortholog_gene_id:
            continue
        group = group_for_tax(row.get("tax_id", ""), taxonomy_groups)
        universe[(gene_id, strategy, "all")].add(ortholog_gene_id)
        universe[(gene_id, strategy, group)].add(ortholog_gene_id)
    return universe


def collect_features(
    variants: list[Variant],
    segments_tsv: Path,
    events_tsv: Path,
    summaries_tsv: Path | None,
    taxonomy_presets_tsv: Path | None,
    feature_coverage_tsv: Path | None,
    strategy_filter: set[str] | None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    taxonomy_groups = load_taxonomy_groups(taxonomy_presets_tsv)
    feature_coverage = load_feature_coverage(feature_coverage_tsv)
    gene_index = variant_index_by_gene(variants)
    universe = load_ortholog_universe(summaries_tsv, taxonomy_groups, strategy_filter)

    coverage: dict[tuple[int, str, str], set[str]] = defaultdict(set)
    alt_support: dict[tuple[int, str, str], set[str]] = defaultdict(set)
    other_support: dict[tuple[int, str, str], set[str]] = defaultdict(set)
    indel_support: dict[tuple[int, str, str], set[str]] = defaultdict(set)
    observed_strategies: set[str] = set()

    for row in iter_tsv(segments_tsv):
        strategy = row.get("strategy", "")
        if strategy_filter is not None and strategy not in strategy_filter:
            continue
        gene_id = row.get("gene_id", "")
        ortholog_gene_id = row.get("ortholog_gene_id", "")
        if not gene_id or not strategy or not ortholog_gene_id:
            continue
        start0 = to_int(row.get("target_start0"))
        end0 = to_int(row.get("target_end0"))
        if start0 is None or end0 is None or end0 <= start0:
            continue
        group = group_for_tax(row.get("tax_id", ""), taxonomy_groups)
        observed_strategies.add(strategy)
        for variant in overlapping_variants(gene_index, gene_id, start0, end0):
            add_grouped(coverage, variant, strategy, group, ortholog_gene_id)
            universe[(gene_id, strategy, "all")].add(ortholog_gene_id)
            universe[(gene_id, strategy, group)].add(ortholog_gene_id)

    for row in iter_tsv(events_tsv):
        strategy = row.get("strategy", "")
        if strategy_filter is not None and strategy not in strategy_filter:
            continue
        gene_id = row.get("gene_id", "")
        ortholog_gene_id = row.get("ortholog_gene_id", "")
        if not gene_id or not strategy or not ortholog_gene_id:
            continue
        start0 = to_int(row.get("target_start0"))
        end0 = to_int(row.get("target_end0"))
        if start0 is None:
            continue
        if end0 is None:
            end0 = start0 + max(1, len(row.get("ref", "")))
        event_ref = clean_allele(row.get("ref", ""))
        event_alt = clean_allele(row.get("alt", ""))
        event_type = row.get("event_type", "")
        group = group_for_tax(row.get("tax_id", ""), taxonomy_groups)
        observed_strategies.add(strategy)
        for variant in overlapping_variants(gene_index, gene_id, start0, max(end0, start0 + 1)):
            # Event-only strategies still prove that this ortholog was callable
            # at the event locus, even when they cannot distinguish REF support.
            add_grouped(coverage, variant, strategy, group, ortholog_gene_id)
            exact = (
                variant.target_start0 == start0
                and variant.target_end0 == end0
                and variant.ref == event_ref
                and variant.alt == event_alt
            )
            if exact:
                add_grouped(alt_support, variant, strategy, group, ortholog_gene_id)
            elif event_type == "snv" and variant.event_type == "snv" and variant.target_start0 == start0:
                add_grouped(other_support, variant, strategy, group, ortholog_gene_id)
            else:
                add_grouped(indel_support, variant, strategy, group, ortholog_gene_id)
            universe[(gene_id, strategy, "all")].add(ortholog_gene_id)
            universe[(gene_id, strategy, group)].add(ortholog_gene_id)

    strategies = sorted(strategy_filter or observed_strategies)
    rows = []
    for variant in variants:
        for strategy in strategies:
            row: dict[str, object] = {
                "variant_id": variant.variant_id,
                "label": variant.label,
                "gene_id": variant.gene_id,
                "genomic_accession": variant.genomic_accession,
                "genomic_start1": variant.genomic_start1,
                "ref": variant.ref,
                "alt": variant.alt,
                "target_start0": variant.target_start0,
                "target_end0": variant.target_end0,
                "event_type": variant.event_type,
                "target_feature_types": variant.target_feature_types,
                "strategy": strategy,
            }
            row.update(feature_coverage_for_variant(feature_coverage, variant, strategy))
            for group in GROUPS:
                covered = coverage.get((variant.index, strategy, group), set())
                alt = alt_support.get((variant.index, strategy, group), set())
                other = other_support.get((variant.index, strategy, group), set())
                indel = indel_support.get((variant.index, strategy, group), set())
                non_ref = alt | other | indel
                ref = covered - non_ref
                universe_count = len(universe.get((variant.gene_id, strategy, group), set()))
                depth = len(covered)
                alt_count = len(alt)
                other_count = len(other)
                indel_count = len(indel)
                ref_count = len(ref)
                no_call = max(0, universe_count - depth)
                prefix = f"gaph_{group}"
                row[f"{prefix}_ortholog_universe"] = universe_count
                row[f"{prefix}_depth"] = depth
                row[f"{prefix}_ref_count"] = ref_count
                row[f"{prefix}_alt_count"] = alt_count
                row[f"{prefix}_other_count"] = other_count
                row[f"{prefix}_indel_count"] = indel_count
                row[f"{prefix}_no_call_count"] = no_call
                row[f"{prefix}_ref_fraction"] = safe_fraction(ref_count, depth)
                row[f"{prefix}_alt_fraction"] = safe_fraction(alt_count, depth)
                row[f"{prefix}_other_fraction"] = safe_fraction(other_count, depth)
                row[f"{prefix}_indel_fraction"] = safe_fraction(indel_count, depth)
                row[f"{prefix}_callable_fraction"] = safe_fraction(depth, universe_count)
                row[f"{prefix}_entropy"] = entropy([ref_count, alt_count, other_count, indel_count])
            rows.append(row)

    summary = {
        "variant_count": len(variants),
        "strategy_count": len(strategies),
        "strategies": strategies,
        "feature_row_count": len(rows),
    }
    return rows, summary


def safe_fraction(numerator: int, denominator: int) -> float:
    if denominator <= 0:
        return 0.0
    return round(numerator / denominator, 6)


def output_fields() -> list[str]:
    fields = PASSTHROUGH_COLUMNS + ["strategy"]
    fields.extend(
        [
            "gaph_feature_context_type",
            "gaph_feature_context_id",
            "gaph_feature_context_coverage_breadth",
            "gaph_feature_context_mean_depth",
            "gaph_feature_context_orthologs_covered",
        ]
    )
    for group in GROUPS:
        prefix = f"gaph_{group}"
        fields.extend(
            [
                f"{prefix}_ortholog_universe",
                f"{prefix}_depth",
                f"{prefix}_ref_count",
                f"{prefix}_alt_count",
                f"{prefix}_other_count",
                f"{prefix}_indel_count",
                f"{prefix}_no_call_count",
                f"{prefix}_ref_fraction",
                f"{prefix}_alt_fraction",
                f"{prefix}_other_fraction",
                f"{prefix}_indel_fraction",
                f"{prefix}_callable_fraction",
                f"{prefix}_entropy",
            ]
        )
    return fields


def main() -> None:
    args = parse_args()
    strategy_filter = parse_strategy_filter(args.strategies)
    variants, skipped = load_variants(args.variants_tsv, args.target_features_tsv)
    rows, summary = collect_features(
        variants=variants,
        segments_tsv=args.segments_tsv,
        events_tsv=args.events_tsv,
        summaries_tsv=args.summaries_tsv,
        taxonomy_presets_tsv=args.taxonomy_presets_tsv,
        feature_coverage_tsv=args.feature_coverage_tsv,
        strategy_filter=strategy_filter,
    )
    summary["skipped_variants"] = skipped
    write_tsv(args.out_tsv, rows, output_fields())
    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
