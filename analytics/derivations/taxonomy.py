"""Taxonomic scopes and grouping semantics for ortholog evidence."""

from __future__ import annotations

import csv
import gzip
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from genomics.taxonomy import TAXONOMY_FIELDS


SCOPE_ANCESTORS = {
    "all": "",
    "eukaryota": "2759",
    "metazoa": "33208",
    "vertebrata": "7742",
    "tetrapoda": "32523",
    "amniota": "32524",
    "mammalia": "40674",
    "primates": "9443",
}
SCOPE_ORDER = tuple(SCOPE_ANCESTORS)
UNIT_ORDER = ("ortholog", "species", "genus", "family", "order")
COUNT_KEYS = tuple(
    f"{scope}__{unit}"
    for scope in SCOPE_ORDER
    for unit in UNIT_ORDER
)
TAXONOMY_SUMMARY_FIELDS = [
    "taxonomic_scope",
    "evidence_unit",
    "gene_count",
    "ortholog_count",
    "taxon_count",
    "unit_count",
    "orthologs_per_gene_min",
    "orthologs_per_gene_median",
    "orthologs_per_gene_mean",
    "orthologs_per_gene_max",
    "units_per_gene_min",
    "units_per_gene_median",
    "units_per_gene_mean",
    "units_per_gene_max",
]


@dataclass(frozen=True)
class TaxonomyProfile:
    tax_id: str
    ancestor_ids: frozenset[str]
    species_id: str
    genus_id: str
    family_id: str
    order_id: str

    def scopes(self) -> tuple[str, ...]:
        return tuple(
            scope
            for scope, ancestor in SCOPE_ANCESTORS.items()
            if not ancestor or ancestor in self.ancestor_ids
        )

    def unit_id(self, unit: str, ortholog_gene_id: str) -> str:
        if unit == "ortholog":
            return f"ortholog:{ortholog_gene_id}"
        value = getattr(self, f"{unit}_id")
        return f"{unit}:{value}" if value else f"{unit}:taxon:{self.tax_id}"


def _split_ids(value: str) -> frozenset[str]:
    return frozenset(
        item.strip()
        for item in value.replace(";", ",").split(",")
        if item.strip()
    )


def load_taxonomy_profiles(path: Path) -> dict[str, TaxonomyProfile]:
    profiles: dict[str, TaxonomyProfile] = {}
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if fieldnames != TAXONOMY_FIELDS:
            raise ValueError(
                f"Taxonomy table {path} must use exact canonical fields: "
                + ", ".join(TAXONOMY_FIELDS)
            )
        for row_number, row in enumerate(reader, start=2):
            tax_id = str(row.get("tax_id") or "").strip()
            if not tax_id:
                raise ValueError(f"Taxonomy table {path} has empty tax_id at row {row_number}")
            if tax_id in profiles:
                raise ValueError(f"Duplicate taxonomy tax_id in {path}: {tax_id}")
            status = str(row.get("taxonomy_status") or "").strip()
            if status not in {"resolved", "not_returned"}:
                raise ValueError(
                    f"Taxonomy table {path} has invalid taxonomy_status for tax_id {tax_id}: "
                    f"{status!r}"
                )
            lineage_value = str(row.get("lineage_tax_ids") or "")
            profiles[tax_id] = TaxonomyProfile(
                tax_id=tax_id,
                ancestor_ids=frozenset(
                    {*_split_ids(lineage_value), tax_id}
                ),
                species_id=str(row.get("species_id") or ""),
                genus_id=str(row.get("genus_id") or ""),
                family_id=str(row.get("family_id") or ""),
                order_id=str(row.get("order_id") or ""),
            )
    return profiles


def member_group_keys(
    ortholog_gene_id: str,
    tax_id: str,
    profiles: dict[str, TaxonomyProfile],
) -> tuple[str, ...]:
    if not ortholog_gene_id:
        return ()
    profile = profiles.get(tax_id)
    if profile is None:
        raise ValueError(
            "Ortholog evidence references tax_id absent from canonical taxonomy: "
            f"ortholog_gene_id={ortholog_gene_id!r}, tax_id={tax_id!r}"
        )
    return tuple(
        f"{scope}__{unit}={profile.unit_id(unit, ortholog_gene_id)}"
        for scope in profile.scopes()
        for unit in UNIT_ORDER
    )


def count_member_groups(
    members: Iterable[tuple[str, str]],
    profiles: dict[str, TaxonomyProfile],
) -> dict[str, int]:
    groups = {key: set() for key in COUNT_KEYS}
    for ortholog_gene_id, tax_id in members:
        for item in member_group_keys(ortholog_gene_id, tax_id, profiles):
            key, group_id = item.split("=", 1)
            groups[key].add(group_id)
    return {key: len(values) for key, values in groups.items()}


def _format_stat(value: float) -> str:
    return f"{value:.1f}"


def build_taxonomy_summary_rows(
    ortholog_rows: Iterable[dict[str, str]],
    profiles: dict[str, TaxonomyProfile],
) -> list[dict[str, object]]:
    members_by_gene: dict[str, list[tuple[str, str]]] = {}
    selected_tax_ids: set[str] = set()
    for row_number, row in enumerate(ortholog_rows, start=2):
        gene_id = str(row.get("query_gene_id") or "")
        ortholog_gene_id = str(row.get("ortholog_gene_id") or "")
        tax_id = str(row.get("tax_id") or "").strip()
        if not tax_id:
            raise ValueError(
                "Selected ortholog table has empty tax_id at row "
                f"{row_number}: ortholog_gene_id={ortholog_gene_id!r}"
            )
        selected_tax_ids.add(tax_id)
        if gene_id and ortholog_gene_id:
            members_by_gene.setdefault(gene_id, []).append((ortholog_gene_id, tax_id))

    taxonomy_tax_ids = set(profiles)
    missing = selected_tax_ids - taxonomy_tax_ids
    unexpected = taxonomy_tax_ids - selected_tax_ids
    if missing or unexpected:
        details = []
        if missing:
            details.append("missing taxonomy tax_id(s): " + ", ".join(sorted(missing)))
        if unexpected:
            details.append("unexpected taxonomy tax_id(s): " + ", ".join(sorted(unexpected)))
        raise ValueError("Canonical taxonomy coverage mismatch; " + "; ".join(details))

    rows: list[dict[str, object]] = []
    for scope in SCOPE_ORDER:
        for unit in UNIT_ORDER:
            ortholog_counts: list[int] = []
            unit_counts: list[int] = []
            all_taxa: set[str] = set()
            all_units: set[str] = set()
            for members in members_by_gene.values():
                selected_orthologs: set[str] = set()
                selected_units: set[str] = set()
                for ortholog_gene_id, tax_id in members:
                    profile = profiles[tax_id]
                    if scope not in profile.scopes():
                        continue
                    selected_orthologs.add(ortholog_gene_id)
                    if tax_id:
                        all_taxa.add(tax_id)
                    group_id = profile.unit_id(unit, ortholog_gene_id)
                    selected_units.add(group_id)
                    all_units.add(group_id)
                ortholog_counts.append(len(selected_orthologs))
                unit_counts.append(len(selected_units))
            if not ortholog_counts:
                continue
            rows.append(
                {
                    "taxonomic_scope": scope,
                    "evidence_unit": unit,
                    "gene_count": len(ortholog_counts),
                    "ortholog_count": sum(ortholog_counts),
                    "taxon_count": len(all_taxa),
                    "unit_count": len(all_units),
                    "orthologs_per_gene_min": min(ortholog_counts),
                    "orthologs_per_gene_median": _format_stat(statistics.median(ortholog_counts)),
                    "orthologs_per_gene_mean": _format_stat(statistics.mean(ortholog_counts)),
                    "orthologs_per_gene_max": max(ortholog_counts),
                    "units_per_gene_min": min(unit_counts),
                    "units_per_gene_median": _format_stat(statistics.median(unit_counts)),
                    "units_per_gene_mean": _format_stat(statistics.mean(unit_counts)),
                    "units_per_gene_max": max(unit_counts),
                }
            )
    return rows


def write_taxonomy_summary(path: Path, rows: Iterable[dict[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TAXONOMY_SUMMARY_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count
