"""Comparable allele-level condition distributions from one ClinVar release."""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter
from itertools import zip_longest
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

from analytics.io.artifacts import path_metadata, write_json_atomic
from genomics.variants import normalize_chrom, variant_type
from .clinvar_validation import clinvar_label, info_value

CONDITION_COLUMNS = [
    "condition_key",
    "condition",
    "variant_count",
    "total_variant_count",
    "named_variant_count",
]


def parse_conditions(names_text: str, ids_text: str) -> list[tuple[str, str, str]]:
    """Use supplied disease identifiers; do not infer disease categories or synonyms."""
    conditions = {}
    for names, identifiers in zip_longest(names_text.split(";"), ids_text.split(";"), fillvalue=""):
        for name, ids in zip_longest(names.split("|"), identifiers.split("|"), fillvalue=""):
            name = unquote(name).strip().replace("_", " ")
            if name.casefold() in {"", ".", "not provided", "not specified", "not applicable"}:
                continue
            ids = "" if ids == "." else unquote(ids)
            values = sorted({item.strip() for item in ids.split(",") if item.strip()})
            key = next(
                (
                    item
                    for prefix in ("MedGen:", "MONDO:", "OMIM:")
                    for item in values
                    if item.startswith(prefix)
                ),
                None,
            )
            key = key or ("ids:" + ",".join(values) if values else "name:" + name.casefold())
            conditions.setdefault(key, (key, name, ",".join(values)))
    return list(conditions.values())


def condition_rows(source: pd.DataFrame) -> pd.DataFrame:
    rows = [
        (str(row.variant_key), key, name, ids)
        for row in source.itertuples(index=False)
        for key, name, ids in parse_conditions(
            str(row.clinvar_disease_names), str(row.clinvar_disease_ids)
        )
    ]
    return pd.DataFrame(
        rows, columns=["variant_key", "condition_key", "condition", "condition_ids"]
    ).drop_duplicates(["variant_key", "condition_key"])


def condition_distribution(keys: set[str], memberships: pd.DataFrame) -> pd.DataFrame:
    selected = memberships[memberships["variant_key"].isin(keys)]
    named = selected["variant_key"].nunique()
    counts = (
        selected.groupby("condition_key")
        .agg(condition=("condition", "min"), variant_count=("variant_key", "nunique"))
        .reset_index()
    )
    if counts.empty:
        counts = pd.DataFrame([{"condition_key": "", "condition": "", "variant_count": 0}])
    return counts.assign(total_variant_count=len(keys), named_variant_count=int(named))[
        CONDITION_COLUMNS
    ]


def global_condition_distribution(clinvar_vcf: Path, cache_dir: Path) -> pd.DataFrame:
    """Stream the release once; keep duplicate-allele state for one chromosome only."""
    inputs = {"schema_version": 1, "vcf": path_metadata(clinvar_vcf)}
    fingerprint = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()[:24]
    path = cache_dir / f"{fingerprint}.json"
    if path.exists():
        cached = json.loads(path.read_text())
        if cached.get("inputs") == inputs and cached.get("status") == "complete":
            return pd.DataFrame(cached["rows"], columns=["variant_type", *CONDITION_COLUMNS])

    totals = Counter()
    named = Counter()
    counts = Counter()
    names = {}
    records = {}

    def finish_chromosome() -> None:
        for (position, ref, alt), (labels, conditions) in records.items():
            if labels != {"pathogenic"}:
                continue
            kind = variant_type(ref, alt)
            totals[kind] += 1
            named[kind] += bool(conditions)
            for key, name in conditions.items():
                counts[(kind, key)] += 1
                names[key] = min(names.get(key, name), name)
        records.clear()

    chromosome = None
    completed = set()
    header_seen = False
    record_count = 0
    opener = gzip.open if clinvar_vcf.suffix == ".gz" else open
    with opener(clinvar_vcf, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                header_seen |= line.startswith("#CHROM\t")
                continue
            fields = line.rstrip("\n").split("\t")
            if not header_seen or len(fields) < 8:
                raise ValueError(f"Invalid ClinVar VCF record at {clinvar_vcf}:{line_number}")
            chrom, pos, _id, ref, alts, _qual, _filter, info = fields[:8]
            chrom = normalize_chrom(chrom)
            if chrom is None:
                continue
            record_count += 1
            if chrom != chromosome:
                finish_chromosome()
                if chrom in completed:
                    raise ValueError(
                        "ClinVar condition background requires chromosome-grouped VCF records"
                    )
                completed.add(chrom)
                chromosome = chrom
            label = clinvar_label(info_value(info, "CLNSIG"))
            if label not in {"benign", "pathogenic"}:
                continue
            conditions = parse_conditions(info_value(info, "CLNDN"), info_value(info, "CLNDISDB"))
            for alt in alts.split(","):
                if variant_type(ref, alt) not in {"snv", "indel"}:
                    continue
                labels, diseases = records.setdefault(
                    (int(pos), ref.upper(), alt.upper()), (set(), {})
                )
                labels.add(label)
                for key, name, _ids in conditions:
                    diseases[key] = min(diseases.get(key, name), name)
    finish_chromosome()
    if not header_seen or record_count == 0:
        raise ValueError(f"ClinVar condition background has no genomic records: {clinvar_vcf}")
    rows = []
    for kind in ("snv", "indel"):
        keys = sorted(key for count_kind, key in counts if count_kind == kind) or [""]
        rows.extend(
            {
                "variant_type": kind,
                "condition_key": key,
                "condition": names.get(key, ""),
                "variant_count": counts[(kind, key)],
                "total_variant_count": totals[kind],
                "named_variant_count": named[kind],
            }
            for key in keys
        )
    cache_dir.mkdir(parents=True, exist_ok=True)
    write_json_atomic(path, {"status": "complete", "inputs": inputs, "rows": rows})
    return pd.DataFrame(rows, columns=["variant_type", *CONDITION_COLUMNS])


def compare_condition_distributions(
    variants: pd.DataFrame,
    universe: pd.DataFrame,
    eligible_gene_ids_by_strategy: dict[str, set[str]],
    global_counts: pd.DataFrame,
) -> pd.DataFrame:
    memberships = condition_rows(universe)
    universe = universe[universe["label_class"].eq("pathogenic")]
    candidate = variants.copy()
    candidate["variant_type"] = candidate["event_type"].map(
        lambda value: "snv" if value == "snv" else "indel"
    )
    candidate["strategy"] = candidate["strategies"].str.split(",")
    candidate = candidate.explode("strategy")
    candidate["strategy"] = candidate["strategy"].str.strip()
    global_all = global_counts.groupby(["condition_key", "condition"], as_index=False)[
        "variant_count"
    ].sum()
    denominators = (
        global_counts.groupby("variant_type")[["total_variant_count", "named_variant_count"]]
        .first()
        .sum()
    )
    global_all = global_all.assign(variant_type="all", **denominators.to_dict())
    global_counts = pd.concat([global_counts, global_all], ignore_index=True)
    parts = [global_counts.assign(strategy="", cohort="global")]
    for strategy, genes in eligible_gene_ids_by_strategy.items():
        targets = universe[
            universe["gene_ids"].map(lambda value: bool(genes.intersection(str(value).split("|"))))
        ]
        for kind in ("all", "snv", "indel"):
            observed = set(
                candidate.loc[
                    candidate["strategy"].eq(strategy)
                    & (candidate["variant_type"].eq(kind) | (kind == "all")),
                    "variant_key",
                ]
            )
            background = set(
                targets.loc[(targets["variant_type"].eq(kind) | (kind == "all")), "variant_key"]
            )
            if not observed <= background:
                raise ValueError(
                    f"P/LP candidates absent from the matching ClinVar background for {strategy}/{kind}"
                )
            for cohort, keys in (("gaph", observed), ("target", background)):
                parts.append(
                    condition_distribution(keys, memberships).assign(
                        strategy=strategy, variant_type=kind, cohort=cohort
                    )
                )
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
