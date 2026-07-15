"""Variant key normalization against GAPH target loci."""

from __future__ import annotations

import bisect
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


DNA_BASES = set("ACGT")


def open_text(path: Path):
    return gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else path.open(newline="")


def read_fasta_sequence(path: Path) -> str:
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            chunks.append(line.strip())
    return "".join(chunks).upper()


def normalize_chrom(value: str | None) -> str | None:
    if not value:
        return None
    chrom = str(value).strip()
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        return "MT"
    return chrom


def refseq_accession_to_chrom(value: str | None) -> str | None:
    chrom = normalize_chrom(value)
    if not chrom:
        return None
    if chrom in {"X", "Y", "MT"}:
        return chrom
    if chrom.isdigit():
        num = int(chrom)
        if num == 23:
            return "X"
        if num == 24:
            return "Y"
        return str(num)

    match = re.search(r"NC_0+(\d+)\.", chrom)
    if not match:
        return None
    num = int(match.group(1))
    if num == 23:
        return "X"
    if num == 24:
        return "Y"
    if num in {12920, 1807}:
        return "MT"
    return str(num)


def variant_key_text(key: tuple[str, int, str, str] | None) -> str:
    if not key:
        return ""
    chrom, pos, ref, alt = key
    return f"{chrom}:{pos}:{ref}>{alt}"


def variant_type(ref: str, alt: str) -> str:
    ref = str(ref or "").upper()
    alt = str(alt or "").upper()
    if not ref or not alt or not set(ref) <= DNA_BASES or not set(alt) <= DNA_BASES:
        return "unsupported"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) != len(alt):
        return "indel"
    return "complex"


def load_target_contexts(genes_tsv: Path, target_sequences_dir: Path) -> dict[str, dict]:
    if not genes_tsv.exists():
        raise FileNotFoundError(f"Target genes table not found: {genes_tsv}")
    if not target_sequences_dir.exists():
        raise FileNotFoundError(f"Target sequences directory not found: {target_sequences_dir}")

    contexts = {}
    with open_text(genes_tsv) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", "genomic_accession", "chromosome", "begin", "end"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Target genes table missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            gene_id = str(row["gene_id"])
            fasta_path = target_sequences_dir / f"{gene_id}.fa.gz"
            if not fasta_path.exists():
                raise FileNotFoundError(f"Target FASTA not found for gene {gene_id}: {fasta_path}")
            contexts[gene_id] = {
                "gene_id": gene_id,
                "accession": row["genomic_accession"],
                "chrom": normalize_chrom(row["chromosome"]) or refseq_accession_to_chrom(row["genomic_accession"]),
                "begin": int(row["begin"]),
                "end": int(row["end"]),
                "fasta_path": fasta_path,
            }
    return contexts


def context_sequence(context: dict) -> str:
    seq = context.get("seq")
    if seq is None:
        seq = read_fasta_sequence(context["fasta_path"])
        context["seq"] = seq
    return seq


def build_context_index(contexts: dict[str, dict]) -> dict[str, tuple[list[dict], list[int]]]:
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for context in contexts.values():
        chrom = context.get("chrom")
        if chrom:
            by_chrom[str(chrom)].append(context)

    index = {}
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda row: (int(row["begin"]), int(row["end"]), str(row["gene_id"])))
        index[chrom] = (rows, [int(row["begin"]) for row in rows])
    return index


def contexts_for_variant(
    context_index: dict[str, tuple[list[dict], list[int]]],
    chrom: str,
    pos: int,
) -> list[dict]:
    rows, starts = context_index.get(chrom, ([], []))
    limit = bisect.bisect_right(starts, pos)
    return [context for context in rows[:limit] if int(context["end"]) >= pos]


def normalize_vcf_key_for_context(
    context: dict,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> tuple[tuple[str, int, str, str] | None, str]:
    ref = str(ref or "").upper()
    alt = str(alt or "").upper()
    if variant_type(ref, alt) == "unsupported":
        return None, "unsupported_allele"

    seq = context_sequence(context)
    pos0 = pos - int(context["begin"])
    if pos0 < 0 or pos0 + len(ref) > len(seq):
        return None, "out_of_target"
    if seq[pos0 : pos0 + len(ref)] != ref:
        return None, "ref_mismatch"

    if len(ref) != len(alt):
        while pos0 > 0 and ref[-1] == alt[-1]:
            prev = seq[pos0 - 1]
            ref = prev + ref[:-1]
            alt = prev + alt[:-1]
            pos0 -= 1

    return (chrom, int(context["begin"]) + pos0, ref, alt), "ok"
