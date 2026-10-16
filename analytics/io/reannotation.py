"""Read-only application of explicitly selected gnomAD supplement snapshots."""

from __future__ import annotations

import csv
import json
from dataclasses import replace
from pathlib import Path

import duckdb

from genomics.gnomad_cache import tiles_for_region
from genomics.variants import (build_context_index, contexts_for_variant, load_target_contexts,
                               normalize_vcf_key_for_context)
from reannotation.store import atomic_json, digest, load_snapshot, tile_key
from analytics.io.duckdb import configure_duckdb_memory

VIEW_SCHEMA = "gaph_gnomad_repaired_view_v1"


def prepare_repaired_source(source, snapshot: Path, cache_root: Path):
    value = load_snapshot(snapshot)
    bindings = [v for v in value["sources"].values()
                if v["root_manifest_sha256"] == digest(source.root_manifest_json)
                and v["evidence_tree_sha256"] == source.evidence_inventory["tree_sha256"]]
    if len(bindings) != 1:
        raise ValueError(f"Supplement does not uniquely bind source {source.run_dir.name}")
    binding = bindings[0]
    output = cache_root / snapshot.stem / source.source_id
    manifest = output / "manifest.json"
    expected = {"schema": VIEW_SCHEMA, "snapshot": snapshot.stem, "source_id": source.source_id,
                "base": str(source.variant_annotations_source.resolve())}
    if manifest.exists():
        old = json.loads(manifest.read_text())
        if all(old.get(k) == v for k, v in expected.items()) and all(
            digest(output / old[name]["path"]) == old[name]["sha256"] for name in ["lookup", "coverage"]
        ):
            return replace(source, variant_annotations_source=manifest)
    output.mkdir(parents=True, exist_ok=True)
    covered = set(value["tiles"]) & set(binding["tiles"])
    objects = sorted({value["tiles"][k]["object"] for k in covered})
    contexts = load_target_contexts(source.genes_tsv, source.target_sequences_dir)
    context_index = build_context_index(contexts)
    temporary = output / "aliases.tsv"
    with temporary.open("w") as handle, duckdb.connect() as connection:
        configure_duckdb_memory(connection, 1)
        connection.execute("SET threads=1")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["variant_key", "gnomad_af", "gnomad_af_source", "gnomad_csq"])
        for sha in objects:
            path = snapshot.parent.parent / value["objects"][sha]["path"]
            active = sorted(k for k in covered if value["tiles"][k]["object"] == sha)
            cursor = connection.execute("SELECT chrom,pos,ref,alt,gnomad_af,gnomad_af_source,gnomad_csq "
                                        "FROM read_parquet(?) WHERE tile IN (SELECT unnest(?))", [str(path), active])
            while rows := cursor.fetchmany(10000):
                for chrom, pos, ref, alt, af, af_source, csq in rows:
                    keys = {(str(chrom), int(pos), ref, alt)}
                    for context in contexts_for_variant(context_index, str(chrom), int(pos)):
                        normalized, status = normalize_vcf_key_for_context(context, str(chrom), int(pos), ref, alt)
                        if status == "ok":
                            keys.add(normalized)
                    for c, p, r, a in keys:
                        writer.writerow([f"{c}:{p}:{r}>{a}", "" if af is None else str(af), af_source or "", csq or ""])
    with duckdb.connect() as connection:
        configure_duckdb_memory(connection, 1)
        connection.execute("SET temp_directory=?", [str(output / "spill")])
        connection.read_csv(str(temporary), delimiter="\t", header=True, all_varchar=True).create_view("aliases")
        conflict = connection.execute("SELECT variant_key FROM aliases GROUP BY variant_key "
                                      "HAVING count(DISTINCT gnomad_af)>1 OR count(DISTINCT gnomad_af_source)>1 LIMIT 1").fetchone()
        if conflict:
            raise ValueError(f"Supplement observations disagree for {conflict[0]}")
        connection.sql("SELECT variant_key, max(coalesce(gnomad_af,'')) AS gnomad_af, "
                       "max(coalesce(gnomad_af_source,'')) AS gnomad_af_source, max(coalesce(gnomad_csq,'')) AS gnomad_csq "
                       "FROM aliases GROUP BY variant_key ORDER BY variant_key").write_parquet(str(output / "lookup.parquet"), compression="zstd")
        connection.execute("CREATE TABLE coverage(chrom VARCHAR,start1 BIGINT,end1 BIGINT,complete BOOLEAN)")
        rows = [(chrom, start, end, all(tile_key(t) in covered for t in tiles_for_region(chrom, max(1,start-100),end+100)))
                for chrom, start, end in binding["regions"]]
        if rows:
            connection.executemany("INSERT INTO coverage VALUES (?,?,?,?)", rows)
        connection.table("coverage").write_parquet(str(output / "coverage.parquet"), compression="zstd")
    temporary.unlink()
    expected.update({name: {"path": f"{name}.parquet", "sha256": digest(output / f"{name}.parquet")}
                     for name in ["lookup", "coverage"]})
    atomic_json(manifest, expected)
    return replace(source, variant_annotations_source=manifest)


def repaired_sql(base_sql: str, repair: dict, columns: tuple[str, ...]) -> str:
    from analytics.io.variant_source import sql_string

    lookup = sql_string(repair["lookup"])
    coverage = sql_string(repair["coverage"])
    scope = (f"EXISTS (SELECT 1 FROM read_parquet({coverage}) c WHERE c.chrom=split_part(b.variant_key,':',1) "
             "AND try_cast(split_part(b.variant_key,':',2) AS BIGINT) BETWEEN c.start1 AND c.end1")
    eligible = f"b.lookup_status='ok' AND nullif(b.gnomad_af,'') IS NULL AND {scope})"
    complete = scope + " AND c.complete)"
    replacements = []
    for name in ["gnomad_af", "gnomad_af_source", "gnomad_csq"]:
        if name in columns:
            replacements.append(f"CASE WHEN ({eligible}) AND r.variant_key IS NOT NULL THEN coalesce(r.{name},'') ELSE b.{name} END AS {name}")
    status = (f"CASE WHEN ({eligible}) AND r.variant_key IS NOT NULL THEN 'found' "
              f"WHEN ({eligible}) AND {complete} THEN 'not_found' ELSE '' END AS gnomad_lookup_status")
    return (f"(SELECT b.* REPLACE ({','.join(replacements)}), {status} FROM {base_sql} b "
            f"LEFT JOIN read_parquet({lookup}) r USING(variant_key))")
