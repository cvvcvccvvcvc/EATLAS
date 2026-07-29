from __future__ import annotations

from pathlib import Path

import pandas as pd

from analytics.core import clinvar_validation


def test_clinvar_vep_artifact_aggregates_genes_and_resumes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    universe_path = tmp_path / "clinvar_universe.snv_indel.tsv.gz"
    universe = pd.DataFrame(
        [
            {
                "variant_key": "1:10:A>G",
                "gene_ids": "1|2",
                "clinvar_mc_terms": "synonymous_variant",
            },
            {
                "variant_key": "1:20:C>T",
                "gene_ids": "1",
                "clinvar_mc_terms": "intron_variant",
            },
        ]
    )
    universe.to_csv(universe_path, sep="\t", index=False, compression="gzip")

    def fake_annotate(requests, _cache_path, **_kwargs):
        assert len(requests) == 3
        return (
            pd.DataFrame(
                [
                    {
                        "variant_key": "1:10:A>G",
                        "gene_id": "1",
                        "status": "ok",
                        "primary_consequence": "missense_variant",
                        "consequence_terms": "missense_variant",
                    },
                    {
                        "variant_key": "1:10:A>G",
                        "gene_id": "2",
                        "status": "ok",
                        "primary_consequence": "splice_region_variant",
                        "consequence_terms": "splice_region_variant",
                    },
                    {
                        "variant_key": "1:20:C>T",
                        "gene_id": "1",
                        "status": "no_target_gene",
                        "primary_consequence": "",
                        "consequence_terms": "",
                    },
                ]
            ),
            {"release": "116", "requested": 3},
        )

    monkeypatch.setattr(clinvar_validation, "annotate_vep_consequences", fake_annotate)
    output_path, manifest_path, enriched, manifest = clinvar_validation.build_or_load_vep_universe(
        universe=universe,
        universe_path=universe_path,
        analytics_dir=tmp_path,
        backend="local",
        release="116",
        vep_executable="vep",
        vep_cache_dir=tmp_path,
        vep_forks=2,
    )

    assert output_path.exists()
    assert manifest_path.exists()
    assert not manifest["cache_hit"]
    first = enriched.set_index("variant_key").loc["1:10:A>G"]
    assert first["vep_primary_consequence"] == "missense_variant"
    assert first["vep_consequence_terms"] == "missense_variant|splice_region_variant"
    assert enriched.set_index("variant_key").loc["1:20:C>T", "vep_status"] == "no_target_gene"

    monkeypatch.setattr(
        clinvar_validation,
        "annotate_vep_consequences",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("cache was not reused")),
    )
    _path, _manifest_path, cached, cached_manifest = clinvar_validation.build_or_load_vep_universe(
        universe=universe,
        universe_path=universe_path,
        analytics_dir=tmp_path,
        backend="local",
        release="116",
        vep_executable="vep",
        vep_cache_dir=tmp_path,
        vep_forks=4,
    )

    assert cached_manifest["cache_hit"]
    assert cached.equals(enriched)
