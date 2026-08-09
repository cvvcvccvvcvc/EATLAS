from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from analytics.analyses import external_evidence
from analytics.analyses.external_evidence import _annotate_gnomad, categorize_clinvar_sig
from genomics.gnomad_cache import GnomadRegionCache
from genomics.gnomad_index import GnomadAlleleIndex


def test_clinvar_categories_keep_missing_classification_separate() -> None:
    assert categorize_clinvar_sig("") == ""
    assert categorize_clinvar_sig("Benign") == "B/LB"
    assert categorize_clinvar_sig("Likely_pathogenic") == "P/LP"
    assert categorize_clinvar_sig("Uncertain_significance") == "VUS"
    assert categorize_clinvar_sig("Conflicting_classifications_of_pathogenicity") == "Other"


def test_clinvar_annotation_queries_exact_alleles_with_temporary_bed(
    tmp_path: Path,
    monkeypatch,
) -> None:
    clinvar_vcf = tmp_path / "clinvar.vcf.gz"
    clinvar_vcf.write_bytes(b"vcf")
    Path(f"{clinvar_vcf}.tbi").write_bytes(b"index")
    variants = pd.DataFrame(
        [
            {
                "variant_key": "1:100:A>G",
                "chrom": "1",
                "pos": 100,
                "ref": "A",
                "alt": "G",
            }
        ]
    )

    monkeypatch.setattr(external_evidence.shutil, "which", lambda _name: "/usr/bin/tabix")

    def fake_run(command, **_kwargs):
        regions_path = Path(command[2])
        assert regions_path.read_text() == "1\t99\t100\n"
        return SimpleNamespace(
            returncode=0,
            stdout="1\t100\t.\tA\tG\t.\t.\tCLNSIG=Benign\n",
            stderr="",
        )

    monkeypatch.setattr(external_evidence.subprocess, "run", fake_run)
    result = external_evidence._annotate_clinvar(variants, clinvar_vcf)

    assert result.to_dict(orient="records") == [
        {
            "variant_key": "1:100:A>G",
            "clinvar_found": True,
            "clinvar_classified": True,
            "clinvar_class": "B/LB",
        }
    ]


def test_gnomad_annotation_distinguishes_absence_from_failed_lookup(
    tmp_path: Path,
    monkeypatch,
) -> None:
    variants = pd.DataFrame(
        [
            {"variant_key": "1:100:A>G", "chrom": "1", "pos": 100, "ref": "A", "alt": "G"},
            {"variant_key": "2:200:C>T", "chrom": "2", "pos": 200, "ref": "C", "alt": "T"},
        ]
    )

    def fake_fetch(chrom, _start, _end, **_kwargs):
        if chrom == "2":
            raise RuntimeError("network unavailable")
        return [
            {
                "chrom": "1",
                "pos": 100,
                "ref": "A",
                "alt": "G",
                "joint": {"an": 100, "ac": 2},
            }
        ]

    monkeypatch.setattr(external_evidence, "fetch_region_variants_recursive", fake_fetch)

    evidence, summary = _annotate_gnomad(
        variants,
        gnomad_cache_dir=tmp_path / "gnomad_cache",
    )
    evidence = evidence.set_index("variant_key")

    assert evidence.loc["1:100:A>G", "gnomad_status"] == "ok"
    assert bool(evidence.loc["1:100:A>G", "gnomad_found"])
    assert evidence.loc["1:100:A>G", "gnomad_af"] == 0.02
    assert evidence.loc["2:200:C>T", "gnomad_status"] == "error"
    assert not bool(evidence.loc["2:200:C>T", "gnomad_found"])
    assert summary["failed_region_count"] == 1
    assert summary["shared_cache"]["enabled"] is True
    assert summary["shared_cache"]["tile_write_count"] == 1


def test_gnomad_region_results_are_consumed_with_bounded_futures(
    monkeypatch,
) -> None:
    class TrackedRecords(list):
        alive = 0
        peak = 0

        def __init__(self):
            super().__init__()
            type(self).alive += 1
            type(self).peak = max(type(self).peak, type(self).alive)

        def __del__(self):
            type(self).alive -= 1

    class FakeCache:
        def __init__(self, *_args, **_kwargs):
            pass

        def fetch_region(self, _chrom, _start, _end):
            return TrackedRecords()

        def snapshot(self):
            return {"enabled": True}

    monkeypatch.setattr(external_evidence, "GNOMAD_CLUSTER_GAP_BP", 0)
    monkeypatch.setattr(external_evidence, "GnomadRegionCache", FakeCache)
    variants = pd.DataFrame(
        [
            {
                "variant_key": f"1:{position}:A>G",
                "chrom": "1",
                "pos": position,
                "ref": "A",
                "alt": "G",
            }
            for position in range(100, 140)
        ]
    )

    evidence, summary = _annotate_gnomad(variants)

    assert len(evidence) == len(variants)
    assert summary["region_count"] == len(variants)
    assert TrackedRecords.peak <= external_evidence.GNOMAD_WORKERS * 2


def test_warm_allele_index_is_equivalent_and_avoids_network(
    tmp_path: Path,
    monkeypatch,
) -> None:
    variants = pd.DataFrame(
        [
            {
                "variant_key": "1:100:A>G",
                "chrom": "1",
                "pos": 100,
                "ref": "A",
                "alt": "G",
            },
            {
                "variant_key": "1:101:A>T",
                "chrom": "1",
                "pos": 101,
                "ref": "A",
                "alt": "T",
            },
        ]
    )
    calls = []

    def fake_fetch(chrom, start, end, **_kwargs):
        calls.append((chrom, start, end))
        return [
            {
                "chrom": "1",
                "pos": 100,
                "ref": "A",
                "alt": "G",
                "joint": {"an": 200, "ac": 5},
            }
        ]

    monkeypatch.setattr(
        external_evidence,
        "fetch_region_variants_recursive",
        fake_fetch,
    )
    first, _first_summary = _annotate_gnomad(
        variants,
        gnomad_cache_dir=tmp_path / "gnomad_cache",
    )
    assert calls

    def unexpected_fetch(*_args, **_kwargs):
        raise AssertionError("warm allele index must not use the network")

    monkeypatch.setattr(
        external_evidence,
        "fetch_region_variants_recursive",
        unexpected_fetch,
    )
    second, second_summary = _annotate_gnomad(
        variants,
        gnomad_cache_dir=tmp_path / "gnomad_cache",
    )

    pd.testing.assert_frame_equal(first, second, check_dtype=False)
    assert second_summary["region_count"] == 0
    assert second_summary["allele_index"]["resolved_count"] == len(variants)
    assert second_summary["shared_cache"]["fetch_batch_count"] == 0


def test_indexed_allele_is_not_failed_by_missing_neighbor_tile(
    tmp_path: Path,
    monkeypatch,
) -> None:
    cache_dir = tmp_path / "gnomad_cache"
    variants = pd.DataFrame(
        [
            {
                "variant_key": "1:25000:A>G",
                "chrom": "1",
                "pos": 25_000,
                "ref": "A",
                "alt": "G",
            }
        ]
    )
    record = {
        "chrom": "1",
        "pos": 25_000,
        "ref": "A",
        "alt": "G",
        "joint": {"an": 100, "ac": 2},
    }
    region_cache = GnomadRegionCache(
        cache_dir,
        fetcher=lambda *_args, **_kwargs: [record],
    )
    region_cache.fetch_region("1", 25_000, 25_000)
    allele_index = GnomadAlleleIndex(cache_dir, region_cache=region_cache)
    indexed, unresolved, _summary = allele_index.lookup(variants)
    assert unresolved.empty
    assert indexed["gnomad_found"].tolist() == [True]

    monkeypatch.setattr(
        external_evidence,
        "fetch_region_variants_recursive",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("an indexed allele must not require a neighboring tile")
        ),
    )
    evidence, summary = _annotate_gnomad(variants, gnomad_cache_dir=cache_dir)

    assert evidence["gnomad_status"].tolist() == ["ok"]
    assert evidence["gnomad_found"].tolist() == [True]
    assert evidence["gnomad_af"].tolist() == [0.02]
    assert summary["region_count"] == 0
    assert summary["allele_index"]["resolved_count"] == 1
    assert summary["shared_cache"]["fetch_batch_count"] == 0


def test_regional_failure_does_not_mask_an_already_fetched_allele(
    tmp_path: Path,
    monkeypatch,
) -> None:
    cache_dir = tmp_path / "gnomad_cache"
    variants = pd.DataFrame(
        [
            {
                "variant_key": "1:100:A>G",
                "chrom": "1",
                "pos": 100,
                "ref": "A",
                "alt": "G",
            },
            {
                "variant_key": "1:60000:C>T",
                "chrom": "1",
                "pos": 60_000,
                "ref": "C",
                "alt": "T",
            },
        ]
    )
    GnomadRegionCache(
        cache_dir,
        fetcher=lambda *_args, **_kwargs: [],
    ).fetch_region("1", 30_000, 30_000)

    def partial_fetch(_chrom, start, _end, **_kwargs):
        if start == 1:
            return [
                {
                    "chrom": "1",
                    "pos": 100,
                    "ref": "A",
                    "alt": "G",
                    "joint": {"an": 100, "ac": 2},
                }
            ]
        raise RuntimeError("later tile failed")

    monkeypatch.setattr(
        external_evidence,
        "fetch_region_variants_recursive",
        partial_fetch,
    )
    evidence, summary = _annotate_gnomad(variants, gnomad_cache_dir=cache_dir)
    evidence = evidence.set_index("variant_key")

    assert evidence.loc["1:100:A>G", "gnomad_status"] == "ok"
    assert bool(evidence.loc["1:100:A>G", "gnomad_found"])
    assert evidence.loc["1:60000:C>T", "gnomad_status"] == "error"
    assert not bool(evidence.loc["1:60000:C>T", "gnomad_found"])
    assert summary["failed_region_count"] == 1
    assert summary["allele_index"]["post_fetch_recovered_count"] == 1
    assert summary["allele_index"]["post_fetch_unresolved_count"] == 1


def test_incomplete_external_evidence_retries_only_failed_alleles(
    tmp_path: Path,
    monkeypatch,
) -> None:
    matched = pd.DataFrame(
        [
            {"variant_key": "1:100:A>G", "chrom": "1", "pos": 100, "ref": "A", "alt": "G"},
            {
                "variant_key": "1:500000:C>T",
                "chrom": "1",
                "pos": 500000,
                "ref": "C",
                "alt": "T",
            },
        ]
    )
    matched_path = tmp_path / "matched.tsv.gz"
    matched.to_csv(matched_path, sep="\t", index=False, compression="gzip")
    clinvar_vcf = tmp_path / "clinvar.vcf.gz"
    clinvar_vcf.write_bytes(b"vcf")
    Path(f"{clinvar_vcf}.tbi").write_bytes(b"index")
    output_path = tmp_path / "evidence.tsv.gz"
    manifest_path = tmp_path / "evidence.manifest.json"

    monkeypatch.setattr(
        external_evidence,
        "_annotate_clinvar",
        lambda variants, _path: pd.DataFrame(
            {
                "variant_key": variants["variant_key"],
                "clinvar_found": False,
                "clinvar_classified": False,
                "clinvar_class": "",
            }
        ),
    )
    calls = []

    def fake_gnomad(variants: pd.DataFrame, gnomad_cache_dir=None):
        calls.append(variants["variant_key"].tolist())
        if len(calls) == 1:
            return pd.DataFrame(
                [
                    {
                        "variant_key": "1:100:A>G",
                        "gnomad_status": "ok",
                        "gnomad_found": True,
                        "gnomad_af": 0.1,
                    },
                    {
                        "variant_key": "1:500000:C>T",
                        "gnomad_status": "error",
                        "gnomad_found": False,
                        "gnomad_af": None,
                    },
                ]
            ), {"failed_region_count": 1, "errors": [{"chrom": "1"}]}
        return pd.DataFrame(
            [
                {
                    "variant_key": "1:500000:C>T",
                    "gnomad_status": "ok",
                    "gnomad_found": False,
                    "gnomad_af": None,
                }
            ]
        ), {"failed_region_count": 0, "errors": []}

    monkeypatch.setattr(external_evidence, "_annotate_gnomad", fake_gnomad)

    _first, first_manifest = external_evidence.build_external_evidence(
        matched=matched,
        matched_path=matched_path,
        clinvar_vcf=clinvar_vcf,
        output_path=output_path,
        manifest_path=manifest_path,
    )
    second, second_manifest = external_evidence.build_external_evidence(
        matched=matched,
        matched_path=matched_path,
        clinvar_vcf=clinvar_vcf,
        output_path=output_path,
        manifest_path=manifest_path,
    )

    assert first_manifest["complete"] is False
    assert second_manifest["complete"] is True
    assert calls == [["1:100:A>G", "1:500000:C>T"], ["1:500000:C>T"]]
    assert second_manifest["gnomad"]["cached_ok_allele_count"] == 1
    assert second_manifest["gnomad"]["queried_allele_count"] == 1
    preserved = second.set_index("variant_key").loc["1:100:A>G"]
    assert bool(preserved["gnomad_found"])
    assert preserved["gnomad_af"] == 0.1
