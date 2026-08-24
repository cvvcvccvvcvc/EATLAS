from __future__ import annotations

import stat
import sys
from contextlib import contextmanager
from pathlib import Path

import pandas as pd
import pytest

from analytics.analyses import clinvar_validation
from analytics.analyses.observed_variant_store import (
    build_or_load_observed_variant_store,
)


def test_clinvar_universe_writer_preserves_schema_and_shared_permissions(
    tmp_path: Path,
) -> None:
    universe_path = tmp_path / "clinvar_universe.tsv.gz"
    row = dict.fromkeys(clinvar_validation.UNIVERSE_FIELDS, "")
    row.update({"variant_key": "1:10:A>G", "variant_type": "snv", "pos": 10})

    clinvar_validation.write_universe(universe_path, [row])

    observed = pd.read_csv(
        universe_path,
        sep="\t",
        compression="gzip",
        keep_default_na=False,
    )
    assert observed.columns.tolist() == clinvar_validation.UNIVERSE_FIELDS
    assert observed.columns.tolist()[:-1] == [
        "variant_key",
        "variant_type",
        "chrom",
        "pos",
        "ref",
        "alt",
        "label_class",
        "clinvar_ids",
        "clinvar_sigs",
        "clinvar_mc_so_ids",
        "clinvar_mc_terms",
        "gene_ids",
    ]
    assert observed.columns[-1] == "clinvar_disease_ids"
    assert observed.loc[0, "variant_key"] == "1:10:A>G"
    assert stat.S_IMODE(universe_path.stat().st_mode) == 0o644


def test_clinvar_universe_preserves_clndisdb_source_structure(
    tmp_path: Path,
    monkeypatch,
) -> None:
    lines = [
        "1\t10\tVCV2\tA\tG\t.\t.\t"
        "CLNSIG=Pathogenic;MC=SO:0001583|missense_variant;"
        "CLNDISDB=MedGen:C1|MONDO:MONDO:1,OMIM:1\n",
        "1\t10\tVCV1\tA\tG\t.\t.\t"
        "CLNSIG=Pathogenic;MC=SO:0001583|missense_variant;"
        "CLNDISDB=.|MedGen:C2\n",
        "1\t10\tVCV3\tA\tG\t.\t.\t"
        "CLNSIG=Pathogenic;MC=SO:0001583|missense_variant;CLNDISDB=.\n",
    ]

    @contextmanager
    def fake_tabix_output_lines(*_args, **_kwargs):
        yield iter(lines)

    monkeypatch.setattr(clinvar_validation.shutil, "which", lambda _name: "tabix")
    monkeypatch.setattr(
        clinvar_validation,
        "tabix_output_lines",
        fake_tabix_output_lines,
    )
    monkeypatch.setattr(
        clinvar_validation,
        "normalize_clinvar_allele_for_targets",
        lambda *_args: [(("1", 10, "A", "G"), "snv", {"gene_id": "1"})],
    )

    rows, _counts = clinvar_validation.query_clinvar_variant_universe(
        tmp_path / "clinvar.vcf.gz",
        tmp_path / "regions.bed",
        {},
    )

    assert len(rows) == 1
    assert rows[0]["clinvar_disease_ids"] == (
        ".|MedGen:C2;MedGen:C1|MONDO:MONDO:1,OMIM:1"
    )


def test_tabix_output_is_streamed_and_errors_are_reported(tmp_path: Path) -> None:
    tabix = tmp_path / "tabix"
    tabix.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        "if any(arg.endswith('fail.vcf.gz') for arg in sys.argv):\n"
        "    print('query failed', file=sys.stderr)\n"
        "    raise SystemExit(2)\n"
        "print('first')\n"
        "print('second')\n"
    )
    tabix.chmod(0o755)

    with clinvar_validation.tabix_output_lines(
        str(tabix),
        tmp_path / "ok.vcf.gz",
        tmp_path / "regions.bed",
    ) as lines:
        assert list(lines) == ["first\n", "second\n"]

    with pytest.raises(RuntimeError, match="query failed"):
        with clinvar_validation.tabix_output_lines(
            str(tabix),
            tmp_path / "fail.vcf.gz",
            tmp_path / "regions.bed",
        ) as lines:
            list(lines)


def test_observed_clinvar_memberships_are_cached_and_reused(
    tmp_path: Path,
    monkeypatch,
) -> None:
    universe_path = tmp_path / "clinvar_universe.tsv.gz"
    annotations_path = tmp_path / "variant_annotations.tsv.gz"
    universe = pd.DataFrame(
        [
            {"variant_key": "1:10:A>G", "variant_type": "snv"},
            {"variant_key": "1:20:A>AT", "variant_type": "indel"},
        ]
    )
    universe.to_csv(universe_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(
        [
            ["1:10:A>G", "gene_a", "snv", "A", "G", "ok", "s1"],
            ["1:10:A>G", "gene_b", "snv", "A", "G", "ok", "s2"],
            ["1:20:A>AT", "gene_a", "ins", "A", "AT", "ok", "s2"],
            ["1:30:C>T", "gene_a", "snv", "C", "T", "ok", "s1"],
        ],
        columns=[
            "variant_key",
            "gene_id",
            "event_type",
            "ref",
            "alt",
            "lookup_status",
            "strategies",
        ],
    ).to_csv(annotations_path, sep="\t", index=False, compression="gzip")
    observed_store = build_or_load_observed_variant_store(
        variant_annotations_source=annotations_path,
        analytics_dir=tmp_path / "analytics",
        strategies=["s1", "s2"],
    )

    observed, manifest, output_path, manifest_path = (
        clinvar_validation.build_or_load_observed_keys_by_strategy_type(
            universe=universe,
            universe_path=universe_path,
            observed_store=observed_store,
            strategies=["s1", "s2"],
            analytics_dir=tmp_path / "analytics",
        )
    )

    assert not manifest["cache_hit"]
    assert manifest["membership_count"] == 3
    assert observed[("s1", "snv")] == {"1:10:A>G"}
    assert observed[("s1", "indel")] == set()
    assert observed[("s2", "snv")] == {"1:10:A>G"}
    assert observed[("s2", "indel")] == {"1:20:A>AT"}
    assert "observed_store" in manifest["inputs"]
    assert output_path.exists()
    assert manifest_path.exists()

    monkeypatch.setattr(
        clinvar_validation,
        "collect_observed_keys_by_strategy_type",
        lambda **_kwargs: (_ for _ in ()).throw(
            AssertionError("annotation scan was not cached")
        ),
    )
    cached, cached_manifest, _output_path, _manifest_path = (
        clinvar_validation.build_or_load_observed_keys_by_strategy_type(
            universe=universe,
            universe_path=universe_path,
            observed_store=observed_store,
            strategies=["s1", "s2"],
            analytics_dir=tmp_path / "analytics",
        )
    )

    assert cached_manifest["cache_hit"]
    assert cached == observed
