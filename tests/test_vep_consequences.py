from __future__ import annotations

from pathlib import Path

import pandas as pd

from analytics.core import vep_consequences as vep


def test_parse_record_selects_target_gene_and_most_severe_term() -> None:
    row = {"variant_key": "1:10:A>G", "gene_id": "25"}
    record = {
        "variant_class": "SNV",
        "transcript_consequences": [
            {
                "gene_id": "999",
                "transcript_id": "NM_other",
                "consequence_terms": ["stop_gained"],
            },
            {
                "gene_id": "25",
                "transcript_id": "NM_target",
                "consequence_terms": ["splice_region_variant", "missense_variant"],
                "impact": "MODERATE",
                "canonical": 1,
                "mane_select": "NM_target.1",
            },
        ],
    }

    result = vep._parse_record(row, record)

    assert result["status"] == "ok"
    assert result["primary_consequence"] == "missense_variant"
    assert result["consequence_terms"] == "missense_variant&splice_region_variant"
    assert result["transcript_id"] == "NM_target"


def test_sqlite_cache_reuses_completed_annotations(tmp_path: Path, monkeypatch) -> None:
    rows = pd.DataFrame(
        [
            {"variant_key": "1:10:A>G", "gene_id": "25", "chrom": "1", "pos": 10, "ref": "A", "alt": "G"},
            {"variant_key": "1:10:A>G", "gene_id": "25", "chrom": "1", "pos": 10, "ref": "A", "alt": "G"},
        ]
    )
    calls = []

    def fake_request(batch, *_args):
        calls.append(batch)
        return [
            {
                "variant_key": item["variant_key"],
                "gene_id": item["gene_id"],
                "status": "ok",
                "primary_consequence": "missense_variant",
                "consequence_terms": "missense_variant",
                "transcript_id": "NM_1",
                "mane_select": "NM_1",
                "canonical": True,
                "impact": "MODERATE",
                "variant_class": "SNV",
            }
            for item in batch
        ]

    monkeypatch.setattr(vep, "_request_batch", fake_request)
    cache = tmp_path / "vep.sqlite"
    first, first_summary = vep.annotate_vep_consequences(rows, cache, release="116")
    second, second_summary = vep.annotate_vep_consequences(rows, cache, release="116")

    assert len(calls) == 1
    assert len(calls[0]) == 1
    assert first_summary["queried"] == 1
    assert second_summary["cached"] == 1
    pd.testing.assert_frame_equal(first, second)
