from __future__ import annotations

import csv
import gzip
from pathlib import Path

from analytics.core.variant_summary import VARIANT_USECOLS, build_variant_summary


def test_variant_summary_accepts_compact_annotation_schema(tmp_path: Path) -> None:
    annotations = tmp_path / "variant_annotations.tsv.gz"
    rows = [
        {
            "variant_key": "1:100:A>G",
            "gene_id": "1",
            "event_type": "snv",
            "ref": "A",
            "alt": "G",
            "strategies": "s1,s2",
            "support_row_count": "4",
            "support_ortholog_count": "3",
            "clinvar_id": "VCV1",
            "clinvar_sig": "Benign",
            "clinvar_review_stars": "2",
            "gnomad_af": "0.01",
            "gnomad_csq": "synonymous_variant",
        },
        {
            "variant_key": "1:200:C>A",
            "gene_id": "1",
            "event_type": "snv",
            "ref": "C",
            "alt": "A",
            "strategies": "s1",
            "support_row_count": "1",
            "support_ortholog_count": "1",
            "clinvar_id": "VCV2",
            "clinvar_sig": "Pathogenic",
            "clinvar_review_stars": "3",
            "gnomad_af": "",
            "gnomad_csq": "missense_variant",
        },
        {
            "variant_key": "1:300:AT>A",
            "gene_id": "1",
            "event_type": "del",
            "ref": "AT",
            "alt": "A",
            "strategies": "s2",
            "support_row_count": "2",
            "support_ortholog_count": "2",
            "clinvar_id": "",
            "clinvar_sig": "",
            "clinvar_review_stars": "",
            "gnomad_af": "",
            "gnomad_csq": "",
        },
    ]
    with gzip.open(annotations, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=VARIANT_USECOLS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in VARIANT_USECOLS})

    summary = build_variant_summary(
        annotations,
        tmp_path / "analytics",
        strategy_label=lambda value: value,
    )

    assert summary.input_row_count == 3
    assert summary.unique_variant_count == 3
    assert summary.strategy_record_count == 4
    assert summary.strategies == ["s1", "s2"]
    by_strategy = summary.strategy_stats.set_index("Strategy")
    assert by_strategy.loc["s1", "Ti/Tv"] == 1.0
    assert by_strategy.loc["s2", "Ti/Tv"] == float("inf")
    assert summary.clinvar_found == 2
    assert summary.gnomad_found == 1
