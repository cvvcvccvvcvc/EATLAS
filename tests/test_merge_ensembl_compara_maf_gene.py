from __future__ import annotations

import csv
import gzip
from pathlib import Path

from bin.merge_ensembl_compara_maf_gene import merge_tsv_gz


def test_merge_tsv_gz_accepts_large_fields(tmp_path: Path) -> None:
    fields = ["gene_id", "alt"]
    large_alt = "A" * 165_969
    source = tmp_path / "source.tsv.gz"
    output = tmp_path / "output.tsv.gz"

    with gzip.open(source, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({"gene_id": "3492", "alt": large_alt})

    assert merge_tsv_gz([source], output, fields) == 1

    with gzip.open(output, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert rows == [{"gene_id": "3492", "alt": large_alt}]
