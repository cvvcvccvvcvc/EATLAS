from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


BASE = Path(__file__).resolve().parent
ROOT = BASE.parent


class SelectClinvarGenesTest(unittest.TestCase):
    def test_selects_balanced_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp) / "clinvar.vcf"
            vcf.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                        "NC_000001.11\t101\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=GENE1:1",
                        "NC_000001.11\t102\t2\tC\tT\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=GENE1:1",
                        "NC_000001.11\t103\t3\tG\tA\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=GENE2:2",
                    ]
                )
                + "\n"
            )
            genes = Path(tmp) / "genes.txt"
            summary = Path(tmp) / "genes.tsv"
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "cadd_validation.select_clinvar_genes",
                    "--clinvar-vcf",
                    str(vcf),
                    "--out-genes",
                    str(genes),
                    "--out-gene-summary",
                    str(summary),
                    "--min-per-label",
                    "1",
                    "--max-genes",
                    "5",
                ],
                check=True,
                cwd=ROOT,
                env={"PYTHONPATH": str(ROOT / "src")},
            )
            self.assertEqual(genes.read_text(), "1\n")
            with summary.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["pathogenic_count"], "1")
            self.assertEqual(rows[0]["benign_count"], "1")


if __name__ == "__main__":
    unittest.main()

