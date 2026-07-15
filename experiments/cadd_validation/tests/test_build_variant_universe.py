from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


BASE = Path(__file__).resolve().parent
FIXTURES = BASE / "fixtures"
ROOT = BASE.parent


class BuildVariantUniverseTest(unittest.TestCase):
    def test_clinvar_snv_universe(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "variants.tsv"
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "cadd_validation.build_variant_universe",
                    "--clinvar-vcf",
                    str(FIXTURES / "clinvar.vcf"),
                    "--target-features-tsv",
                    str(FIXTURES / "target_features.tsv"),
                    "--out-tsv",
                    str(out),
                ],
                check=True,
                cwd=ROOT,
                env={"PYTHONPATH": str(ROOT / "src")},
            )
            with out.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["label"], "pathogenic")
        self.assertEqual(rows[0]["target_start0"], "0")
        self.assertEqual(rows[1]["label"], "benign")


if __name__ == "__main__":
    unittest.main()

