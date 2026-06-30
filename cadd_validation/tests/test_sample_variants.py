from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


BASE = Path(__file__).resolve().parent
ROOT = BASE.parent


class SampleVariantsTest(unittest.TestCase):
    def test_limits_each_group(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "variants.tsv"
            out = Path(tmp) / "sampled.tsv"
            with source.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["variant_id", "gene_id", "label"],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for gene in ["g1", "g2"]:
                    for label in ["benign", "pathogenic"]:
                        for index in range(5):
                            writer.writerow({"variant_id": f"{gene}_{label}_{index}", "gene_id": gene, "label": label})
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "cadd_validation.sample_variants",
                    "--variants-tsv",
                    str(source),
                    "--out-tsv",
                    str(out),
                    "--max-per-group",
                    "2",
                ],
                check=True,
                cwd=ROOT,
                env={"PYTHONPATH": str(ROOT / "src")},
            )
            with out.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 8)


if __name__ == "__main__":
    unittest.main()

