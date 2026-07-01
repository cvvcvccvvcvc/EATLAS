from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


BASE = Path(__file__).resolve().parent
ROOT = BASE.parent


class EvaluateAblationTest(unittest.TestCase):
    def test_ablation_writes_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dataset = Path(tmp) / "dataset.tsv"
            outdir = Path(tmp) / "eval"
            with dataset.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["variant_id", "gene_id", "label", "CADD_RAW", "gaph_all_alt_fraction"],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for index in range(12):
                    label = "pathogenic" if index % 2 else "benign"
                    writer.writerow(
                        {
                            "variant_id": f"v{index}",
                            "gene_id": f"g{index}",
                            "label": label,
                            "CADD_RAW": 2.0 if label == "pathogenic" else -1.0,
                            "gaph_all_alt_fraction": 0.8 if label == "pathogenic" else 0.1,
                        }
                    )
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "cadd_validation.evaluate_ablation",
                    "--dataset-tsv",
                    str(dataset),
                    "--outdir",
                    str(outdir),
                    "--baseline-features",
                    "CADD_RAW",
                    "--group-column",
                    "gene_id",
                    "--folds",
                    "3",
                    "--bootstrap-iterations",
                    "10",
                ],
                check=True,
                cwd=ROOT,
                env={"PYTHONPATH": str(ROOT / "src")},
            )
            with (outdir / "metrics.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["feature_set"] for row in rows}, {
            "baseline",
            "gaph",
            "baseline_plus_gaph",
            "baseline_plus_shuffled_gaph",
        })
        self.assertIn("auprc_oof", rows[0])


if __name__ == "__main__":
    unittest.main()

