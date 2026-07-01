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


class JoinBaselineTest(unittest.TestCase):
    def test_join_adds_baseline_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            features = Path(tmp) / "features.tsv"
            joined = Path(tmp) / "joined.tsv"
            env = {"PYTHONPATH": str(ROOT / "src")}
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "cadd_validation.build_features",
                    "--variants-tsv",
                    str(FIXTURES / "variants.tsv"),
                    "--segments-tsv",
                    str(FIXTURES / "segments.tsv"),
                    "--events-tsv",
                    str(FIXTURES / "events.tsv"),
                    "--summaries-tsv",
                    str(FIXTURES / "summaries.tsv"),
                    "--taxonomy-presets-tsv",
                    str(FIXTURES / "taxonomy_presets.tsv"),
                    "--strategies",
                    "minimap2_asm20",
                    "--out-tsv",
                    str(features),
                ],
                check=True,
                cwd=ROOT,
                env=env,
            )
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "cadd_validation.join_baseline",
                    "--gaph-features-tsv",
                    str(features),
                    "--baseline-tsv",
                    str(FIXTURES / "baseline.tsv"),
                    "--out-tsv",
                    str(joined),
                ],
                check=True,
                cwd=ROOT,
                env=env,
            )
            with joined.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 2)
        self.assertIn("CADD_RAW", rows[0])
        self.assertIn("baseline_label", rows[0])


if __name__ == "__main__":
    unittest.main()

