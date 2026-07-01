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


class BuildFeaturesTest(unittest.TestCase):
    def test_counts_ref_alt_other_and_taxonomy_groups(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "features.tsv"
            summary = Path(tmp) / "summary.json"
            cmd = [
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
                str(out),
                "--summary-json",
                str(summary),
            ]
            subprocess.run(cmd, check=True, cwd=ROOT, env={"PYTHONPATH": str(ROOT / "src")})
            with out.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 2)
        alt = rows[0]
        self.assertEqual(alt["variant_id"], "var_alt")
        self.assertEqual(alt["gaph_all_depth"], "3")
        self.assertEqual(alt["gaph_all_alt_count"], "1")
        self.assertEqual(alt["gaph_all_other_count"], "1")
        self.assertEqual(alt["gaph_all_ref_count"], "1")
        self.assertEqual(alt["gaph_all_no_call_count"], "1")
        self.assertEqual(alt["gaph_primates_alt_count"], "1")
        self.assertEqual(alt["gaph_other_mammals_other_count"], "1")
        self.assertEqual(alt["gaph_non_mammal_vertebrates_ref_count"], "1")


if __name__ == "__main__":
    unittest.main()

