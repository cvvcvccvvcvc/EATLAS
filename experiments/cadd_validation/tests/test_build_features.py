from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


BASE = Path(__file__).resolve().parent
FIXTURES = BASE / "fixtures"
ROOT = BASE.parent
PROJECT_ROOT = ROOT.parents[1]
TEST_PYTHONPATH = os.pathsep.join((str(PROJECT_ROOT), str(ROOT / "src")))
sys.path.insert(0, str(ROOT / "src"))

from cadd_validation.build_features import group_for_tax, load_taxonomy_groups  # noqa: E402


class BuildFeaturesTest(unittest.TestCase):
    def test_canonical_taxonomy_lineage_defines_groups(self) -> None:
        groups = load_taxonomy_groups(FIXTURES / "taxonomy.tsv")

        self.assertEqual(groups["9443"], "primates")
        self.assertEqual(groups["10090"], "other_mammals")
        self.assertEqual(groups["7955"], "non_mammal_vertebrates")

    def test_noncanonical_taxonomy_columns_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            taxonomy = Path(tmp) / "taxonomy.tsv"
            taxonomy.write_text(
                "tax_id\tis_primate\tis_mammal\tis_vertebrate\tparent_tax_ids\n"
                "9598\ttrue\ttrue\ttrue\tnone\n"
            )

            with self.assertRaisesRegex(ValueError, "lineage_tax_ids"):
                load_taxonomy_groups(taxonomy)

    def test_tax_id_absent_from_canonical_taxonomy_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "absent from canonical taxonomy"):
            group_for_tax("999999", {"9598": "primates"})

    def test_taxonomy_argument_is_required(self) -> None:
        completed = subprocess.run(
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
                "--out-tsv",
                "unused.tsv",
            ],
            cwd=ROOT,
            env={"PYTHONPATH": TEST_PYTHONPATH},
            capture_output=True,
            text=True,
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("--taxonomy-tsv", completed.stderr)

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
                "--taxonomy-tsv",
                str(FIXTURES / "taxonomy.tsv"),
                "--strategies",
                "minimap2_asm20",
                "--out-tsv",
                str(out),
                "--summary-json",
                str(summary),
            ]
            subprocess.run(cmd, check=True, cwd=ROOT, env={"PYTHONPATH": TEST_PYTHONPATH})
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
