from __future__ import annotations

import unittest

from cadd_validation.fetch_cadd_scores import cadd_chrom, normalize_response


class FetchCaddScoresTest(unittest.TestCase):
    def test_ncbi_accession_mapping(self) -> None:
        self.assertEqual(cadd_chrom("NC_000001.11"), "1")
        self.assertEqual(cadd_chrom("NC_000007.14"), "7")
        self.assertEqual(cadd_chrom("NC_000019.10"), "19")
        self.assertEqual(cadd_chrom("NC_000023.11"), "X")
        self.assertEqual(cadd_chrom("chr2"), "2")

    def test_json_response_parser(self) -> None:
        raw, phred, status = normalize_response('[{"RawScore": "1.25", "PHRED": "12.3"}]')
        self.assertEqual(status, "ok")
        self.assertEqual(raw, "1.25")
        self.assertEqual(phred, "12.3")

    def test_tsv_response_parser(self) -> None:
        raw, phred, status = normalize_response("Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n1\t10\tA\tG\t0.4\t8.2\n")
        self.assertEqual(status, "ok")
        self.assertEqual(raw, "0.4")
        self.assertEqual(phred, "8.2")


if __name__ == "__main__":
    unittest.main()
