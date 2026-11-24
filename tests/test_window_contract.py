"""구형 cache 정규화와 좌표 변조 차단을 검사한다."""
import csv
import tempfile
import unittest
from pathlib import Path

from integrate_wgs_rna.window_contract import ensure_window_contract, read_windows


class WindowContractTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name) / "windows.tsv"
        self.header = "contig\thaplotype\twindow_index\traw_start\tactual_length\tkept\n"

    def tearDown(self):
        self.temp.cleanup()

    def test_upgrade_preserves_original_cells_and_is_idempotent(self):
        text = self.header + "chr1\tHAP1\t0\t10\t100\t1\nchr1\tHAP1\t1\t10\t100\t0\n"
        self.path.write_text(text)
        self.assertTrue(ensure_window_contract(self.path))
        self.assertFalse(ensure_window_contract(self.path))
        backups = list(self.path.parent.glob("*.legacy.*"))
        self.assertEqual(len(backups), 1)
        self.assertEqual(backups[0].read_text(), text)
        original = list(csv.DictReader(text.splitlines(), delimiter="\t"))
        current = read_windows(self.path)
        self.assertEqual([{key: row[key] for key in original[0]} for row in current], original)
        self.assertEqual(current[0]["personal_end"], "110")

    def test_mismatching_coordinate_representations_are_rejected(self):
        self.path.write_text(self.header.rstrip("\n") + "\tpersonal_start\tpersonal_end\n"
                             + "chr1\tHAP1\t0\t10\t100\t1\t11\t110\n")
        with self.assertRaisesRegex(ValueError, "불일치"):
            read_windows(self.path)

    def test_kept_overlap_duplicate_and_invalid_bounds_are_rejected(self):
        for rows in (
            "chr1\tHAP1\t0\t10\t100\t1\nchr1\tHAP1\t1\t100\t100\t1\n",
            "chr1\tHAP1\t0\t10\t100\t1\nchr1\tHAP1\t0\t110\t100\t1\n",
            "chr1\tHAP1\t0\t-1\t100\t1\n",
            "chr1\tHAP1\t0\t0\t0\t1\n",
        ):
            with self.subTest(rows=rows):
                self.path.write_text(self.header + rows)
                with self.assertRaises(ValueError):
                    read_windows(self.path)


if __name__ == "__main__":
    unittest.main()
