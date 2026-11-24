from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from integrate_wgs_rna import full_preprocess as full


def _result(row: dict[str, str]) -> dict[str, object]:
    sample = row["sample_id"]
    return {
        "manifest_order": int(row["_manifest_order"]),
        "sample_id": sample,
        "status": "complete",
        "elapsed_seconds": "1.000",
        "elapsed_hms": "00:00:01",
        "output_sample_root": f"/output/{sample}",
        "size_bytes": 1,
        "size_human": "1.00 B",
        "error": "",
    }


class FullManifestResumeTest(unittest.TestCase):
    def test_interrupted_run_reuses_only_completed_rows(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            manifest = root / "manifest.tsv"
            manifest.write_text("sample_id\nS1\nS2\n", encoding="utf-8")
            status = root / "status.tsv"
            config = full.FullManifestConfig(manifest=manifest, status_tsv=status)
            first_calls: list[str] = []

            def interrupted(row: dict[str, str], _: full.FullManifestConfig):
                first_calls.append(row["sample_id"])
                if row["sample_id"] == "S2":
                    raise RuntimeError("intentional interruption")
                return _result(row)

            with patch.object(full, "run_full_manifest_row", side_effect=interrupted):
                with self.assertRaisesRegex(RuntimeError, "intentional interruption"):
                    full.run_full_manifest(config)
            self.assertEqual(first_calls, ["S1", "S2"])
            journal = Path(f"{status}.progress.sqlite3")
            self.assertTrue(journal.exists())
            self.assertFalse(status.exists())

            resumed_calls: list[str] = []

            def resumed(row: dict[str, str], _: full.FullManifestConfig):
                resumed_calls.append(row["sample_id"])
                return _result(row)

            with patch.object(full, "run_full_manifest_row", side_effect=resumed):
                results = full.run_full_manifest(config)
            self.assertEqual(resumed_calls, ["S2"])
            self.assertEqual([item["sample_id"] for item in results], ["S1", "S2"])
            self.assertFalse(journal.exists())
            with status.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["sample_id"] for row in rows], ["S1", "S2"])


if __name__ == "__main__":
    unittest.main()
