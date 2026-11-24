from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from integrate_wgs_rna.pipeline.manifest_progress import (
    ManifestProgressJournal,
    manifest_progress_signature,
)


def _result(order: int, sample: str, status: str = "complete") -> dict[str, object]:
    return {
        "manifest_order": order,
        "sample_id": sample,
        "status": status,
        "output_sample_root": f"/output/{sample}",
    }


class ManifestProgressJournalTest(unittest.TestCase):
    def test_signature_is_stable_and_covers_rows_and_config(self) -> None:
        rows = [{"sample_id": "S1", "_manifest_order": "0"}]
        config = {"jobs": 4, "status_tsv": Path("/tmp/status.tsv")}
        signature = manifest_progress_signature(rows, config)
        self.assertEqual(signature, manifest_progress_signature(rows, config))
        self.assertNotEqual(
            signature,
            manifest_progress_signature(rows + [{"sample_id": "S2"}], config),
        )
        self.assertNotEqual(
            signature,
            manifest_progress_signature(rows, config | {"jobs": 8}),
        )

    def test_resume_order_retry_reset_and_cleanup(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "status.progress.sqlite3"
            journal = ManifestProgressJournal(path, "signature-a")
            journal.record(_result(2, "S3"))
            journal.record(_result(0, "S1"))
            journal.record(_result(1, "S2", "failed"))
            self.assertEqual(
                [item["sample_id"] for item in journal.results()],
                ["S1", "S2", "S3"],
            )
            self.assertEqual(
                [item["sample_id"] for item in journal.results(completion_order=True)],
                ["S3", "S1", "S2"],
            )
            self.assertEqual(journal.completed_sample_ids(), {"S1", "S3"})
            journal.close()

            resumed = ManifestProgressJournal(path, "signature-a")
            self.assertEqual(resumed.completed_sample_ids(), {"S1", "S3"})
            resumed.record(_result(1, "S2"))
            self.assertEqual(
                [item["sample_id"] for item in resumed.results(completion_order=True)],
                ["S3", "S1", "S2"],
            )
            resumed.close()

            reset = ManifestProgressJournal(path, "signature-b")
            self.assertEqual(reset.results(), [])
            reset.record(_result(0, "NEW"))
            reset.cleanup()
            self.assertFalse(path.exists())
            self.assertFalse(Path(f"{path}-wal").exists())
            self.assertFalse(Path(f"{path}-shm").exists())


if __name__ == "__main__":
    unittest.main()
