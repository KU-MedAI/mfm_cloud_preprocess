from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from integrate_wgs_rna.cohort_priority import (
    PAIR,
    WGS_ONLY,
    REQUIRED_COLUMNS,
    default_targets,
    select_candidates,
    write_selection,
)


class CohortPriorityTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def row(
        self,
        sample: str,
        dataset_role: str,
        *,
        platform: str = "ont",
        rna_index: bool = True,
    ) -> dict[str, str]:
        stem = f"{sample}.{dataset_role}"
        alignment = self.root / f"{stem}.cram"
        alignment.touch()
        alignment.with_suffix(".crai").touch()

        is_pair = dataset_role == "post_training"
        rna = self.root / f"{stem}.rna.bam"
        if is_pair:
            rna.touch()
            if rna_index:
                rna.with_suffix(".bai").touch()

        return {
            "sample_id": sample,
            "dataset_role": dataset_role,
            "platform": platform,
            "alignment_path": str(alignment),
            "alignment_exists": "1",
            "alignment_has_index": "1",
            "rna_path": str(rna) if is_pair else "-",
            "rna_exists": "1" if is_pair else "0",
            "rna_has_index": "1" if is_pair and rna_index else "0",
        }

    @staticmethod
    def targets(pair: int = 1, wgs: int = 1) -> dict[str, dict[str, int]]:
        return {
            PAIR: {"ont": pair, "illumina": 0, "hifi": 0},
            WGS_ONLY: {"ont": wgs, "illumina": 0, "hifi": 0},
        }

    def test_production_defaults_prioritise_hifi_and_illumina(self) -> None:
        targets = default_targets()
        self.assertEqual(targets[PAIR], {"ont": 0, "illumina": 5, "hifi": 5})
        self.assertEqual(targets[WGS_ONLY], {"ont": 0, "illumina": 5, "hifi": 5})

    def test_pair_donor_is_reserved_before_wgs_only(self) -> None:
        rows = [
            self.row("SAME", "post_training"),
            self.row("SAME", "pre_training"),
            self.row("WGS_ONLY", "pre_training"),
        ]
        result = select_candidates(rows, self.targets(), seed=42)
        selected = [(row["selection_role"], row["sample_id"]) for row in result["selected"]]
        self.assertEqual(selected, [(PAIR, "SAME"), (WGS_ONLY, "WGS_ONLY")])
        self.assertEqual(result["qc"]["selected_unique_donors"], 2)
        self.assertIn("donor_already_selected", result["qc"]["rejection_counts"])

    def test_pair_without_rna_index_is_rejected(self) -> None:
        result = select_candidates(
            [self.row("NO_INDEX", "post_training", rna_index=False)],
            self.targets(pair=1, wgs=0),
            seed=42,
        )
        self.assertEqual(result["selected"], [])
        self.assertEqual(result["qc"]["summary"][0]["missing"], 1)
        self.assertEqual(result["qc"]["rejection_counts"], {"rna_index_not_verified": 1})

    def test_seeded_selection_is_reproducible(self) -> None:
        rows = [self.row(f"P{i}", "post_training") for i in range(4)]
        first = select_candidates(rows, self.targets(pair=2, wgs=0), seed=123)
        second = select_candidates(rows, self.targets(pair=2, wgs=0), seed=123)
        self.assertEqual(
            [row["sample_id"] for row in first["selected"]],
            [row["sample_id"] for row in second["selected"]],
        )

    def test_result_writes_separate_manifests_and_qc(self) -> None:
        rows = [self.row("PAIR", "post_training"), self.row("WGS", "pre_training")]
        result = select_candidates(rows, self.targets(), seed=42)
        paths = write_selection(self.root / "out", list(REQUIRED_COLUMNS), result)
        self.assertTrue(paths["pair_manifest"].is_file())
        self.assertTrue(paths["wgs_only_manifest"].is_file())
        self.assertTrue(paths["qc"].is_file())
        with paths["pair_manifest"].open(newline="") as handle:
            pair_rows = list(csv.DictReader(handle, delimiter="\t"))
        with paths["wgs_only_manifest"].open(newline="") as handle:
            wgs_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["sample_id"] for row in pair_rows], ["PAIR"])
        self.assertEqual([row["sample_id"] for row in wgs_rows], ["WGS"])


if __name__ == "__main__":
    unittest.main()
