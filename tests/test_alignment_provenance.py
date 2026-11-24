from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from integrate_wgs_rna.personal_realign import alignment_provenance
from integrate_wgs_rna.realign.config import PersonalRealignConfig


class AlignmentProvenanceTest(unittest.TestCase):
    def config(self, **kwargs) -> PersonalRealignConfig:
        values = {
            "cohort": "test",
            "output_root": Path("/data/project"),
            "sample_list": None,
            "samples": ["S1"],
            "chroms": ["chr1"],
            "sequence_variant_policy": "snv_indel",
            "dataset_role": "test",
        }
        values.update(kwargs)
        return PersonalRealignConfig(**values)

    def test_explicit_alignment_is_reported_as_actual_input(self):
        path = Path("/data/project/S1/haplotagged/S1.illumina.bam")
        self.assertEqual(
            alignment_provenance(self.config(alignment=path)),
            {
                "alignment_source": "explicit_path",
                "alignment_pattern": "",
                "alignment_path": str(path),
            },
        )

    def test_pattern_is_reported_when_no_explicit_alignment_is_given(self):
        config = self.config()
        self.assertEqual(
            alignment_provenance(config),
            {
                "alignment_source": "pattern",
                "alignment_pattern": config.alignment_pattern,
                "alignment_path": "",
            },
        )


if __name__ == "__main__":
    unittest.main()
