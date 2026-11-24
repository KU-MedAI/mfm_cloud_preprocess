#!/usr/bin/env python3
"""현행·구번호 RNA entrypoint가 공통 helper를 공유하되 정책 경계는 유지하는지 확인한다."""
from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import catalog_common
import rna_allelic_common
import rna_ase_stats
import rna_interval_common
import rna_metadata_common


def load_script(name: str):
    path = ROOT / name
    spec = importlib.util.spec_from_file_location(f"test_{path.stem}", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestSharedRnaHelpers(unittest.TestCase):
    def test_allelic_helpers_share_but_stage21_policy_stays_separate(self):
        old16 = load_script("16_allelic_counts.py")
        old20 = load_script("20_allelic_counts.py")
        current21 = load_script("21_allelic_counts.py")
        for module in (old16, old20, current21):
            self.assertIs(module.split_bucket, rna_allelic_common.split_bucket)
            self.assertIs(module.load_exons, rna_allelic_common.load_exons)
            self.assertIs(module.genes_for_indexed, rna_allelic_common.genes_for_indexed)
        self.assertIs(old16.assign_fragments, rna_allelic_common.assign_fragments)
        self.assertIs(old20.assign_fragments, rna_allelic_common.assign_fragments)
        self.assertIsNot(current21.assign_fragments, rna_allelic_common.assign_fragments)

    def test_ase_helpers_preserve_legacy_fit_mode(self):
        old18 = load_script("18_ase_qc.py")
        old22 = load_script("22_ase_qc.py")
        current23 = load_script("23_ase_qc.py")
        self.assertIs(old22.binom_two_sided, rna_ase_stats.binom_two_sided)
        self.assertIs(current23.betabinom_two_sided, rna_ase_stats.betabinom_two_sided)
        rows = [(1, 2), (2, 4), (3, 6), (12, 20)] * 80
        self.assertEqual(
            old18.fit_overdispersion(rows),
            rna_ase_stats.fit_overdispersion(rows, include_low_depth_null=True),
        )
        self.assertEqual(old22.fit_overdispersion(rows), rna_ase_stats.fit_overdispersion(rows))

    def test_other_exact_helpers_are_reused(self):
        self.assertIs(load_script("17_cnv_filter.py").merge, rna_interval_common.merge)
        self.assertIs(load_script("21_cnv_filter.py").merge, rna_interval_common.merge)
        self.assertIs(load_script("22_cnv_filter.py").merge, rna_interval_common.merge)
        for name in ("19_stamp_rna_metadata.py", "23_stamp_rna_metadata.py", "24_stamp_rna_metadata.py"):
            self.assertTrue(callable(load_script(name).stamp))
        for name in ("21_make_catalog.py", "27_make_catalog.py"):
            self.assertIs(load_script(name).read_json, catalog_common.read_json)
            self.assertIs(load_script(name).first_existing, catalog_common.first_existing)


if __name__ == "__main__":
    unittest.main()
