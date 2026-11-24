from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from integrate_wgs_rna import personal_realign as MODULE


class PlatformPresetTest(unittest.TestCase):
    def test_clr_uses_map_pb_and_pacbio_read_group(self) -> None:
        self.assertEqual(MODULE.platform_minimap2_preset("clr", "auto"), "map-pb")
        self.assertEqual(MODULE.platform_minimap2_preset("pacbio_clr", "auto"), "map-pb")
        self.assertEqual(MODULE.sam_rg_platform("clr"), "PACBIO")

    def test_hifi_and_ont_are_unchanged(self) -> None:
        self.assertEqual(MODULE.platform_minimap2_preset("hifi", "auto"), "map-hifi")
        self.assertEqual(MODULE.platform_minimap2_preset("ont", "auto"), "map-ont")


if __name__ == "__main__":
    unittest.main()
