from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "13_lift_annotation.py"
SPEC = importlib.util.spec_from_file_location("lift_annotation", SCRIPT)
assert SPEC and SPEC.loader
LIFT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LIFT)


def _fixture(root: Path) -> tuple[Path, Path]:
    """참조 6–9 가 결실 구간에 걸리는 chain + exon 하나, 그리고 대상 외 contig 하나."""
    chains = root / "chains"
    chains.mkdir()
    (chains / "chr1.chain").write_text(
        "chain 100 chr1 10 + 0 10 chr1 11 + 0 11 1\n5 1 2\n4\n"
    )
    gtf = root / "input.gtf"
    gtf.write_text(
        "chr1\ttest\texon\t6\t9\t.\t+\t.\tgene_id \"G1\";\n"
        "chrX\ttest\texon\t1\t3\t.\t+\t.\tgene_id \"GX\";\n"
    )
    return chains, gtf


def _run(chains: Path, gtf: Path, out_gtf: Path, out_exons: Path, qc_json: Path,
         *extra: str) -> None:
    before = sys.argv[:]
    try:
        sys.argv = [
            str(SCRIPT), "--chain-dir", str(chains), "--gtf", str(gtf),
            "--out-gtf", str(out_gtf), "--out-exons", str(out_exons),
            "--qc-json", str(qc_json), "--contigs", "chr1", *extra,
        ]
        assert LIFT.main() == 0
    finally:
        sys.argv = before


class LiftAnnotationQcTest(unittest.TestCase):
    def test_default_drops_unlifted_features(self) -> None:
        """기본값 drop — AlleleSeq2 는 lift 되지 않은 좌표를 버린다.

        스냅은 HAP1/HAP2 에서 비대칭으로 일어날 수 있어 ASE 에서 위험하다.
        """
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            chains, gtf = _fixture(root)
            out_gtf, out_exons, qc_json = root / "o.gtf", root / "e.tsv", root / "q.json"
            _run(chains, gtf, out_gtf, out_exons, qc_json)

            qc = json.loads(qc_json.read_text())
            self.assertEqual(qc["schema"], "mfm_gtf_lift_qc_v2")
            self.assertEqual(qc["on_unlifted"], "drop")
            self.assertEqual(qc["features_input_total"], 2)
            self.assertEqual(qc["features_input_target_contigs"], 1)   # chrX 는 대상 외
            self.assertEqual(qc["features_output"], 0)
            self.assertEqual(qc["features_dropped_unlifted"], 1)
            self.assertEqual(qc["features_snapped"], 0)

    def test_snap_mode_keeps_feature_and_flags_it(self) -> None:
        """--on-unlifted snap 은 이전 동작 — 버리지 않고 플래그만 남긴다."""
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            chains, gtf = _fixture(root)
            out_gtf, out_exons, qc_json = root / "o.gtf", root / "e.tsv", root / "q.json"
            _run(chains, gtf, out_gtf, out_exons, qc_json, "--on-unlifted", "snap")

            qc = json.loads(qc_json.read_text())
            self.assertEqual(qc["on_unlifted"], "snap")
            self.assertEqual(qc["features_output"], 1)
            self.assertEqual(qc["features_snapped"], 1)
            self.assertEqual(qc["features_dropped_unlifted"], 0)
            self.assertEqual(qc["snap_fraction_of_output"], 1.0)
            self.assertIn('lift_snapped "1"', out_gtf.read_text())
            self.assertTrue(out_exons.is_file())


if __name__ == "__main__":
    unittest.main()
