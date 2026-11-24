from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from integrate_wgs_rna.realign.config import PersonalRealignConfig
from integrate_wgs_rna.realign.qc import scan_realigned_output


class ParallelQcTest(unittest.TestCase):
    def test_indexed_bam_and_cram_match_serial_scan(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            reference = root / "reference.fa"
            reference.write_text(">chr1\n" + "A" * 1000 + "\n>chr2\n" + "C" * 1000 + "\n")
            pysam.faidx(str(reference))
            header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [
                {"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000},
            ]}
            records = [
                ("mapped", 0, 0, 10, 1, False, False),
                ("duplicate", 1024, 0, 20, 2, False, True),
                ("secondary", 256, 0, 30, 1, False, False),
                ("supplementary", 2048, 0, 40, 2, False, False),
                ("placed_unmapped", 4, 0, 50, 0, False, False),
                ("qcfail", 512, 1, 10, 0, False, False),
                ("rescue", 0, 1, 20, 0, True, False),
                ("unplaced", 4, -1, -1, 0, True, False),
            ]
            for suffix, mode in (("bam", "wb"), ("cram", "wc")):
                output = root / f"output.{suffix}"
                kwargs = {"reference_filename": str(reference)} if suffix == "cram" else {}
                with pysam.AlignmentFile(str(output), mode, header=header, **kwargs) as handle:
                    for name, flag, tid, pos, hp, rescue, od in records:
                        read = pysam.AlignedSegment(handle.header)
                        read.query_name = name
                        read.query_sequence = ("C" if tid == 1 else "A") * 10
                        read.query_qualities = pysam.qualitystring_to_array("I" * 10)
                        read.flag = flag
                        read.reference_id = tid
                        read.reference_start = pos
                        read.mapping_quality = 0 if flag & 4 else 60
                        if not flag & 4:
                            read.cigar = ((0, 10),)
                        read.set_tag("HP", hp)
                        read.set_tag("NM", 1)
                        if rescue:
                            read.set_tag("XU", 1)
                        if od:
                            read.set_tag("OD", 1)
                        handle.write(read)
                config = PersonalRealignConfig(
                    cohort="test", output_root=root, sample_list=None, samples=["S"],
                    chroms=["chr1", "chr2"], sequence_variant_policy="snv_indel",
                    dataset_role="test", output_format=suffix, sort_threads=2,
                )
                serial = scan_realigned_output(config, output, reference, 1)
                self.assertEqual(serial["primary"], 6)
                self.assertEqual(serial["secondary"], 1)
                self.assertEqual(serial["supplementary"], 1)
                self.assertEqual(serial["rescue_primary"], 2)
                # Without an index, requesting parallel workers retains the serial path.
                self.assertEqual(scan_realigned_output(config, output, reference, 1, workers=3), serial)
                pysam.index(str(output))
                self.assertEqual(scan_realigned_output(config, output, reference, 1, workers=3), serial)


if __name__ == "__main__":
    unittest.main()
