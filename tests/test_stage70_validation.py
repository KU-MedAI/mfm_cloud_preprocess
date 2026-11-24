from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from integrate_wgs_rna import direct_compact_writer as writer_module
from integrate_wgs_rna import cli as cli_module


def write_fasta(path: Path, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">chr1\n{sequence}\n", encoding="ascii")
    pysam.faidx(str(path))


def write_cram(path: Path, reference: Path, sequence: str) -> None:
    md5 = hashlib.md5(sequence.encode("ascii")).hexdigest()
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": len(sequence), "M5": md5}]}
    read = pysam.AlignedSegment()
    read.query_name = "r1"
    read.query_sequence = sequence[:4]
    read.query_qualities = pysam.qualitystring_to_array("IIII")
    read.flag = 0
    read.reference_id = 0
    read.reference_start = 0
    read.mapping_quality = 60
    read.cigar = ((0, 4),)
    with pysam.AlignmentFile(str(path), "wc", header=header, reference_filename=str(reference)) as output:
        output.write(read)
    pysam.index(str(path))


class Stage70ValidationTest(unittest.TestCase):
    def test_cli_exports_compact_validator(self) -> None:
        self.assertEqual(cli_module.CompactValidationConfig.__name__, "CompactValidationConfig")
        self.assertTrue(callable(cli_module.run_compact_validation))

    def _personal_job(self, root: Path, hap2_target_sequence: str = "ACGTACGT") -> dict:
        hap1 = root / "individual" / "HAP1" / "chr1.fa"
        hap2 = root / "individual" / "HAP2" / "chr1.fa"
        write_fasta(hap1, "ACGTACGT")
        write_fasta(hap2, "ACGTACGT")

        alignments: dict[int, Path] = {}
        for hap, target_sequence in ((1, "ACGTACGT"), (2, hap2_target_sequence)):
            out_dir = root / "realigned" / f"HAP{hap}"
            target = out_dir / "target.fa"
            write_fasta(target, target_sequence)
            cram = out_dir / f"HAP{hap}.cram"
            write_cram(cram, target, target_sequence)
            alignments[hap] = cram

        windows = root / "windows"
        windows.mkdir()
        (windows / "windows.tsv").write_text(
            "contig\thaplotype\twindow_index\traw_start\tactual_length\nchr1\tHAP1\t0\t0\t8\n",
            encoding="utf-8",
        )
        return {
            "contig": "chr1",
            "read_coordinate_mode": "personal_realigned",
            "validation_level": "full",
            "compact_compressor": "zstd",
            "compact_compression_level": 5,
            "config": {"window_size": 8, "offset_unit": 1, "curriculum_unit": 8, "seed": 42,
                       "min_mapping_quality": 30, "exclude_unmapped": True, "exclude_secondary": True,
                       "exclude_supplementary": False, "exclude_qcfail": True, "exclude_duplicate": False,
                       "require_proper_pair": False, "max_n_ratio": 0.0, "min_coordinate_coverage": 0.0},
            "window_index_dir": str(windows),
            "hap1_alignment": str(alignments[1]),
            "hap2_alignment": str(alignments[2]),
            "hap1_fasta": str(hap1),
            "hap2_fasta": str(hap2),
        }

    def test_personal_provenance_requires_same_hap_target_and_cram_axis(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            job = self._personal_job(Path(raw))
            provenance = writer_module._personal_input_provenance(job, {})
            self.assertEqual(provenance["schema"], writer_module.PERSONAL_INPUT_PROVENANCE_SCHEMA)
            self.assertEqual(len(provenance["fingerprint_sha256"]), 64)

        with tempfile.TemporaryDirectory() as raw:
            job = self._personal_job(Path(raw), hap2_target_sequence="ACGTTCGT")
            with self.assertRaisesRegex(ValueError, "HAP FASTA != CRAM target FASTA"):
                writer_module._personal_input_provenance(job, {})

    def test_resume_rejects_missing_or_changed_input_fingerprint(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            contig = Path(raw) / "chr1"
            contig.mkdir()
            expected = {"fingerprint_sha256": "a" * 64}
            (contig / "contig_metadata.json").write_text(json.dumps({"input_provenance": expected}), encoding="utf-8")
            writer_module._assert_resume_provenance(contig, expected)
            with self.assertRaisesRegex(ValueError, "fingerprint differs"):
                writer_module._assert_resume_provenance(contig, {"fingerprint_sha256": "b" * 64})

    def test_full_validation_recalculates_cigar_and_depth(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw) / "compact"
            compact = writer_module.DirectCompactChrWriter(
                root,
                final_target=root,
                sample="TEST",
                contig="chr1",
                hap_window_counts={"HAP1": 1},
                window_size=8,
                source_alignment=Path(raw) / "source.cram",
                reference_fasta=None,
                sample_metadata={},
                attrs={"read_coordinate_mode": "personal_realigned"},
            )
            compact.append_hap_chunk(
                "HAP1",
                {
                    "windows": [SimpleNamespace(sequence="ACGT" + "P" * 4, raw_start=0, actual_length=4)],
                    "records": [[{
                        "read_key": ("r1", 0),
                        "shared_read_key": ("r1", 0, 4, "4M", False, False),
                        "read_sequence": "ACGT",
                        "read_haplotype": 1,
                        "read_phase_set": 1,
                        "query_start": 0,
                        "query_end": 4,
                        "read_start": 0,
                        "read_end": 4,
                        "read_strand": 1,
                        "reference_start": 0,
                        "reference_end": 4,
                        "cigar_ops_packed": [4 << 4],
                        "segment_type": 0,
                    }]],
                    "attrs": {},
                    "counts": {},
                },
            )
            metadata = compact.finalize(validation_level="full")
            self.assertEqual(metadata["validation"]["checks"]["semantic"]["depth_windows"], 1)
            original_packed = int(compact.hap_arrays["HAP1"]["bases_packed"][0])
            compact.hap_arrays["HAP1"]["bases_packed"][0] = 0x11
            with self.assertRaisesRegex(ValueError, "HAP window sequence mismatch"):
                compact._validate_full_semantics()
            compact.hap_arrays["HAP1"]["bases_packed"][0] = original_packed
            compact.depth_arrays["HAP1"]["depth"][0] = 0
            with self.assertRaisesRegex(ValueError, "CIGAR-derived depth mismatch"):
                compact._validate_full_semantics()
            compact.depth_arrays["HAP1"]["depth"][0] = 1
            compact.cigar_ops[0] = 3 << 4
            with self.assertRaisesRegex(ValueError, "CIGAR/query mismatch"):
                compact._validate_full_semantics()
            compact.cigar_ops[0] = 4 << 4
            original_read_base = int(compact.read_arrays["bases_packed"][0])
            compact.read_arrays["bases_packed"][0] = 0xF1
            with self.assertRaisesRegex(ValueError, "invalid 4-bit read code"):
                compact._validate_full_semantics()
            compact.read_arrays["bases_packed"][0] = original_read_base
            compact.read_arrays["ref_end"][0] = 3
            with self.assertRaisesRegex(ValueError, "CIGAR/reference mismatch"):
                compact._validate_full_semantics()


    def test_long_read_preset_separates_primary_and_supplementary_depth(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw) / "compact"
            policy = {
                "window_size": 8, "offset_unit": 1, "curriculum_unit": 8, "seed": 42,
                "min_mapping_quality": 30, "exclude_unmapped": True, "exclude_secondary": True,
                "exclude_supplementary": False, "exclude_qcfail": True, "exclude_duplicate": False,
                "require_proper_pair": False, "max_n_ratio": 0.0, "min_coordinate_coverage": 0.0,
                "filter_preset": "long_read_primary_supp_evidence_v1",
                "depth_segment_types": (0,), "supplementary_depth_segment_types": (1,), "model_segment_types": (0,),
            }
            compact = writer_module.DirectCompactChrWriter(
                root, final_target=root, sample="TEST", contig="chr1", hap_window_counts={"HAP1": 1},
                window_size=8, source_alignment=Path(raw) / "source.cram", reference_fasta=None,
                sample_metadata={"read_type": "long_read"},
                attrs={"read_coordinate_mode": "personal_realigned", "config": policy},
            )
            records = []
            for index, segment_type in enumerate((0, 1)):
                records.append({
                    "read_key": (f"r{index}", 0),
                    "shared_read_key": (f"r{index}", 0, 4, "4M", False, bool(segment_type)),
                    "read_sequence": "ACGT", "read_haplotype": 1, "read_phase_set": 1,
                    "query_start": 0, "query_end": 4, "read_start": 0, "read_end": 4, "read_strand": 1,
                    "reference_start": 0, "reference_end": 4, "cigar_ops_packed": [4 << 4],
                    "segment_type": segment_type,
                })
            compact.append_hap_chunk(
                "HAP1",
                {"windows": [SimpleNamespace(sequence="ACGT" + "P" * 4, raw_start=0, actual_length=4)],
                 "records": [records], "attrs": {}, "counts": {}},
            )
            # One record per validation chunk forces a metadata/CIGAR chunk
            # boundary while retaining the exact semantic checks.
            with patch.object(writer_module, "VALIDATION_RECORD_CHUNK", 1):
                metadata = compact.finalize(validation_level="full")
            self.assertEqual(compact.depth_arrays["HAP1"]["depth"][:4].tolist(), [1, 1, 1, 1])
            self.assertEqual(compact.supplementary_depth_arrays["HAP1"]["depth"][:4].tolist(), [1, 1, 1, 1])
            self.assertEqual(metadata["filter_policy"]["model_segment_types"], [0])
            self.assertEqual(metadata["depth_summary"]["supplementary_evidence"]["depth_sum"], 4)
            self.assertEqual(metadata["validation"]["checks"]["semantic"]["supplementary_depth_windows"], 1)
if __name__ == "__main__":
    unittest.main()
