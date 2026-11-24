import csv
import json
import tempfile
import unittest
from pathlib import Path

from integrate_wgs_rna.full_preprocess import FullPreprocessConfig, validate_personal_fasta_reference_contract
from integrate_wgs_rna.personal_fasta import parse_overlap_skip_count, reference_fingerprint
from integrate_wgs_rna.personal_realign import (
    PersonalRealignConfig,
    previous_matches_qc_contract,
    verify_realign_counts,
    write_realign_counts,
)


class PilotContractTest(unittest.TestCase):
    def test_relaxed_pilot_is_not_reused_by_production(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            config = PersonalRealignConfig(
                cohort="test", output_root=Path(tmp), sample_list=None, samples=["S1"],
                chroms=["chr22"], sequence_variant_policy="snv_indel",
                dataset_role="pre_training", min_mapped_frac=0.90,
                ratio_qc_policy="strict",
            )
            relaxed = {"qc_contract": {"min_mapped_frac": 0.85}, "routed": {"mapped_frac": 0.881227}}
            strict = {"qc_contract": {"min_mapped_frac": 0.90}, "routed": {"mapped_frac": 0.95}}
            self.assertFalse(previous_matches_qc_contract(relaxed, config))
            self.assertTrue(previous_matches_qc_contract(strict, config))
            self.assertFalse(previous_matches_qc_contract({"routed": {"mapped_frac": 0.95}}, config))

            record_only = PersonalRealignConfig(
                cohort="test", output_root=Path(tmp), sample_list=None, samples=["S1"],
                chroms=["chr22"], sequence_variant_policy="snv_indel",
                dataset_role="pre_training", ratio_qc_policy="record_only",
            )
            observed_low = {
                "qc_contract": {
                    "min_mapped_frac": 0.90,
                    "ratio_qc_policy": "record_only",
                },
                "routed": {"mapped_frac": 0.20},
            }
            self.assertTrue(previous_matches_qc_contract(observed_low, record_only))
            self.assertFalse(previous_matches_qc_contract(observed_low, config))

    def test_realign_metric_keys_are_unique(self) -> None:
        scan = {
            "primary": 10, "secondary": 1, "supplementary": 2,
            "paired_primary": 10, "proper_pair_primary": 8, "paired_primary_frac": 1.0,
            "duplicate_primary": 1, "od_tag_primary": 1,
            "routed_primary": 10, "routed_mapped": 9,
            "routed_hp1": 4, "routed_hp2": 3, "routed_untagged": 3,
            "rescue_primary": 0, "rescue_mapped": 0,
            "routed_mapped_frac": 0.9, "rescue_mapped_frac": 0.0,
            "nm_sum": 2, "mapped_bp": 900, "nm_per_base_all_mapped": 0.002,
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "counts.tsv"
            write_realign_counts(path, scan, 1)
            with path.open() as handle:
                keys = [row[0] for row in list(csv.reader(handle, delimiter="\t"))[1:]]
        self.assertEqual(len(keys), len(set(keys)))
        self.assertEqual(keys.count("routed_hp1"), 1)

    def test_low_mapped_fraction_is_record_only_unless_strict(self) -> None:
        scan = {
            "routed_hp1": 1,
            "routed_mapped_frac": 0.20,
            "paired_primary_frac": 0.20,
        }
        with tempfile.TemporaryDirectory() as tmp:
            common = dict(
                cohort="test", output_root=Path(tmp), sample_list=None, samples=["S1"],
                chroms=["chr22"], sequence_variant_policy="snv_indel",
                dataset_role="pre_training", verify_counts=False,
            )
            record_only = PersonalRealignConfig(**common, ratio_qc_policy="record_only")
            strict = PersonalRealignConfig(**common, ratio_qc_policy="strict")
            verify_realign_counts(record_only, 1, {}, scan, None)
            with self.assertRaisesRegex(RuntimeError, "routed mapped fraction"):
                verify_realign_counts(strict, 1, {}, scan, None)

    def test_overlap_parser_and_reference_fingerprint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "ref.fa"
            fasta.write_text(">chr22\\nACGT\\n")
            fai = Path(str(fasta) + ".fai")
            fai.write_text("chr22\\t4\\t7\\t4\\t5\\n")
            first = reference_fingerprint(fasta)
            self.assertEqual(parse_overlap_skip_count(
                "The site overlaps with another variant, skipping...\\n"
                "overlaps with another variant, skipping...\\n"
            ), 2)
            fai.write_text("chr22\\t5\\t7\\t5\\t6\\n")
            self.assertNotEqual(reference_fingerprint(fasta)["fai_sha256"], first["fai_sha256"])

    def test_stage12_personal_reference_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = root / "ref.fa"
            reference.write_text(">chr22\\nACGT\\n")
            Path(str(reference) + ".fai").write_text("chr22\\t4\\t7\\t4\\t5\\n")
            alignment = root / "source.bam"
            alignment.touch()
            hap1 = root / "individual_fasta" / "HAP1"
            hap2 = root / "individual_fasta" / "HAP2"
            hap1.mkdir(parents=True)
            hap2.mkdir(parents=True)
            qc_path = hap1.parent / "personal_fasta_qc.json"
            rows = [{"contig": "chr22", "chain_qc": {"status": "pass"}}]
            qc_path.write_text(json.dumps({
                "schema": "mfm_personal_fasta_sample_qc_v1", "status": "pass", "sample": "S1",
                "reference_fasta": str(reference.resolve()),
                "reference_fingerprint": reference_fingerprint(reference),
                "hap_results": {"HAP1": rows, "HAP2": rows},
            }))
            config = FullPreprocessConfig(
                sample_id="S1", dataset_role="pre_training", read_type="long_read", project_id="test",
                hap1_fasta_dir=hap1, hap1_chain_dir=hap1, hap2_fasta_dir=hap2, hap2_chain_dir=hap2,
                alignment=alignment, contigs=["chr22"], sequence_reference_fasta=reference,
                personal_fasta_qc=qc_path, sequence_variant_policy="snv_indel",
            )
            result = validate_personal_fasta_reference_contract(config, ["chr22"])
            self.assertEqual(result["status"], "pass")
            payload = json.loads(qc_path.read_text())
            payload["reference_fingerprint"]["fai_sha256"] = "0" * 64
            qc_path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "fingerprint"):
                validate_personal_fasta_reference_contract(config, ["chr22"])


if __name__ == "__main__":
    unittest.main()
