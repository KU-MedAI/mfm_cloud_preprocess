from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pysam

from integrate_wgs_rna.full_preprocess import FullPreprocessConfig, validate_personal_fasta_reference_contract
from integrate_wgs_rna.personal_realign import PersonalRealignConfig, validate_config
from integrate_wgs_rna.provided_assembly import AssemblyDirectConfig, run_assembly_direct
from integrate_wgs_rna.window_index import WindowIndexConfig, run_window_index_generation


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n{sequence}\n")
    pysam.faidx(str(path))


class ProvidedAssemblyTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.reference = self.root / "ref.fa"
        write_fasta(self.reference, {"chr1": "ACGT" * 20, "chr2": "TGCA" * 20})
        self.h1 = self.root / "h1.fa"
        self.h2 = self.root / "h2.fa"
        write_fasta(self.h1, {"S1#1#A1": "ACGT" * 21, "S1#1#A2": "TGCA" * 19, "S1#1#U": "N" * 10})
        write_fasta(self.h2, {"S1#2#B1": "ACGT" * 19, "S1#2#B2": "TGCA" * 21})
        self.a1 = self.root / "h1.alias"
        self.a2 = self.root / "h2.alias"
        self.a1.write_text("# assembly\tucsc\tgenbank\nS1#1#A1\tchr1\tA1\nS1#1#A2\tchr2\tA2\nS1#1#U\tchrUn_U\tU\n")
        self.a2.write_text("# assembly\tucsc\tgenbank\nS1#2#B1\tchr1\tB1\nS1#2#B2\tchr2\tB2\n")

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def config(self, **kwargs) -> AssemblyDirectConfig:
        values = dict(
            cohort="TEST", sample="S1", output_root=self.root / "out",
            reference_fasta=self.reference, hap1_fasta=self.h1, hap1_chrom_alias=self.a1,
            hap2_fasta=self.h2, hap2_chrom_alias=self.a2, chroms=["chr1-chr2"],
        )
        values.update(kwargs)
        return AssemblyDirectConfig(**values)

    def test_extract_and_qc(self) -> None:
        result = run_assembly_direct(self.config())
        self.assertEqual(result["status"], "pass")
        self.assertEqual(result["phase_label_status"], "unverified")
        self.assertEqual(result["hap_summaries"]["HAP1"]["omitted_contigs"], 1)
        out = self.root / "out/S1/individual_fasta/HAP1/chr1.fa"
        with pysam.FastaFile(str(out)) as fasta:
            self.assertEqual(fasta.references, ["chr1"])
            self.assertEqual(fasta.fetch("chr1"), "ACGT" * 21)
        reused = run_assembly_direct(self.config())
        self.assertEqual(reused["run_status"], "SKIP")

    def test_fragmented_chromosome_fails(self) -> None:
        self.a1.write_text("# assembly\tucsc\tgenbank\nS1#1#A1\tchr1\tA1\nS1#1#A2\tchr2_A2_random\tA2\n")
        with self.assertRaisesRegex(ValueError, "exactly one chromosome-level contig"):
            run_assembly_direct(self.config())

    def test_window_blacklist_gate_and_contract(self) -> None:
        run_assembly_direct(self.config())
        with self.assertRaisesRegex(ValueError, "requires disabled GRCh38 blacklist"):
            run_window_index_generation(WindowIndexConfig(
                cohort="TEST", output_root=self.root / "out", samples=["S1"], chroms=["chr1-chr2"],
                sequence_variant_policy="snv_indel", sequence_source="assembly_direct",
                blacklist_bed=self.root / "blacklist.bed",
            ))
        result = run_window_index_generation(WindowIndexConfig(
            cohort="TEST", output_root=self.root / "out", samples=["S1"], chroms=["chr1-chr2"],
            sequence_variant_policy="snv_indel", sequence_source="assembly_direct", blacklist_bed=None,
            window_size=16, offset_unit=1,
        ))
        self.assertEqual(result["fail_count"], 0)
        meta = json.loads((self.root / "out/S1/windows/windows.json").read_text())
        self.assertEqual(meta["blacklist_policy"], "disabled_for_assembly_direct_no_multichain_projection")

    def test_phase_routing_gate(self) -> None:
        config = PersonalRealignConfig(
            cohort="TEST", output_root=self.root, sample_list=None, samples=["S1"],
            chroms=["chr1"], sequence_variant_policy="snv_indel", dataset_role="pre_training",
            sequence_source="assembly_direct", hap_read_mode="hap_specific",
        )
        with self.assertRaisesRegex(ValueError, "assembly_direct defaults to all_to_both"):
            validate_config(config)


if __name__ == "__main__":
    unittest.main()
