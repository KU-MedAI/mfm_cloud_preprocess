#!/usr/bin/env python3
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pysam

from integrate_wgs_rna.sample_vcf import PerSampleVcfConfig, validate_tag_vcf


class SampleTagVcfPhaseTest(unittest.TestCase):
    def make_vcf(self, root: Path, genotypes: list[str]) -> Path:
        plain = root / "tag.vcf"
        with plain.open("w") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("##contig=<ID=chr1,length=1000>\n")
            handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
            for index, gt in enumerate(genotypes, 1):
                handle.write(f"chr1\t{index}\t.\tA\tC\t.\tPASS\t.\tGT\t{gt}\n")
        compressed = root / "tag.vcf.gz"
        pysam.tabix_compress(str(plain), str(compressed), force=True)
        pysam.tabix_index(str(compressed), preset="vcf", force=True)
        return compressed

    def config(self, root: Path, *, fast: bool = False) -> PerSampleVcfConfig:
        return PerSampleVcfConfig(
            cohort="test", common_vcf_dir=root, output_dir=root,
            chroms=["chr1"], sequence_variant_policy="snv_indel",
            bcftools_bin="bcftools", fast_validate=fast,
        )

    def test_complete_phased_heterozygotes_pass_even_in_fast_mode(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            result = validate_tag_vcf(self.config(root, fast=True), self.make_vcf(root, ["0|1", "1|0"]))
            self.assertEqual(result["complete_phased_het_count"], 2)
            self.assertEqual(result["unphased_het_count"], 0)
            self.assertEqual(result["missing_gt_count"], 0)

    def test_unphased_heterozygote_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            with self.assertRaisesRegex(RuntimeError, "complete phased"):
                validate_tag_vcf(self.config(root), self.make_vcf(root, ["0|1", "0/1"]))

    def test_partial_missing_genotype_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            with self.assertRaisesRegex(RuntimeError, "missing=1"):
                validate_tag_vcf(self.config(root), self.make_vcf(root, ["0|1", "1|."]))


if __name__ == "__main__":
    unittest.main()
