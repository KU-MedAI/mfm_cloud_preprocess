from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pysam


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from integrate_wgs_rna import personal_realign as MODULE
SAMTOOLS = "samtools"


def make_read(name: str, flag: int, hp: int | None, unmapped: bool = False) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = "A" * 50
    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
    read.flag = flag
    if unmapped:
        read.reference_id = -1
        read.reference_start = -1
    else:
        read.reference_id = 0
        read.reference_start = 100
        read.mapping_quality = 60
        read.cigar = ((0, 50),)
    if hp is not None:
        read.set_tag("HP", hp, value_type="i")
    return read


class HapReadModeTest(unittest.TestCase):
    def test_mode_resolver(self) -> None:
        config = SimpleNamespace(
            hap_read_mode="hap_specific",
            read_routing_policy="all_to_both",
            unknown_read_policy="include_both",
            include_unmapped=True,
        )
        self.assertEqual(MODULE.resolve_hap_read_mode(config), "hap_specific")
        self.assertEqual(
            (config.read_routing_policy, config.unknown_read_policy, config.include_unmapped),
            ("hp_routed", "drop_unknown", False),
        )

    def test_fragment_classifier(self) -> None:
        self.assertEqual(
            MODULE.classify_illumina_fragment(
                [make_read("a", 99, 1), make_read("a", 147, None)]
            ),
            ("hap1", 1),
        )
        self.assertEqual(
            MODULE.classify_illumina_fragment(
                [make_read("b", 99, 1), make_read("b", 147, 2)]
            ),
            ("conflict", 0),
        )
        self.assertEqual(
            MODULE.classify_illumina_fragment(
                [make_read("c", 99, None), make_read("c", 147, None)]
            ),
            ("unknown", 0),
        )

    def test_strict_illumina_fragment_routing_is_disjoint(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            source = root / "source.bam"
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
            with pysam.AlignmentFile(source, "wb", header=header) as out:
                out.write(make_read("h1", 99, 1))
                out.write(make_read("h1", 147, 1))
                out.write(make_read("h2", 99, 2))
                out.write(make_read("h2", 147, 2))
                out.write(make_read("h1_one", 99, 1))
                out.write(make_read("h1_one", 133, None, unmapped=True))
                out.write(make_read("unknown", 99, None))
                out.write(make_read("unknown", 147, None))
                out.write(make_read("conflict", 99, 1))
                out.write(make_read("conflict", 147, 2))
                out.write(make_read("unmapped_unknown", 77, None, unmapped=True))
                out.write(make_read("unmapped_unknown", 141, None, unmapped=True))

            config = MODULE.PersonalRealignConfig(
                cohort="test",
                output_root=root,
                sample_list=None,
                samples=["TEST"],
                chroms=["chr1"],
                sequence_variant_policy="snv_indel",
                dataset_role="pre_training",
                platform="illumina",
                samtools_bin=SAMTOOLS,
                sort_threads=1,
                sort_mem="64M",
            )
            prepared, qc = MODULE.prepare_illumina_hap_specific_inputs(
                config, "TEST", source, ["chr1"], root
            )
            with pysam.AlignmentFile(prepared[1].path, "rb") as bam:
                h1_names = {read.query_name for read in bam.fetch(until_eof=True)}
            with pysam.AlignmentFile(prepared[2].path, "rb") as bam:
                h2_names = {read.query_name for read in bam.fetch(until_eof=True)}

            self.assertEqual(h1_names, {"h1", "h1_one"})
            self.assertEqual(h2_names, {"h2"})
            self.assertFalse(h1_names & h2_names)
            self.assertEqual(qc["qname_intersection_hap1_hap2"], 0)
            self.assertEqual(qc["source_fragments"], 6)
            self.assertEqual(qc["hap1_fragments"], 2)
            self.assertEqual(qc["hap2_fragments"], 1)
            self.assertEqual(qc["unknown_fragments"], 2)
            self.assertEqual(qc["conflict_fragments"], 1)
            self.assertEqual(qc["one_mate_hp_evidence_fragments"], 1)
            self.assertEqual(qc["unmapped_without_hp_evidence_fragments"], 1)
            self.assertEqual(prepared[1].stats["rescue_unmapped_primary"], 1)
            self.assertTrue(prepared[1].name_collated)
            self.assertFalse((root / "TEST.illumina.hap_specific.name_sorted.bam").exists())

    def test_hap_specific_dry_run_writes_mode_and_qc(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            source = root / "source.bam"
            header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
            with pysam.AlignmentFile(source, "wb", header=header) as out:
                out.write(make_read("h1", 99, 1))
                out.write(make_read("h1", 147, 1))
                out.write(make_read("h2", 99, 2))
                out.write(make_read("h2", 147, 2))
                out.write(make_read("h1_one", 99, 1))
                out.write(make_read("h1_one", 133, None, unmapped=True))
            pysam.index(str(source))
            hap_dirs = {}
            for hap in (1, 2):
                hap_dir = root / f"HAP{hap}"
                hap_dir.mkdir()
                fasta = hap_dir / "chr1.fa"
                fasta.write_text(">chr1\n" + "A" * 1000 + "\n")
                MODULE.subprocess.run([SAMTOOLS, "faidx", str(fasta)], check=True)
                (hap_dir / "chr1.chain").write_text(
                    "chain 0 chr1 1000 + 0 1000 chr1 1000 + 0 1000 1\n1000\n"
                )
                hap_dirs[hap] = hap_dir
            reference = root / "reference.fa"
            reference.write_text(">chr1\n" + "A" * 1000 + "\n")
            MODULE.subprocess.run([SAMTOOLS, "faidx", str(reference)], check=True)

            config = MODULE.PersonalRealignConfig(
                cohort="test",
                output_root=root,
                sample_list=None,
                samples=["TEST"],
                chroms=["chr1"],
                sequence_variant_policy="snv_indel",
                dataset_role="pre_training",
                alignment=source,
                hap1_fasta_dir=hap_dirs[1],
                hap2_fasta_dir=hap_dirs[2],
                reference_fasta=reference,
                platform="illumina",
                hap_read_mode="hap_specific",
                samtools_bin=SAMTOOLS,
                minimap2_bin="minimap2",
                sort_threads=1,
                sort_mem="64M",
                dry_run=True,
            )
            result = MODULE.run_personal_realign_generation(config)
            self.assertEqual(result["ok_count"], 1)
            qc = root / "TEST" / "qc" / "TEST.hap_read_routing.json"
            self.assertEqual(qc.is_file(), True)
            qc_payload = MODULE.json.loads(qc.read_text())
            self.assertEqual(qc_payload["hap_read_mode"], "hap_specific")
            self.assertEqual(qc_payload["input_record_counts"]["unmapped"], 1)
            self.assertEqual(
                qc_payload["illumina_fragment_routing"]["per_haplotype"]["HAP1"]["rescue_unmapped_primary"],
                1,
            )
            h1_json = root / "TEST" / "realigned" / "HAP1" / "TEST.realign.json"
            dry_info = MODULE.json.loads(h1_json.read_text())
            self.assertEqual(dry_info["status"], "dry_run")
            self.assertEqual(dry_info["hap_read_mode"], "hap_specific")

    def test_strict_prepared_source_skips_second_collate(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            config = MODULE.PersonalRealignConfig(
                cohort="test",
                output_root=root,
                sample_list=None,
                samples=["TEST"],
                chroms=["chr1"],
                sequence_variant_policy="snv_indel",
                dataset_role="pre_training",
                platform="illumina",
                samtools_bin=SAMTOOLS,
            )
            prepared = MODULE.PreparedPairedInput(
                path=root / "strict_source.bam",
                route_label="hap1_specific",
                stats={},
                name_collated=True,
            )
            script = MODULE.write_illumina_paired_realign_script(
                config, "TEST", 1, prepared, root / "target.fa", root / "out.bam", root, root,
                {"index_prefix": str(root / "index")},
            )
            content = script.read_text()
            self.assertIn("-T HP,PS,XU,OD", content)
            self.assertIn("-0 /dev/null -s /dev/null", content)
            self.assertIn('"$BWA_MEM2" mem -p -Y -C', content)
            self.assertNotIn('"$SAMTOOLS" collate -u -O', content)


    def test_personal_bam_output_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            config = MODULE.PersonalRealignConfig(
                cohort="test",
                output_root=Path(raw),
                sample_list=None,
                samples=["TEST"],
                chroms=["chr1"],
                sequence_variant_policy="snv_indel",
                dataset_role="pre_training",
                platform="illumina",
                output_format="bam",
            )
            with self.assertRaisesRegex(ValueError, "must be cram"):
                MODULE.validate_config(config)

if __name__ == "__main__":
    unittest.main()
