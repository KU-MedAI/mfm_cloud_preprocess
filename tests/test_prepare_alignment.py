from __future__ import annotations

import gzip
import importlib.util
import json
import random
import subprocess
import tempfile
import unittest
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "prepare_alignment", ROOT / "src/integrate_wgs_rna/prepare_alignment.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)
SAMTOOLS = "samtools"
MINIMAP2 = "minimap2"
BWA_MEM2 = "bwa-mem2"


class PrepareAlignmentTest(unittest.TestCase):
    def test_platform_presets_distinguish_hifi_and_clr(self) -> None:
        self.assertEqual(MODULE.platform_alignment_parameters("hifi"), ("PACBIO", "map-hifi"))
        self.assertEqual(MODULE.platform_alignment_parameters("clr"), ("PACBIO", "map-pb"))
        self.assertEqual(MODULE.platform_alignment_parameters("ont"), ("ONT", "map-ont"))

    def test_bwa_mem2_command_keeps_supplementary_semantics(self) -> None:
        command = MODULE.bwa_mem2_command(
            sample="ILM1", run_id="runA", read_group_id="rgA", reference=Path("ref.fa"),
            bwa_mem2=BWA_MEM2, align_threads=12,
        )
        self.assertEqual(command[:5], [BWA_MEM2, "mem", "-t", "12", "-Y"])
        self.assertNotIn("-M", command)
        self.assertIn("PL:ILLUMINA", command[command.index("-R") + 1])
        policy = MODULE.alignment_policy(
            states=["unaligned"], platform="illumina", illumina_aligner="minimap2")
        self.assertEqual(policy["aligner"], "minimap2")
        self.assertEqual(policy["preset"], "sr")

    def test_duplicate_plan_is_platform_and_provenance_aware(self) -> None:
        trusted = [{"trusted_picard_marking": True}]
        untrusted = [{"trusted_picard_marking": False}]
        preserved = MODULE.resolve_duplicate_plan(
            platform="illumina", states=["aligned"], requested_policy="auto", provenance=trusted)
        self.assertEqual(preserved["action"], "preserved_source_picard")
        marked = MODULE.resolve_duplicate_plan(
            platform="illumina", states=["aligned"], requested_policy="auto", provenance=untrusted)
        self.assertEqual(marked["action"], "picard_marked")
        long_read = MODULE.resolve_duplicate_plan(
            platform="hifi", states=["unaligned"], requested_policy="auto", provenance=[])
        self.assertEqual(long_read["action"], "not_applicable")
        with self.assertRaises(MODULE.PrepareAlignmentError):
            MODULE.resolve_duplicate_plan(
                platform="illumina", states=["unaligned"], requested_policy="preserve", provenance=[])

    def test_picard_marks_duplicate_pair_without_removing_records(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            source = root / "source.bam"
            marked = root / "marked.bam"
            metrics = root / "metrics.txt"
            header = {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "chr1", "LN": 2000}],
                "RG": [{"ID": "rg1", "SM": "ILM", "PL": "ILLUMINA", "LB": "lib1"}],
            }
            records = []
            for pair_index, qname in enumerate(("A:1:FC:1:1101:100:100", "A:1:FC:1:1101:300:300")):
                for is_read2, start, mate_start, flag, tlen in (
                    (False, 100, 300, 99, 250), (True, 300, 100, 147, -250),
                ):
                    read = pysam.AlignedSegment()
                    read.query_name = qname
                    read.query_sequence = "A" * 50
                    quality = "I" if pair_index == 0 else "H"
                    read.query_qualities = pysam.qualitystring_to_array(quality * 50)
                    read.flag = flag
                    read.reference_id = 0
                    read.reference_start = start
                    read.mapping_quality = 60
                    read.cigar = ((0, 50),)
                    read.next_reference_id = 0
                    read.next_reference_start = mate_start
                    read.template_length = tlen
                    read.set_tag("RG", "rg1")
                    records.append(read)
            records.sort(key=lambda read: (read.reference_start, read.query_name))
            with pysam.AlignmentFile(source, "wb", header=header) as out:
                for read in records:
                    out.write(read)
            with (root / "picard.log").open("wb") as log:
                MODULE.mark_illumina_duplicates(
                    source, marked, metrics, picard="picard", heap_gb=1, temp_dir=root, log=log)
            stats = MODULE.tag_record_stats(marked)
            self.assertEqual(stats["records"], 4)
            self.assertEqual(stats["duplicate_records"], 2)
            self.assertTrue(metrics.is_file())

    def test_clr_fastq_uses_map_pb_and_becomes_analysis_bam(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            rng = random.Random(31)
            sequence = "".join(rng.choice("ACGT") for _ in range(5000))
            reference = root / "ref.fa"
            reference.write_text(">chr1\n" + sequence + "\n")
            subprocess.run([SAMTOOLS, "faidx", str(reference)], check=True)
            fastq = root / "clr.fastq.gz"
            read = sequence[1000:3000]
            quality = "".join("HI"[i % 2] for i in range(len(read)))
            with gzip.open(fastq, "wt") as out:
                out.write(f"@clr_read\n{read}\n+\n{quality}\n")
            contract = {
                "sample_id": "CLR1", "platform": "clr", "input_mode": "read_files_manifest",
                "source": {"files": [{"sample_id": "CLR1", "run_id": "run1",
                    "file_role": "long_reads", "format": "fastq_gz", "alignment_state": "unaligned",
                    "path": str(fastq), "read_group_id": "clr_rg"}],
                    "file_count": 1, "alignment_states": ["unaligned"]},
                "references": {"analysis": str(reference)}, "prepared_alignment": None,
            }
            output = root / "CLR1.analysis.bam"
            result = MODULE.prepare(
                contract, output=output, qc_json=root / "qc.json", minimap2=MINIMAP2,
                samtools=SAMTOOLS, align_threads=2, sort_threads=1, sort_memory="64M",
                extract_threads=1, overwrite=False, keep_run_bams=False,
            )
            self.assertEqual(result["preset"], "map-pb")
            self.assertEqual(result["metrics"]["mapped_reads"], 1)
            with pysam.AlignmentFile(output, "rb") as alignment:
                self.assertEqual(alignment.header.to_dict()["RG"][0]["PL"], "PACBIO")

    def test_ubam_and_fastq_runs_become_one_analysis_bam(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            rng = random.Random(17)
            reference_sequence = "".join(rng.choice("ACGT") for _ in range(6000))
            reference = root / "ref.fa"
            reference.write_text(">chr1\n" + reference_sequence + "\n")
            subprocess.run([SAMTOOLS, "faidx", str(reference)], check=True)
            sequences = [reference_sequence[500:1500], reference_sequence[3000:4000]]

            ubam = root / "run1.ubam"
            with pysam.AlignmentFile(ubam, "wb", header={"HD": {"VN": "1.6"}}) as out:
                read = pysam.AlignedSegment()
                read.query_name = "ubam_read"
                read.query_sequence = sequences[0]
                read.query_qualities = pysam.qualitystring_to_array("I" * len(sequences[0]))
                read.flag = 4
                out.write(read)
            fastq = root / "run2.fastq.gz"
            with gzip.open(fastq, "wt") as out:
                out.write(f"@fastq_read\n{sequences[1]}\n+\n{'H' * len(sequences[1])}\n")

            files = [
                {"sample_id": "S1", "run_id": "run1", "file_role": "alignment",
                 "format": "ubam", "alignment_state": "unaligned", "path": str(ubam),
                 "read_group_id": "rg1"},
                {"sample_id": "S1", "run_id": "run2", "file_role": "long_reads",
                 "format": "fastq_gz", "alignment_state": "unaligned", "path": str(fastq),
                 "read_group_id": "rg2"},
            ]
            contract = {
                "sample_id": "S1", "platform": "hifi", "input_mode": "read_files_manifest",
                "source": {"files": files, "file_count": 2,
                           "alignment_states": ["unaligned"]},
                "references": {"analysis": str(reference)}, "prepared_alignment": None,
            }
            output = root / "S1.analysis.bam"
            qc = root / "S1.stage04.json"
            try:
                result = MODULE.prepare(
                    contract, output=output, qc_json=qc, minimap2=MINIMAP2,
                    samtools=SAMTOOLS, align_threads=2, sort_threads=1, sort_memory="64M",
                    extract_threads=1, overwrite=False, keep_run_bams=False,
                )
            except Exception:
                log = Path(str(output) + ".stage04.log")
                print(log.read_text() if log.is_file() else "Stage 04 log 없음")
                raise
            self.assertEqual(result["status"], "created_raw_longread")
            self.assertEqual(result["metrics"]["mapped_reads"], 2)
            self.assertEqual(result["metrics"]["base_quality"]["distinct_phred_values"], 2)
            self.assertTrue(Path(str(output) + ".bai").is_file())
            with pysam.AlignmentFile(output, "rb") as alignment:
                header = alignment.header.to_dict()
            self.assertEqual({rg["ID"] for rg in header["RG"]}, {"rg1", "rg2"})
            self.assertEqual({rg["SM"] for rg in header["RG"]}, {"S1"})
            second = MODULE.prepare(
                contract, output=output, qc_json=qc, minimap2=MINIMAP2,
                samtools=SAMTOOLS, align_threads=2, sort_threads=1, sort_memory="64M",
                extract_threads=1, overwrite=False, keep_run_bams=False,
            )
            self.assertEqual(second["status"], "reused_canonical_bam")
            self.assertEqual(json.loads(qc.read_text())["normalization_reference"]["path"],
                             str(reference.resolve()))


    def test_illumina_paired_fastq_becomes_analysis_bam(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            rng = random.Random(79)
            sequence = "".join(rng.choice("ACGT") for _ in range(7000))
            reference = root / "ref.fa"
            reference.write_text(">chr1\n" + sequence + "\n")
            subprocess.run([SAMTOOLS, "faidx", str(reference)], check=True)
            r1, r2 = root / "R1.fastq.gz", root / "R2.fastq.gz"
            with gzip.open(r1, "wt") as out:
                out.write(f"@pair1\n{sequence[1000:1150]}\n+\n{'I' * 150}\n")
            with gzip.open(r2, "wt") as out:
                out.write(f"@pair1\n{sequence[1200:1350]}\n+\n{'H' * 150}\n")
            files = [
                {"sample_id": "ILM1", "run_id": "runA", "file_role": "read1",
                 "format": "fastq_gz", "alignment_state": "unaligned", "path": str(r1),
                 "read_group_id": "rgA"},
                {"sample_id": "ILM1", "run_id": "runA", "file_role": "read2",
                 "format": "fastq_gz", "alignment_state": "unaligned", "path": str(r2),
                 "read_group_id": "rgA"},
            ]
            contract = {
                "sample_id": "ILM1", "platform": "illumina", "input_mode": "read_files_manifest",
                "source": {"files": files, "file_count": 2, "alignment_states": ["unaligned"]},
                "references": {"analysis": str(reference), "source_decode": str(reference)},
                "prepared_alignment": None,
            }
            output = root / "ILM1.analysis.bam"
            result = MODULE.prepare(
                contract, output=output, qc_json=root / "qc.json", minimap2=MINIMAP2,
                bwa_mem2=BWA_MEM2,
                samtools=SAMTOOLS, align_threads=2, sort_threads=1, sort_memory="64M",
                extract_threads=1, overwrite=False, keep_run_bams=False,
            )
            self.assertEqual(result["status"], "created_raw_illumina")
            self.assertEqual(result["preset"], "mem")
            self.assertEqual(result["alignment_policy"]["aligner"], "bwa-mem2")
            self.assertIn("bwa_mem2", result["tools"])
            self.assertRegex(result["tools"]["bwa_mem2"], r"^bwa-mem2 \d+\.\d+")
            self.assertEqual(result["bwa_mem2_index"]["status"], "created")
            self.assertTrue(all(Path(path).is_file() for path in result["bwa_mem2_index"]["files"]))
            with (root / "cache_reuse.log").open("ab") as log:
                reused = MODULE.ensure_bwa_mem2_index(
                    reference, bwa_mem2=BWA_MEM2,
                    index_root=root / "_reference_cache" / "bwa_mem2", log=log)
            self.assertEqual(reused["status"], "reused")
            self.assertEqual(result["output_record_tag_counts"]["records"], 2)
            self.assertEqual(result["duplicate_marking"]["action"], "picard_marked")
            self.assertFalse(result["duplicate_marking"]["remove_duplicates"])
            self.assertTrue(Path(result["duplicate_marking"]["metrics_path"]).is_file())
            with pysam.AlignmentFile(output, "rb") as alignment:
                programs = alignment.header.to_dict().get("PG", [])
            self.assertTrue(any("markduplicates" in " ".join(map(str, program.values())).lower()
                                for program in programs))
            second = MODULE.prepare(
                contract, output=output, qc_json=root / "qc.json", minimap2=MINIMAP2,
                bwa_mem2=BWA_MEM2,
                samtools=SAMTOOLS, align_threads=2, sort_threads=1, sort_memory="64M",
                extract_threads=1, overwrite=False, keep_run_bams=False,
            )
            self.assertEqual(second["status"], "reused_canonical_bam")
            self.assertEqual(second["tools"], result["tools"])
            self.assertEqual(second["bwa_mem2_index"], result["bwa_mem2_index"])
            with pysam.AlignmentFile(output, "rb") as alignment:
                self.assertEqual(alignment.header.to_dict()["RG"][0]["PL"], "ILLUMINA")

    def test_aligned_cram_is_transcoded_to_bam_with_hp_ps_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            rng = random.Random(83)
            sequence = "".join(rng.choice("ACGT") for _ in range(5000))
            reference = root / "ref.fa"
            reference.write_text(">chr1\n" + sequence + "\n")
            subprocess.run([SAMTOOLS, "faidx", str(reference)], check=True)
            source_bam = root / "source.bam"
            header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                      "SQ": [{"SN": "chr1", "LN": len(sequence)}],
                      "RG": [{"ID": "rg1", "SM": "S1", "PL": "ONT"}]}
            with pysam.AlignmentFile(source_bam, "wb", header=header) as out:
                read = pysam.AlignedSegment()
                read.query_name = "tagged_read"
                read.query_sequence = sequence[1000:1150]
                read.query_qualities = pysam.qualitystring_to_array("IH" * 75)
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 1000
                read.mapping_quality = 60
                read.cigar = ((0, 150),)
                read.set_tag("RG", "rg1")
                read.set_tag("HP", 1, value_type="i")
                read.set_tag("PS", 100, value_type="i")
                out.write(read)
            subprocess.run([SAMTOOLS, "index", str(source_bam)], check=True)
            source_cram = root / "source.cram"
            subprocess.run([SAMTOOLS, "view", "-C", "-T", str(reference), "-o", str(source_cram), str(source_bam)],
                           check=True)
            subprocess.run([SAMTOOLS, "index", str(source_cram)], check=True)
            contract = {
                "sample_id": "S1", "platform": "ont", "input_mode": "legacy_alignment",
                "source": {"alignment": str(source_cram), "files": [], "file_count": 1,
                           "alignment_states": ["aligned"]},
                "references": {"analysis": str(reference), "source_decode": str(reference)},
                "prepared_alignment": None,
            }
            output = root / "S1.analysis.bam"
            result = MODULE.prepare(
                contract, output=output, qc_json=root / "qc.json", minimap2=MINIMAP2,
                samtools=SAMTOOLS, align_threads=2, sort_threads=1, sort_memory="64M",
                extract_threads=1, overwrite=False, keep_run_bams=False,
            )
            self.assertEqual(result["status"], "converted_aligned")
            self.assertEqual(result["source_record_tag_counts"][0],
                             {"records": 1, "mapped_records": 1, "primary_records": 1,
                              "secondary_records": 0, "supplementary_records": 0,
                              "duplicate_records": 0, "primary_duplicate_records": 0,
                              "hp_tagged_records": 1, "ps_tagged_records": 1})
            self.assertEqual(result["output_record_tag_counts"], result["source_record_tag_counts"][0])
            with pysam.AlignmentFile(output, "rb") as alignment:
                read = next(alignment.fetch(until_eof=True))
                self.assertEqual(read.get_tag("HP"), 1)
                self.assertEqual(read.get_tag("PS"), 100)

            bam_contract = {
                **contract,
                "source": {"alignment": str(source_bam), "files": [], "file_count": 1,
                           "alignment_states": ["aligned"]},
            }
            bam_output = root / "S1.from_bam.analysis.bam"
            bam_result = MODULE.prepare(
                bam_contract, output=bam_output, qc_json=root / "bam_qc.json", minimap2=MINIMAP2,
                samtools=SAMTOOLS, align_threads=2, sort_threads=1, sort_memory="64M",
                extract_threads=1, overwrite=False, keep_run_bams=False,
            )
            self.assertEqual(bam_result["status"], "converted_aligned")
            self.assertEqual(bam_result["output_record_tag_counts"],
                             {"records": 1, "mapped_records": 1, "primary_records": 1,
                              "secondary_records": 0, "supplementary_records": 0,
                              "duplicate_records": 0, "primary_duplicate_records": 0,
                              "hp_tagged_records": 1, "ps_tagged_records": 1})

if __name__ == "__main__":
    unittest.main()
