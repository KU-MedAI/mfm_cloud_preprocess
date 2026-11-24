"""실제 WhatsHap 직렬/병렬 결과 동등성과 실패 시 원본 보존을 검증한다."""
import gzip
import json
import os
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

from integrate_wgs_rna.haplotag import HaplotagConfig, run_whatshap
from integrate_wgs_rna.haplotag_parallel import ChromosomeBudget


class ParallelHaplotagTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.ref = self.root / "ref.fa"
        names = ["chr2", "chr1", "chrX", "chrMissing"]
        self.ref.write_text("".join(f">{c}\n" + "A" * 100 + "\n" for c in names))
        pysam.faidx(str(self.ref))
        self.bam = self.root / "input.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": c, "LN": 100} for c in names],
                  "RG": [{"ID": "rg", "SM": "ALIAS"}]}
        with pysam.AlignmentFile(str(self.bam), "wb", header=header) as out:
            for tid in range(4):
                for flag, start in [(0, 0), (256, 1), (2048, 20)]:
                    read = pysam.AlignedSegment(out.header)
                    read.query_name = f"read{tid}"
                    read.query_sequence = "AAAAAGAAAA"
                    read.query_qualities = pysam.qualitystring_to_array("IIIIIIIIII")
                    read.reference_id = tid
                    read.reference_start = start
                    read.flag = flag
                    read.mapping_quality = 60
                    read.cigarstring = "10M"
                    read.set_tag("RG", "rg")
                    read.set_tag("HP", 7)
                    read.set_tag("PS", 9)
                    read.set_tag("PC", 3)
                    out.write(read)
            read = pysam.AlignedSegment(out.header)
            read.query_name = "unmapped"
            read.query_sequence = "AAAAA"
            read.query_qualities = pysam.qualitystring_to_array("IIIII")
            read.reference_id = read.reference_start = -1
            read.flag = 4
            read.set_tag("RG", "rg")
            read.set_tag("HP", 7)
            out.write(read)
        pysam.index(str(self.bam))
        self.vcf = self.root / "tag.vcf.gz"
        header = pysam.VariantHeader()
        for name in names[:-1]:
            header.contigs.add(name, length=100)
        header.formats.add("GT", 1, "String", "Genotype")
        header.formats.add("PS", 1, "Integer", "Phase set")
        header.add_sample("S")
        with pysam.VariantFile(str(self.vcf), "wz", header=header) as out:
            for chrom in names[:2]:
                rec = out.new_record(contig=chrom, start=5, stop=6, alleles=("A", "G"))
                rec.samples["S"]["GT"] = (0, 1)
                rec.samples["S"].phased = True
                rec.samples["S"]["PS"] = 1
                out.write(rec)
        pysam.tabix_index(str(self.vcf), preset="vcf", force=True)
        self.config = HaplotagConfig(cohort="test", raw_root=self.root,
            tag_vcf_dir=self.root, output_root=self.root, reference_fasta=self.ref,
            platform="illumina", source_sample_id="ALIAS", tag_supplementary=True,
            output_threads=1)

    def tearDown(self):
        self.tmp.cleanup()

    def compare(self, source):
        outputs = []
        for cpus, label in [(1, "serial"), (3, "parallel")]:
            bam = self.root / f"{label}.bam"
            tags = self.root / f"{label}.tags.tsv.gz"
            with patch.dict(os.environ, MFM_HAPLOTAG_CPU_BUDGET=str(cpus)):
                run_whatshap(self.config, "S", source, self.vcf, self.ref, bam, tags)
            pysam.index(str(bam))
            with pysam.AlignmentFile(str(bam)) as handle:
                records = [read.to_string() for read in handle]
                rg = handle.header.to_dict()["RG"]
            with gzip.open(tags, "rb") as handle:
                tag_bytes = handle.read()
            outputs.append((records, rg, tag_bytes))
        self.assertEqual(outputs[0], outputs[1])
        self.assertEqual(10, len(outputs[1][0]))
        self.assertIn("unmapped", outputs[1][0][-1])
        reports = list((self.root / "parallel_haplotag").glob("*/execution.json"))
        report = json.loads(reports[-1].read_text())
        self.assertEqual(["chr2", "chr1", "chrX"], report["contigs"])
        self.assertEqual(3, report["jobs"])
        self.assertEqual("success", report["status"])

    def test_real_bam_records_tags_order_and_unmapped_equal_serial(self):
        self.compare(self.bam)

    def test_real_cram_records_tags_order_and_unmapped_equal_serial(self):
        cram = self.root / "input.cram"
        with pysam.AlignmentFile(str(self.bam)) as src:
            with pysam.AlignmentFile(str(cram), "wc", header=src.header,
                                      reference_filename=str(self.ref)) as dst:
                for read in src:
                    dst.write(read)
        pysam.index(str(cram))
        self.compare(cram)

    def test_failure_preserves_partial_outputs_and_does_not_publish(self):
        output = self.root / "output.bam"
        part = Path(str(output) + ".part")
        part.write_bytes(b"previous incomplete BAM")
        tags = self.root / "tags.tsv.gz"
        tags.write_bytes(b"previous incomplete tags")
        self.config.whatshap_bin = "/bin/false"
        with patch.dict(os.environ, MFM_HAPLOTAG_CPU_BUDGET="3"):
            with self.assertRaises(Exception):
                run_whatshap(self.config, "S", self.bam, self.vcf, self.ref, output, tags)
        self.assertFalse(output.exists())
        self.assertEqual(b"previous incomplete BAM", part.read_bytes())
        self.assertEqual(b"previous incomplete tags", tags.read_bytes())

    def test_free_cpu_growth_respects_total_and_tail_returns_cpus(self):
        pool = self.root / "pool"
        (pool / "leases").mkdir(parents=True)
        (pool / "waiters").mkdir()
        ours = pool / "leases/ours.json"
        other = pool / "leases/other.json"
        ours.write_text(json.dumps({"cpus": 2, "total_cpus": 5, "manager_pid": os.getpid()}))
        other.write_text(json.dumps({"cpus": 3, "total_cpus": 5, "manager_pid": os.getpid()}))
        with patch.dict(os.environ, MFM_CPU_POOL_LEASE_FILE=str(ours)):
            budget = ChromosomeBudget(2, 4)
        budget.acquire()
        budget.acquire()
        acquired = threading.Event()
        thread = threading.Thread(target=lambda: (budget.acquire(), acquired.set()))
        thread.start()
        self.assertFalse(acquired.wait(0.2))
        # 다른 단계 종료를 모사한다.
        other.unlink()
        self.assertTrue(acquired.wait(3))
        thread.join()
        self.assertEqual(4, json.loads(ours.read_text())["cpus"])
        self.assertLessEqual(json.loads(ours.read_text())["cpus"], 5)
        budget.release()
        self.assertEqual(3, json.loads(ours.read_text())["cpus"])
        budget.release()
        self.assertEqual(2, json.loads(ours.read_text())["cpus"])
        budget.release()
        self.assertEqual(1, json.loads(ours.read_text())["cpus"])

    def test_success_preserves_previous_part_and_taglist_copy(self):
        output = self.root / "output.bam"
        part = Path(str(output) + ".part")
        part.write_bytes(b"previous partial BAM")
        tags = self.root / "tags.tsv.gz"
        tags.write_bytes(b"previous partial tags")
        with patch.dict(os.environ, MFM_HAPLOTAG_CPU_BUDGET="3"):
            run_whatshap(self.config, "S", self.bam, self.vcf, self.ref, output, tags)
        self.assertEqual(b"previous partial BAM", part.read_bytes())
        backup = list((self.root / "parallel_haplotag").glob("*/previous.tags.tsv.gz"))
        self.assertEqual(b"previous partial tags", backup[0].read_bytes())
        pysam.quickcheck(str(output))


    def test_native_tail_reacquires_index_and_preserves_records(self):
        from dataclasses import replace
        from integrate_wgs_rna.haplotag import process_one_sample
        baseline=self.root/'baseline.bam';baseline_tags=self.root/'baseline.tags.tsv.gz'
        config=replace(self.config,jobs=1,alignment=self.bam,allow_longread_vcf_haplotag=True,output_threads=4)
        with patch.dict(os.environ,MFM_HAPLOTAG_CPU_BUDGET='1'):
            run_whatshap(config,'S',self.bam,self.vcf,self.ref,baseline,baseline_tags)
        for suffix in ['', '.tbi']:os.link(Path(str(self.vcf)+suffix),self.root/('S.tag.vcf.gz'+suffix))
        pool=self.root/'pool';(pool/'leases').mkdir(parents=True);(pool/'waiters').mkdir()
        lease=pool/'leases/native.json'
        lease.write_text(json.dumps({'branch':'WGS','stage':'06','sample':'S','cpus':3,'total_cpus':8,'memory_bytes':128*1024**3,'manager_pid':os.getpid()}))
        (pool/'policy.json').write_text(json.dumps({'memory_limit_bytes':256*1024**3,'cpu_pool_limit':8,'live_limits':True}))
        env={'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_HAPLOTAG_CPU_BUDGET':'3'}
        with patch.dict(os.environ,env):
            result=process_one_sample(config,'S',self.root,self.root,self.root/'native',self.ref)
        self.assertNotEqual(result.status,'FAIL',result.message)
        with pysam.AlignmentFile(str(baseline)) as handle:expected=[r.to_string() for r in handle]
        with pysam.AlignmentFile(result.haplotagged_bam) as handle:actual=[r.to_string() for r in handle]
        self.assertEqual(expected,actual)
        pysam.quickcheck(result.haplotagged_bam)
        self.assertTrue(Path(result.haplotagged_bam+'.bai').exists())
        row=json.loads(lease.read_text());events={v['phase']:v for v in row['haplotag_resource_events']}
        self.assertEqual(events['combine_and_tags']['memory_bytes'],4*1024**3)
        self.assertEqual(events['output_index']['cpus'],5)
        self.assertEqual(events['inspect_tags']['cpus'],2)
        self.assertEqual(row['cpus'],1);self.assertEqual(row['memory_bytes'],4*1024**3)
        self.assertEqual(row['haplotag_initial_memory_bytes'],128*1024**3)



if __name__ == "__main__":
    unittest.main()
