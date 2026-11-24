from __future__ import annotations

import json
import os
import threading
import time
import random
import sys
import tempfile
import unittest
from unittest.mock import patch
from dataclasses import replace
from pathlib import Path

import pysam


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from integrate_wgs_rna import personal_realign as module

SAMTOOLS = "samtools"
BWA_MEM2 = "bwa-mem2"


def make_record(name: str, flag: int, start: int, mate_start: int, sequence: str, duplicate: bool = False):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = sequence
    read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    read.flag = flag | (0x400 if duplicate else 0)
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigar = ((0, len(sequence)),)
    if read.is_paired:
        read.next_reference_id = 0
        read.next_reference_start = mate_start
    read.set_tag("HP", 1, value_type="i")
    read.set_tag("PS", 100, value_type="i")
    return read


class Stage22BwaIntegrationTest(unittest.TestCase):
    def test_pair_singleton_duplicate_and_chain_contract(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            rng = random.Random(17)
            sequence = "".join(rng.choice("ACGT") for _ in range(3000))
            reference = root / "reference.fa"
            reference.write_text(">chr1\n" + sequence + "\n")
            module.subprocess.run([SAMTOOLS, "faidx", str(reference)], check=True)

            source = root / "source.bam"
            header = {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "chr1", "LN": len(sequence)}],
                "RG": [{"ID": "rg1", "SM": "TEST", "PL": "ILLUMINA"}],
            }
            with pysam.AlignmentFile(source, "wb", header=header) as out:
                rows = [
                    make_record("dup", 99, 100, 300, sequence[100:200], True),
                    make_record("dup", 147, 300, 100, sequence[300:400], True),
                    make_record("normal", 99, 500, 700, sequence[500:600]),
                    make_record("normal", 147, 700, 500, sequence[700:800]),
                    make_record("singleton", 73, 900, 900, sequence[900:1000]),
                    make_record("other", 0, 1100, 1100, sequence[1100:1200]),
                ]
                for read in rows:
                    read.set_tag("RG", "rg1")
                    out.write(read)
            pysam.index(str(source))

            hap_dirs = {}
            for hap in (1, 2):
                hap_dir = root / f"HAP{hap}"
                hap_dir.mkdir()
                fasta = hap_dir / "chr1.fa"
                fasta.write_text(">chr1\n" + sequence + "\n")
                module.subprocess.run([SAMTOOLS, "faidx", str(fasta)], check=True)
                (hap_dir / "chr1.chain").write_text(
                    f"chain 0 chr1 {len(sequence)} + 0 {len(sequence)} chr1 {len(sequence)} + 0 {len(sequence)} 1\n{len(sequence)}\n"
                )
                hap_dirs[hap] = hap_dir

            config = module.PersonalRealignConfig(
                cohort="test", output_root=root, sample_list=None, samples=["TEST"],
                chroms=["chr1"], sequence_variant_policy="snv_indel",
                dataset_role="pre_training", alignment=source,
                hap1_fasta_dir=hap_dirs[1], hap2_fasta_dir=hap_dirs[2],
                reference_fasta=reference, platform="illumina", hap_read_mode="all_to_both",
                include_unmapped=False, min_mapped_frac=0.5, threads=2,
                sort_threads=1, sort_mem="64M", samtools_bin=SAMTOOLS,
                bwa_mem2_bin=BWA_MEM2,
                illumina_aligner="bwa-mem2", illumina_duplicate_mode="preserve_source",
                inspect_reads=10, hap_jobs=2, route_jobs=4,
            )
            counter = module.add_alignment_input_metrics
            def force_parallel_metrics(*args, **kwargs):
                return counter(*args, workers=2, chunk_bytes=257, **kwargs)
            with patch.object(module, "add_alignment_input_metrics", side_effect=force_parallel_metrics):
                result = module.run_personal_realign_generation(config)
            self.assertEqual(result["ok_count"], 1)

            routing = json.loads((root / "TEST/qc/TEST.hap_read_routing.json").read_text())
            for hap in ("HAP1", "HAP2"):
                stats = routing["prepared_inputs"][hap]
                self.assertEqual(stats["paired_fragments"], 2)
                self.assertEqual(stats["input_metrics_execution"]["mode"], "parallel_sam_stream")
                self.assertEqual(stats["alignment_input_primary"], 4)
                self.assertEqual(stats["singleton_discard_records"], 1)
                self.assertEqual(stats["read_other_discard_records"], 1)
                self.assertEqual(stats["alignment_input_duplicate"], 2)

            realign_json = root / "TEST/realigned/HAP1/TEST.realign.json"
            qc = json.loads(realign_json.read_text())
            self.assertTrue(qc["aligner"].startswith("bwa-mem2"))
            self.assertEqual(qc["alignment_contract"]["args"], "mem -p -Y -C")
            self.assertEqual(qc["records"]["duplicate_primary"], 2)
            self.assertEqual(qc["records"]["od_tag_primary"], 0)
            self.assertEqual(qc["chain_provenance"]["contigs"][0]["status"], "pass")
            self.assertEqual(len(qc["chain_provenance"]["contigs"][0]["md5"]), 32)
            timing = json.loads((root / "TEST/qc/TEST.stage11_timing.json").read_text())
            self.assertEqual(timing["hap_jobs"], 2)
            self.assertEqual(timing["route_jobs"], 4)
            self.assertEqual(timing["status"], "ok")
            self.assertTrue(timing["ephemeral_bwa_index_removed"])
            self.assertEqual(timing['bwa_index_jobs'], 2)
            execution = timing['bwa_index_execution']
            overlap = min(row['command_finished_at'] for row in execution) - max(row['command_started_at'] for row in execution)
            self.assertGreater(overlap, 0)
            self.assertTrue(all(row['actual_cpu_seconds'] > 0 and row['max_rss_bytes'] > 0 for row in execution))
            self.assertFalse((root / "TEST/work/personal_realign/bwa_mem2_index").exists())
            self.assertEqual(timing["input_count_mode"], "fused_with_prepare")
            self.assertEqual(routing["input_record_counts"]["all_primary_mapped"], 6)
            self.assertEqual(routing["prepared_inputs"]["HAP1"]["route_jobs"], 2)
            self.assertIn("align_restore_sort", timing["haplotypes"]["HAP1"])
            self.assertIn("align_restore_sort", timing["haplotypes"]["HAP2"])
            for hap in (1, 2):
                resource_qc = json.loads(
                    (root / f"TEST/logs/personal_realign/TEST.HAP{hap}.resource.json").read_text()
                )
                self.assertEqual(resource_qc["exit_code"], 0)
                self.assertGreater(resource_qc["max_rss_kib"], 0)

            # 같은 입력을 HAP 직렬로 다시 실행해 병렬화가 record/QC 계약을
            # 바꾸지 않는지 확인한다. 경로·시간·cache created/reused 값은 제외한다.
            serial_root = root / "serial"
            serial_config = replace(config, output_root=serial_root, hap_jobs=1, route_jobs=1)
            with patch.dict(os.environ, {'MFM_REALIGN_PHASE_RESOURCES': '0'}):
                serial_result = module.run_personal_realign_generation(serial_config)
            self.assertEqual(serial_result["ok_count"], 1)
            for hap in (1, 2):
                parallel_dir = root / f"TEST/realigned/HAP{hap}"
                serial_dir = serial_root / f"TEST/realigned/HAP{hap}"
                parallel_cram = parallel_dir / f"TEST.personal.HAP{hap}.cram"
                serial_cram = serial_dir / f"TEST.personal.HAP{hap}.cram"
                def records(path: Path, ref: Path):
                    with pysam.AlignmentFile(path, "rc", reference_filename=str(ref)) as handle:
                        return [
                            (read.query_name, read.flag, read.reference_id, read.reference_start,
                             read.mapping_quality, read.cigarstring, read.query_sequence,
                             read.qual, tuple(read.get_tags()))
                            for read in handle.fetch(until_eof=True)
                        ]
                self.assertEqual(
                    records(parallel_cram, parallel_dir / "target.fa"),
                    records(serial_cram, serial_dir / "target.fa"),
                )
                parallel_qc = json.loads((parallel_dir / "TEST.realign.json").read_text())
                serial_qc = json.loads((serial_dir / "TEST.realign.json").read_text())
                for key in ("records", "routed", "rescue_unmapped", "qc_contract", "alignment_contract"):
                    self.assertEqual(parallel_qc[key], serial_qc[key])
            print('[REALIGN-VALIDATION] ' + json.dumps({'index_command_overlap_seconds': overlap,
                'old_new_alignment_records_equal': True, 'qc_contract_equal': True,
                'index_execution': execution, 'phases': timing['resource_phases']}, sort_keys=True))

            # 한 index가 실패해도 다른 index가 끝나기 전 큰 예약을 반환하지 않음.
            pool = root / 'failure_pool'
            (pool / 'leases').mkdir(parents=True)
            lease = pool / 'leases/test.json'
            lease.write_text(json.dumps({'cpus': 72, 'total_cpus': 72, 'memory_bytes': 192 * 1024 ** 3}))
            other_finished = threading.Event()
            original_index = module.ensure_bwa_mem2_index
            def fail_one(reference, **kwargs):
                if reference.parent.name == 'HAP1':
                    raise RuntimeError('의도한 index 실패')
                time.sleep(.05)
                try:
                    return original_index(reference, **kwargs)
                finally:
                    other_finished.set()
            with patch.dict(os.environ, {'MFM_CPU_POOL_LEASE_FILE': str(lease), 'MFM_REALIGN_PHASE_RESOURCES': '1'}):
                with patch.object(module, 'ensure_bwa_mem2_index', side_effect=fail_one):
                    with self.assertRaisesRegex(RuntimeError, 'personal realign failed'):
                        module.run_personal_realign_generation(replace(config, output_root=root / 'failed'))
            failed = json.loads((root / 'failed/test.personal_realign.metadata.json').read_text())
            self.assertEqual(failed['fail_count'], 1)
            self.assertTrue(other_finished.is_set())
            self.assertEqual(json.loads(lease.read_text())['cpus'], 1)
            self.assertEqual(json.loads(lease.read_text())['memory_bytes'], 4 * 1024 ** 3)


if __name__ == "__main__":
    unittest.main()
