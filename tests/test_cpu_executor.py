"""실제 프로세스 병렬 계산의 CPU 제한, 증액·반환과 예외를 검증한다."""
import json
import hashlib
import os
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from integrate_wgs_rna.cpu_executor import WorkerBudget, cpu_results


def task(value, delay=0.1):
    started = time.monotonic()
    time.sleep(delay)
    if value < 0:
        raise ValueError("실패 작업")
    return value, started, time.monotonic(), os.getpid()


class CPUExecutorTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / "leases").mkdir()
        (self.root / "waiters").mkdir()
        self.lease = self.root / "leases/own.json"

    def tearDown(self):
        self.tmp.cleanup()

    def configure(self, cpus, total):
        self.lease.write_text(json.dumps({"cpus": cpus, "total_cpus": total}))
        return patch.dict(os.environ, MFM_CPU_POOL_LEASE_FILE=str(self.lease))

    def test_prime_grant_returns_unused_capacity(self):
        for grant in (11, 61):
            with self.configure(grant, 72):
                budget = WorkerBudget(grant, 22, True)
                self.assertEqual(min(grant, 22), budget.capacity(22, grow=False))
                self.assertEqual(2, budget.capacity(2, grow=False))
                self.assertEqual(0 if budget.memory_enabled else 1, budget.capacity(0, grow=False))
                self.assertEqual(1,json.loads(self.lease.read_text())["cpus"])

    def test_multithread_worker_returns_all_but_metadata_cpu(self):
        with self.configure(12, 12):
            budget = WorkerBudget(4, 4, True, cpus_per_worker=3)
            self.assertEqual(4, budget.capacity(4, grow=False))
            self.assertEqual(12, json.loads(self.lease.read_text())["cpus"])
            budget.capacity(0, grow=False)
            self.assertEqual(1, json.loads(self.lease.read_text())["cpus"])

    def test_waiting_stage_has_priority_and_explicit_limit_is_preserved(self):
        with self.configure(2, 8):
            waiter = self.root / "waiters/new.json"
            waiter.write_text("{}")
            budget = WorkerBudget(2, 6, True)
            self.assertEqual(2, budget.capacity(10))
            waiter.unlink()
            self.assertEqual(6, budget.capacity(10))
            explicit = WorkerBudget(2, 2, False)
            self.assertEqual(2, explicit.capacity(10))

    def test_unfittable_ram_waiter_does_not_starve_two_cpu_worker(self):
        with self.configure(1, 8), patch.dict(os.environ, MFM_RAM_LIMIT_BYTES=str(64 * 1024**3)):
            (self.root / "waiters/large.json").write_text(json.dumps(
                {"min_cpus": 2, "memory_bytes": 128 * 1024**3}))
            budget = WorkerBudget(1, 2, True, 2 * 1024**3, cpus_per_worker=2)
            self.assertEqual(2, budget.capacity(2))
            self.assertEqual(4, json.loads(self.lease.read_text())["cpus"])

    def test_fitting_waiter_keeps_priority_during_two_cpu_growth(self):
        with self.configure(1, 8):
            (self.root / "waiters/small.json").write_text(json.dumps(
                {"min_cpus": 2, "memory_bytes": 0}))
            self.assertEqual(0, WorkerBudget(1, 2, True, cpus_per_worker=2).capacity(2))

    def test_idle_loan_does_not_block_regular_worker_growth(self):
        (self.root / "leases/loan.json").write_text(json.dumps(
            {"cpus": 8, "memory_bytes": 0, "cpu_backfill": {"nice": 19}}))
        with self.configure(1, 8):
            self.assertEqual(4, WorkerBudget(1, 4, True).capacity(4))

    def test_fixed_base_memory_waits_without_overreserving(self):
        with self.configure(1, 8), patch.dict(os.environ, MFM_RAM_LIMIT_BYTES=str(64 * 1024**3)):
            own = json.loads(self.lease.read_text());own['memory_bytes'] = 4 * 1024**3
            self.lease.write_text(json.dumps(own))
            (self.root / 'leases/other.json').write_text(json.dumps(
                {'cpus': 1, 'memory_bytes': 48 * 1024**3}))
            budget = WorkerBudget(1, 8, True);budget.base_memory = 32 * 1024**3
            self.assertEqual(0, budget.capacity(8))
            self.assertEqual(4 * 1024**3, json.loads(self.lease.read_text())['memory_bytes'])
            (self.root / 'leases/other.json').unlink()
            self.assertEqual(8, budget.capacity(8))
            self.assertEqual(32 * 1024**3, json.loads(self.lease.read_text())['memory_bytes'])

    def test_actual_workers_grow_when_other_stage_releases(self):
        other = self.root / "leases/other.json"
        other.write_text(json.dumps({"cpus": 4}))
        with self.configure(2, 6):
            release = threading.Timer(0.18, other.unlink)
            release.start()
            try:
                rows = list(cpu_results(task, [(x, 0.14) for x in range(12)], 2,
                                        elastic=True, max_workers=6))
            finally:
                release.join()
        values = [value for _, value, error in rows if error is None]
        self.assertEqual(list(range(12)), sorted(row[0] for row in values))
        peak = max(sum(start <= row[1] < end for _, start, end, _ in values) for row in values)
        self.assertGreaterEqual(peak, 3)
        self.assertLessEqual(peak, 6)
        self.assertEqual(1, json.loads(self.lease.read_text())["cpus"])

    def test_failure_is_reported_and_lease_returned(self):
        with self.configure(8, 8):
            rows = list(cpu_results(task, [(1,), (-1,), (2,)], 8, elastic=True))
        self.assertEqual(1, sum(isinstance(error, ValueError) for _, _, error in rows))
        self.assertEqual(1, json.loads(self.lease.read_text())["cpus"])

    def test_real_compact_parallel_matches_serial_and_resumes(self):
        import pysam
        import zarr
        from pre_training.preprocess_raw_readset.config import ExportConfig
        from integrate_wgs_rna.direct_compact_writer import run_direct_compact_jobs
        contigs = ["chr1", "chr2"]
        sequence = "ACGT" * 16
        rows = [{"sample": "TEST", "contig": chrom} for chrom in contigs]
        for hap in (1, 2):
            target = self.root / f"HAP{hap}/target.fa"
            target.parent.mkdir()
            target.write_text("".join(f">{chrom}\n{sequence}\n" for chrom in contigs))
            pysam.faidx(str(target))
            cram = target.parent / "source.cram"
            header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [
                {"SN": chrom, "LN": len(sequence), "M5": hashlib.md5(sequence.encode()).hexdigest()}
                for chrom in contigs]}
            with pysam.AlignmentFile(str(cram), "wc", header=header, reference_filename=str(target)) as out:
                for tid, chrom in enumerate(contigs):
                    read = pysam.AlignedSegment(out.header)
                    read.query_name = f"read{tid}"
                    read.query_sequence = "ACGT"
                    read.query_qualities = pysam.qualitystring_to_array("IIII")
                    read.reference_id, read.reference_start = tid, 0
                    read.flag, read.mapping_quality, read.cigarstring = 0, 60, "4M"
                    out.write(read)
                    crossing = pysam.AlignedSegment(out.header)
                    crossing.query_name = f"cross{tid}"
                    crossing.query_sequence = "ACGTACGTACGT"
                    crossing.query_qualities = pysam.qualitystring_to_array("I" * 12)
                    crossing.reference_id, crossing.reference_start = tid, 4
                    crossing.flag, crossing.mapping_quality, crossing.cigarstring = 0, 60, "12M"
                    crossing.set_tag("HP", hap)
                    crossing.set_tag("PS", 17)
                    out.write(crossing)
                    fasta = target.parent / f"{chrom}.fa"
                    fasta.write_text(f">{chrom}\n{sequence}\n")
                    pysam.faidx(str(fasta))
                    rows[tid].update({f"hap{hap}_fasta": str(fasta), f"hap{hap}_chain": "",
                                      f"hap{hap}_alignment": str(cram), "bam": str(cram)})
            pysam.index(str(cram))
        record_tsv = self.root / "records.tsv"
        record_tsv.write_text("sample\tcontig\nTEST\tchr1\nTEST\tchr2\n")
        windows = self.root / "windows"
        windows.mkdir()
        (windows / "windows.tsv").write_text(
            "contig\thaplotype\traw_start\tactual_length\tkept\n" + "".join(
                f"{chrom}\tHAP{hap}\t{start}\t8\t1\n"
                for chrom in contigs for hap in (1, 2) for start in (0, 8)))
        config = ExportConfig(window_size=8, offset_unit=1, curriculum_unit=8, seed=42,
                              min_mapping_quality=30, min_coordinate_coverage=0.0)
        def run(output, workers, batch=1):
            return run_direct_compact_jobs(records=rows, record_tsv=record_tsv, output_root=output,
                windows_per_chunk=batch, jobs=workers, schedule="longest-first", max_windows=2,
                validation_level="full", overwrite=False, resume=True, reference_fasta=None,
                export_config=config, sample_metadata={}, read_coordinate_mode="personal_realigned",
                window_index_dir=str(windows))
        def digests(output):
            result = {}
            for group_path in output.glob("chr*/*.zarr"):
                group = zarr.open_group(str(group_path), mode="r")
                def visit(name, obj):
                    if isinstance(obj, zarr.Array):
                        result[str(group_path.relative_to(output)) + "/" + name] = hashlib.sha256(obj[:].tobytes()).hexdigest()
                group.visititems(visit)
            return result
        serial, parallel = self.root / "serial", self.root / "parallel"
        with patch.dict(os.environ, MFM_COMPACT_INDEX_MODE="full"):
            run(serial, 1, batch=32)
        with self.configure(2, 4), patch.dict(os.environ, MFM_CPU_UTILIZATION="maximum",
                                             MFM_COMPACT_JOBS_REQUESTED="auto", MFM_BLOSC_THREADS="1"):
            run(parallel, 2)
            before = digests(parallel)
            self.assertTrue(before)
            self.assertEqual(digests(serial), before)
            resumed = run(parallel, 2)
            self.assertTrue(all(row["status"] == "skip_existing" for row in resumed["records"].values()))
            self.assertEqual(before, digests(parallel))
            self.assertEqual(1, json.loads(self.lease.read_text())["cpus"])


if __name__ == "__main__":
    unittest.main()
