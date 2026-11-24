"""QNAME 경계·기존 집계와의 일치·자원 예산·실패 시 미반영을 검증함."""
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from integrate_wgs_rna.realign.input_metrics import METRIC_KEYS, complete_qname_batches, count_sam_batch, metric_plan
from integrate_wgs_rna.realign.qc import add_alignment_input_metrics, _add_alignment_input_metrics_serial


def record(name, flag, tags=()):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.flag = flag
    read.query_sequence = "A" * 30
    read.query_qualities = pysam.qualitystring_to_array("I" * 30)
    read.reference_id = 0
    read.reference_start = 10
    read.mapping_quality = 60
    read.cigarstring = "30M"
    for key, value, kind in tags:
        read.set_tag(key, value, value_type=kind)
    return read


class InputMetricsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.path = self.root / "collated.bam"
        self.rows = []
        for i in range(30):
            self.rows.extend([
                record(f"normal{i}", 64, [("HP", 1, "i"), ("OD", 1, "i")]),
                record(f"normal{i}", 128, [("HP", 2, "i")]),
                record(f"rescue{i}", 64, [("XU", 1, "i"), ("OD", "1", "Z")]),
                record(f"rescue{i}", 128, [("HP", 1.9, "f")]),
                record(f"extra{i}", 64), record(f"extra{i}", 64, [("HP", 2, "i")]),
                record(f"secondary{i}", 256 | 64),
                record(f"extra{i}", 128, [("HP", "bad", "Z")]),
                record(f"extra{i}", 192),
                record(f"singleton{i}", 128), record(f"other{i}", 0),
                record(f"supplementary{i}", 2048 | 128),
            ])
        self.write(self.rows)

    def tearDown(self):
        self.temp.cleanup()

    def write(self, rows):
        with pysam.AlignmentFile(str(self.path), "wb", header={"HD": {"VN": "1.6", "SO": "unsorted"},
                                 "SQ": [{"SN": "chr1", "LN": 1000}]}) as out:
            for read in rows:
                out.write(read)

    def test_parallel_matches_serial_at_many_pair_boundaries(self):
        baseline = {}
        _add_alignment_input_metrics_serial(self.path, baseline, threads=0)
        actual = {}
        add_alignment_input_metrics(self.path, actual, workers=2, chunk_bytes=257)
        self.assertEqual({k: actual[k] for k in METRIC_KEYS}, baseline)
        execution = actual["input_metrics_execution"]
        self.assertEqual(execution["mode"], "parallel_sam_stream")
        self.assertGreater(execution["workers_used"], 0)
        self.assertLessEqual(execution["peak_inflight"], 4)
        self.assertEqual({p.name for p in self.root.iterdir()}, {"collated.bam"})

    def test_chunker_never_splits_a_qname_group(self):
        data = b"".join(f"q{i}\t{flag}\tx\n".encode() for i in range(20) for flag in [64, 128, 0])
        batches = list(complete_qname_batches(io.BytesIO(data), 17))
        self.assertEqual(b"".join(batches), data)
        seen = set()
        for batch in batches:
            names = {line.split(b"\t", 1)[0] for line in batch.splitlines()}
            self.assertFalse(names & seen)
            seen.update(names)
        with self.assertRaisesRegex(ValueError, "buffer"):
            list(complete_qname_batches(io.BytesIO(b"q\t64\tx\n" * 40), 16, max_group_bytes=32))

    def test_resource_plan_honors_grants_memory_and_actual_lease(self):
        with patch.dict(os.environ, {"MFM_CPU_POOL_GRANTED": "8", "MFM_RAM_LIMIT_BYTES": "1000000000000"}, clear=True):
            plan = metric_plan(self.path, 4, workers=100)
            self.assertLessEqual(plan["workers_requested"] + plan["samtools_aux_threads"] + 2, 8)
            self.assertLessEqual(plan["estimated_peak_rss_bytes"], plan["ram_budget_bytes"])
            lease = self.root / "lease.json"
            lease.write_text(json.dumps({"cpus": 12, "memory_bytes": 8 * 1024 ** 3}))
            os.environ["MFM_CPU_POOL_LEASE_FILE"] = str(lease)
            plan = metric_plan(self.path, 4, workers=100)
            self.assertEqual(plan["cpu_budget"], min(12, len(os.sched_getaffinity(0))))
            self.assertLessEqual(plan["ram_budget_bytes"], 8 * 1024 ** 3)
        with patch.dict(os.environ, {"MFM_CPU_POOL_GRANTED": "1", "SLURM_MEM_PER_NODE": "0", "SLURM_MEM_PER_CPU": "0"}, clear=True):
            plan = metric_plan(self.path, 4, workers=8)
            self.assertEqual((plan["mode"], plan["samtools_aux_threads"]), ("serial_pysam", 0))

    def test_failed_child_does_not_commit_partial_metrics(self):
        self.write([record("bad", 64, [("XU", "bad", "Z")]), record("bad", 128)])
        actual = {"preserved": 17}
        with self.assertRaises(ValueError):
            add_alignment_input_metrics(self.path, actual, workers=2, chunk_bytes=257)
        self.assertEqual(actual, {"preserved": 17})

    def test_empty_and_failed_samtools(self):
        self.write([])
        actual = {}
        add_alignment_input_metrics(self.path, actual, workers=2, chunk_bytes=257)
        self.assertTrue(all(actual[k] == 0 for k in METRIC_KEYS))
        actual = {}
        with self.assertRaisesRegex(RuntimeError, "samtools"):
            add_alignment_input_metrics(self.path, actual, workers=2, samtools_bin="/bin/false", chunk_bytes=257)
        self.assertEqual(actual, {})


if __name__ == "__main__":
    unittest.main()
