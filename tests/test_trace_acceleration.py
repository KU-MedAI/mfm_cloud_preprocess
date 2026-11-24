"""독립 검증 가속의 경계·무작위 oracle·변조 회귀. writer helper를 oracle로 쓰지 않는다."""

from __future__ import annotations

from collections import Counter
import hashlib
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from integrate_wgs_rna import stage70_trace_validator as trace
from integrate_wgs_rna.trace_acceleration import DepthAccumulator, WindowLookup


def make_read(start, cigar, flag=0, name="r"):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = "A" * sum(length for op, length in cigar if op in (0, 1, 4, 7, 8))
    read.query_qualities = pysam.qualitystring_to_array("I" * len(read.query_sequence))
    read.reference_id, read.reference_start, read.mapping_quality = 0, start, 60
    read.flag, read.cigartuples = flag, cigar
    read.set_tag("HP", 1)
    read.set_tag("PS", 7)
    return read


class MemoryBag:
    def __init__(self):
        self.rows = {"src": Counter(), "out": Counter()}

    def add(self, side, kind, digest):
        self.rows[side][(kind, digest)] += 1


class TraceAccelerationTest(unittest.TestCase):
    def test_window_candidates_match_exhaustive_projection(self):
        rng = random.Random(261001)
        layouts = [
            [(i * 10, 10, 12) for i in range(30)],
            [(i * 10, 6, 8) for i in range(30) if i % 3 != 1],
            [(i * 7, 13, 16) for i in range(30)],
            [(0, 30, 30), (10, 5, 5), (20, 20, 20)],  # end 비단조: 전수 fallback
        ]
        reads = [
            make_read(0, [(1, 3), (0, 10), (1, 4)]),
            make_read(10, [(5, 2), (4, 3), (0, 10), (1, 4), (4, 2), (5, 1)]),
            make_read(0, [(0, 5), (2, 10), (1, 3), (3, 10), (7, 3), (6, 2), (8, 2)]),
            make_read(10, [(1, 4)]),
        ]
        for _ in range(500):
            cigar = [(4, 2)] + [(rng.choice((0, 1, 2, 3, 7, 8)), rng.randint(1, 12)) for _ in range(6)] + [(0, 2), (1, 3), (4, 2)]
            reads.append(make_read(rng.randrange(250), cigar))
        for windows in layouts:
            lookup = WindowLookup(windows)
            for read in reads:
                old = [(i, trace._project(read, a, a + b)) for i, (a, b, _) in enumerate(windows)]
                old = [(i, p) for i, p in old if p is not None]
                new = [(i, trace._project(read, windows[i][0], windows[i][0] + windows[i][1])) for i in lookup.candidates(read)]
                new = [(i, p) for i, p in new if p is not None]
                self.assertEqual(old, new, (windows, read.cigarstring, read.reference_start))

    def test_invalid_cigar_not_hidden_by_empty_candidates(self):
        lookup = WindowLookup([(0, 10, 10)])
        for op, length, message in ((0, 0, "invalid source CIGAR length"), (9, 1, "unsupported source CIGAR op: 9")):
            with self.assertRaisesRegex(ValueError, message):
                lookup.candidates(make_read(100, [(0, 1), (op, length)]))

    def test_depth_endpoints_equal_per_base_oracle_and_wrap(self):
        windows = [(0, 100, 110), (105, 90, 100)]
        accumulator = DepthAccumulator(windows)
        expected = [np.zeros(stored, np.uint32) for _, _, stored in windows]
        rng = random.Random(261001)
        for _ in range(1000):
            index = rng.randrange(2)
            start = rng.randrange(windows[index][1])
            end = rng.randrange(start + 1, windows[index][1] + 1)
            accumulator.add(index, start, end)
            expected[index][start:end] += 1
        for i in range(2):
            np.testing.assert_array_equal(accumulator.finish_window(i), expected[i])
        wrap = DepthAccumulator([(0, 3, 3)])
        wrap.buffers[0][0] = np.uint32(0xFFFFFFFF)
        wrap.add(0, 0, 3)
        np.testing.assert_array_equal(wrap.finish_window(0), np.zeros(3, np.uint32))

    def test_source_audit_digests_segments_equal_exhaustive(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            target = root / "target.fa"
            target.write_text(">chr1\n" + "A" * 1000 + "\n")
            pysam.faidx(str(target))
            bam = root / "reads.bam"
            reads = []
            for i in range(90):
                flag = (0, 1024, 2048, 256, 512, 16)[i % 6]
                read = make_read(i * 7, [(4, 2), (0, 10), (2, 5), (1, 3), (7, 15), (3, 6), (8, 10), (1, 4), (4, 2)], flag, f"r{i}")
                reads.extend((read, read) if i % 9 == 0 else (read,))
            with pysam.AlignmentFile(str(bam), "wb", header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}) as out:
                for read in reads:
                    out.write(read)
            pysam.index(str(bam))
            windows = [(i * 50, 40, 50) for i in range(18) if i % 3 != 1]
            for duplicate, supplementary in ((False, False), (True, False), (True, True)):
                config = {"min_mapping_quality": 30, "exclude_duplicate": duplicate, "exclude_supplementary": supplementary}
                observed = trace._source(bam, target, "chr1", windows, config)
                with patch.object(WindowLookup, "candidates", lambda self, read: range(len(self.starts))):
                    oracle = trace._source(bam, target, "chr1", windows, config)
                self.assertEqual(observed, oracle)
                from integrate_wgs_rna.compact.trace_store import ReplayRecords
                records = ReplayRecords(root / "replay.records")
                try:
                    replay, audit = trace._source(bam, target, "chr1", windows, config, records)
                    self.assertEqual(observed, (list(replay), audit))
                    self.assertEqual(list(replay), list(replay))
                finally:
                    records.close()

    def test_disk_replay_dedup_survives_pruning_and_disjoint_intervals(self):
        from integrate_wgs_rna.compact.trace_store import ReplayRecords
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw); target = root / "target.fa"
            target.write_text(">chr1\n" + "A" * 20000 + "\n")
            pysam.faidx(str(target))
            bam = root / "reads.bam"
            with pysam.AlignmentFile(str(bam), "wb", header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 20000}]}) as out:
                # 한 read가 떨어진 두 fetch interval에 걸쳐도 한 번만 기록한다.
                out.write(make_read(0, [(0, 1), (3, 12999), (0, 1)], name="crossing"))
                for i in range(12020):
                    read = make_read(i, [(0, 12)], name=f"r{i}")
                    out.write(read)
                    if i % 1000 == 0:
                        out.write(read)
            pysam.index(str(bam))
            windows = [(0, 11000, 11000), (12000, 2000, 2000)]
            config = {"min_mapping_quality": 30}
            expected = trace._source(bam, target, "chr1", windows, config)
            records = ReplayRecords(root / "source.records")
            try:
                source, audit = trace._source(bam, target, "chr1", windows, config, records)
                self.assertEqual(expected, (list(source), audit))
                self.assertEqual(1, sum(record["ref_start"] == 0 and len(record["cigar"]) == 3 for record in source))
            finally:
                records.close()

    def test_projection_primary_supplementary_and_mutation(self):
        windows = [(0, 10, 12), (10, 10, 12)]
        source = []
        for i, read in enumerate((make_read(0, [(0, 5), (2, 2), (0, 8), (1, 3)]), make_read(5, [(0, 2), (3, 4), (0, 8)], 2048))):
            sequence, clip, cigar = trace._sequence_and_cigar(read)
            segments = []
            for index, (start, actual, _) in enumerate(windows):
                p = trace._project(read, start, start + actual)
                if p:
                    q0, q1, h0, h1 = p
                    segments.append((index, q0 - clip, q1 - clip, h0 - start, h1 - start))
            source.append({"ref_start": read.reference_start, "cigar": cigar, "segment_type": int(read.is_supplementary), "digest": hashlib.sha256(str(i).encode()).digest(), "strand": 1, "segments": segments})
        columns = {k: [] for k in ("read_id", "query_start", "query_end", "hap_start", "hap_end", "strand")}
        offsets = [0]
        for index in range(len(windows)):
            for read_id, record in enumerate(source):
                for wi, q0, q1, h0, h1 in record["segments"]:
                    if wi == index:
                        for name, value in zip(columns, (read_id, q0, q1, h0, h1, 1)):
                            columns[name].append(value)
            offsets.append(len(columns["read_id"]))
        segments = {name: np.asarray(values) for name, values in columns.items()}
        segments["window_read_offsets"] = np.asarray(offsets)
        ids = [r["digest"] for r in source]
        for allowed in ((0,), (1,), (0, 1)):
            oracle = [np.zeros(stored, np.uint32) for _, _, stored in windows]
            for record in source:
                if record["segment_type"] in allowed:
                    for index, *_ in record["segments"]:
                        start, actual, _ = windows[index]
                        for a, b in trace._match_intervals(record["ref_start"], record["cigar"]):
                            lo, hi = max(a, start) - start, min(b, start + actual) - start
                            if hi > lo:
                                oracle[index][lo:hi] += 1
            depth = {"offsets": np.asarray([0, 12, 24]), "depth": np.concatenate(oracle)}
            bag = MemoryBag()
            result = trace._check_projection(segments, depth, windows, source, ids, bag, 2, allowed)
            self.assertEqual(result["depth_mismatch_bases"], 0)
            self.assertEqual(bag.rows["src"], bag.rows["out"])
            unused = MemoryBag()
            self.assertEqual(result, trace._check_projection(segments, depth, windows, source, ids, unused, 2, allowed, check_segments=False))
            self.assertFalse(unused.rows["src"] or unused.rows["out"])
            depth["depth"][0] += 1
            failed = trace._check_projection(segments, depth, windows, source, ids, MemoryBag(), 2, allowed)
            self.assertEqual(failed["depth_mismatch_bases"], 1)
            self.assertEqual(failed["examples"][0]["raw_position"], 0)
        segments["query_start"][0] += 1
        changed = MemoryBag()
        trace._check_projection(segments, depth, windows, source, ids, changed, 2, (0,))
        self.assertNotEqual(changed.rows["src"], changed.rows["out"])

    def test_sqlite_without_rowid_counts_and_order(self):
        rng = random.Random(261001)
        expected = {"src": Counter(), "out": Counter()}
        with tempfile.TemporaryDirectory() as raw:
            bag = trace._Bag(Path(raw) / "bag.sqlite3", batch_size=7)
            try:
                sql = bag.db.execute("SELECT sql FROM sqlite_master WHERE name='x'").fetchone()[0]
                self.assertIn("WITHOUT ROWID", sql)
                self.assertEqual(bag.db.execute("PRAGMA cache_size").fetchone()[0], -1024 * trace.SQLITE_CACHE_MIB)
                for _ in range(3000):
                    side, kind = rng.choice(("src", "out")), rng.choice(("read", "segment"))
                    digest = hashlib.sha256(str(rng.randrange(80)).encode()).digest()
                    expected[side][(kind, digest)] += 1
                    bag.add(side, kind, digest)
                bag.result(3)
                observed = list(bag.db.execute("SELECT kind,digest,src,out FROM x ORDER BY kind,digest"))
                keys = sorted(expected["src"].keys() | expected["out"].keys())
                self.assertEqual(observed, [(kind, digest, expected["src"][(kind, digest)], expected["out"][(kind, digest)]) for kind, digest in keys])
                result = bag.result(3)
                self.assertFalse(result["pass"])
                for kind, item in result["mismatches"].items():
                    mismatched = [(digest, expected["src"][(k, digest)], expected["out"][(k, digest)]) for k, digest in keys if k == kind and expected["src"][(k, digest)] != expected["out"][(k, digest)]]
                    self.assertEqual(item["identity_count"], len(mismatched))
                    self.assertEqual(item["examples"], [{"sha256": d.hex(), "source_count": s, "zarr_count": o} for d, s, o in mismatched[:3]])
            finally:
                bag.close()


if __name__ == "__main__":
    unittest.main()
