"""분할정렬과 SQLite의정체성/중복/오염판정이동일함을검증함."""
from pathlib import Path
from collections import Counter
import hashlib,json,os,random,sys,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from integrate_wgs_rna import stage70_trace_validator as trace
from integrate_wgs_rna.compact.trace_multiset import PartitionedBag
from integrate_wgs_rna.compact.trace_decode import decode_reads,read_digest,segment_digest


class MultisetTest(unittest.TestCase):
    def test_random_counts_mismatch_examples_and_partition_limit(self):
        rng=random.Random(261020)
        with tempfile.TemporaryDirectory() as raw:
            r=Path(raw);old=trace._Bag(r/'old.db',batch_size=19);new=PartitionedBag(r/'new',leaf_bytes=1024,batch_records=19)
            try:
                for _ in range(7000):
                    row=(rng.choice(['src','out']),rng.choice(['read','segment']),hashlib.sha256(str(rng.randrange(170)).encode()).digest())
                    old.add(*row);new.add(*row)
                self.assertEqual(old.result(5),new.result(5))
                self.assertEqual(old.result(2),new.result(2))
                self.assertLessEqual(new.stats['leaf_peak_raw_bytes'],1024)
                self.assertGreater(new.stats['partition_passes'],0)
                row=('src','other',b'x'*32);old.add(*row);new.add(*row)
                self.assertEqual(old.result(3),new.result(3))
            finally:old.close();new.close()

    def test_extreme_duplicate_prefix_and_zero_bytes(self):
        with tempfile.TemporaryDirectory() as raw:
            r=Path(raw);old=trace._Bag(r/'old.db');new=PartitionedBag(r/'new',leaf_bytes=64)
            try:
                for i in range(220):
                    digest=bytes([i%4])+bytes(31)
                    for side in ['src','out']:old.add(side,'read',digest);new.add(side,'read',digest)
                old.add('out','read',bytes(32));new.add('out','read',bytes(32))
                self.assertEqual(old.result(10),new.result(10))
                self.assertLessEqual(new.stats['leaf_peak_raw_bytes'],64)
                self.assertGreaterEqual(new.stats['partition_passes'],32)
            finally:old.close();new.close()

    def test_digest_and_batch_decode_match_original_canonical(self):
        rng=random.Random(42)
        for _ in range(500):
            seq=''.join(rng.choices('ACGTN',k=rng.randrange(300)))
            values=[rng.randint(-100,100000) for _ in range(5)];cigar=tuple(rng.randint(0,2**32-1) for _ in range(rng.randrange(8)))
            self.assertEqual(read_digest(seq,*values,cigar),trace._json_hash((seq,*values,cigar)))
            fields=[rng.randint(-100,100000) for _ in range(6)];digest=hashlib.sha256(seq.encode()).digest()
            self.assertEqual(segment_digest(digest,*fields),trace._json_hash((digest.hex(),*fields)))
        # record별홀수nibble·invalid padding을기존과같이무시하고실제invalid는검출함.
        packed=np.array([0x01,0x4f,0x23],dtype=np.uint8);offsets=np.array([10,12,13]);lengths=np.array([3,2])
        self.assertEqual(decode_reads(packed,offsets,lengths),['ACN','GT'])
        with self.assertRaisesRegex(ValueError,'invalid packed read'):decode_reads(packed,offsets,[4,2])


    def test_packed_batch_adapts_without_dropping_long_reads(self):
        from unittest.mock import patch
        n, length = 9, 600000
        reads = {'lengths': np.full(n, length), 'offsets': np.arange(n+1)*((length+1)//2),
                 'bases_packed': np.zeros(n*((length+1)//2), dtype=np.uint8)}
        for name in ['hp','ps','ref_start','ref_end','segment_type','cigar_offset','cigar_len']:
            reads[name] = np.zeros(n, dtype=np.int64)
        root = {'reads': reads, 'cigar': {'ops_packed': np.empty(0, dtype=np.uint32)}}
        class Bag:
            def __init__(self): self.rows = []
            def add(self, side, kind, digest): self.rows.append((side,kind,digest))
        bag = Bag()
        with patch.dict(os.environ, {'MFM_TRACE_PREFETCH_PACKED_MIB':'1'}):
            ids = trace._output_reads(root, bag, prefetch_records=65536)
        expected = trace._json_hash(('A'*length,0,0,0,0,0,()))
        self.assertEqual(len(bag.rows), n)
        self.assertTrue(all(bytes(value)==expected for value in ids))
        reads['offsets'] = np.arange(n+1)*(2<<20)
        with patch.dict(os.environ, {'MFM_TRACE_PREFETCH_PACKED_MIB':'1'}):
            with self.assertRaisesRegex(ValueError,'single read exceeds'):
                trace._output_reads(root, Bag())

    def test_scratch_admission_and_sqlite_fallback(self):
        from unittest.mock import patch
        from types import SimpleNamespace
        from integrate_wgs_rna.compact import trace_resources as resources
        estimate = resources.scratch_estimate(100,200,300)
        self.assertGreater(estimate['partitioned_estimate_bytes'], estimate['sqlite_estimate_bytes'])
        with tempfile.TemporaryDirectory() as directory:
            for free, expected in [(1000000,'partitioned'),(100000,'sqlite')]:
                with patch.object(resources.shutil,'disk_usage',return_value=SimpleNamespace(free=free)):
                    plan = resources.select_backend(directory,[estimate],'auto')
                    self.assertEqual(plan['backend'], expected)
                    self.assertLessEqual(plan['estimated_scratch_bytes'], free*9//10)
            with patch.object(resources.shutil,'disk_usage',return_value=SimpleNamespace(free=100000)):
                with self.assertRaisesRegex(ValueError,'scratch'):
                    resources.select_backend(directory,[estimate],'partitioned')
            with patch.object(resources.shutil,'disk_usage',return_value=SimpleNamespace(free=1)):
                with self.assertRaisesRegex(ValueError,'scratch'):
                    resources.select_backend(directory,[estimate],'auto')
            with patch.object(Path,'read_text',return_value='1 0 0:1 / / rw - ext4 /dev/root rw\n2 1 0:2 / /tmp rw - tmpfs tmpfs rw\n'):
                self.assertTrue(resources.is_memory_filesystem('/tmp/a'))
                self.assertFalse(resources.is_memory_filesystem('/NHNHOME/a'))


if __name__=='__main__':unittest.main()
