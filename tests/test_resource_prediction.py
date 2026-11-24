"""RAM에 맞는 계산 진입과 RSS 예약·주소공간·캐시의 독립성을 검증함."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import json
import os
import resource
import tempfile
import unittest
import numpy as np
from integrate_wgs_rna.memory_policy import GIB
from integrate_wgs_rna.realign.align_resources import alignment_plan
from integrate_wgs_rna.compact.memory_preflight import estimate_peak, forecast_compact
from integrate_wgs_rna.cpu_executor import cpu_results, WorkerBudget
from integrate_wgs_rna.compact.validation import ReadValidationChunkCache
from integrate_wgs_rna.direct_compact_writer import _ReadValidationChunkCache
from integrate_wgs_rna.compact.runtime import nonnegative_env_int

AUTO = {'MFM_WGS_BACKEND':'cpu','WGS_AUTO_THREADS':'1',
        'MFM_WGS_ALIGN_REQUESTED':'auto','MFM_WGS_SORT_REQUESTED':'auto',
        'MFM_WGS_HAP_REQUESTED':'auto','MFM_CPU_UTILIZATION':'maximum','MFM_WGS_CPU_BALANCE':'legacy'}

def address_limit():
    return resource.getrlimit(resource.RLIMIT_AS)[0]

class ResourcePredictionTests(unittest.TestCase):
    def config(self):
        return SimpleNamespace(platform='hifi',sort_mem='2G',threads=30,sort_threads=4,
                               hap_jobs=2,input_decode_threads=0)

    def summary(self):
        return dict(contig_lengths={'chr1':249_086_379},total_reads=10_157_288,
                    total_reference_bases=2_876_676_610,read_length=31_888,
                    count_source='정렬 QC의 primary 읽기 수',read_length_samples=256)

    def test_memory_limited_cpu_can_start_instead_of_waiting(self):
        lengths={1:2_876_676_610,2:2_876_520_048}
        p=alignment_plan(self.config(),lengths,68,120*GIB,AUTO)
        self.assertIsNotNone(p)
        self.assertEqual(p['reserved_cpus'],68)
        self.assertLessEqual(p['reserved_memory_bytes'],120*GIB)
        self.assertTrue(p['allocation_diagnostics']['memory_rebalanced_sort'])
        self.assertFalse(p['allocation_diagnostics']['memory_limited_cpu'])
        self.assertGreater(p['threads'],40)
        # RAM이 더 작아 전체 CPU가 불가능할 때만 CPU를 낮춘다.
        smaller=alignment_plan(self.config(),lengths,68,100*GIB,AUTO)
        self.assertIsNotNone(smaller)
        self.assertLess(smaller['reserved_cpus'],68)
        self.assertLessEqual(smaller['reserved_memory_bytes'],100*GIB)
        ample=alignment_plan(self.config(),lengths,68,1485000000000,AUTO)
        self.assertEqual((ample['hap_jobs'],ample['threads'],ample['sort_threads']),(2,17,8))
        self.assertFalse(ample['allocation_diagnostics']['memory_limited_cpu'])

    def test_explicit_threads_are_never_lowered_to_fit_ram(self):
        env=dict(AUTO,MFM_WGS_ALIGN_REQUESTED='40',MFM_WGS_SORT_REQUESTED='19',
                 MFM_WGS_HAP_REQUESTED='1',MFM_MINIMAP_DECODE_THREADS='4')
        self.assertIsNone(alignment_plan(self.config(),{1:2_876_676_610},68,120*GIB,env))

    def test_input_reservation_keeps_address_space_and_unknown_input_fallback(self):
        summary=self.summary()
        with patch.dict(os.environ,{},clear=True):
            p=estimate_peak(summary,'chr1',2,1<<20)
            self.assertGreaterEqual(p['memory_reservation_bytes'],p['estimated_peak_bytes'])
            self.assertLess(p['memory_reservation_bytes'],p['address_space_limit_bytes'])
            self.assertEqual(p['address_space_limit_bytes'],32*GIB)
            unknown=dict(summary,count_source='압축 파일 크기 추정')
            q=estimate_peak(unknown,'chr1',2,1<<20)
            self.assertEqual(q['memory_reservation_bytes'],q['address_space_limit_bytes'])
        with patch.dict(os.environ,MFM_COMPACT_RESERVATION_MODEL='legacy'):
            q=estimate_peak(summary,'chr1',2,1<<20)
            self.assertEqual(q['memory_reservation_bytes'],32*GIB)
        with patch.dict(os.environ,MFM_COMPACT_RESERVATION_MODEL='invalid'):
            with self.assertRaises(ValueError):estimate_peak(summary,'chr1',2,1<<20)

    def test_two_hap_estimates_use_larger_peak_even_when_vms_floors_tie(self):
        first=self.summary();second=dict(first,read_length=50000)
        job=dict(contig='chr1',sample='S',read_coordinate_mode='personal_realigned',
                 hap1_alignment='/tmp/hap1',hap2_alignment='/tmp/hap2')
        with patch.dict(os.environ,{},clear=True), \
             patch('integrate_wgs_rna.compact.memory_preflight.alignment_summary',side_effect=[first,second]), \
             patch('integrate_wgs_rna.compact.memory_preflight.memory_limit',return_value=64*GIB):
            p=forecast_compact([job],2,1<<20)
            expected=estimate_peak(second,'chr1',2,1<<20)
        self.assertEqual(p['estimated_peak_per_worker_bytes'],expected['estimated_peak_bytes'])
        self.assertGreater(p['single_stage_memory_worker_cap'],1)

    def test_memory_worker_cap_replans_decode_threads_within_cpu_and_ram(self):
        from integrate_wgs_rna.compact.resources import fit_memory
        forecast=dict(worker_ceiling_bytes=32*GIB,worker_reservation_bytes=8*GIB,
                      estimated_peak_per_worker_bytes=4*GIB,memory_limit_bytes=128*GIB)
        with patch.dict(os.environ,{'MFM_CPU_BUDGET':'68','MFM_BLOSC_THREADS':'1'},clear=True), \
             patch('os.sched_getaffinity',return_value=set(range(72))), \
             patch('integrate_wgs_rna.compact.resources.cgroup_cpu_limit',return_value=68):
            execution=fit_memory(22,22,forecast,hap_split=True)
        self.assertLess(execution['workers'],17)
        self.assertGreater(execution['alignment_threads'],2)
        self.assertLessEqual(execution['workers']*execution['cpu_per_worker'],68)
        self.assertLessEqual(execution['workers']*forecast['worker_reservation_bytes']+4*GIB,128*GIB*9//10)
        self.assertGreater(forecast['worker_address_space_bytes'],forecast['worker_reservation_bytes'])

    def test_actual_child_gets_vms_limit_without_reserving_that_amount(self):
        with patch.dict(os.environ,{'MFM_RAM_LIMIT_BYTES':str(64*GIB)},clear=True):
            rows=list(cpu_results(address_limit,[()],1,memory_per_worker_bytes=256<<20,
                                  address_space_per_worker_bytes=8*GIB))
        self.assertIsNone(rows[0][2]);self.assertEqual(rows[0][1],8*GIB)
        with patch.dict(os.environ,{},clear=True):
            rows=list(cpu_results(address_limit,[()],1,address_space_per_worker_bytes=8*GIB))
        self.assertIsNone(rows[0][2]);self.assertEqual(rows[0][1],8*GIB)
        with self.assertRaises(ValueError):
            list(cpu_results(address_limit,[()],1,memory_per_worker_bytes=GIB,
                             address_space_per_worker_bytes=GIB-1))

    def test_worker_ledger_reports_memory_limit_and_keeps_running_reservation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'leases').mkdir();(root/'waiters').mkdir()
            own=root/'leases/a.json';own.write_text(json.dumps(dict(cpus=8,total_cpus=8,memory_bytes=4*GIB)))
            other=root/'leases/b.json';other.write_text(json.dumps(dict(cpus=0,memory_bytes=48*GIB)))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(own),'MFM_RAM_LIMIT_BYTES':str(64*GIB)},clear=True):
                budget=WorkerBudget(8,8,True,8*GIB)
                self.assertEqual(budget.capacity(8),0)
                d=json.loads(own.read_text())
                self.assertIn('RAM 예약',d['worker_allocation']['limiting_factors'])
                self.assertEqual(budget.capacity(8,running=2),2)
                self.assertGreaterEqual(json.loads(own.read_text())['memory_bytes'],20*GIB)

    def test_memory_snapshot_distinguishes_current_and_lifetime_peaks(self):
        from integrate_wgs_rna.memory_policy import process_memory_snapshot
        def read(path):
            if str(path).endswith('/status'):
                return 'VmRSS: 10 kB\nVmSize: 100 kB\nVmHWM: 20 kB\nVmPeak: 200 kB\n'
            return 'Pss: 8 kB\n'
        with patch.object(Path,'read_text',autospec=True,side_effect=read):
            p=process_memory_snapshot()
        self.assertEqual((p['rss_bytes'],p['pss_bytes'],p['vms_bytes']),(10240,8192,102400))
        self.assertEqual(p['process_lifetime_peak_vms_bytes'],204800)
        with patch.object(Path,'read_text',side_effect=PermissionError):
            self.assertIsNone(process_memory_snapshot()['rss_bytes'])

    def test_zero_cache_preserves_validation_result_for_both_call_paths(self):
        arrays={k:np.array([v],dtype=d) for k,v,d in [
            ('cigar_offset',0,'i8'),('cigar_len',1,'i4'),('ref_start',0,'i8'),
            ('ref_end',3,'i8'),('lengths',3,'i4'),('segment_type',0,'i1')]}
        cigar=np.array([3<<4],dtype='u4')
        for cls in [ReadValidationChunkCache,_ReadValidationChunkCache]:
            cached=cls(arrays,cigar,1,max_chunks=4);uncached=cls(arrays,cigar,1,max_chunks=0)
            self.assertEqual(cached.decode(0),uncached.decode(0))
            self.assertEqual(len(uncached._cache),0)
            with self.assertRaises(ValueError):cls(arrays,cigar,1,max_chunks=-1)
        with patch.dict(os.environ,{'CACHE':'0'}):self.assertEqual(nonnegative_env_int('CACHE',64),0)
        with patch.dict(os.environ,{'CACHE':'-1'}):
            with self.assertRaises(ValueError):nonnegative_env_int('CACHE',64)

if __name__=='__main__':unittest.main()
