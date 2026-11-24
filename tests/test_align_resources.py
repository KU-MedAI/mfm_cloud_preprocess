"""가용 CPU 변화·명시 요청·RAM 제한·원자 예약을 검증함."""
import json
import os
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch
from integrate_wgs_rna.memory_policy import GIB
from integrate_wgs_rna.realign.config import PersonalRealignConfig
from integrate_wgs_rna.realign.align_resources import alignment_plan
from integrate_wgs_rna.realign.phases import RealignPhases, _atomic

AUTO = {'MFM_WGS_ALIGN_REQUESTED': 'auto', 'MFM_WGS_SORT_REQUESTED': 'auto',
        'MFM_WGS_HAP_REQUESTED': 'auto'}


class Plans(unittest.TestCase):
    def config(self, root=Path('/tmp')):
        return PersonalRealignConfig(cohort='test', output_root=root, sample_list=None,
            samples=['S'], chroms=['chr1'], sequence_variant_policy='snv_indel',
            dataset_role='pre_training', threads=26, sort_threads=4, hap_jobs=2)

    def test_실제BWA에서_자동thread변경후에도_read와QC계약동일(self):
        from test_stage22_bwa import Stage22BwaIntegrationTest, module
        measured=[]
        original=module.run_personal_realign_generation
        def observe(config):
            result=original(config)
            p=config.output_root / 'TEST/qc/TEST.stage11_timing.json'
            if p.exists():measured.append(json.loads(p.read_text()))
            return result
        with patch.dict(os.environ, dict(AUTO, MFM_CPU_BUDGET='12')), \
             patch.object(module, 'run_personal_realign_generation', side_effect=observe):
            Stage22BwaIntegrationTest('test_pair_singleton_duplicate_and_chain_contract').debug()
        self.assertEqual(measured[0]['effective_threads']['hap_jobs'],1)
        self.assertEqual(measured[0]['effective_threads']['align_per_hap'],7)
        self.assertEqual(measured[0]['resource_phases'][-1]['reserved_cpus'],1)

    def test_58cpu계획도_12cpu로즉시시작(self):
        p = alignment_plan(self.config(), {1: 3_000_000_000, 2: 3_000_000_000}, 12, 900 * GIB, AUTO)
        self.assertEqual(p['hap_jobs'], 1)
        self.assertEqual(p['reserved_cpus'], 12)
        self.assertEqual(p['threads'], 7)
        self.assertEqual(p['sort_threads'], 2)

    def test_단계시작후_큰여유면_전체가용예산활용(self):
        c = self.config(); c.threads = 2; c.hap_jobs = 1
        p = alignment_plan(c, {1: 3_000_000_000, 2: 3_000_000_000}, 68, 900 * GIB, AUTO)
        self.assertEqual(p['hap_jobs'], 2)
        self.assertEqual(p['reserved_cpus'], 68)
        self.assertGreater(p['threads'], c.threads)

    def test_숫자명시값보존(self):
        env = dict(AUTO, MFM_WGS_ALIGN_REQUESTED='6', MFM_WGS_SORT_REQUESTED='2', MFM_WGS_HAP_REQUESTED='2')
        p = alignment_plan(self.config(), {1: 3_000_000_000, 2: 3_000_000_000}, 68, 900 * GIB, env)
        self.assertEqual((p['threads'], p['sort_threads'], p['hap_jobs']), (6, 2, 2))
        self.assertEqual(p['reserved_cpus'], 22)
        self.assertIsNone(alignment_plan(self.config(), {1: 3_000_000_000, 2: 3_000_000_000}, 12, 900 * GIB, env))

    def test_자동thread비활성화존중(self):
        self.assertIsNone(alignment_plan(self.config(), {1: 3_000_000_000, 2: 3_000_000_000}, 12, 900 * GIB,
                                        dict(AUTO, WGS_AUTO_THREADS='0')))

    def test_RAM부족이면_HAP1개로재선택(self):
        p = alignment_plan(self.config(), {1: 4_000_000_000, 2: 4_000_000_000}, 68, 211_000_000_000, AUTO)
        self.assertEqual(p['hap_jobs'], 1)
        self.assertLessEqual(p['reserved_memory_bytes'], 211_000_000_000)
        smaller = alignment_plan(self.config(), {1: 4_000_000_000}, 68, 100 * GIB, AUTO)
        self.assertIsNotNone(smaller)
        self.assertLess(smaller['reserved_cpus'], 68)
        self.assertLessEqual(smaller['reserved_memory_bytes'], 100 * GIB)
        self.assertIsNone(alignment_plan(self.config(), {1: 4_000_000_000}, 68, 90 * GIB, AUTO))

    def test_완료HAP은_동시자원에서제외(self):
        p = alignment_plan(self.config(), {2: 3_000_000_000}, 68, 900 * GIB, AUTO)
        self.assertEqual(p['hap_jobs'], 1)

    def test_최소4CPU미만일때_명시대기(self):
        self.assertIsNone(alignment_plan(self.config(), {1: 1000}, 3, 900 * GIB, AUTO))

    def test_HAP병렬거부존중(self):
        p = alignment_plan(self.config(), {1: 1000, 2: 1000}, 68, 900 * GIB,
                           dict(AUTO, MFM_HAP_PARALLEL_ALLOWED='0'))
        self.assertEqual(p['hap_jobs'], 1)

    def test_작은실제lease에서_기존큰요청축소(self):
        with tempfile.TemporaryDirectory() as directory:
            r = Path(directory);(r / 'leases').mkdir();(r / 'waiters').mkdir()
            owner = r / 'leases/owner.json'
            _atomic(owner, {'cpus': 1, 'total_cpus': 68, 'memory_bytes': 4 * GIB, 'manager_pid': os.getpid()})
            _atomic(r / 'leases/other.json', {'cpus': 56, 'memory_bytes': 500 * GIB})
            with patch.dict(os.environ, dict(AUTO, MFM_CPU_POOL_LEASE_FILE=str(owner)), clear=True):
                phases = RealignPhases()
                result = phases.align_config(self.config(r), {1: 3_000_000_000, 2: 3_000_000_000})
            self.assertEqual((result.hap_jobs, result.threads), (1, 7))
            row = json.loads(owner.read_text())
            self.assertEqual(row['cpus'], 12)
            self.assertLessEqual(row['memory_bytes'] + 500 * GIB, 900_000_000_000)

    def test_다른예약반환후_대기작업재계산(self):
        with tempfile.TemporaryDirectory() as directory:
            r = Path(directory);(r / 'leases').mkdir();(r / 'waiters').mkdir()
            owner=r / 'leases/owner.json';other=r / 'leases/other.json'
            _atomic(owner, {'cpus': 1, 'total_cpus': 68, 'memory_bytes': 4 * GIB, 'manager_pid': os.getpid()})
            _atomic(other, {'cpus': 65, 'memory_bytes': 100 * GIB})
            results=[];errors=[]
            with patch.dict(os.environ, dict(AUTO, MFM_CPU_POOL_LEASE_FILE=str(owner)), clear=True):
                phase=RealignPhases()
                def run():
                    try:results.append(phase.align_config(self.config(r), {1: 1000, 2: 1000}))
                    except Exception as e:errors.append(e)
                thread=threading.Thread(target=run);thread.start()
                try:
                    end=time.monotonic()+2
                    while time.monotonic()<end:
                        if json.loads(owner.read_text()).get('realign_phase','').endswith(':waiting'):break
                        time.sleep(.01)
                    self.assertTrue(thread.is_alive())
                    self.assertEqual(json.loads(owner.read_text())['memory_bytes'],4*GIB)
                finally:
                    other.unlink();thread.join(3)
                self.assertFalse(thread.is_alive());self.assertFalse(errors)
                self.assertEqual(json.loads(owner.read_text())['cpus'],68)
                self.assertEqual(results[0].hap_jobs,2)

if __name__ == '__main__':
    unittest.main()
