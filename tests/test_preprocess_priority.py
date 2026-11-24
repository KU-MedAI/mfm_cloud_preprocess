"""전처리 일반우선순위와 준비입장/실제정렬 예약 분리를 검사한다."""
from pathlib import Path
import importlib.util
import json
import os
import tempfile
import unittest
from unittest.mock import patch
import test_cpu_backfill as fixtures
from integrate_wgs_rna import pool_backfill as module

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('priority_cpu_pool',ROOT/'scripts/cpu_pool.py')
POOL=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(POOL)


class PreprocessPriorityTests(unittest.TestCase):
    def test_명시정책이면새전처리를nice0으로실행(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=fixtures.BackfillTests().fixture(raw)
            (root/'policy.json').write_text(json.dumps({'cpu_backfill':{'enabled':True,'nice':0}}))
            with patch.object(module,'eligible_waiter',return_value=True),patch.object(module,'in_cpu_group',return_value=True),patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=1000),patch.object(os,'setpriority') as priority,patch.object(os,'getpriority',return_value=0):
                result=module.lend_idle_cpus(root,8,obs,observer,guard)
                self.assertEqual(result['launched'][0]['nice'],0)
                priority.assert_called_with(os.PRIO_PROCESS,os.getpid(),0)
            self.assertEqual(json.loads((root/'leases/a.json').read_text())['cpu_backfill']['nice'],0)

    def test_기존nice19는무리하게상향하지않고계산진행(self):
        with tempfile.TemporaryDirectory() as raw:
            root,obs,observer,guard=fixtures.BackfillTests().fixture(raw)
            (root/'policy.json').write_text(json.dumps({'cpu_backfill':{'enabled':True,'nice':0}}))
            with patch.object(module,'eligible_waiter',return_value=True),patch.object(module,'in_cpu_group',return_value=True),patch('integrate_wgs_rna.memory_policy.usable_memory',return_value=1000),patch.object(os,'setpriority'),patch.object(os,'getpriority',return_value=19):
                self.assertEqual(module.lend_idle_cpus(root,8,obs,observer,guard)['launched'][0]['nice'],19)

    def test_동적HiFi준비만2CPU이며실제정렬예산은별도(self):
        plan={'adaptive_realign':True,'stage_estimates_bytes':{'personal_realign_entry':4*1024**3}}
        env={'MFM_CPU_UTILIZATION':'maximum','MFM_CPU_POOL_PROFILE':'throughput','MFM_CPU_POOL_LAYOUT':'original','MFM_LONGREAD_MEMORY_PLAN':json.dumps(plan),'MFM_REALIGN_PHASE_RESOURCES':'auto'}
        with patch.dict(os.environ,env):
            self.assertEqual(POOL.stage_cap('WGS','10',68),(2,68))
            with patch.dict(os.environ,{'MFM_REALIGN_PHASE_RESOURCES':'0'}):
                self.assertEqual(POOL.stage_cap('WGS','10',68),(8,68))
            with patch.dict(os.environ,{'MFM_LONGREAD_MEMORY_PLAN':'{}'}):
                self.assertEqual(POOL.stage_cap('WGS','10',68),(8,68))
            with patch.dict(os.environ,{'MFM_LONGREAD_MEMORY_PLAN':''}):
                self.assertEqual(POOL.stage_cap('WGS','10',68),(8,68))


if __name__=='__main__':unittest.main()
