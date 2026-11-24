"""전송 병렬수의 명시 설정 우선순위와 RAM 예약 일치를 검증함."""
import json,os,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from handoff.download_budget import memory_plan

class RangePolicyTests(unittest.TestCase):
    def test_policy_and_explicit_override_use_matching_memory(self):
        with tempfile.TemporaryDirectory() as raw, patch.dict(os.environ,{},clear=True):
            root=Path(raw);(root/'policy.json').write_text(json.dumps({'download_range_workers':4}))
            c={'pipeline_env':{'MFM_CPU_POOL_DIR':str(root)}}
            plan=memory_plan(c)
            self.assertEqual(plan['range_workers'],4)
            self.assertEqual(plan['reservation_bytes'],((512<<20)+4*3*(32<<20))*12//10)
            c['pipeline_env']['MFM_HTTP_RANGE_WORKERS']='2'
            self.assertEqual(memory_plan(c)['range_workers'],2)
            c['http_range_workers']=8
            self.assertEqual(memory_plan(c)['range_workers'],8)
            c['http_range_workers']=True
            with self.assertRaises(ValueError):memory_plan(c)
    def test_invalid_policy_is_rejected_and_default_stays_one(self):
        with tempfile.TemporaryDirectory() as raw, patch.dict(os.environ,{},clear=True):
            self.assertEqual(memory_plan({})['range_workers'],1)
            root=Path(raw);(root/'policy.json').write_text(json.dumps({'download_range_workers':17}))
            with self.assertRaises(ValueError):memory_plan({'pipeline_env':{'MFM_CPU_POOL_DIR':str(root)}})
