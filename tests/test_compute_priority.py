"""디스크 포화 재발 방지와 RAM·CPU 동시 입장 제어를 검증한다."""
from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch
from integrate_wgs_rna.pool_leases import fitting_waiters
from handoff.queue_resources import download_admit
from handoff.local_queue import Adapter,QueueStopped


class ComputePriorityTests(unittest.TestCase):
    def test_중지요청뒤에는새프로세스를투입하지않음(self):
        adapter=Adapter.__new__(Adapter);adapter.stop_requested=True
        with patch('subprocess.Popen') as spawn:
            with self.assertRaises(QueueStopped):adapter.start('worker','tag',None,{})
            spawn.assert_not_called()

    def test_큰BAM준비는한개만선택하고다른단계는계속투입(self):
        rows=[('a',dict(branch='WGS',stage='01',min_cpus=4,memory_bytes=60)),
              ('b',dict(branch='WGS',stage='01',min_cpus=4,memory_bytes=60)),
              ('c',dict(branch='RNA',stage='16',min_cpus=8,memory_bytes=100))]
        self.assertEqual([p for p,_ in fitting_waiters(rows,68,500,stage_limits={'WGS:01':1})],['a','c'])
        self.assertEqual([p for p,_ in fitting_waiters(rows,68,500,stage_limits={'WGS:01':1},active=[{'branch':'WGS','stage':'01'}])],['c'])

    def test_메모리안맞는앞샘플은stage슬롯을점유하지않음(self):
        rows=[('large',dict(branch='WGS',stage='01',min_cpus=4,memory_bytes=500)),
              ('small',dict(branch='WGS',stage='01',min_cpus=4,memory_bytes=50))]
        self.assertEqual([p for p,_ in fitting_waiters(rows,68,100,stage_limits={'WGS:01':1})],['small'])

    def test_디스크부족과전송동시성은계산중단없이다운로드만대기(self):
        with tempfile.TemporaryDirectory() as raw:
            pool=Path(raw);(pool/'policy.json').write_text(json.dumps({'max_download_jobs':1,'download_min_free_bytes':500}))
            config={'data_root':raw,'memory_limit_bytes':1_650_000_000_000}
            with patch('handoff.queue_resources.live_leases',return_value={'A':[{'branch':'DOWNLOAD','memory_bytes':4_000_000_000}]}):
                ok,reason=download_admit(config,pool);self.assertFalse(ok);self.assertIn('동시성',reason)
            class Stat:
                f_bavail=100
                f_frsize=1
            with patch('handoff.queue_resources.live_leases',return_value={}),patch('os.statvfs',return_value=Stat()):
                ok,reason=download_admit(config,pool);self.assertFalse(ok);self.assertIn('디스크',reason)


if __name__=='__main__':unittest.main()
