"""읽기큐의잠금대기제거와기존RAM예약판정을확인한다."""
from pathlib import Path
from unittest.mock import patch
import fcntl
import json
import os
import sys
import tempfile
import threading
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from handoff import queue_resources as module
GIB=1024**3


class SnapshotTests(unittest.TestCase):
    def fixture(self,root):
        (root/'leases').mkdir()
        (root/'leases/one.json').write_text(json.dumps({'sample':'S','branch':'WGS','stage':'10',
            'cpus':4,'memory_bytes':48*GIB,'manager_pid':os.getpid()}))

    def test_read_snapshot_returns_while_allocator_lock_is_held(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);self.fixture(root);results=[]
            with (root/'pool.lock').open('a+') as handle:
                fcntl.flock(handle,fcntl.LOCK_EX)
                thread=threading.Thread(target=lambda:results.append(module.live_leases(root)))
                thread.start()
                try:
                    thread.join(.5)
                    self.assertFalse(thread.is_alive(),'읽기스냅샷이할당기잠금에대기함')
                finally:
                    fcntl.flock(handle,fcntl.LOCK_UN)
                    thread.join(2)
            self.assertEqual(results[0]['S'][0]['memory_bytes'],48*GIB)

    def test_corrupt_and_dead_records_keep_existing_skip_contract(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);self.fixture(root)
            (root/'leases/bad.json').write_text('{broken')
            (root/'leases/dead.json').write_text(json.dumps({'sample':'D','manager_pid':12345678}))
            self.assertEqual(set(module.live_leases(root)),{'S'})

    def test_download_still_refuses_shared_ram_over_budget(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);self.fixture(root)
            with patch.object(module,'cgroup_limits',return_value=(0,100*GIB,68)),patch.object(module,'effective_memory',return_value=50*GIB):
                allowed,reason=module.download_admit({'data_root':str(root)},root)
            self.assertFalse(allowed)
            self.assertIn('공유 RAM 예약',reason)


if __name__=='__main__':unittest.main()
