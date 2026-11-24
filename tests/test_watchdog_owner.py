"""pool.lock의원자소유권으로고아감시기락·중복/PID재사용을검증한다."""
from pathlib import Path
from contextlib import ExitStack
from unittest.mock import patch
import fcntl
import importlib.util
import json
import os
import socket
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]
spec=importlib.util.spec_from_file_location('watchdog_owner_test',ROOT/'scripts/resource_watchdog.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class OwnerTests(unittest.TestCase):
    def fixture(self,root,prior):
        for name in ['leases','waiters']:(root/name).mkdir()
        (root/'watchdog.json').write_text(json.dumps(prior))

    def tick(self,root,*,exclusive=True):
        with ExitStack() as stack:
            for name,value in [('pool_cpu_limit',68),('refresh_storage_capacity',None),
                               ('allocate',None),('protect_pool_memory',{'rss_bytes':0,'action':'continue'}),
                               ('memory_limit',1650000000000),('usable_memory',1485000000000),
                               ('lend_idle_cpus',{'ready':False})]:
                stack.enter_context(patch.object(module,name,return_value=value))
            stack.enter_context(patch.object(module.OBSERVER,'observe',return_value={'ready':False}))
            return module.tick(root,68,os.getpid(),exclusive=exclusive)

    def test_dead_owner_reclaims_even_with_unavailable_legacy_lock(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);self.fixture(root,{'pid':12345678})
            with (root/'watchdog.lock').open('a+') as old_lock:
                fcntl.flock(old_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                state=self.tick(root)
            self.assertEqual(state['watchdog_owner']['pid'],os.getpid())
            self.assertEqual(state['watchdog_owner']['started'],module.process_start(os.getpid()))

    def test_live_owner_blocks_duplicate_before_allocation(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);self.fixture(root,{'watchdog_owner':{'host':socket.gethostname(),'pid':12345678,'started':99}})
            with patch.object(module,'process_alive',return_value=True),patch.object(module,'process_start',return_value=99),patch.object(module,'allocate') as allocate:
                self.assertIsNone(module.tick(root,68,0,exclusive=True))
                allocate.assert_not_called()

    def test_reused_pid_does_not_block_new_owner(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);self.fixture(root,{'watchdog_owner':{'host':socket.gethostname(),'pid':12345678,'started':1}})
            actual_start=module.process_start
            with patch.object(module,'process_alive',return_value=True),patch.object(module,'process_start',
                    side_effect=lambda pid:2 if pid==12345678 else actual_start(pid)):
                self.assertIsNotNone(self.tick(root))

    def test_other_host_owner_is_preserved_even_if_stale(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);self.fixture(root,{'watchdog_owner':{'host':'another-host','pid':12345678,'started':1},'updated':0})
            self.assertIsNone(self.tick(root))

    def test_manual_tick_preserves_daemon_identity(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);owner={'host':socket.gethostname(),'pid':12345678,'started':1}
            self.fixture(root,{'watchdog_owner':owner})
            state=self.tick(root,exclusive=False)
            self.assertEqual(state['watchdog_owner'],owner)


if __name__=='__main__':unittest.main()
