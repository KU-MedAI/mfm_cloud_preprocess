"""종료 PID·zombie·생존 분리 그룹의 예약 회수를 검증함."""
import json
import os
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch
from integrate_wgs_rna.pool_leases import process_alive, lease_alive, prune_leases


class Leases(unittest.TestCase):
    def test_수동예약은임의삭제안함(self):
        self.assertTrue(lease_alive({'cpus': 4}))

    def test_분리된자식그룹살아있으면보존(self):
        with patch('integrate_wgs_rna.pool_leases.process_alive', return_value=False), \
             patch('integrate_wgs_rna.pool_leases.group_alive', return_value=True):
            self.assertTrue(lease_alive({'manager_pid':123,'child_pid':124}))

    def test_모두종료된예약만회수(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'leases').mkdir()
            dead=root/'leases/dead.json';manual=root/'leases/manual.json'
            dead.write_text(json.dumps({'manager_pid':123,'child_pid':124,'cpus':1}))
            manual.write_text(json.dumps({'cpus':4}))
            with patch('integrate_wgs_rna.pool_leases.process_alive', return_value=False), \
                 patch('integrate_wgs_rna.pool_leases.group_alive', return_value=False):
                self.assertEqual(prune_leases(root), ['dead.json'])
            self.assertFalse(dead.exists());self.assertTrue(manual.exists())

    def test_실제zombie_PID는생존계산으로보지않음(self):
        child=os.fork()
        if child==0:os._exit(0)
        try:
            deadline=time.monotonic()+2
            while time.monotonic()<deadline:
                raw=Path(f'/proc/{child}/stat').read_text()
                if raw[raw.rfind(')')+2:].split()[0]=='Z':break
                time.sleep(.01)
            self.assertFalse(process_alive(child))
        finally:os.waitpid(child,0)
        self.assertTrue(process_alive(os.getpid()))

if __name__ == '__main__':
    unittest.main()
