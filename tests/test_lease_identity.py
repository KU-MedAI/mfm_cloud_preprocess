"""별도 GPU 비교가 CPU 재개를 막거나 RAM을 중복 집계하지 않는지 검증함."""
from pathlib import Path
from unittest.mock import patch
import unittest
from handoff.lease_identity import lease_matches_config
from handoff.local_queue import Adapter
from handoff.queue_resources import active_ram,GIB


class Tests(unittest.TestCase):
    config={'output_root':'/tmp/production','output_layout':'run_root','run_id':'main'}

    def test_foreign_output_is_not_our_detached_stage(self):
        lease={'output_root':'/tmp/gpu_compare','run_id':'gpu'}
        self.assertFalse(lease_matches_config(lease,self.config))

    def test_imported_old_run_in_same_output_is_ours(self):
        self.assertTrue(lease_matches_config({'output_root':'/tmp/production','run_id':'old_batch'},self.config))

    def test_legacy_resolved_output_subdirectories_are_distinct(self):
        self.assertFalse(lease_matches_config({'output_root':'/tmp/production/other','run_id':'other'},dict(self.config,output_root='/tmp/production/main',output_layout='run_id')))

    def test_same_resolved_output_is_protected_even_with_different_run_id(self):
        self.assertTrue(lease_matches_config({'output_root':'/tmp/production','run_id':'other'},dict(self.config,output_layout='run_id')))

    def test_unknown_old_owner_is_preserved_conservatively(self):
        self.assertTrue(lease_matches_config({},self.config))

    def test_child_environment_is_used_when_manager_is_gone(self):
        with patch('handoff.lease_identity._process_environment',side_effect=[{}, {'OUTPUT_ROOT':'/tmp/gpu_compare','RUN_ID':'gpu'}]):
            self.assertFalse(lease_matches_config({'manager_pid':1,'child_pid':2},self.config))

    def test_download_lease_is_not_a_compute_stage(self):
        adapter=Adapter.__new__(Adapter);adapter.config=self.config;adapter.SHARED_CPU_POOL=Path('/tmp/fake');adapter.sample_state=lambda sample:Path('/tmp/fake')
        with patch('handoff.local_queue.locked',return_value=False),patch('handoff.local_queue.live_leases',return_value={'S':[{'branch':'DOWNLOAD'}]}):
            self.assertFalse(adapter.detached_running({'sample_id':'S'}))

    def test_gpu_lease_does_not_block_cpu_resume(self):
        adapter=Adapter.__new__(Adapter);adapter.config=self.config;adapter.SHARED_CPU_POOL=Path('/tmp/fake');adapter.sample_state=lambda sample:Path('/tmp/fake')
        leases={'S':[{'branch':'WGS','output_root':'/tmp/gpu_compare','memory_bytes':40*GIB}]}
        with patch('handoff.local_queue.locked',return_value=False),patch('handoff.local_queue.live_leases',return_value=leases):
            self.assertFalse(adapter.detached_running({'sample_id':'S'}))

    def test_cpu_and_gpu_ram_are_not_attributed_twice_to_cpu_sample(self):
        leases={'S':[{'branch':'WGS','output_root':'/tmp/production','memory_bytes':10*GIB},
                     {'branch':'WGS','output_root':'/tmp/gpu_compare','memory_bytes':40*GIB}]}
        row={'sample_id':'S','platform':'hifi'}
        with patch('handoff.queue_resources.live_leases',return_value=leases),patch('handoff.queue_resources.full_ram',return_value=80*GIB),patch('handoff.queue_resources.starting_memory',return_value=64*GIB):
            self.assertEqual(26*GIB,active_ram('/tmp/fake',[row],config=self.config))


if __name__=='__main__':unittest.main()
