from pathlib import Path
import json
import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from integrate_wgs_rna.stage_storage import wait_for_storage
from integrate_wgs_rna.memory_policy import GIB


class StageStorageTests(unittest.TestCase):
    def test_대기중변경된동시성정책을다시읽음(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,row=self.fixture(raw);policy=json.loads((root/'policy.json').read_text());policy['stage_concurrency']['WGS:01']=0;(root/'policy.json').write_text(json.dumps(policy))
            def advance(_):
                policy['stage_concurrency']['WGS:01']=2;(root/'policy.json').write_text(json.dumps(policy))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)},clear=True),patch('os.statvfs',return_value=SimpleNamespace(f_bavail=400_000_000_000,f_frsize=1)),patch('integrate_wgs_rna.stage_storage.time.sleep',side_effect=advance) as sleep:
                self.assertEqual(wait_for_storage(),8)
            sleep.assert_called_once()

    def test_기존8CPU전부대신가용6CPU로계산을재개함(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,row=self.fixture(raw);(root/'leases/other.json').write_text(json.dumps({'branch':'RNA','stage':'16','cpus':62,'memory_bytes':64*GIB}));free=[0]
            def advance(_):free[0]=400_000_000_000
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)},clear=True),patch('os.statvfs',side_effect=lambda _:SimpleNamespace(f_bavail=free[0],f_frsize=1)),patch('integrate_wgs_rna.stage_storage.time.sleep',side_effect=advance):
                self.assertEqual(wait_for_storage(),6)
            self.assertEqual(json.loads(path.read_text())['cpus'],6)

    def fixture(self,raw):
        root=Path(raw);(root/'leases').mkdir();path=root/'leases/a.json'
        (root/'policy.json').write_text(json.dumps({'prepare_min_free_bytes':350_000_000_000,'storage_path':raw,'stage_concurrency':{'WGS:01':1}}))
        row={'sample':'S','branch':'WGS','stage':'01','cpus':8,'total_cpus':68,'memory_bytes':64*GIB}
        path.write_text(json.dumps(row));return root,path,row

    def test_공간대기중CPU와RAM반환후원래예약복구(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,row=self.fixture(raw);free=[0]
            def advance(_):
                paused=json.loads(path.read_text());self.assertEqual(paused['cpus'],0);self.assertEqual(paused['memory_bytes'],4*GIB)
                free[0]=400_000_000_000
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)},clear=True),patch('os.statvfs',side_effect=lambda _:SimpleNamespace(f_bavail=free[0],f_frsize=1)),patch('integrate_wgs_rna.stage_storage.time.sleep',side_effect=advance):
                wait_for_storage()
            current=json.loads(path.read_text());self.assertEqual(current['cpus'],8);self.assertEqual(current['memory_bytes'],64*GIB);self.assertFalse(current['stage_storage_wait'])

    def test_공간이충분하면예약변경없음(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,row=self.fixture(raw)
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)},clear=True),patch('os.statvfs',return_value=SimpleNamespace(f_bavail=400_000_000_000,f_frsize=1)):
                wait_for_storage()
            current=json.loads(path.read_text());self.assertEqual({key:current[key] for key in row},row)
            self.assertEqual(current['prepare_phase'],'writing_bam')

    def test_실제입력두배와50GB도확보(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,row=self.fixture(raw);source=root/'input.bam'
            with source.open('wb') as handle:handle.truncate(200_000_000_000)
            contract=root/'contract.json';contract.write_text(json.dumps({'source':{'files':[{'path':str(source)}]}}));free=[400_000_000_000]
            def advance(_):free[0]=500_000_000_000
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)},clear=True),patch('os.statvfs',side_effect=lambda _:SimpleNamespace(f_bavail=free[0],f_frsize=1)),patch('integrate_wgs_rna.stage_storage.time.sleep',side_effect=advance) as sleep:
                wait_for_storage(contract)
            sleep.assert_called_once()


if __name__=='__main__':unittest.main()
