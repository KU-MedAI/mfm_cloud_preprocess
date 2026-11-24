from pathlib import Path
from types import SimpleNamespace
import json
import os
import tempfile
import unittest
from unittest.mock import patch
from integrate_wgs_rna.prepare_phases import prepare_entry,reserve_prepare_write,mark_read_only_prepare
from integrate_wgs_rna.pool_leases import fitting_waiters
from integrate_wgs_rna.memory_policy import GIB


class PreparePhasesTests(unittest.TestCase):
    def fixture(self,raw):
        root=Path(raw);(root/'leases').mkdir();path=root/'leases/a.json'
        policy={'prepare_phased':True,'prepare_min_free_bytes':350_000_000_000,
                'storage_path':raw,'stage_concurrency':{'WGS:01':1}}
        (root/'policy.json').write_text(json.dumps(policy))
        row={'branch':'WGS','stage':'01','cpus':8,'total_cpus':68,'memory_bytes':64*GIB}
        path.write_text(json.dumps(row));source=root/'input.bam';source.write_bytes(b'input')
        contract=root/'contract.json';contract.write_text(json.dumps({'source':{'alignment_states':['aligned'],
            'files':[{'format':'bam','path':str(source)}]}}));return root,path,row,contract

    def test_읽기QC는디스크쓰기슬롯없이CPU_RAM예산으로입장(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,row,contract=self.fixture(raw)
            (root/'leases/writing.json').write_text(json.dumps({'branch':'WGS','stage':'01','cpus':8,'memory_bytes':64*GIB}))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)},clear=True),patch('os.statvfs',return_value=SimpleNamespace(f_bavail=0,f_frsize=1)):
                self.assertEqual(prepare_entry(contract),8)
            current=json.loads(path.read_text());self.assertEqual(current['prepare_phase'],'read_only_qc')
            self.assertEqual(current['memory_bytes'],64*GIB)

    def test_큰BAM쓰기직전에는다시공간확인하고대기예약반환(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,row,contract=self.fixture(raw);free=[0]
            def advance(_):
                paused=json.loads(path.read_text());self.assertTrue(paused['stage_storage_wait']);self.assertEqual(paused['cpus'],0)
                free[0]=400_000_000_000
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path),'MFM_PREPARE_INPUT_CONTRACT':str(contract)},clear=True),patch('os.statvfs',side_effect=lambda _:SimpleNamespace(f_bavail=free[0],f_frsize=1)),patch('integrate_wgs_rna.stage_storage.time.sleep',side_effect=advance):
                mark_read_only_prepare();self.assertEqual(reserve_prepare_write(),8)
                self.assertEqual(os.environ['MFM_CPU_POOL_GRANTED'],'8')
            self.assertEqual(json.loads(path.read_text())['prepare_phase'],'writing_bam')

    def test_raw입력은기존쓰기보호를우회하지않음(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,row,contract=self.fixture(raw);data=json.loads(contract.read_text());data['source']['alignment_states']=['unaligned'];data['source']['files'][0]['format']='fastq';contract.write_text(json.dumps(data));free=[0]
            def advance(_):free[0]=400_000_000_000
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)},clear=True),patch('os.statvfs',side_effect=lambda _:SimpleNamespace(f_bavail=free[0],f_frsize=1)),patch('integrate_wgs_rna.stage_storage.time.sleep',side_effect=advance) as sleep:
                self.assertEqual(prepare_entry(contract),8)
            sleep.assert_called_once();self.assertEqual(json.loads(path.read_text())['prepare_phase'],'writing_bam')

    def test_마지막색인은CPU유지하고계산한RAM만예약(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,row,contract=self.fixture(raw)
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)},clear=True):
                mark_read_only_prepare(index_threads=7)
            current=json.loads(path.read_text())
            self.assertEqual(current['cpus'],8)
            self.assertEqual(current['memory_bytes'],4*GIB+8*64*1024**2)
            self.assertEqual(current['prepare_finish_plan']['prior_memory_bytes'],64*GIB)
            self.assertEqual(current['prepare_phase'],'read_only_qc')

    def test_색인예약은기존보다늘리지않고일반읽기예약유지(self):
        with tempfile.TemporaryDirectory() as raw:
            root,path,row,contract=self.fixture(raw)
            row['memory_bytes']=4*GIB;path.write_text(json.dumps(row))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(path)},clear=True):
                mark_read_only_prepare(index_threads=7)
            self.assertEqual(json.loads(path.read_text())['memory_bytes'],4*GIB)

    def test_읽기입장은쓰기동시성과분리하되CPU_RAM은지킴(self):
        rows=[(str(n),{'branch':'WGS','stage':'01','min_cpus':4,'memory_bytes':64}) for n in range(4)]
        active=[{'branch':'WGS','stage':'01','prepare_phase':'writing_bam'}]
        selected=fitting_waiters(rows,12,192,stage_limits={'WGS:01':1},active=active,prepare_read_only=True)
        self.assertEqual(len(selected),3)
        self.assertFalse(fitting_waiters(rows,12,192,stage_limits={'WGS:01':1},active=active))


if __name__=='__main__':unittest.main()
