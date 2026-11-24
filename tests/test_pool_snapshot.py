"""종료된 다른 임대와의 읽기 경합은 허용하고 본인 예약 누락은 차단함."""
from pathlib import Path
from unittest.mock import patch
import json,os,tempfile,unittest
from integrate_wgs_rna.pool_leases import read_pool_records
from integrate_wgs_rna.realign.phases import RealignPhases


class Tests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        (self.root/'leases').mkdir();(self.root/'waiters').mkdir()
        self.own=self.root/'leases/own.json';self.other=self.root/'leases/other.json'
        self.own.write_text(json.dumps({'cpus':1,'total_cpus':8,'memory_bytes':4*1024**3}))
        self.other.write_text(json.dumps({'cpus':7,'memory_bytes':0}))

    def tearDown(self):self.temp.cleanup()

    def vanish(self):
        original=Path.read_text;other=self.other
        def read(path,*args,**kwargs):
            if path==other:
                path.unlink(missing_ok=True);raise FileNotFoundError(str(path))
            return original(path,*args,**kwargs)
        return patch.object(Path,'read_text',read)

    def test_snapshot_skips_just_returned_other_lease(self):
        with self.vanish():records=list(read_pool_records(self.root/'leases'))
        self.assertEqual([self.own],[path for path,_ in records])

    def test_phase_survives_other_lease_return(self):
        with patch.dict(os.environ,MFM_CPU_POOL_LEASE_FILE=str(self.own)):
            phases=RealignPhases()
            with patch('integrate_wgs_rna.realign.phases.prune_leases'),self.vanish():
                self.assertEqual(2,phases.reserve('test',max_units=2))
        self.assertEqual(2,json.loads(self.own.read_text())['cpus'])

    def test_missing_own_lease_still_fails(self):
        with patch.dict(os.environ,MFM_CPU_POOL_LEASE_FILE=str(self.own)):
            phases=RealignPhases();self.own.unlink()
            with self.assertRaises(FileNotFoundError):phases.reserve('test',max_units=2)

    def test_corrupt_other_record_is_not_silently_ignored(self):
        self.other.write_text('{broken')
        with self.assertRaises(json.JSONDecodeError):list(read_pool_records(self.root/'leases'))


if __name__=='__main__':unittest.main()
