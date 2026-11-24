from types import SimpleNamespace
import unittest
from unittest.mock import patch
from integrate_wgs_rna.storage_admission import refresh_storage_capacity
from integrate_wgs_rna.pool_leases import fitting_waiters


class StorageAdmissionTests(unittest.TestCase):
    def test_공간증가와감소에따라과거한개상한을자동갱신(self):
        policy={'stage_concurrency':{'WGS:01':1},'adaptive_prepare_storage':{
            'enabled':True,'path':'/tmp','peak_bytes_per_job':350_000_000_000,'reserve_bytes':150_000_000_000}}
        for free,slots in [(1_840_000_000_000,4),(400_000_000_000,0),(2_250_000_000_000,6)]:
            with patch('os.statvfs',return_value=SimpleNamespace(f_bavail=free,f_frsize=1)):
                result=refresh_storage_capacity(policy,68)
            self.assertEqual(result['slots'],slots)
            self.assertEqual(policy['stage_concurrency']['WGS:01'],slots)

    def test_공간반환대기와일시정지는활성슬롯에포함하지않음(self):
        waiting=[('new',{'branch':'WGS','stage':'01','min_cpus':4,'memory_bytes':64})]
        active=[{'branch':'WGS','stage':'01','cpus':0,'stage_storage_wait':True},
                {'branch':'WGS','stage':'01','storage_paused':True}]
        self.assertEqual(len(fitting_waiters(waiting,8,128,stage_limits={'WGS:01':1},active=active)),1)

    def test_정책변경이없으면매주기파일갱신하지않음(self):
        policy={'stage_concurrency':{'WGS:01':1},'adaptive_prepare_storage':{'enabled':True,'path':'/tmp'}}
        with patch('os.statvfs',return_value=SimpleNamespace(f_bavail=600_000_000_000,f_frsize=1)):
            self.assertTrue(refresh_storage_capacity(policy,68)['policy_changed'])
            self.assertFalse(refresh_storage_capacity(policy,68)['policy_changed'])


if __name__=='__main__':unittest.main()
