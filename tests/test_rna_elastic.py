"""RNA 입장 상한·작업경계 재확보·부모 및 실행 중 CPU/RAM 예약을 확인함."""
from pathlib import Path
from unittest.mock import patch
import json,os,tempfile,unittest,time,threading,hashlib
from integrate_wgs_rna.cpu_executor import WorkerBudget,cpu_results
from integrate_wgs_rna.rna_wgs_native import cpu_upper,plan
GIB=1024**3


def burn_task(value):
    begin=time.monotonic();cpu=time.process_time();data=b'ACGT'*1024
    while time.process_time()-cpu<.08:data=hashlib.sha256(data).digest()
    return value,begin,time.monotonic(),time.process_time()-cpu


class ElasticRnaTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);(self.root/'leases').mkdir();(self.root/'waiters').mkdir()
        self.path=self.root/'leases/own.json'
        self.row=dict(branch='RNA',stage='20',cpus=17,total_cpus=68,memory_bytes=96*GIB)
        self.write()
        self.env=dict(MFM_CPU_POOL_LEASE_FILE=str(self.path),MFM_CPU_BUDGET='68',MFM_SAMPLE_CPU_BUDGET='17',
            MFM_CPU_POOL_STAGE='20',MFM_CPU_POOL_STAGE_MAX='17',MFM_RNA_DAG_DYNAMIC_CPU_MAX='68',MFM_RNA_DAG_ROLE='variants',
            MFM_CPU_POOL_LAYOUT='original',MFM_DYNAMIC_CPU_POOL_ACTIVE='1',MFM_RAM_LIMIT_BYTES=str(1650000000000))
    def write(self):self.path.write_text(json.dumps(self.row))

    def test_본인RNAlease에서만상한을재확보함(self):
        with patch.dict(os.environ,self.env,clear=True):
            self.assertEqual(cpu_upper(17),68)
            p=plan(22);self.assertEqual(p['workers_initial'],5);self.assertEqual(p['workers_upper'],22)
            self.assertEqual(p['parent_cpus'],1)
            self.assertLessEqual(p['workers_upper']*p['cpu_per_worker']+p['parent_cpus'],68)
            for key,value in [('branch','WGS'),('stage','16'),('cpu_backfill',{'cpus':17})]:
                old=dict(self.row);self.row[key]=value;self.write();self.assertEqual(cpu_upper(17),17)
                self.row=old;self.write()
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':''}):self.assertEqual(cpu_upper(17),17)
            with patch.dict(os.environ,{'MFM_RNA_DAG_DYNAMIC_CPU_MAX':'8'}):self.assertEqual(cpu_upper(17),8)
            with patch.dict(os.environ,{'MFM_CPU_UTILIZATION':'balanced'}):self.assertEqual(cpu_upper(17),17)
            with patch.dict(os.environ,{'MFM_CPU_POOL_PROFILE':'balanced'}):self.assertEqual(cpu_upper(17),17)

    def test_전달번호및명시CPU제한을지킴(self):
        self.row['stage']='18';self.write()
        with patch.dict(os.environ,{**self.env,'MFM_CPU_POOL_LAYOUT':'legacy','MFM_CPU_POOL_STAGE':'18'},clear=True):
            self.assertEqual(cpu_upper(17),68)
            with patch.dict(os.environ,{'MFM_CPU_BUDGET':'12'}):self.assertEqual(cpu_upper(17),12)
            self.assertEqual(plan(22,cpu_request=4)['cpu_upper'],68)
            with patch.dict(os.environ,{'MFM_DYNAMIC_CPU_POOL_ACTIVE':'0'}):self.assertEqual(plan(22,cpu_request=4)['cpu_upper'],4)

    def test_입장17에서남은공유풀을재확보하되대기자를우선함(self):
        other=self.root/'leases/peer.json';other.write_text(json.dumps({'cpus':51,'memory_bytes':0}))
        with patch.dict(os.environ,self.env,clear=True):
            b=WorkerBudget(5,22,True,4*GIB,cpus_per_worker=3,parent_cpus=1)
            self.assertEqual(b.capacity(22),5)
            self.assertEqual(json.loads(self.path.read_text())['cpus'],16)
            other.unlink()
            waiter=self.root/'waiters/next.json';waiter.write_text(json.dumps({'min_cpus':8,'memory_bytes':0}))
            self.assertEqual(b.capacity(22),5)
            waiter.unlink();self.assertEqual(b.capacity(22),22)
            self.assertEqual(json.loads(self.path.read_text())['cpus'],67)
            b.capacity(0,grow=False);self.assertEqual(json.loads(self.path.read_text())['cpus'],1)

    def test_실제프로세스가반환된CPU를쓰고부모1CPU를보존함(self):
        self.row.update(cpus=3,total_cpus=7);self.write()
        peer=self.root/'leases/peer.json';peer.write_text(json.dumps({'cpus':4,'memory_bytes':0}))
        env={**self.env,'MFM_CPU_BUDGET':'7','MFM_SAMPLE_CPU_BUDGET':'3'}
        timer=threading.Timer(.12,peer.unlink)
        with patch.dict(os.environ,env,clear=True):
            timer.start()
            try:rows=list(cpu_results(burn_task,[(i,) for i in range(12)],2,elastic=True,max_workers=6,parent_cpus=1))
            finally:timer.join()
        self.assertTrue(all(error is None for _,_,error in rows))
        values=[v for _,v,_ in rows];self.assertEqual(sorted(v[0] for v in values),list(range(12)))
        peak=max(sum(a[1]<=b[1]<a[2] for a in values) for b in values)
        self.assertGreaterEqual(peak,3);self.assertLessEqual(peak,6)
        row=json.loads(self.path.read_text());self.assertLessEqual(row['peak_cpus'],7);self.assertEqual(row['cpus'],1)
        self.assertGreater(sum(v[3] for v in values),.5)

    def test_RAM부족으로현재worker예약을조기에반환하지않음(self):
        self.row.update(cpus=16,memory_bytes=24*GIB);self.write()
        other=self.root/'leases/peer.json';other.write_text(json.dumps({'cpus':0,'memory_bytes':63*GIB}))
        with patch.dict(os.environ,{**self.env,'MFM_RAM_LIMIT_BYTES':str(64*GIB)},clear=True):
            b=WorkerBudget(5,22,True,4*GIB,cpus_per_worker=3,parent_cpus=1)
            self.assertEqual(b.capacity(12,running=5),5)
            row=json.loads(self.path.read_text());self.assertEqual(row['cpus'],16);self.assertGreaterEqual(row['memory_bytes'],24*GIB)


if __name__=='__main__':unittest.main()
