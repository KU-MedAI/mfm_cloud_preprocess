"""FASTQ backend 선택·실제 read/mate/제외수·cache 계약을 확인함."""
from pathlib import Path
from unittest.mock import patch
import collections,json,os,shutil,sys,tempfile,unittest
import test_rna_cache_contract as cache
from test_parabricks_bam import ParabricksBamTest,PB


class FastqBackendTest(unittest.TestCase):
    def test_기존자동선택과명시선택및잘못된값(self):
        stage=cache.module('11')
        for bam,requested,expected in [('parabricks','auto','samtools'),('cpu','auto','samtools'),
                                       ('parabricks','samtools','samtools'),('cpu','parabricks','parabricks')]:
            with patch.dict(os.environ,{'MFM_BAM_BACKEND':bam,'RNA_FASTQ_BACKEND':requested}):self.assertEqual(stage.fastq_backend(),expected)
        with patch.dict(os.environ,{'RNA_FASTQ_BACKEND':'wrong'}),self.assertRaises(ValueError):stage.fastq_backend()

    def test_samtools부모_해독_IOthread합이grant이하(self):
        stage=cache.module('11')
        for grant in [2,3,4,5,8,17,68]:
            with patch.dict(os.environ,{'MFM_SAMPLE_CPU_BUDGET':str(grant),'MFM_CPU_BUDGET':'68'}):
                plan=stage.samtools_cpu_plan(100,100)
                self.assertLessEqual(plan['cpu_upper'],grant)
        with patch.dict(os.environ,{'MFM_SAMPLE_CPU_BUDGET':'1'}),self.assertRaises(ValueError):stage.samtools_cpu_plan(8,4)

    def test_실행전고정thread를넘는정상예약만반환함(self):
        from integrate_wgs_rna.rna_fastq_native import reserve_conversion
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir();lease=root/'leases/own.json'
            row=dict(branch='RNA',stage='13',cpus=68,total_cpus=68,memory_bytes=64<<30)
            env={'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_CPU_POOL_LAYOUT':'original','MFM_RAM_LIMIT_BYTES':str(128<<30)}
            plan={'cpu_upper':16,'predicted_peak_ram_bytes':6<<30}
            for branch,loan in [('RNA',False),('BENCH',False),('RNA',True)]:
                original=dict(row,branch=branch)
                if loan:original['cpu_backfill']={'cpus':68}
                lease.write_text(json.dumps(original))
                with patch.dict(os.environ,env):result=reserve_conversion(plan)
                actual=json.loads(lease.read_text())
                if branch=='RNA' and not loan:self.assertEqual(actual['cpus'],16);self.assertEqual(actual['memory_bytes'],6<<30)
                else:self.assertIsNone(result);self.assertEqual(actual,original)

    def test_변환완료후native검사는사용하는CPU와RAM만보유함(self):
        from integrate_wgs_rna import rna_fastq_native as native
        with tempfile.TemporaryDirectory() as raw:
            folder=Path(raw);(folder/'leases').mkdir();(folder/'waiters').mkdir();lease=folder/'leases/self.json'
            row={'cpus':20,'total_cpus':68,'branch':'RNA','stage':'13','memory_bytes':64<<30};lease.write_text(json.dumps(row))
            observed=[]
            def process(*_):
                observed.append(json.loads(lease.read_text())['cpus'])
                return {'records':1,'max_read_length':2}
            env={'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_CPU_POOL_LAYOUT':'original','MFM_RAM_LIMIT_BYTES':str(128<<30)}
            with patch.dict(os.environ,env),patch.object(native,'process',side_effect=process):
                values=native.batch([(folder/'a',None,None),(folder/'b',None,None)],20)
            self.assertEqual(observed,[3,3]);self.assertEqual(len(values),2)
            after=json.loads(lease.read_text());self.assertEqual(after['cpus'],1);self.assertEqual(after['memory_bytes'],4<<30)

    @unittest.skipUnless(Path(PB).is_file() and shutil.which('samtools'),'실제 Parabricks/samtools 필요함')
    def test_실제두backend_read_서열_품질_짝_QCfail_제외수와cache(self):
        stage=cache.module('11')
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);ref,bam=ParabricksBamTest.fixture(self,root)
            stat=bam.stat();pre=root/'preflight.json'
            pre.write_text(json.dumps({'schema':'mfm_rna_preflight_v1','status':'pass',
                'rna':{'path':str(bam),'size':stat.st_size,'mtime_ns':stat.st_mtime_ns,'inode':stat.st_ino},
                'checks':{'sm_status':'pass','read_group_samples':['S'],'expected_samples':['S']}}))
            records={};qc={}
            for mode in ['parabricks','samtools']:
                out=root/mode
                args=['--sample','S','--rna',bam,'--out-dir',out,'--log-dir',out/'logs',
                      '--preflight-json',pre,'--samtools',shutil.which('samtools'),'--collate-threads','100','--fastq-threads','100']
                env={'MFM_BAM_BACKEND':'parabricks','RNA_FASTQ_BACKEND':mode,'MFM_PBRUN':PB,'MFM_GPU_RUNTIME':'native',
                     'MFM_CPU_BUDGET':'4','MFM_SAMPLE_CPU_BUDGET':'4','MFM_CPU_POOL_GRANTED':'4','MFM_RAM_LIMIT_BYTES':str(16<<30)}
                with patch.dict(os.environ,env):
                    cache.invoke(stage,args);before=(out/'R1.fq').stat().st_mtime_ns
                    cache.invoke(stage,args);self.assertEqual(before,(out/'R1.fq').stat().st_mtime_ns)
                def read(path):
                    lines=path.read_bytes().splitlines();return collections.Counter(tuple(lines[i:i+4]) for i in range(0,len(lines),4))
                records[mode]=[read(out/n) for n in ['R1.fq','R2.fq']]
                qc[mode]=json.loads((out/'fastq_qc.json').read_text())
                self.assertEqual(qc[mode]['fastq_backend'],mode)
            self.assertEqual(records['parabricks'],records['samtools'])
            for key in ['r1_pairs','r2_pairs','singleton_discard_records','unpaired_discard_records','paired_fraction','max_read_length']:
                self.assertEqual(qc['parabricks'][key],qc['samtools'][key],key)
            self.assertEqual(qc['samtools']['r1_pairs'],4)
            self.assertEqual(qc['samtools']['singleton_discard_records'],1)
            self.assertEqual(qc['samtools']['unpaired_discard_records'],1)
            self.assertEqual(qc['samtools']['backend_execution']['compression'],'uncompressed')
            self.assertNotEqual(qc['parabricks']['cache_fingerprint'],qc['samtools']['cache_fingerprint'])


if __name__=='__main__':unittest.main()
