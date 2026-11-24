"""GPU를 점유하지 않고 HAP 자원·입력 공유·출력 데이터 계약을 검증한다."""
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import gzip
import io
import json
import os
import subprocess
import sys
import tempfile
import threading
import unittest

import numpy as np
import pysam

from integrate_wgs_rna.gpu_fastq import split_fastq, cached_inputs
from integrate_wgs_rna.gpu_stage_resources import StageBudget
from integrate_wgs_rna.gpu_wgs import restore_bam
from integrate_wgs_rna.memory_policy import GIB
from integrate_wgs_rna.realign.phases import _atomic
from integrate_wgs_rna.realign.align_resources import alignment_plan
from integrate_wgs_rna import gpu_arrays

RAW=b'@read/1 HP:i:1 PS:i:123 OD:i:1\nACGT\n+\nIIII\n@read/2 HP:i:1 PS:i:123 OD:i:1\nTGCA\n+\nIIII\n'
AUTO={'WGS_AUTO_THREADS':'1','MFM_WGS_ALIGN_REQUESTED':'auto','MFM_WGS_SORT_REQUESTED':'auto','MFM_WGS_HAP_REQUESTED':'auto'}


class PipelineTests(unittest.TestCase):
    def test_병렬gzip은읽기순서와태그를바꾸지않음(self):
        with tempfile.TemporaryDirectory() as raw:
            p=Path(raw);files=[p/'r1.gz',p/'r2.gz']
            self.assertEqual(split_fastq(io.BytesIO(RAW*1000),files,True,cpus=16,producer_cpus=3),2000)
            self.assertEqual(gzip.open(files[0],'rb').read(),RAW.split(b'@read/2')[0]*1000)
            self.assertEqual(gzip.open(files[1],'rb').read(),(b'@read/2'+RAW.split(b'@read/2')[1])*1000)

    def test_잘린입력과다른mate는완료로처리하지않음(self):
        with tempfile.TemporaryDirectory() as raw:
            files=[Path(raw)/'r1.gz',Path(raw)/'r2.gz']
            for value in [RAW[:-3],RAW.replace(b'@read/2',b'@other/2'),RAW.split(b'@read/2')[0]]:
                with self.assertRaises(ValueError):split_fastq(io.BytesIO(value),files,True,cpus=16)

    def test_공통FASTQ는동시HAP에서도한번만생성됨(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{},clear=True):
            p=Path(raw);data=p/'input.fq';data.write_bytes(RAW)
            script=p/'source.sh';script.write_text('echo run >> '+str(p/'calls')+'\ncat '+str(data)+'\n')
            results=[];errors=[]
            def run():
                try:results.append(cached_inputs(script,'same verified source',p,True,StageBudget(),0))
                except BaseException as e:errors.append(e)
            threads=[threading.Thread(target=run) for _ in range(2)]
            for t in threads:t.start()
            for t in threads:t.join(5)
            self.assertFalse(any(t.is_alive() for t in threads));self.assertEqual(errors,[])
            self.assertEqual(len(results),2);self.assertEqual((p/'calls').read_text(),'run\n')
            self.assertEqual(sorted(r[2] for r in results),[False,True])
            results[0][0][0].write_bytes(b'bad')
            with self.assertRaises(ValueError):cached_inputs(script,'same verified source',p,True,StageBudget(),0)

    def test_실패한source의부분파일은재사용하지않음(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{},clear=True):
            p=Path(raw);s=p/'source.sh';s.write_text('printf "@r\\nA\\n+\\nI\\n"\nexit 3\n')
            with self.assertRaises(RuntimeError):cached_inputs(s,'key',p,False,StageBudget(),0)
            self.assertEqual(list(p.glob('gpu_input_*/complete.json')),[])

    def test_CPU14개에서도가능한두HAP를선택함(self):
        c=SimpleNamespace(platform='illumina',threads=16,sort_threads=4,sort_mem='2G',hap_jobs=2)
        result=alignment_plan(c,{1:3000,2:3000},14,300*GIB,AUTO)
        self.assertEqual(result['hap_jobs'],2);self.assertEqual(result['threads'],4)
        self.assertEqual(result['reserved_cpus'],14)
        self.assertEqual(alignment_plan(c,{1:3000,2:3000},14,300*GIB,dict(AUTO,MFM_WGS_HAP_REQUESTED='1'))['hap_jobs'],1)

    def test_실제CPU14개두HAP정렬과준비단계QC가동일함(self):
        from test_stage22_bwa import Stage22BwaIntegrationTest,module
        original=module.run_personal_realign_generation
        timings=[]
        def run(config):
            result=original(config)
            path=config.output_root/'TEST/qc/TEST.stage11_timing.json'
            if path.exists():timings.append(json.loads(path.read_text()))
            return result
        with patch.dict(os.environ,dict(AUTO,MFM_CPU_BUDGET='14')),patch.object(module,'run_personal_realign_generation',side_effect=run):
            Stage22BwaIntegrationTest('test_pair_singleton_duplicate_and_chain_contract').debug()
        self.assertEqual(timings[0]['effective_threads']['hap_jobs'],2)
        self.assertEqual(timings[0]['effective_threads']['align_per_hap'],4)
        self.assertEqual(timings[0]['resource_phases'][0]['phase'],'input_prepare')
        self.assertTrue(all(r['reserved_cpus']<=14 for r in timings[0]['resource_phases']))

    def test_작은진입RAM은확인된단계예약에서만사용함(self):
        import importlib.util
        root=Path(__file__).resolve().parents[1]
        spec=importlib.util.spec_from_file_location('entry_budget_test',root/'scripts/cpu_pool.py')
        pool=importlib.util.module_from_spec(spec);spec.loader.exec_module(pool)
        env={'MFM_CPU_POOL_LAYOUT':'original','MFM_REALIGN_PREPARE_MAX_CPUS':'4'}
        with patch.dict(os.environ,env,clear=True):
            self.assertEqual(pool.stage_memory('WGS','10'),4*GIB)
            with patch.dict(os.environ,MFM_REALIGN_PHASE_RESOURCES='0'):
                self.assertEqual(pool.stage_memory('WGS','10'),192*GIB)
            with patch.dict(os.environ,WGS_REALIGN_JOBS='2'):
                self.assertEqual(pool.stage_memory('WGS','10'),192*GIB)

    def test_GPU입력예약과실제host피크는구분함(self):
        c=SimpleNamespace(platform='hifi',threads=16,sort_threads=4,sort_mem='2G',hap_jobs=2)
        r=alignment_plan(c,{1:3_000_000_000,2:3_000_000_000},32,20*GIB,dict(AUTO,MFM_WGS_BACKEND='parabricks',MFM_GPU_HAP_JOBS='2'))
        self.assertEqual(r['hap_jobs'],2);self.assertEqual(r['reserved_memory_bytes'],8*GIB)
        self.assertGreater(r['memory_model']['gpu_compute_peak_bytes'],250*GIB)

    def test_GPU대기는CPU_RAM을반환하고확보실패GPU도즉시반환함(self):
        with tempfile.TemporaryDirectory() as raw:
            p=Path(raw);(p/'leases').mkdir();(p/'waiters').mkdir()
            own=p/'leases/own.json';other=p/'leases/other.json'
            _atomic(own,{'cpus':16,'memory_bytes':4*GIB,'total_cpus':68,'effective_alignment':{'hap_jobs':1}})
            _atomic(other,{'cpus':64,'memory_bytes':4*GIB})
            closed=[]
            def lease(*a,**kw):return SimpleNamespace(close=lambda:closed.append(True),device={'uuid':'test'},reserved=GIB)
            def idle(_):
                row=json.loads(own.read_text());self.assertEqual(row['cpus'],1);self.assertEqual(row['memory_bytes'],4*GIB)
                _atomic(other,{'cpus':40,'memory_bytes':4*GIB})
            env={'MFM_CPU_POOL_LEASE_FILE':str(own),'MFM_GPU_HAP_ID':'1','MFM_GPU_STAGE_MEMORY_BYTES':str(150*GIB)}
            with patch.dict(os.environ,env,clear=True):
                budget=StageBudget()
                with patch('integrate_wgs_rna.gpu_pool.acquire',side_effect=lease),patch('integrate_wgs_rna.gpu_stage_resources.time.sleep',side_effect=idle):
                    with budget.gpu(GIB,purpose='test'):
                        r=json.loads(own.read_text());self.assertEqual(r['cpus'],4);self.assertEqual(r['memory_bytes'],150*GIB)
                self.assertEqual(len(closed),1)

    def test_SlurmRAM보다큰GPU요청은입력생성전에거부함(self):
        with patch.dict(os.environ,{'SLURM_MEM_PER_NODE':'65536'},clear=True):
            with self.assertRaisesRegex(ValueError,'실제 전체 예약'):
                StageBudget()

    def test_다른HAP의GPU대기예약을출력QC가덮어쓰지않음(self):
        from integrate_wgs_rna.realign.hap_resources import HaplotypeResources
        with tempfile.TemporaryDirectory() as raw:
            p=Path(raw);(p/'leases').mkdir();own=p/'leases/own.json'
            _atomic(own,{'cpus':32,'memory_bytes':8*GIB,'total_cpus':68,'manager_pid':os.getpid(),'child_pid':os.getpid()})
            phase=SimpleNamespace(path=own,total=68,ceiling=lambda:900*GIB)
            owner=HaplotypeResources(phase,None,[1,2])
            row=json.loads(own.read_text());row['hap_resources']['2']={'cpus':1,'memory_bytes':4*GIB,'phase':'gpu_wait'};_atomic(own,row)
            owner.controller(1).reserve('output_qc',max_units=2)
            row=json.loads(own.read_text());self.assertEqual(row['hap_resources']['2']['cpus'],1)
            owner.release(1);self.assertEqual(json.loads(own.read_text())['cpus'],1)

    def test_쓰기피크는다른HAP와합산하고공간대기CPU를반환함(self):
        with tempfile.TemporaryDirectory() as raw:
            p=Path(raw);(p/'leases').mkdir();own=p/'leases/own.json'
            _atomic(own,{'cpus':16,'memory_bytes':4*GIB,'total_cpus':68,'effective_alignment':{'hap_jobs':1}})
            _atomic(p/'leases/other.json',{'cpus':1,'memory_bytes':4*GIB,'hap_resources':{'2':{
                'disk_peak_bytes':15*GIB,'write_paths':[str(p/'other.bam')]}}})
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(own),'MFM_GPU_HAP_ID':'1'},clear=True), \
                 patch('integrate_wgs_rna.gpu_stage_resources.shutil.disk_usage',return_value=SimpleNamespace(free=150_000_000_000+20*GIB)):
                b=StageBudget()
                self.assertFalse(b.try_phase('gpu_compute',16,150*GIB,write_paths=[p/'out.bam'],disk_peak_bytes=10*GIB))
                self.assertEqual(b.last_wait_reason['resource'],'disk')
                b.wait('storage_wait');self.assertEqual(json.loads(own.read_text())['cpus'],0)
                self.assertTrue(b.try_phase('gpu_compute',16,150*GIB,write_paths=[p/'out.bam'],disk_peak_bytes=4*GIB))
                b.wait('output_wait',preserve_disk=True)
                self.assertEqual(json.loads(own.read_text())['hap_resources']['1']['disk_peak_bytes'],4*GIB)
                b.phase('output_restore_sort',2,4*GIB,preserve_disk=True)
                self.assertEqual(json.loads(own.read_text())['hap_resources']['1']['disk_peak_bytes'],4*GIB)

    def test_BAM직접출력은RG와duplicate만복원하고태그를보존함(self):
        with tempfile.TemporaryDirectory() as raw:
            p=Path(raw);src=p/'in.bam';out=p/'out.bam'
            with pysam.AlignmentFile(src,'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':32}],'RG':[{'ID':'old','SM':'old'}]}) as f:
                r=pysam.AlignedSegment();r.query_name='read';r.query_sequence='ACGT';r.query_qualities=pysam.qualitystring_to_array('IIII')
                r.flag=0;r.reference_id=0;r.reference_start=4;r.cigar=[(0,4)];r.mapping_quality=60
                for k,v in [('RG','old'),('OD',1),('HP',2),('PS',123),('NM',0)]:r.set_tag(k,v)
                f.write(r)
            restore_bam(src,'@RG\\tID:S\\tSM:S','S','preserve_source',out)
            with pysam.AlignmentFile(out,'rb') as f:
                r=next(f);self.assertEqual(r.query_sequence,'ACGT');self.assertEqual(r.reference_start,4)
                self.assertTrue(r.is_duplicate);self.assertFalse(r.has_tag('OD'))
                self.assertEqual(r.get_tag('RG'),'S');self.assertEqual(r.get_tag('HP'),2);self.assertEqual(r.get_tag('PS'),123)
            legacy=p/'legacy.sam'
            restore_bam(src,'@RG\\tID:S\\tSM:S','S',destination=legacy,output_mode='wh')
            with pysam.AlignmentFile(str(legacy),'r') as f:
                r=next(f);self.assertTrue(r.has_tag('OD'));self.assertFalse(r.is_duplicate)
                self.assertEqual(r.get_tag('RG'),'S')

    def test_깊이와QC를한GPU배열에서계산해도CPU와동일함(self):
        @contextmanager
        def session(*a,**kw):
            yield SimpleNamespace(**{k:getattr(np,k) for k in ['zeros','int64','uint32','uint64','asarray','add','cumsum','iinfo','sum','count_nonzero']},asnumpy=np.asarray)
        with patch.object(gpu_arrays,'session',session):
            values,metrics=gpu_arrays.depth_from_intervals([0,2],[4,6],8,return_metrics=True)
        np.testing.assert_array_equal(values,[1,1,2,2,1,1,0,0]);self.assertEqual(metrics['depth_sum'],8)
        self.assertEqual(metrics['depth_nonzero_count'],6);self.assertEqual(metrics['depth_max'],2)

    def test_실제생성쉘은한입력으로두HAP의CRAM과태그를보존함(self):
        self._실제생성쉘_검증(False)

    def test_로컬임시CRAM링크도두HAP_read와태그를그대로보존함(self):
        self._실제생성쉘_검증(True)

    def _실제생성쉘_검증(self, scratch):
        from integrate_wgs_rna.gpu_wgs import write_script
        from integrate_wgs_rna.realign.config import PersonalRealignConfig
        from types import SimpleNamespace
        import shutil
        with tempfile.TemporaryDirectory(prefix='gpu pipeline ') as raw:
            p=Path(raw);ref=p/'ref.fa';ref.write_text('>chr1\n'+'ACGT'*250+'\n');pysam.faidx(str(ref))
            bam=p/'prepared.bam'
            with pysam.AlignmentFile(str(bam),'wb',header={'SQ':[{'SN':'chr1','LN':1000}]}) as f:
                for i,seq in enumerate(['ACGT','TGCA']):
                    r=pysam.AlignedSegment();r.query_name='read';r.flag=65 if i==0 else 129
                    r.reference_id=0;r.reference_start=10+i;r.next_reference_id=0;r.next_reference_start=11-i
                    r.query_sequence=seq;r.query_qualities=pysam.qualitystring_to_array('IIII');r.cigar=[(0,4)]
                    r.set_tag('HP',1);r.set_tag('PS',123);r.set_tag('OD',1);f.write(r)
            # 정렬 도구만 대역으로 바꾸고 실제 source shell/압축/복원/samtools CRAM을 실행한다.
            (p/'sitecustomize.py').write_text('''import gzip
from types import SimpleNamespace
import pysam
from integrate_wgs_rna import gpu_pool
from integrate_wgs_rna.parabricks_runtime import Runtime
gpu_pool.acquire=lambda *a,**k:SimpleNamespace(device={'uuid':'test'},reserved=1,close=lambda:None)
def run(self,args,gpu,paths=(),**kwargs):
    fq=args[args.index('--in-fq')+1:args.index('--out-bam')]
    out=args[args.index('--out-bam')+1]
    with pysam.AlignmentFile(out,'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}]}) as f:
        for i,path in enumerate(fq):
            with gzip.open(path,'rt') as s:
                head=s.readline();seq=s.readline().strip();s.readline();qual=s.readline().strip()
            r=pysam.AlignedSegment();r.query_name=head.split()[0][1:];r.flag=65 if i==0 else 129
            r.reference_id=0;r.reference_start=10+i;r.cigar=[(0,len(seq))];r.query_sequence=seq
            r.query_qualities=pysam.qualitystring_to_array(qual)
            for field in head.split()[1:]:
                key,kind,value=field.split(':',2);r.set_tag(key,int(value) if kind=='i' else value)
            f.write(r)
    return {'gpu_uuid':'test'}
Runtime.run=run
''')
            cfg=PersonalRealignConfig(cohort='test',output_root=p,sample_list=None,samples=['S'],chroms=['chr1'],
                sequence_variant_policy='snv_indel',dataset_role='pre_training',platform='illumina',preset='sr',
                threads=4,sort_threads=0,reference_fasta=ref,hap_read_mode='all_to_both',read_routing_policy='all_to_both')
            root=Path(__file__).resolve().parents[1]
            bwa=Path(os.environ.get('MFM_GPU_BWA',str(root/'tools/gpu/bwa')))
            if not bwa.is_file():self.skipTest('실제 BWA 인덱스 생성 도구가 필요한 검사임')
            env=dict(os.environ,PYTHONPATH=str(p)+':'+str(root/'src'),MFM_GPU_RUNTIME='native',MFM_PBRUN='/bin/true',
                     MFM_GPU_BWA=str(bwa),MFM_WGS_BACKEND='parabricks')
            env.pop('MFM_CPU_POOL_LEASE_FILE',None)
            if scratch:env.update(MFM_GPU_SCRATCH_ROOT=str(p/'scratch'),MFM_GPU_ALIGNMENT_GPUS='3')
            for hap in [1,2]:
                out=p/f'hap{hap}.cram'
                with patch.dict(os.environ,env):script=write_script(cfg,'S',hap,bam,ref,out,['chr1'],p,p,SimpleNamespace(path=bam,name_collated=True))
                if hap==2:(p/'S.HAP2.input.sh').write_text('exit 99\n')
                result=subprocess.run(['bash',str(script)],env=env,capture_output=True,timeout=20)
                self.assertEqual(result.returncode,0,result.stderr.decode(errors='replace'))
                self.assertEqual(out.is_symlink(),scratch)
                with pysam.AlignmentFile(str(out),'rc',reference_filename=str(ref)) as f:
                    reads=list(f)
                self.assertEqual([r.query_sequence for r in reads],['ACGT','TGCA'])
                self.assertTrue(all(r.is_duplicate and r.get_tag('RG')=='S' and r.get_tag('PS')==123 and not r.has_tag('OD') for r in reads))
            from integrate_wgs_rna.gpu_wgs import scratch_directory
            with patch.dict(os.environ,env):cache=scratch_directory(p)
            self.assertEqual(len(list(cache.glob('gpu_input_*/complete.json'))),1)
            self.assertTrue(json.loads((p/'S.HAP2.gpu.json').read_text())['input_reused'])


if __name__=='__main__':unittest.main()
