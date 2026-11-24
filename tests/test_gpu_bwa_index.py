"""실제 BWA의 서열/도구 지문·동시 생성·무결성 및 stream 중첩 연결을 검증함."""
from concurrent.futures import ThreadPoolExecutor
from threading import Event
from contextlib import nullcontext
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import json,os,random,shutil,subprocess,tempfile,unittest
from integrate_wgs_rna import gpu_bwa_index as cache
from integrate_wgs_rna import gpu_wgs

ROOT=Path(__file__).resolve().parents[1]

class IndexTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        binary=Path(os.environ.get('MFM_GPU_BWA',str(ROOT/'tools/gpu/bwa')))
        if not binary.is_file():self.skipTest('실제 BWA 설치가 필요한 검사임')
        self.binary=binary
        self.env=patch.dict(os.environ,{'MFM_GPU_BWA':str(binary),'MFM_GPU_BWA_INDEX_ROOT':str(self.root/'cache')})
        self.env.start();self.addCleanup(self.env.stop)
        self.ref=self.root/'a.fa';rng=random.Random(17)
        self.sequence=''.join(rng.choices('ACGT',k=80000));self.ref.write_text('>chr1\n'+self.sequence+'\n')

    def test_경로_mtime가달라도같은서열은생성1회_실제정렬이기준과같음(self):
        first={};a=cache.ensure_index(self.ref,self.root/'work1',metrics=first)
        other=self.root/'other.fa';shutil.copyfile(self.ref,other);os.utime(other,(10,10))
        second={};b=cache.ensure_index(other,self.root/'work2',metrics=second)
        self.assertEqual(a,b);self.assertFalse(first['reused']);self.assertTrue(second['reused'])
        self.assertEqual(first['index_build_resources'],second['index_build_resources'])
        self.assertEqual(cache.lookup(other,self.root/'work3'),a)
        reads=self.root/'reads.fq';reads.write_text(''.join(f'@r{i}\n{self.sequence[i:i+100]}\n+\n'+('I'*100)+'\n' for i in range(20,220,10)))
        subprocess.run([str(self.binary),'index',str(self.ref)],check=True,capture_output=True)
        alignments=[]
        for ref in [self.ref,a]:
            result=subprocess.run([str(self.binary),'mem','-t','1',str(ref),str(reads)],check=True,capture_output=True,text=True)
            alignments.append([line for line in result.stdout.splitlines() if not line.startswith('@')])
        self.assertEqual(*alignments);self.assertTrue(alignments[0])

    def test_서열과도구가달라지면캐시를공유하지않음(self):
        original=cache.ensure_index(self.ref,self.root)
        st=self.ref.stat();self.ref.write_text(self.ref.read_text().replace('A','T',1));os.utime(self.ref,ns=(st.st_atime_ns,st.st_mtime_ns))
        self.assertIsNone(cache.lookup(self.ref,self.root))
        different=cache.ensure_index(self.ref,self.root);self.assertNotEqual(original,different)
        binary=self.root/'bwa_changed';shutil.copy2(self.binary,binary)
        with binary.open('ab') as stream:stream.write(b'changed-tool-identity')
        with patch.dict(os.environ,{'MFM_GPU_BWA':str(binary)}):self.assertIsNone(cache.lookup(self.ref,self.root))

    def test_동일서열_동시요청은1회만생성함(self):
        real=subprocess.run
        with patch.object(cache.subprocess,'run',wraps=real) as run:
            with ThreadPoolExecutor(2) as pool:results=list(pool.map(lambda _:cache.ensure_index(self.ref,self.root),range(2)))
        self.assertEqual(results[0],results[1]);self.assertEqual(run.call_count,1)

    def test_불완전하거나손상된캐시는재사용하지않음(self):
        ref=cache.ensure_index(self.ref,self.root);bwt=ref.with_suffix('.fa.bwt')
        bwt.write_bytes(bwt.read_bytes()[:-1])
        with self.assertRaises(ValueError):cache.ensure_index(self.ref,self.root)
        with self.assertRaises(ValueError):cache.lookup(self.ref,self.root)
        (ref.parent/'complete.json').unlink()
        with self.assertRaises(ValueError):cache.ensure_index(self.ref,self.root)

    def test_생성실패나입력변경은완료캐시를남기지않음(self):
        fail=self.root/'fail_bwa';fail.write_text('#!/bin/sh\nexit 9\n');fail.chmod(0o755)
        with patch.dict(os.environ,{'MFM_GPU_BWA':str(fail)}):
            with self.assertRaises(subprocess.CalledProcessError):cache.ensure_index(self.ref,self.root)
        self.assertFalse(list((self.root/'cache').glob('*/complete.json')))
        real=subprocess.run
        def mutate(*a,**kw):
            result=real(*a,**kw);self.ref.write_text(self.ref.read_text()+'\n');return result
        with patch.object(cache.subprocess,'run',side_effect=mutate):
            with self.assertRaises(ValueError):cache.ensure_index(self.ref,self.root)
        self.assertFalse(list((self.root/'cache').glob('*/complete.json')))
        self.assertFalse(list((self.root/'cache').glob('.bwa_build_*')))

    def test_부족한RAM에서인덱스를시작하지않음(self):
        budget=SimpleNamespace(memory=cache.GIB,phase=lambda *a,**k:None)
        with patch.object(cache.subprocess,'run') as run:
            with self.assertRaises(ValueError):cache.ensure_index(self.ref,self.root,budget)
        run.assert_not_called()

    def test_공유경로를명시하지않으면작업폴더안에둠(self):
        with patch.dict(os.environ,{},clear=True):self.assertEqual(cache.cache_root(self.root/'work'),self.root/'work/bwa_indexes')

class PreparationBudgetTests(unittest.TestCase):
    def test_입력종료뒤인덱스가남아도CPU는1개만유지함(self):
        from integrate_wgs_rna.gpu_prepare import prepare_short
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);ref=root/'ref';ref.write_text('>chr1\nACGT\n')
            allocated=Event();released=Event();events=[]
            def phase(label,cpus,memory,**kwargs):
                events.append((label,cpus,memory))
                if label=='bwa_index_and_input':allocated.set()
                if label=='input_cache_wait_index_active' and allocated.is_set():released.set()
                return True
            budget=SimpleNamespace(maximum_cpus=16,available_cpus=lambda:16,memory=128*cache.GIB,
                phase=phase,try_phase=phase)
            def indexer(reference,work):
                if not released.wait(5):raise AssertionError('입력 종료 뒤 CPU가 반환되지 않음')
                return reference
            def inputs(source,key,work,paired,share,producer,expected):
                share.phase('input',15,4*cache.GIB)
                return [root/'r1',root/'r2'],2,False
            with patch('integrate_wgs_rna.gpu_prepare.cached_inputs',side_effect=inputs):
                result=prepare_short(ref,root,'source','key',root,budget,3,100,indexer)
            self.assertEqual(events[-1][1],1)
            self.assertEqual(result[-1]['post_input_reserved_cpus'],1)

    def test_공유인덱스경로도다음단계설정에전달함(self):
        from integrate_wgs_rna.stage_gpu_config import apply_stage_environment
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);manifest=root/'sample_manifests/S/input.tsv';manifest.parent.mkdir(parents=True)
            (root/'effective_config.json').write_text(json.dumps(dict(run_id='R',data_root=str(root/'raw'),output_root=str(root/'out'),
                pipeline_env={'MFM_GPU_BWA_INDEX_ROOT':str(root/'cache')})))
            env=dict(RUN_ID='R',INPUT_ROOT=str(root/'raw'),OUTPUT_ROOT=str(root/'out'),PIPELINE_MANIFEST=str(manifest),MFM_WGS_BACKEND='parabricks')
            apply_stage_environment(env);self.assertEqual(env['MFM_GPU_BWA_INDEX_ROOT'],str(root/'cache'))

class StreamTests(unittest.TestCase):
    def run_case(self,warm,overlap_result=True):
        with tempfile.TemporaryDirectory() as tmp,patch.dict(os.environ,{'MFM_GPU_INPUT_MODE':'stream','MFM_GPU_PREPARE_OVERLAP':'1'},clear=True):
            root=Path(tmp);ref=root/'ref';ref.write_text('>chr1\nACGT\n');fq=[root/'r1',root/'r2']
            for p in fq:p.write_bytes(b'input')
            budget=SimpleNamespace(cpus=16,maximum_cpus=16,events=[],memory=128*cache.GIB,
                phase=lambda *a,**k:None,gpu=lambda *a,**k:nullcontext(None),output=lambda:None)
            def run(command,*a,**kw):Path(command[command.index('--out-bam')+1]).write_bytes(b'BAM');return {}
            runtime=SimpleNamespace(validate=lambda:None,run=run)
            stream=SimpleNamespace(paths=fq,reused=False,result=(fq,2),cpus=3,lengths=[],check=lambda:None,finish=lambda:(fq,2,False))
            args=['--reference',str(ref),'--platform','illumina','--preset','sr','--sample','S','--read-group','@RG\\tID:S',
                '--work',str(root),'--metadata',str(root/'meta.json'),'--threads','16','--input-source',str(root/'source.sh'),
                '--input-cache-key',json.dumps({'alignment':['unused',100]})]
            prepared=(ref,fq,2,False,{'index_seconds':3,'input_seconds':2,'elapsed_seconds':3}) if overlap_result else None
            with patch.object(gpu_wgs,'Runtime',return_value=runtime),patch.object(gpu_wgs,'StageBudget',return_value=budget),\
                 patch('integrate_wgs_rna.gpu_bwa_index.lookup',return_value=ref if warm else None),\
                 patch('integrate_wgs_rna.gpu_prepare.prepare_short',return_value=prepared) as overlap,\
                 patch('integrate_wgs_rna.gpu_stream_inputs.StreamInputs',return_value=nullcontext(stream)) as streaming,\
                 patch.object(gpu_wgs,'bwa_reference',return_value=ref) as index,patch.object(gpu_wgs,'restore_bam',return_value={}):
                self.assertEqual(gpu_wgs.main(args),0)
            meta=json.loads((root/'meta.json').read_text())
            self.assertEqual(overlap.call_count,0 if warm else 1)
            self.assertEqual(index.call_count,0 if warm or overlap_result else 1)
            self.assertEqual(streaming.call_count,int(warm or not overlap_result))
            self.assertEqual(meta['input_streamed'],warm or not overlap_result)
            self.assertEqual(meta['input_records'],2)
    def test_cold_stream은먼저인덱스와입력을겹쳐준비함(self):self.run_case(False)
    def test_warm_stream은중간FASTQ생성을추가하지않음(self):self.run_case(True)
    def test_중첩자원이부족하면기존순차준비를유지함(self):self.run_case(False,False)

if __name__=='__main__':unittest.main()
