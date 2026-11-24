"""GPU 예약/명령/수치·Zarr 형식의 CPU 동등성을 검증함. GPU 실행은 명시할 때만 수행함."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import gzip,io,json,os,subprocess,sys,tempfile,unittest
import numpy as np
import pysam,zarr

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src'),str(ROOT/'vendor/MFM-DNARNAFM')]
from handoff import acceleration
from integrate_wgs_rna import gpu_pool,gpu_wgs,gpu_rna,gpu_arrays
from integrate_wgs_rna.parabricks_runtime import Runtime

class SelectionTests(unittest.TestCase):
    def test_CPU설정은기존설정을보존함(self):
        original={'pipeline_env':{'WINDOW_STORAGE_POLICY':'annotated_acgtn_v1'}}
        self.assertEqual(acceleration.apply(original),original)
        c=acceleration.apply(original,{'wgs_backend':'parabricks','coverage_backend':'cupy','gpu_devices':'0,4'})
        self.assertEqual(c['pipeline_env']['MFM_WGS_BACKEND'],'parabricks')
        self.assertEqual(c['pipeline_env']['MFM_GPU_DEVICES'],'0,4')
        self.assertEqual(original,{'pipeline_env':{'WINDOW_STORAGE_POLICY':'annotated_acgtn_v1'}})

    def test_할당범위보다GPU목록을넓히지못함(self):
        devices=[{'index':i,'uuid':'GPU-'+str(i)} for i in range(4)]
        result=gpu_pool.allowed_devices(devices,'0,1,2',{'CUDA_VISIBLE_DEVICES':'1,2',
                         'NVIDIA_VISIBLE_DEVICES':'GPU-2,GPU-3'})
        self.assertEqual([d['index'] for d in result],[2])
        self.assertEqual(gpu_pool.allowed_devices(devices,'auto',{'CUDA_VISIBLE_DEVICES':''}),[])

    def test_다른작업GPU는피하고공유계산은분산함(self):
        with tempfile.TemporaryDirectory() as tmp,patch.dict(os.environ,{'MFM_GPU_POOL_DIR':tmp,'MFM_GPU_DEVICES':'auto'}):
            devices=[dict(index=i,uuid='GPU-'+str(i),pci=f'0000:0{i}:00.0',total=32*gpu_pool.GIB,
                          used=(4*gpu_pool.GIB if i==0 else 0),utilization=(99 if i==0 else 0)) for i in range(3)]
            with patch.object(gpu_pool,'inventory',return_value=devices):
                one=gpu_pool.acquire(2*gpu_pool.GIB,purpose='test',timeout=0)
                two=gpu_pool.acquire(2*gpu_pool.GIB,purpose='test',timeout=0)
                self.assertEqual(one.device['index'],1);self.assertEqual(two.device['index'],2)
                with self.assertRaises(TimeoutError):gpu_pool.acquire(2*gpu_pool.GIB,purpose='align',exclusive=True,timeout=0)
                two.close()
                exclusive=gpu_pool.acquire(2*gpu_pool.GIB,purpose='align',exclusive=True,timeout=0)
                self.assertEqual(exclusive.device['index'],2);self.assertEqual(exclusive.reserved,32*gpu_pool.GIB)
                with self.assertRaises(TimeoutError):gpu_pool.acquire(2*gpu_pool.GIB,purpose='busy',exclusive=True,timeout=0)
                exclusive.close();one.close();two.close()

    def test_WGS입력주석과paired정합성유지함(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths=[Path(tmp)/'r1.fq.gz',Path(tmp)/'r2.fq.gz']
            raw=b'@r/1 HP:i:1 OD:i:1\nACGT\n+\nIIII\n@r/2 HP:i:1 OD:i:1\nTGCA\n+\nIIII\n'
            self.assertEqual(gpu_wgs.split_fastq(io.BytesIO(raw),paths,True),2)
            self.assertIn(b'OD:i:1',gzip.open(paths[0],'rb').read())
            with self.assertRaisesRegex(ValueError,'mate 누락'):gpu_wgs.split_fastq(io.BytesIO(raw.split(b'@r/2')[0]),paths,True)
        command=gpu_wgs.commands('/ref.fa',['/r1.fq','/r2.fq'],'/out.bam','illumina','sr',8,'S')
        self.assertIn('--no-markdups',command);self.assertIn('--bwa-options=-Y -C',command)
        self.assertEqual(command[command.index('--read-group-id-prefix')+1],'S')
        self.assertNotIn('--knownSites',command)
        command=gpu_wgs.commands('/ref.fa',['/r.fq'],'/out.bam','hifi','map-hifi',8,'S')
        self.assertEqual(command[command.index('--preset')+1],'map-hifi');self.assertIn('--copy-comment',command)

    def test_GPU_STAR는실제차이를기록하고알수없는옵션을거부함(self):
        import rna_profile
        values=gpu_rna.parse(['--runMode','alignReads','--runThreadN','8','--genomeDir','/index',
                  '--readFilesIn','/r1','/r2','--outFileNamePrefix','/output/',
                  '--readFilesCommand','cat','--outSAMattrRGline','ID:S','SM:S',*rna_profile.STAR_ALIGN_ARGS])
        with self.assertRaisesRegex(RuntimeError,'BySJout'):gpu_rna.align_arguments(values,'/ref',Path('/out/a.bam'),8)
        with patch.dict(os.environ,{'MFM_GPU_RNA_POLICY':'experimental'}):command=gpu_rna.align_arguments(values,'/ref.fa',Path('/output/a.bam'),8,74)
        self.assertIn('--two-pass-mode',command);self.assertIn('--no-markdups',command)
        self.assertNotIn('--max-collapsed-junctions',command)
        self.assertEqual(command[command.index('--sjdb-overhang')+1],'74')
        self.assertEqual(gpu_rna.effective_contract()['differences']['--limitSjdbInsertNsj']['gpu'],['1000000'])
        self.assertEqual(command[command.index('--max-out-filter-multimap')+1],'20')
        self.assertFalse(gpu_rna.effective_contract()['profile_equivalent'])
        self.assertEqual(gpu_rna.effective_contract()['differences']['--outFilterType']['gpu'],['Normal'])
        with patch.dict(os.environ,{'MFM_GPU_RNA_POLICY':'experimental'}),self.assertRaisesRegex(ValueError,'미지원'):gpu_rna.align_arguments({**values,'--unknown':['1']},'/ref',Path('/output/a.bam'),8)

    def test_컨테이너에는예약GPU만공개하고명령은shell문자열로합치지않음(self):
        with tempfile.TemporaryDirectory(prefix='gpu paths ') as tmp,patch.dict(os.environ,{'MFM_GPU_RUNTIME':'docker'}),patch('shutil.which',return_value='/bin/docker'):
            root=Path(tmp);ref=root/'ref.fa';ref.touch();out=root/'out.bam'
            command=Runtime().command(['minimap2','--ref',str(ref),'--out-bam',str(out)],'GPU-example',[ref,out])
            self.assertIn('device=GPU-example',command)
            self.assertIn(str(ref),command);self.assertNotIn('--gpus all',' '.join(command))
            self.assertIn('--memory',command);self.assertIn('--cpus',command)
            self.assertIn('HOME=/tmp',command);self.assertIn('TMPDIR=/tmp',command)

    def test_GPU는두HAP준비출력을겹치고장치단계만직렬화함(self):
        from integrate_wgs_rna.realign.align_resources import alignment_plan
        config=SimpleNamespace(platform='hifi',threads=24,sort_threads=4,sort_mem='2G',hap_jobs=2,input_decode_threads=0)
        env={'MFM_WGS_BACKEND':'parabricks','WGS_AUTO_THREADS':'1','MFM_WGS_ALIGN_REQUESTED':'auto',
             'MFM_WGS_SORT_REQUESTED':'auto','MFM_WGS_HAP_REQUESTED':'auto'}
        plan=alignment_plan(config,{1:3_000_000_000,2:3_000_000_000},68,900_000_000_000,env)
        self.assertEqual(plan['hap_jobs'],2);self.assertLessEqual(plan['reserved_cpus'],68)
        self.assertEqual(plan['reserved_memory_bytes'],8*1024**3)
        self.assertGreaterEqual(plan['memory_model']['gpu_compute_peak_bytes'],128*1024**3)

    def test_Slurm_GPU요청은worker에만추가함(self):
        from handoff.slurm import commands
        c=acceleration.apply({'state_root':'/tmp/state','slurm':{'cpus':16,'memory':'192G',
            'time':'12:00:00','concurrency':2}}, {'wgs_backend':'parabricks'})
        prepare,array,finish=commands(c,ROOT,Path('/tmp/c.json'),Path('/tmp/selection.tsv'),2)
        self.assertIn('--gpus-per-task=1',array)
        self.assertFalse(any('--gpus' in x for x in prepare+finish))

    def test_WGS생성스크립트는준비된paired입력과기존태그복원을사용함(self):
        from integrate_wgs_rna.realign.config import PersonalRealignConfig
        with tempfile.TemporaryDirectory(prefix='gpu script ') as tmp:
            root=Path(tmp)
            config=PersonalRealignConfig(cohort='test',output_root=root,sample_list=None,samples=['S'],
                     chroms=['chr1'],sequence_variant_policy='snv_indel',dataset_role='pre_training',
                     platform='illumina',preset='sr',threads=4,sort_threads=1)
            prepared=SimpleNamespace(path=root/'prepared.bam',name_collated=True)
            script=gpu_wgs.write_script(config,'S',1,root/'original.bam',root/'target.fa',root/'out.cram',
                    ['chr1'],root,root,prepared)
            subprocess.run(['bash','-n',str(script)],check=True)
            text=script.read_text()
            self.assertIn('prepared.bam',text);self.assertNotIn('INPUT='+str(root/'original.bam'),text)
            self.assertIn('--duplicate-mode preserve_source',text)
            self.assertIn('HP,PS,XU,OD',(root/'S.HAP1.input.sh').read_text())
            self.assertIn('integrate_wgs_rna.gpu_wgs',text)

    def test_명시한두HAP는대기자수가많아도실행가능한최소CPU를받음(self):
        code = """
import json,os,tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from integrate_wgs_rna.realign.phases import RealignPhases
with tempfile.TemporaryDirectory() as raw:
    root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir()
    path=root/'leases/own.json'
    path.write_text(json.dumps({'cpus':1,'total_cpus':68,'memory_bytes':4*1024**3,'manager_pid':os.getpid(),'child_pid':os.getpid()}))
    (root/'leases/other.json').write_text(json.dumps({'cpus':48,'memory_bytes':1024**3,'manager_pid':os.getppid(),'child_pid':os.getppid()}))
    env={'MFM_CPU_POOL_LEASE_FILE':str(path),'MFM_CPU_BUDGET':'68','MFM_RAM_LIMIT_BYTES':'1000000000000',
         'MFM_WGS_BACKEND':'parabricks','MFM_GPU_HAP_JOBS':'2','MFM_GPU_HOST_MEMORY_BYTES':str(32*1024**3),
         'MFM_WGS_HAP_REQUESTED':'2','MFM_WGS_ALIGN_REQUESTED':'auto','MFM_WGS_SORT_REQUESTED':'auto','WGS_AUTO_THREADS':'1','MFM_WGS_MIN_ALIGN_THREADS':'1'}
    from integrate_wgs_rna.realign.config import PersonalRealignConfig
    config=PersonalRealignConfig(cohort='test',output_root=root,sample_list=None,samples=['S'],chroms=['chr1'],sequence_variant_policy='snv_indel',dataset_role='pre_training',platform='illumina',threads=2,sort_threads=1,sort_mem='64M',hap_jobs=2,input_decode_threads=0)
    with patch.dict(os.environ,env),patch('integrate_wgs_rna.realign.phases._living_waiters',return_value=100):
        phase=RealignPhases();actual=phase.align_config(config,{1:3000,2:3000})
        assert actual.hap_jobs==2
        assert phase.resource_limits['cpus']==8
"""
        subprocess.run([sys.executable,'-c',code],check=True,timeout=5)

    def test_GPU도구의hostRSS상한초과를실제로중지함(self):
        with tempfile.TemporaryDirectory() as tmp,patch.dict(os.environ,{'MFM_GPU_RUNTIME':'native','MFM_PBRUN':sys.executable}):
            runtime=Runtime();gpu=SimpleNamespace(device={'uuid':'GPU-test'},reserved=24*1024**3)
            with self.assertRaises(MemoryError):
                runtime.run(['-c','import time;time.sleep(2)'],gpu,host_limit=1,log=Path(tmp)/'tool.log')


@unittest.skipUnless(os.environ.get('MFM_RUN_GPU_TESTS')=='1','실제 GPU 검증은 MFM_RUN_GPU_TESTS=1로 선택함')
class NumericGpuTests(unittest.TestCase):
    def test_정수depth산술및통계일치함(self):
        random=np.random.default_rng(42);starts=random.integers(0,100000,10000);ends=np.minimum(starts+random.integers(0,1000,10000),100000)
        delta=np.zeros(100001,np.int64);np.add.at(delta,starts,1);np.add.at(delta,ends,-1)
        expected=np.cumsum(delta[:-1]).astype(np.uint32)
        observed=gpu_arrays.depth_from_intervals(starts,ends,100000)
        np.testing.assert_array_equal(expected,observed)
        total=int(expected.sum(dtype=np.uint64));nz=int(np.count_nonzero(expected))
        self.assertEqual(gpu_arrays.depth_metrics(observed,len(observed)),{'depth_sum':total,'depth_mean':total/len(expected),
                     'depth_zero_fraction':(len(expected)-nz)/len(expected),'depth_max':int(expected.max())})

    def test_4bit바이트와오류위치가CPU와일치함(self):
        from integrate_wgs_rna.legacy_shared_hap_converter import pack_4bit_acgtnp_chunk,_pack_4bit_acgtnp_chunk_legacy
        sets=[[],[''],['A','NNN','acgt','P','Rr?',''],['ACGTN'*1001]*1200]
        with patch.dict(os.environ,{'MFM_ZARR_BACKEND':'cupy'}):
            for sequences in sets:
                cpu=_pack_4bit_acgtnp_chunk_legacy(sequences,frozenset('ACGTNP'),17,'to_n')
                gpu=pack_4bit_acgtnp_chunk(sequences,frozenset('ACGTNP'),17,'to_n')
                for a,b in zip(cpu[:3],gpu[:3]):np.testing.assert_array_equal(a,b)
                self.assertEqual(cpu[3:],gpu[3:])
            for function in [_pack_4bit_acgtnp_chunk_legacy,pack_4bit_acgtnp_chunk]:
                with self.assertRaisesRegex(ValueError,"record_index=18 position=1"):
                    function(['A','T?'],frozenset('ACGTN'),17,'fail')

    def test_CIGAR_GPU검증은오류를숨기지않음(self):
        from integrate_wgs_rna.compact.chunk_validation import validate_chunk
        arrays={'cigar_offset':np.asarray([0,3]),'cigar_len':np.asarray([3,1]),
                'lengths':np.asarray([5,2]),'ref_start':np.asarray([0,10]),'ref_end':np.asarray([6,12])}
        ops=np.asarray([(2<<4)|0,(1<<4)|2,(3<<4)|0,(2<<4)|0],np.uint32)
        with patch.dict(os.environ,{'MFM_QC_BACKEND':'cupy'}):
            self.assertEqual(validate_chunk(arrays,ops,0,2),2)
            arrays['lengths'][1]=3
            with self.assertRaisesRegex(ValueError,'CIGAR/query mismatch'):validate_chunk(arrays,ops,0,2)

    def test_RNA_GPUdepth는samtools_J중복과스플라이스규칙까지일치함(self):
        from integrate_wgs_rna.gpu_rna_coverage import contig_depth
        import rna_profile
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);bam=root/'in.bam'
            specs=[(0,255,[(0,5),(2,2),(0,3)]),(1024,255,[(0,3),(3,5),(0,4)]),
                   (0,254,[(0,15)]),(256,255,[(0,15)]),(512,255,[(0,15)]),(2048,255,[(0,15)])]
            with pysam.AlignmentFile(bam,'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':32}]}) as out:
                for i,(flag,mapq,cigar) in enumerate(specs):
                    read=pysam.AlignedSegment();read.query_name='r'+str(i);read.flag=flag
                    read.reference_id=0;read.reference_start=0;read.mapping_quality=mapq;read.cigar=cigar
                    n=sum(length for op,length in cigar if op in {0,1,4,7,8})
                    read.query_sequence='A'*n;read.query_qualities=pysam.qualitystring_to_array('I'*n);out.write(read)
            pysam.index(str(bam))
            expected=np.zeros(32,np.uint32)
            result=subprocess.run(['samtools','depth',*rna_profile.COVERAGE_DEPTH_ARGS,'-r','chr1',str(bam)],text=True,capture_output=True,check=True)
            for line in result.stdout.splitlines():
                _,pos,depth=line.split('\t');expected[int(pos)-1]=int(depth)
            store=root/'coverage.zarr';group=zarr.open_group(str(store),mode='w');group.zeros('chr1',shape=(32,),dtype='u4',chunks=(8,))
            observed=contig_depth(('chr1',32,str(bam),'samtools',str(store)))
            np.testing.assert_array_equal(expected,group['chr1'][:]);self.assertEqual(observed[2],int(expected.sum()))

    def test_GPU생성Zarr는CPUfullQC와변조탐지를통과함(self):
        # 생성·패킹·CIGAR QC는 GPU 선택, 독립적인 semantic depth/sequence 비교는 CPU 구현임.
        from test_stage70_validation import Stage70ValidationTest
        case=Stage70ValidationTest('test_full_validation_recalculates_cigar_and_depth')
        result=unittest.TestResult()
        with patch.dict(os.environ,{'MFM_COVERAGE_BACKEND':'cupy','MFM_QC_BACKEND':'cupy','MFM_ZARR_BACKEND':'cupy'}):case.run(result)
        self.assertTrue(result.wasSuccessful(),str(result.errors)+str(result.failures))

if __name__=='__main__':unittest.main()
