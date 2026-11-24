"""실제 BAM/BWA 대조와 단일 lease의 합산·반환·실패 정리를 검증함."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from dataclasses import replace
import json,os,sys,tempfile,unittest
import pysam
import test_prepared_cache as fixture
from integrate_wgs_rna import personal_realign as personal
from integrate_wgs_rna.realign import gpu_prepare_overlap as module
from integrate_wgs_rna.realign.phases import RealignPhases,_atomic
from integrate_wgs_rna.gpu_bwa_index import lookup
from integrate_wgs_rna.memory_policy import GIB

ROOT=Path(__file__).resolve().parents[1]

class OverlapTests(unittest.TestCase):
    def setUp(self):
        self.fixture=fixture.PreparedCacheTests();self.fixture.setUp();self.addCleanup(self.fixture.doCleanups)
        self.root=self.fixture.root;self.cfg=replace(self.fixture.config,route_jobs=2)
        self.work=self.root/'overlap';self.work.mkdir();self.realigned=self.root/'realigned';self.logs=self.root/'logs';self.logs.mkdir()
        self.haps={h:self.root/f'HAP{h}' for h in [1,2]}
        for h,d in self.haps.items():d.mkdir();(d/'chr1.fa').write_text('>chr1\n'+('ACGT' if h==1 else 'TGCA')*500+'\n');pysam.faidx(str(d/'chr1.fa'))
        self.pool=self.root/'pool';(self.pool/'leases').mkdir(parents=True);(self.pool/'waiters').mkdir()
        self.lease=self.pool/'leases/own.json';_atomic(self.lease,dict(manager_pid=os.getpid(),child_pid=None,cpus=8,total_cpus=8,memory_bytes=4*GIB,branch='WGS',stage='10',sample='S'))
        _atomic(self.pool/'policy.json',dict(cpu_pool_limit=8,memory_limit_bytes=64*GIB,cpu_backfill={'enabled':False}))
        self.binary=Path(os.environ.get('MFM_GPU_BWA',str(ROOT/'tools/gpu/bwa')))
        if not self.binary.is_file():self.skipTest('실제 BWA가 필요한 검사임')
        self.env=patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(self.lease),'MFM_CPU_BUDGET':'8','MFM_CPU_POOL_GRANTED':'8',
            'MFM_RAM_LIMIT_BYTES':str(64*GIB),'MFM_WGS_BACKEND':'parabricks','MFM_GPU_BWA':str(self.binary),
            'PYTHONPATH':os.pathsep.join([str(ROOT),str(ROOT/'src'),str(ROOT/'tests')]),'MFM_GPU_COMPARE_UNRESERVED':'0',
            'MFM_GPU_SCRATCH_ROOT':'','MFM_GPU_BWA_INDEX_ROOT':'','MFM_GPU_BAM_INDEX_OVERLAP':'1',
            'MFM_WGS_HAP_REQUESTED':'1','MFM_WGS_ALIGN_REQUESTED':'auto','MFM_WGS_SORT_REQUESTED':'auto','MFM_WGS_MIN_ALIGN_THREADS':'1'})
        self.env.start();self.addCleanup(self.env.stop)

    def run_prepare(self,metrics,phase=None):
        # KB 규모 fixture만350GB 디스크 계획 하한을 우회하며 RAM/CPU lease는 실제 검사함.
        with patch.object(module,'_disk_fits',return_value=True):
            return module.prepare(self.cfg,'S',self.fixture.source,['chr1'],self.haps,self.realigned,self.work,self.logs,[1,2],phase or RealignPhases(),metrics)

    def test_실제BAM의모든read_태그_통계가순차와같고인덱스가재사용됨(self):
        baseline=self.fixture.prepare(self.cfg);metrics={};result=self.run_prepare(metrics)
        def reads(path):
            with pysam.AlignmentFile(str(path),'rb') as stream:return sorted(r.to_string() for r in stream)
        self.assertEqual(reads(baseline.path),reads(result.path));self.assertEqual(baseline.input_counts,result.input_counts)
        for key in baseline.stats:
            if key not in {'input_metrics_execution','prepared_cache'}:self.assertEqual(baseline.stats[key],result.stats[key],key)
        self.assertEqual(metrics['status'],'pass');self.assertEqual(metrics['initial_reserved_cpus'],8)
        self.assertEqual(metrics['input_cpus'],5);self.assertEqual(metrics['index_cpus'],2)
        self.assertLessEqual(metrics['initial_reserved_memory_bytes'],64*GIB*9//10)
        for hap in [1,2]:self.assertIsNotNone(lookup(self.realigned/f'HAP{hap}/target.fa',self.work))
        row=json.loads(self.lease.read_text());self.assertEqual(row['cpus'],1);self.assertEqual(row['memory_bytes'],4*GIB);self.assertEqual(row['prepare_components'],{})
        for record in metrics['runs'].values():self.assertEqual(record['exit_code'],0);self.assertFalse(module._group_members(record['pid']))
        warm={};self.assertIsNone(self.run_prepare(warm));self.assertEqual(warm['mode'],'indexes_already_cached')
        self.assertEqual(json.loads(self.lease.read_text())['cpus'],1)

    def test_BAM이먼저끝나면인덱스2개와제어1CPU만유지함(self):
        # /tmp는 noexec일 수 있어 실행 대역만 패키지의 임시 폴더에 둠.
        executable=tempfile.TemporaryDirectory(prefix='.overlap_exec_',dir=ROOT);self.addCleanup(executable.cleanup)
        binary=Path(executable.name)/'wait_bwa';meta=self.logs/'S.prepare_overlap.json'
        code='#!'+sys.executable+'\nfrom pathlib import Path\nimport os,sys,json,time\n'
        code+='until=time.monotonic()+15\nwhile time.monotonic()<until:\n'
        code+='    r=json.loads(Path('+repr(str(meta))+').read_text())\n'
        code+="    if 'input' in r.get('runs',{}):break\n    time.sleep(.02)\nelse:raise SystemExit(8)\n"
        code+='os.execv('+repr(str(self.binary))+',['+repr(str(self.binary))+',*sys.argv[1:]])\n'
        binary.write_text(code);binary.chmod(0o755);phases=RealignPhases();metrics={}
        with patch.dict(os.environ,{'MFM_GPU_BWA':str(binary)}):self.run_prepare(metrics,phases)
        states=[r for r in phases.events if set(r.get('active_components',[]))=={'index_1','index_2'}]
        self.assertTrue(states);self.assertTrue(all(r['reserved_cpus']==3 for r in states))
        self.assertTrue(metrics['overlap_observed'])

    def test_다른샘플의RAM때문에동시피크가안맞으면순차경로로돌아감(self):
        for h in [1,2]:personal.build_target_fasta(self.cfg,self.haps[h],self.realigned/f'HAP{h}/target.fa',['chr1'],'chr1')
        _atomic(self.pool/'leases/other.json',dict(cpus=1,total_cpus=8,memory_bytes=50*GIB,manager_pid=os.getppid(),child_pid=None))
        metrics={}
        with patch.object(module.subprocess,'Popen') as popen:
            result=module.prepare(self.cfg,'S',self.fixture.source,['chr1'],self.haps,self.realigned,self.work,self.logs,[1,2],RealignPhases(),metrics)
        self.assertIsNone(result);popen.assert_not_called();self.assertEqual(metrics['mode'],'memory_fallback')
        self.assertEqual(json.loads(self.lease.read_text())['cpus'],1)

    def test_실패시다른자식도회수한뒤자원을반환함(self):
        fail=self.root/'bwa_fail';fail.write_text('#!/bin/sh\nexit 7\n');fail.chmod(0o755)
        metrics={}
        with patch.dict(os.environ,{'MFM_GPU_BWA':str(fail)}):
            with self.assertRaises(RuntimeError):self.run_prepare(metrics)
        self.assertEqual(metrics['status'],'failed')
        row=json.loads(self.lease.read_text());self.assertEqual(row['cpus'],1);self.assertEqual(row['memory_bytes'],4*GIB)
        for record in metrics['runs'].values():self.assertFalse(module._group_members(record['pid']))
        self.assertFalse(list(self.work.glob('bwa_indexes/*/complete.json')))

    def test_디스크부족은실행없이자원반환후순차경로로돌아감(self):
        metrics={}
        for h in [1,2]:personal.build_target_fasta(self.cfg,self.haps[h],self.realigned/f'HAP{h}/target.fa',['chr1'],'chr1')
        with patch.object(module,'_disk_fits',return_value=False),patch.object(module.subprocess,'Popen') as popen:
            result=module.prepare(self.cfg,'S',self.fixture.source,['chr1'],self.haps,self.realigned,self.work,self.logs,[1,2],RealignPhases(),metrics)
        self.assertIsNone(result);popen.assert_not_called();self.assertEqual(metrics['mode'],'disk_fallback')
        self.assertEqual(json.loads(self.lease.read_text())['cpus'],1)

    def test_작은CPU풀에서는기존경로를유지함(self):
        phase=RealignPhases();phase.total=6
        with patch.object(module.subprocess,'Popen') as popen:
            self.assertIsNone(module.prepare(self.cfg,'S',self.fixture.source,['chr1'],self.haps,self.realigned,self.work,self.logs,[1,2],phase,{}))
        popen.assert_not_called()

    def test_호출부는중첩준비BAM을재생성하지않고그대로사용함(self):
        prepared=self.fixture.prepare(self.cfg);sample=self.root/'S';(sample/'haplotagged').mkdir(parents=True)
        expected=sample/'haplotagged/S.ont.hg38.haplotagged.bam';expected.symlink_to(self.fixture.source)
        Path(str(expected)+'.bai').symlink_to(Path(str(self.fixture.source)+'.bai'))
        def early(*args):args[-1].update(status='pass');return prepared
        cfg=replace(self.cfg,hap1_fasta_dir=self.haps[1],hap2_fasta_dir=self.haps[2])
        with (patch.object(module,'prepare',side_effect=early),patch.object(personal,'prepare_illumina_paired_input',side_effect=AssertionError('BAM 재생성')),
              patch.object(personal,'needs_illumina_paired_realign',return_value=True),patch.object(personal,'run_one_haplotype',return_value={'status':'OK','output':'mock','timing_sec':{}})):
            row=personal.process_one_sample(cfg,'S',['chr1'])
        self.assertEqual(row['status'],'OK',row.get('message'))

    def test_명시적으로끄면기존BAM준비경로를사용함(self):
        sample=self.root/'S';(sample/'haplotagged').mkdir(parents=True)
        expected=sample/'haplotagged/S.ont.hg38.haplotagged.bam';expected.symlink_to(self.fixture.source)
        Path(str(expected)+'.bai').symlink_to(Path(str(self.fixture.source)+'.bai'))
        cfg=replace(self.cfg,hap1_fasta_dir=self.haps[1],hap2_fasta_dir=self.haps[2])
        with (patch.dict(os.environ,{'MFM_GPU_BAM_INDEX_OVERLAP':'0'}),
              patch.object(module,'prepare',side_effect=AssertionError('비활성 경로에서 중첩 호출')),
              patch.object(personal,'prepare_illumina_paired_input',wraps=personal.prepare_illumina_paired_input) as original,
              patch.object(personal,'needs_illumina_paired_realign',return_value=True),
              patch.object(personal,'run_one_haplotype',return_value={'status':'OK','output':'mock','timing_sec':{}})):
            row=personal.process_one_sample(cfg,'S',['chr1'])
        self.assertEqual(row['status'],'OK',row.get('message'));self.assertEqual(original.call_count,1)

class BudgetTests(unittest.TestCase):
    def test_겹친시간은빈틈을제외하고중복합산하지않음(self):
        row=lambda a,b:{'_work_started_at':a,'_work_finished_at':b}
        self.assertEqual(module.overlap_seconds(row(1,10),[row(0,3),row(7,12)]),5)
        self.assertEqual(module.overlap_seconds(row(1,10),[row(0,7),row(4,12)]),9)
        self.assertEqual(module.overlap_seconds(row(1,10),[]),0)

    def test_실행중인요소는보존하고끝난요소만반환함(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir();path=root/'leases/own.json'
            phase=SimpleNamespace(path=path,events=[],unreserved=False)
            _atomic(path,dict(cpus=8,memory_bytes=20*GIB,realign_phase='gpu_bam_index_prepare'))
            index=dict(cpus=1,memory_bytes=5*GIB,write_paths=[],disk_peak_bytes=0)
            module.commit_components(phase,{'index_1':index,'index_2':index})
            row=json.loads(path.read_text());self.assertEqual(row['cpus'],3);self.assertEqual(row['memory_bytes'],14*GIB)
            with self.assertRaises(ValueError):module.commit_components(phase,{'input':dict(cpus=8,memory_bytes=16*GIB)})
            module.commit_components(phase,{})
            self.assertEqual(json.loads(path.read_text())['cpus'],1)

    def test_디스크점유는다른샘플의중첩예약과GPU예약을합산함(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);part=dict(write_paths=[str(root/'pending')],disk_peak_bytes=100)
            with patch.object(module.os,'statvfs',return_value=SimpleNamespace(f_bavail=350,f_frsize=1)):
                self.assertTrue(module._disk_fits({'input':part},[],100))
                self.assertFalse(module._disk_fits({'input':part},[{'prepare_components':{'a':part},'hap_resources':{'1':part}}],100))

    def test_고정기초CPU가부모grant에포함됨(self):
        with patch.dict(os.environ,{'MFM_CPU_BUDGET':'8','MFM_RAM_LIMIT_BYTES':str(64*GIB)},clear=True):
            phases=RealignPhases();grant=phases.reserve('test',max_units=5,min_units=4,base_cpus=3,memory_per_unit=2*GIB)
            self.assertEqual(grant,5);self.assertEqual(phases.resource_limits['cpus'],8)
            with self.assertRaises(ValueError):phases.reserve('over',max_units=8,min_units=6,base_cpus=3)

if __name__=='__main__':unittest.main()
