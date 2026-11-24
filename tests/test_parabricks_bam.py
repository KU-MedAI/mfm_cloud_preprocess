"""GPU 모드와 실제 Parabricks FASTQ·GPU duplicate 동등성을 검사함."""
from pathlib import Path
from unittest.mock import patch
import collections
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import pysam
from handoff import acceleration
from integrate_wgs_rna import parabricks_bam as module
from integrate_wgs_rna.parabricks_runtime import Runtime

ROOT=Path(__file__).resolve().parents[1]
PB=os.environ.get('MFM_TEST_PBRUN',str(ROOT/'tools/gpu/parabricks_4.8.0-1/pbrun-native'))


class ParabricksBamTest(unittest.TestCase):
    def fixture(self,root):
        ref=root/'ref.fa';ref.write_text('>chr1\n'+'A'*10000+'\n');pysam.faidx(str(ref))
        bam=root/'input.bam';header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':10000}],
            'RG':[{'ID':'S','SM':'S','LB':'L','PL':'ILLUMINA','PU':'unit'}]}
        records=[]
        for n,start in enumerate([50,50,300,500]):
            for flag,pos,mate in [(99,start,start+100),(147,start+100,start)]:
                read=pysam.AlignedSegment();read.query_name=f'INST:1:FLOW:1:1101:{1000+n}:1000';read.query_sequence='A'*40
                read.query_qualities=pysam.qualitystring_to_array('I'*40);read.flag=flag|(512 if n==3 else 0)
                read.reference_id=0;read.reference_start=pos;read.mapping_quality=255;read.cigar=[(0,40)]
                read.next_reference_id=0;read.next_reference_start=mate;read.template_length=140 if flag==99 else -140
                read.set_tag('RG','S');read.set_tag('HP',1);read.set_tag('PS',5);records.append(read)
        for name,flag,pos in [('orphan',73,800),('single',0,900),('secondary',256,1000),('supplementary',2048,1100)]:
            read=pysam.AlignedSegment();read.query_name=name;read.query_sequence='C'*40
            read.query_qualities=pysam.qualitystring_to_array('J'*40);read.flag=flag;read.reference_id=0;read.reference_start=pos
            read.mapping_quality=30;read.cigar=[(0,40)];read.set_tag('RG','S');records.append(read)
        with pysam.AlignmentFile(str(bam),'wb',header=header) as out:
            for read in sorted(records,key=lambda r:r.reference_start):out.write(read)
        return ref,bam

    def fastq(self,path):
        records=[]
        with Path(path).open('rb') as stream:
            while first:=stream.readline():
                seq=stream.readline();plus=stream.readline();qual=stream.readline()
                self.assertTrue(plus.startswith(b'+'));records.append((first.strip(),seq.strip(),qual.strip()))
        return collections.Counter(records)

    def test_GPU모드와기존GPU설정은BAM도구자동선택_CPU개별선택보존(self):
        config=acceleration.apply({'acceleration':{'mode':'gpu','rna_gpu_policy':'experimental'}})
        self.assertEqual(config['pipeline_env']['MFM_BAM_BACKEND'],'parabricks')
        self.assertEqual(config['pipeline_env']['MFM_WGS_BACKEND'],'parabricks')
        self.assertEqual(config['pipeline_env']['MFM_COVERAGE_BACKEND'],'cupy')
        config=acceleration.apply({'acceleration':{'wgs_backend':'parabricks','bam_backend':'cpu'}})
        self.assertEqual(config['pipeline_env']['MFM_BAM_BACKEND'],'cpu')
        config=acceleration.apply({'acceleration':{'wgs_backend':'parabricks'}})
        self.assertEqual(config['pipeline_env']['MFM_BAM_BACKEND'],'parabricks')
        with patch.dict(os.environ,{'MFM_WGS_BACKEND':'parabricks'},clear=True):self.assertTrue(module.enabled())
        with patch.dict(os.environ,{},clear=True):self.assertFalse(module.enabled())

    def test_CPU제공도구는GPU미요청이며실제containerRAM_CPU상한사용(self):
        with patch.dict(os.environ,{'MFM_GPU_RUNTIME':'docker'}),patch('shutil.which',return_value='/bin/docker'):
            command=Runtime().command(['bam2fq'],None,host_limit=123456789,cpu_limit=2)
        self.assertNotIn('--gpus',command)
        self.assertEqual(command[command.index('--memory')+1],'123456789')
        self.assertEqual(command[command.index('--cpus')+1],'2')
        self.assertIn('NVIDIA_VISIBLE_DEVICES=none',command)

    def test_WGS만선택한GPU모드는사용하지않는RNAstrict검사를요구하지않음(self):
        with patch.dict(os.environ,{},clear=True),patch('integrate_wgs_rna.gpu_pool.inventory',return_value=[]),\
             patch('integrate_wgs_rna.gpu_pool.allowed_devices',return_value=[{}]),\
             patch('importlib.util.find_spec',return_value=True),patch.object(Runtime,'validate'),\
             patch('integrate_wgs_rna.gpu_rna.validate_policy',side_effect=AssertionError('RNA 불필요')):
            acceleration.preflight({'acceleration':{'mode':'gpu'}},[{'data_group':'wgs'}])

    def test_GPU_BAM중복표시는두CPU이상과공유풀전체상한을사용함(self):
        import importlib.util
        spec=importlib.util.spec_from_file_location('bam_pool',ROOT/'scripts/cpu_pool.py')
        pool=importlib.util.module_from_spec(spec);spec.loader.exec_module(pool)
        for layout,stage in [('original','17'),('legacy','15')]:
            with patch.dict(os.environ,{'MFM_CPU_POOL_LAYOUT':layout,'MFM_CPU_POOL_PROFILE':'throughput','MFM_BAM_BACKEND':'parabricks'},clear=True):
                self.assertEqual(pool.stage_cap('RNA',stage,68),(2,68))

    def test_GPU없는실제자식실행도정확한메타데이터기록(self):
        with tempfile.TemporaryDirectory() as raw,patch.dict(os.environ,{'MFM_GPU_RUNTIME':'native','MFM_PBRUN':sys.executable}):
            result=Runtime().run(['-c','import os;assert os.environ["CUDA_VISIBLE_DEVICES"]==""'],None,
                                 host_limit=128*1024**2,cpu_limit=1,log=Path(raw)/'run.log')
        self.assertIsNone(result['gpu_uuid']);self.assertEqual(result['gpu_reserved_bytes'],0)

    def test_FASTQ손상과mate이상은차단함(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);path=root/'reads.gz'
            with gzip.open(path,'wb') as out:out.write(b'@read/1\nAC\n+\nII\n')
            self.assertEqual(module.unpack_fastq(path,root/'out',1),1)
            self.assertEqual((root/'out').read_bytes(),b'@read\nAC\n+\nII\n')
            with self.assertRaisesRegex(ValueError,'suffix'):module.unpack_fastq(path,root/'out',2)
            with gzip.open(path,'wb') as out:out.write(b'@read/1\nAC\n+\nI\n')
            with self.assertRaisesRegex(ValueError,'손상'):module.unpack_fastq(path,root/'out',1)

    def test_GPU검증실패_코드변경_미완료는실제출력생성전차단(self):
        import hashlib
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);marker=root/'validation.json'
            env={'MFM_PARABRICKS_MARKDUP_VALIDATION':str(marker),'MFM_GPU_WAIT_SECONDS':'0'}
            with patch.dict(os.environ,env):
                with self.assertRaises(TimeoutError):module.mark_duplicates(None,None,None,None,None,None)
                for status,digest in [('failed',hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()),('pass','stale')]:
                    marker.write_text(json.dumps({'status':status,'backend_sha256':digest}))
                    with self.assertRaisesRegex(RuntimeError,'검증 미통과'):module.mark_duplicates(None,None,None,None,None,None)

    def test_GPU검증진행표식은완료까지대기후판정함(self):
        import hashlib
        with tempfile.TemporaryDirectory() as raw:
            marker=Path(raw)/'validation.json'
            marker.write_text(json.dumps({'status':'running'}))
            def finish(_):
                marker.write_text(json.dumps({'status':'pass','backend_sha256':hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()}))
            with patch.dict(os.environ,{'MFM_PARABRICKS_MARKDUP_VALIDATION':str(marker),'MFM_GPU_WAIT_SECONDS':'3'}),patch.object(module.time,'sleep',side_effect=finish),patch.object(module,'resource_plan',side_effect=ValueError('검증 후 자원진입')):
                with self.assertRaisesRegex(ValueError,'검증 후 자원진입'):module.mark_duplicates(None,None,None,None,None,None)

    def test_GPUwriter중복HD는read와HP_PS_순서그대로정상화함(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);ref,bam=self.fixture(root);header=root/'bad.sam';bad=root/'bad.bam'
            with pysam.AlignmentFile(str(bam)) as handle:
                header.write_text('@HD\tVN:1.6\tSO:queryname\n'+handle.text)
                before=[r.to_string() for r in handle]
            with bad.open('wb') as out:subprocess.run(['samtools','reheader','-P',str(header),str(bam)],stdout=out,check=True)
            with pysam.AlignmentFile(str(bad)) as handle:
                with self.assertRaises(ValueError):handle.header.to_dict()
            info=module.normalize_markdup_header(bad,root,root/'header.log')
            with pysam.AlignmentFile(str(bad)) as handle:
                self.assertEqual(handle.header.to_dict()['HD']['SO'],'coordinate')
                self.assertEqual([r.to_string() for r in handle],before)
            self.assertEqual(info['duplicate_hd_lines_removed'],1)

    def test_GPU중복표시인자는기존coordinate_baseline과광학거리보존(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);ref,bam=self.fixture(root);output=root/'marked.bam';metrics=root/'metrics.txt';calls=[]
            def run(args,paths,log,plan,*,gpu):
                calls.append((args,gpu))
                if args[0]=='markdup':shutil.copyfile(bam,output);metrics.write_text('LIBRARY\tREAD_PAIR_DUPLICATES\nL\t1\n')
                return {'tool':args[0]}
            with patch.dict(os.environ,{'MFM_CPU_BUDGET':'4','MFM_CPU_POOL_GRANTED':'4'},clear=True),patch.object(module,'run_tool',side_effect=run):
                result=module.mark_duplicates(bam,ref,output,metrics,root,root/'log')
            self.assertEqual([c[0][0] for c in calls],['bamsort','markdup'])
            self.assertTrue(all(c[1] for c in calls))
            self.assertIn('--gpusort',calls[0][0]);self.assertIn('--gpuwrite',calls[1][0])
            self.assertNotIn('--markdups-assume-sortorder-queryname',calls[1][0])
            self.assertEqual(calls[1][0][calls[1][0].index('--optical-duplicate-pixel-distance')+1],'100')
            self.assertFalse(result['remove_duplicates'])

    @unittest.skipUnless(Path(PB).is_file(),'설치된 Parabricks native 필요함')
    def test_실제bam2fq는기존read_이름_서열_품질_제외수를보존(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);ref,bam=self.fixture(root);parts={k:root/(k+'.fq') for k in ['r1','r2','un','sg']}
            env={'MFM_PBRUN':PB,'MFM_GPU_RUNTIME':'native','MFM_CPU_BUDGET':'2','MFM_CPU_POOL_GRANTED':'2'}
            with patch.dict(os.environ,env,clear=True):info=module.rna_fastq(bam,None,parts,root/'pb.log')
            col=subprocess.Popen(['samtools','collate','-u','-O',str(bam)],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
            cpu={k:root/('cpu_'+k+'.fq') for k in parts}
            result=subprocess.run(['samtools','fastq','-1',str(cpu['r1']),'-2',str(cpu['r2']),'-0',str(cpu['un']),'-s',str(cpu['sg']),'-n','-'],stdin=col.stdout,capture_output=True)
            col.stdout.close();self.assertEqual(col.wait(),0);self.assertEqual(result.returncode,0)
            for key in parts:self.assertEqual(self.fastq(parts[key]),self.fastq(cpu[key]),key)
            self.assertFalse(info['gpu_kernels_requested']);self.assertIsNone(info['gpu_uuid'])


@unittest.skipUnless(os.environ.get('MFM_RUN_PARABRICKS_GPU_TESTS')=='1','실제 GPU 반환 후 별도 실행함')
class ActualGpuMarkdupTest(unittest.TestCase):
    fixture=ParabricksBamTest.fixture
    def test_실제GPUmarkdup의read와duplicateflag가Picard와동일(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'tests',prefix='.tmp_pb_markdup_') as raw:
            root=Path(raw);ref,bam=self.fixture(root);gpu=root/'gpu.bam';metrics=root/'gpu.metrics.txt'
            env={'MFM_PBRUN':PB,'MFM_GPU_RUNTIME':'native','MFM_CPU_BUDGET':'4','MFM_CPU_POOL_GRANTED':'4','MFM_RAM_LIMIT_BYTES':str(16*1024**3)}
            try:
                with patch.dict(os.environ,env):module.mark_duplicates(bam,ref,gpu,metrics,root,root/'pb.log')
            except BaseException:
                if (root/'pb.log').exists():print((root/'pb.log').read_text(errors='replace')[-16000:],file=sys.stderr,flush=True)
                raise
            cpu=root/'cpu.bam'
            subprocess.run(['picard','-Xmx1g','-XX:ActiveProcessorCount=2',f'-Djava.io.tmpdir={root}',
                'MarkDuplicates',f'I={bam}',f'O={cpu}',f'M={root}/cpu.metrics.txt',f'TMP_DIR={root}',
                'PROGRAM_RECORD_ID=null','ASSUME_SORT_ORDER=coordinate','TAGGING_POLICY=DontTag','OPTICAL_DUPLICATE_PIXEL_DISTANCE=100'],check=True,capture_output=True)
            def rows(path):
                # PG는 도구 실행 이력이며 Picard PROGRAM_RECORD_ID=null과 Parabricks가 다름.
                # 다른 모든 read 필드·부가태그·duplicate flag는 그대로 전량 비교함.
                with pysam.AlignmentFile(str(path)) as handle:
                    return sorted((r.query_name,r.flag,r.reference_id,r.reference_start,r.mapping_quality,r.cigarstring,
                                   r.next_reference_id,r.next_reference_start,r.template_length,r.query_sequence,r.qual,
                                   tuple(sorted((k,v) for k,v in r.get_tags() if k!='PG'))) for r in handle)
            self.assertEqual(rows(cpu),rows(gpu))


if __name__=='__main__':unittest.main()
