"""헤더 정상·내부 손상, 캐시 갱신, 별도 경로 자동 복구 및 복구 한도를 검증한다."""
import json
import os
import subprocess
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch
from handoff import rna_integrity, cli, downloads

class IntegrityTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.bam = self.root/'original.bam'
        self.bam.write_bytes(b'bad_input')
        self.tool = self.root/'samtools'
        self.tool.write_bytes(b'executable')
        self.report = self.root/'validation.json'
    def tearDown(self):
        self.directory.cleanup()
    def decode(self,command,**kwargs):
        if 'flagstat' in command:
            path = Path(command[-1])
            if path.read_bytes().startswith(b'bad'):
                return subprocess.CompletedProcess(command,1,'','[E::bgzf_uncompress] CRC error')
            return subprocess.CompletedProcess(command,0,json.dumps({'QC-passed reads':{'total':5}}),'')
        return subprocess.CompletedProcess(command,0,'','')
    def contexts(self):
        return patch.object(rna_integrity.shutil,'which',return_value=str(self.tool))
    def test_internal_corruption_after_quickcheck_pass(self):
        with self.contexts(),patch.object(rna_integrity.subprocess,'run',side_effect=self.decode) as run:
            with self.assertRaises(rna_integrity.InputCorruptionError):
                rna_integrity.validate_full(self.bam,self.report,'HG00097')
        self.assertEqual(run.call_count,2)
        self.assertEqual(json.loads(self.report.read_text())['failed_check'],'full_decode')
    def test_cache_reused_only_for_identical_input(self):
        self.bam.write_bytes(b'good_data')
        old = time.time_ns()-5_000_000_000
        os.utime(self.bam,ns=(old,old))
        with self.contexts(),patch.object(rna_integrity.time,'time_ns',return_value=time.time_ns()+5_000_000_000),patch.object(rna_integrity.subprocess,'run',side_effect=self.decode) as run:
            rna_integrity.validate_full(self.bam,self.report,'HG00097')
            rna_integrity.validate_full(self.bam,self.report,'HG00097')
            self.assertEqual(run.call_count,2)
            self.bam.write_bytes(b'bad_input')
            with self.assertRaises(rna_integrity.InputCorruptionError):
                rna_integrity.validate_full(self.bam,self.report,'HG00097')
            self.assertEqual(run.call_count,4)
    def test_tool_or_sample_change_invalidates_cache(self):
        self.bam.write_bytes(b'good_data')
        with self.contexts(),patch.object(rna_integrity.time,'time_ns',return_value=time.time_ns()+5_000_000_000),patch.object(rna_integrity.subprocess,'run',side_effect=self.decode) as run:
            rna_integrity.validate_full(self.bam,self.report,'A')
            rna_integrity.validate_full(self.bam,self.report,'B')
            self.tool.write_bytes(b'new executable')
            rna_integrity.validate_full(self.bam,self.report,'B')
            self.assertEqual(run.call_count,6)
    def test_signal_failure_does_not_trigger_corruption_recovery(self):
        def run(command,**kwargs):
            return subprocess.CompletedProcess(command,-9,'','killed') if 'flagstat' in command else self.decode(command)
        with self.contexts(),patch.object(rna_integrity.subprocess,'run',side_effect=run):
            with self.assertRaises(RuntimeError):rna_integrity.validate_full(self.bam,self.report,'HG00097')
    def test_file_change_during_scan_rejected(self):
        self.bam.write_bytes(b'good_data')
        def run(command,**kwargs):
            result=self.decode(command)
            if 'flagstat' in command:self.bam.write_bytes(b'changed_data')
            return result
        with self.contexts(),patch.object(rna_integrity.subprocess,'run',side_effect=run):
            with self.assertRaises(RuntimeError):rna_integrity.validate_full(self.bam,self.report,'HG00097')
    def configuration(self):
        config={'data_root':str(self.root),'state_root':str(self.root/'state'),'sample_cpus':4}
        row={'sample_id':'HG00097','rna_path':'original.bam'}
        resources=[{'kind':'rna_alignment','path':'original.bam','url':'https://example.test/original.bam','size_bytes':9,'checksum':''},
                   {'kind':'rna_index','path':'original.bam.bai','url':'https://example.test/original.bam.bai','size_bytes':4,'checksum':''}]
        return config,row,resources
    def fetch(self,resource,*args,**kwargs):
        path=downloads.safe_path(self.root,resource['path'])
        if not path.exists():
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(b'good_data' if resource['kind']=='rna_alignment' else b'BAI1')
        return path
    def test_corrupt_original_preserved_replacement_bound_and_reused(self):
        config,row,resources=self.configuration()
        before=self.bam.read_bytes()
        with self.contexts(),patch.object(rna_integrity.subprocess,'run',side_effect=self.decode),patch.object(rna_integrity,'fetch',side_effect=self.fetch) as fetch,patch.object(rna_integrity,'fetch_or_index',side_effect=self.fetch):
            local,replaced=rna_integrity.prepare_rna(config,row,resources)
            relative=local['sample_rna_paths']['HG00097']
            self.assertNotEqual(relative,row['rna_path'])
            self.assertEqual(self.bam.read_bytes(),before)
            self.assertEqual(resources[0]['path'],'original.bam')
            self.assertNotIn('sample_rna_paths',config)
            self.assertEqual(replaced[1]['path'],relative+'.bai')
            again,reused=rna_integrity.prepare_rna(config,row,resources)
            self.assertEqual(again['sample_rna_paths']['HG00097'],relative)
            self.assertEqual(fetch.call_count,1) # 손상 기존 파일과 검증된 대체 파일은 다시 덮어쓰지 않는다.
        refs={'ncbi':'ncbi.fa','autosomes':'auto.fa','broad':'broad.fa','gtf':'genes.gtf','rna_extra':'extra.fa'}
        bound_row=dict(row,platform='hifi',alignment_path='wgs.bam',public_vcf_path='',data_group='wgs_rna_pair',project_id='HPRC_R2')
        local.update(output_root=str(self.root/'out'),run_id='test')
        with patch.object(cli.planner,'read_tsv',return_value=[]),patch.object(cli.bundled_assets,'panel_directory',return_value=self.root/'vcf'):
            _,bound=cli.bind_manifest(local,bound_row,replaced,refs)
        self.assertEqual(bound['rna_path'],str(self.root/relative))
    def test_persistent_corruption_has_bounded_retries(self):
        config,row,resources=self.configuration()
        def fetch_bad(resource,*args,**kwargs):
            path=downloads.safe_path(self.root,resource['path']);path.parent.mkdir(parents=True,exist_ok=True)
            if not path.exists():path.write_bytes(b'bad_input')
            return path
        with self.contexts(),patch.object(rna_integrity.subprocess,'run',side_effect=self.decode),patch.object(rna_integrity,'fetch',side_effect=fetch_bad) as fetch:
            with self.assertRaises(rna_integrity.InputCorruptionError):rna_integrity.prepare_rna(config,row,resources)
            self.assertEqual(fetch.call_count,2)
        self.assertEqual(self.bam.read_bytes(),b'bad_input')
        self.assertEqual(len(list((self.root/'rna_recovery').rglob('*.bam'))),2)
    def test_bad_size_original_is_preserved(self):
        config,row,resources=self.configuration()
        self.bam.write_bytes(b'short')
        with self.contexts(),patch.object(rna_integrity.subprocess,'run',side_effect=self.decode),patch.object(rna_integrity,'fetch',side_effect=self.fetch),patch.object(rna_integrity,'fetch_or_index',side_effect=self.fetch):
            local,_=rna_integrity.prepare_rna(config,row,resources)
        self.assertEqual(self.bam.read_bytes(),b'short')
        self.assertEqual((self.root/local['sample_rna_paths']['HG00097']).read_bytes(),b'good_data')
    def test_memory_error_is_not_redownloaded(self):
        config,row,resources=self.configuration()
        def run(command,**kwargs):
            return subprocess.CompletedProcess(command,1,'','[E::bgzf_read] Cannot allocate memory')
        with self.contexts(),patch.object(rna_integrity.subprocess,'run',side_effect=run),patch.object(rna_integrity,'fetch') as fetch:
            with self.assertRaises(RuntimeError):rna_integrity.prepare_rna(config,row,resources)
            fetch.assert_not_called()
    def test_real_internal_bgzf_corruption_recovers_from_same_source(self):
        import shutil,struct
        if not shutil.which('samtools'):self.skipTest('samtools가 필요합니다.')
        import pysam
        source=self.root/'public.bam'
        with pysam.AlignmentFile(str(source),'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000000}]}) as output:
            for number in range(1000):
                read=pysam.AlignedSegment();read.query_name='read'+str(number)
                read.query_sequence='ACGT'*50;read.flag=0;read.reference_id=0;read.reference_start=number*10
                read.mapping_quality=60;read.cigarstring='200M';read.query_qualities=pysam.qualitystring_to_array('I'*200)
                output.write(read)
        pysam.index(str(source))
        data=bytearray(source.read_bytes());blocks=[];offset=0
        while offset<len(data):
            blocks.append(offset);offset+=struct.unpack_from('<H',data,offset+16)[0]+1
        self.assertGreaterEqual(len(blocks),3)
        data[blocks[1]+20]^=1;self.bam.write_bytes(data)
        quick=subprocess.run(['samtools','quickcheck','-v',str(self.bam)],capture_output=True,text=True)
        self.assertEqual(quick.returncode,0,quick.stderr)
        with self.assertRaises(rna_integrity.InputCorruptionError):
            rna_integrity.validate_full(self.bam,self.report,'HG00097')
        config,row,resources=self.configuration()
        resources[0].update(url=source.as_uri(),size_bytes=source.stat().st_size)
        resources[1].update(url=Path(str(source)+'.bai').as_uri(),size_bytes=Path(str(source)+'.bai').stat().st_size)
        local,replaced=rna_integrity.prepare_rna(config,row,resources)
        self.assertEqual(self.bam.read_bytes(),data)
        new=self.root/local['sample_rna_paths']['HG00097']
        self.assertEqual(new.read_bytes(),source.read_bytes())
        self.assertTrue(Path(str(new)+'.bai').is_file())
        self.assertEqual(json.loads(Path(local['rna_full_decode_report']).read_text())['status'],'pass')

    def test_wgs_only_skips_rna_validation(self):
        config,row,resources=self.configuration();row['rna_path']=''
        with patch.object(rna_integrity,'fetch') as fetch:
            actual=rna_integrity.prepare_rna(config,row,resources)
            self.assertEqual(actual,(config,resources));fetch.assert_not_called()

if __name__=='__main__':unittest.main()
