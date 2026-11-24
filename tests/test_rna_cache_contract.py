"""작은 임시 자료로 RNA 입력 변경·재사용·강제 재생성 계약을 확인한다."""
from __future__ import annotations
import importlib.util, json, os, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
import pysam

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from rna_cache import cache_fields, cache_valid, file_fingerprint

def module(stage):
    stems={'11':'rna_to_fastq','12':'lift_annotation','13':'star_index','14':'star_align','15':'mark_duplicates'}
    stem=stems[stage]
    if (ROOT/'09_personal_fasta.sh').is_file():stage=f'{int(stage)+2:02d}'
    path=ROOT/f'{stage}_{stem}.py'
    spec=importlib.util.spec_from_file_location(f'stage{stage}',path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

def invoke(mod,args):
    with patch.object(sys,'argv',[mod.__file__,*map(str,args)]):return mod.main()

def write_bam(path,base='A'):
    header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}],'RG':[{'ID':'S','SM':'S'}]}
    with pysam.AlignmentFile(str(path),'wb',header=header) as out:
        for flag,pos,mate in [(99,100,200),(147,200,100)]:
            r=pysam.AlignedSegment(out.header);r.query_name='pair';r.query_sequence=base*40;r.flag=flag
            r.reference_id=0;r.reference_start=pos;r.mapping_quality=60;r.cigar=[(0,40)]
            r.next_reference_id=0;r.next_reference_start=mate;r.template_length=140 if flag==99 else -140
            r.query_qualities=pysam.qualitystring_to_array('I'*40);r.set_tag('RG','S');out.write(r)

class CacheContract(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='mfm_cache_test_');self.root=Path(self.temp.name)
        self.hap=self.root/'HAP1';self.hap.mkdir();self.work=self.root/'work';self.work.mkdir()
        (self.hap/'chr1.fa').write_text('>chr1\n'+'A'*100+'\n');pysam.faidx(str(self.hap/'chr1.fa'))
        self.gtf=self.root/'g.gtf';self.gtf.write_text('##genome-build GRCh38\nchr1\tT\texon\t55\t70\t.\t+\t.\tgene_id "G"; transcript_id "T";\n')
        self.qc=self.root/'fastq_qc.json';self.r1=self.root/'R1.fq';self.r2=self.root/'R2.fq'
        for p in [self.r1,self.r2]:p.write_text('@R\nAAAA\n+\nIIII\n')
        self.qc.write_text(json.dumps({'max_read_length':4,'r1_fastq':str(self.r1),'r2_fastq':str(self.r2)}))
        self.index=module('13');self.align=module('14');self.fake_calls=[]
    def tearDown(self):self.temp.cleanup()
    def fake_run(self,cmd,**kwargs):
        self.fake_calls.append(cmd)
        if 'faidx' in cmd:pysam.faidx(cmd[-1]);return
        if 'genomeGenerate' in cmd:
            d=Path(cmd[cmd.index('--genomeDir')+1]);d.mkdir(exist_ok=True)
            for name in ['Genome','SA','SAindex','chrName.txt','chrLength.txt','chrStart.txt','genomeParameters.txt']:(d/name).write_text('완료\n')
        elif 'alignReads' in cmd:
            d=Path(cmd[cmd.index('--outFileNamePrefix')+1]);d.mkdir(exist_ok=True)
            (d/'Aligned.sortedByCoord.out.bam').write_bytes('임시BAM'.encode())
            (d/'Log.final.out').write_text('Uniquely mapped reads % | 90.0%\n')
            (d/'ReadsPerGene.out.tab').write_text('G\t1\t1\t1\n')
        else:raise AssertionError(cmd)
    def index_args(self):
        return ['--sample','S','--hap','1','--hap-fasta-dir',self.hap,'--contigs','chr1','--gtf',self.gtf,'--fastq-qc',self.qc,'--work-dir',self.work,'--log-dir',self.root/'logs','--cache-tag','test','--star',sys.executable,'--samtools',sys.executable]
    def build_index(self,*extra):
        with patch.object(self.index,'run',self.fake_run):invoke(self.index,[*self.index_args(),*extra])
        return json.loads((self.work/'hap1.index.json').read_text())
    def build_align(self,*extra):
        args=['--sample','S','--hap','1','--index-json',self.work/'hap1.index.json','--fastq-qc',self.qc,'--work-dir',self.work,'--out-dir',self.root/'output','--log-dir',self.root/'logs','--star',sys.executable,*extra]
        with patch.object(self.align,'run',self.fake_run):invoke(self.align,args)
        return json.loads((self.work/'hap1.align.json').read_text())
    def test_stat_input_same_size_mtime_restored_invalidates(self):
        record=cache_fields({'source':file_fingerprint(self.r1)},[self.r2]);before=self.r1.stat()
        self.r1.write_text('@R\nCCCC\n+\nIIII\n');os.utime(self.r1,ns=(before.st_atime_ns,before.st_mtime_ns))
        self.assertFalse(cache_valid(record))
    def test_large_file_probe_detects_same_size_edit(self):
        p=self.root/'large.fa';p.write_bytes(b'A'*(2<<20));record=cache_fields({'source':file_fingerprint(p)},[self.r2]);before=p.stat()
        with p.open('r+b') as f:f.write(b'C')
        os.utime(p,ns=(before.st_atime_ns,before.st_mtime_ns))
        self.assertFalse(cache_valid(record))
    def test_lift_reuse_detects_chain_and_output_change(self):
        chain=self.hap/'chr1.chain';head='chain 1 chr1 100 + 0 100 chr1 100 + 0 100 1\n';chain.write_text(head+'100\n')
        lift=module('12');out=self.root/'lift.gtf';qc=self.root/'lift.json'
        args=['--chain-dir',self.hap,'--gtf',self.gtf,'--contigs','chr1','--reference-fai',self.hap/'chr1.fa.fai','--out-gtf',out,'--out-exons',self.root/'exons.tsv','--qc-json',qc,'--reuse']
        invoke(lift,args);mtime=out.stat().st_mtime_ns;invoke(lift,args);self.assertEqual(mtime,out.stat().st_mtime_ns)
        chain.write_text(head+'40 1 0\n19 0 1\n40\n');invoke(lift,args);self.assertIn('\texon\t54\t70\t',out.read_text())
        out.write_text('외부 변경\n');invoke(lift,args);self.assertIn('\texon\t54\t70\t',out.read_text())
    def test_index_reuses_complete_contract(self):
        first=self.build_index();n=len(self.fake_calls);second=self.build_index()
        self.assertEqual(first['index_dir'],second['index_dir']);self.assertEqual(n,len(self.fake_calls))
    def test_index_changed_fasta_keeps_old_generation(self):
        first=self.build_index();p=self.hap/'chr1.fa';before=p.stat();p.write_text('>chr1\n'+'C'*100+'\n');os.utime(p,ns=(before.st_atime_ns,before.st_mtime_ns))
        second=self.build_index();self.assertNotEqual(first['index_dir'],second['index_dir'])
        self.assertTrue(Path(first['index_dir']).exists());self.assertIn('C'*100,Path(second['target_fasta']).read_text())
    def test_index_changed_gtf_and_force_regenerate(self):
        first=self.build_index();self.gtf.write_text(self.gtf.read_text().replace('\t55\t70\t','\t54\t70\t'))
        second=self.build_index();third=self.build_index('--force')
        self.assertNotEqual(first['index_dir'],second['index_dir']);self.assertNotEqual(second['index_dir'],third['index_dir'])
    def test_index_partial_output_is_not_reused(self):
        first=self.build_index();(Path(first['index_dir'])/'SA').unlink();second=self.build_index()
        self.assertNotEqual(first['index_dir'],second['index_dir'])
    def test_alignment_reuse_fastq_change_and_force(self):
        self.build_index();first=self.build_align();self.assertEqual(first['align_dir'],self.build_align()['align_dir'])
        self.r1.write_text('@R\nCCCC\n+\nIIII\n');second=self.build_align();self.assertNotEqual(first['align_dir'],second['align_dir'])
        third=self.build_align('--force');self.assertNotEqual(second['align_dir'],third['align_dir']);self.assertTrue(Path(first['align_dir']).exists())
    def test_alignment_rejects_changed_index(self):
        ix=self.build_index();(Path(ix['index_dir'])/'SAindex').write_text('외부변경\n')
        with self.assertRaises(SystemExit):self.build_align()
    def test_fastq_reuse_refreshes_for_new_preflight_input(self):
        samtools=shutil.which('samtools')
        if not samtools:self.skipTest('samtools 필요')
        stage=module('11');bam=self.root/'rna.bam';pre=self.root/'preflight.json';write_bam(bam,'A')
        def preflight():
            s=bam.stat();pre.write_text(json.dumps({'schema':'mfm_rna_preflight_v1','status':'pass','rna':{'path':str(bam),'size':s.st_size,'mtime_ns':s.st_mtime_ns,'inode':s.st_ino},'checks':{'sm_status':'pass','read_group_samples':['S'],'expected_samples':['S']}}))
        preflight();out=self.root/'fastq';args=['--sample','S','--rna',bam,'--out-dir',out,'--log-dir',self.root/'logs','--preflight-json',pre,'--samtools',samtools,'--decode-reference',self.root/'unused_for_bam.fa']
        invoke(stage,args);mtime=(out/'R1.fq').stat().st_mtime_ns;invoke(stage,args);self.assertEqual(mtime,(out/'R1.fq').stat().st_mtime_ns)
        write_bam(bam,'C');preflight();invoke(stage,args);self.assertIn('C'*40,(out/'R1.fq').read_text())
    def test_markdup_new_star_bam_regenerates(self):
        stage=module('15');ad=self.root/'aln';ad.mkdir();bam=ad/'star.bam';write_bam(bam,'A')
        al=self.root/'align.json';al.write_text(json.dumps({'align_dir':str(ad),'star_bam':str(bam)}));calls=[]
        def fake(cmd,**kwargs):
            if 'MarkDuplicates' in cmd:
                calls.append(cmd);values={x.split('=',1)[0]:x.split('=',1)[1] for x in cmd if '=' in x}
                shutil.copyfile(values['I'],values['O']);Path(values['M']).write_text('READ_PAIRS_EXAMINED\tREAD_PAIR_DUPLICATES\tPERCENT_DUPLICATION\n1\t0\t0\n')
            elif 'index' in cmd:pysam.index(cmd[-1])
            else:raise AssertionError(cmd)
        args=['--sample','S','--hap','1','--align-json',al,'--out-dir',self.root/'final','--work-dir',self.work,'--log-dir',self.root/'logs','--picard',sys.executable,'--samtools',sys.executable]
        with patch.object(stage,'run',fake):
            invoke(stage,args);invoke(stage,args);self.assertEqual(len(calls),1)
            write_bam(bam,'C');invoke(stage,args);self.assertEqual(len(calls),2)

if __name__=='__main__':unittest.main()
