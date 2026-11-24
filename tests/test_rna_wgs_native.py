"""CIGAR 경계/필터/HP와 겹치는 BED의 BAM·CRAM 결과를 독립 경로와 비교함."""
from pathlib import Path
import copy,importlib.util,os,random,subprocess,sys,tempfile,unittest
from unittest.mock import patch
import pysam
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'src')]
from integrate_wgs_rna.record_stats_parallel import native_tool
from integrate_wgs_rna.rna_wgs_native import count_support,plan,select
from integrate_wgs_rna.record_stats_parallel import NativeUnavailable
VARIANT='20_hap_variants.py' if (ROOT/'09_personal_fasta.sh').exists() else '18_hap_variants.py'
spec=importlib.util.spec_from_file_location('variant_reference',ROOT/VARIANT);legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy)


class WgsNativeTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.ref=self.root/'ref.fa';self.ref.write_text('>chr1\n'+'A'*3000+'\n>chr2\n'+'A'*100+'\n');pysam.faidx(str(self.ref))
        self.bam=self.root/'reads.bam';header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':3000},{'SN':'chr2','LN':100}]}
        rng=random.Random(42);reads=[]
        cigars=[[(0,100)],[(4,5),(0,20),(1,3),(7,20),(2,3),(8,20)],[(0,10),(3,30),(0,30)],[(2,5),(0,30),(2,5)],[(3,5),(0,30),(3,5)],[(0,10),(1,2),(0,20)]]
        for i in range(240):
            read=pysam.AlignedSegment();read.query_name='r'+str(i);read.reference_id=0;read.reference_start=rng.randrange(30,1800)
            read.flag=[0,64,128,1024,512,256,2048][i%7];read.mapping_quality=[0,19,20,60,255][i%5];read.cigar=cigars[i%len(cigars)]
            size=sum(n for op,n in read.cigar if op in [0,1,4,7,8]);read.query_sequence=''.join(rng.choice('ACGTN') for _ in range(size));read.query_qualities=pysam.qualitystring_to_array('I'*size)
            if i%5 in [0,1]:read.set_tag('HP',i%5+1)
            elif i%5==2:read.set_tag('HP',3.0,value_type='f')
            elif i%5==3:read.set_tag('HP','1',value_type='Z')
            reads.append(read)
        with pysam.AlignmentFile(str(self.bam),'wb',header=header) as out:
            for r in sorted(reads,key=lambda r:r.reference_start):out.write(r)
        pysam.index(str(self.bam));self.cram=self.root/'reads.cram'
        with pysam.AlignmentFile(str(self.bam),'rb') as source,pysam.AlignmentFile(str(self.cram),'wc',template=source,reference_filename=str(self.ref)) as out:
            for read in source:out.write(read)
        pysam.index(str(self.cram))

    def groups(self):
        result=[]
        for begin,end in [(0,600),(550,1100),(1400,1800),(2000,2500)]:
            variants=[]
            for pos in range(max(0,begin-20),end+10,13):
                for kind,ref,alt in [('SNV','A','T'),('MNV','AA','GT'),('INDEL','AAAA','A'),('INDEL','A','AGGG')]:
                    variants.append(dict(p0=pos,ref=ref,alt=alt,kind=kind,wgs_depth=0,wgs_alt=0,hp1_alt=0,hp1_depth=0,hp2_alt=0,hp2_depth=0))
            result.append(((begin,end),variants))
        return result

    def test_variant_all_counts_and_region_occurrences(self):
        for bam in [self.bam,self.cram]:
            for mapq in [0,20,255]:
                expected=self.groups();actual=copy.deepcopy(expected)
                with pysam.AlignmentFile(str(bam),reference_filename=str(self.ref)) as source:
                    for region,variants in expected:legacy.count_support(source,str(self.ref),'chr1',region,variants,mapq)
                count_support(str(native_tool('rna_variant_support')),bam,self.ref,'chr1',actual,mapq,1)
                self.assertEqual(actual,expected)

    def test_padding_compatibility_and_invalid_float_hp(self):
        bam=self.root/'padding.bam'
        with pysam.AlignmentFile(str(self.bam),'rb') as template,pysam.AlignmentFile(str(bam),'wb',template=template) as out:
            read=pysam.AlignedSegment();read.query_name='padding';read.reference_id=0;read.reference_start=10;read.mapping_quality=60;read.query_sequence='ACGTACGT';read.cigar=[(0,2),(6,2),(0,6)];out.write(read)
        pysam.index(str(bam))
        rows=[dict(p0=p,ref='A',alt='A',kind='SNV',wgs_depth=0,wgs_alt=0,hp1_alt=0,hp1_depth=0,hp2_alt=0,hp2_depth=0) for p in range(10,16)]
        expected=[((0,100),rows)];actual=copy.deepcopy(expected)
        with pysam.AlignmentFile(str(bam),'rb') as stream:legacy.count_support(stream,None,'chr1',(0,100),rows,20)
        count_support(str(native_tool('rna_variant_support')),bam,None,'chr1',actual,20)
        self.assertEqual(actual,expected)
        bad=self.root/'hp_float.bam'
        with pysam.AlignmentFile(str(bam),'rb') as template,pysam.AlignmentFile(str(bad),'wb',template=template) as out:
            read=next(template);read.set_tag('HP',1.0,value_type='f');out.write(read)
        pysam.index(str(bad));rows=copy.deepcopy(expected[0][1])
        with pysam.AlignmentFile(str(bad),'rb') as stream:
            with self.assertRaises(KeyError):legacy.count_support(stream,None,'chr1',(0,100),rows,20)
        with self.assertRaisesRegex(RuntimeError,'non-integer HP'):count_support(str(native_tool('rna_variant_support')),bad,None,'chr1',copy.deepcopy(expected),20)

    def test_cnv_all_intervals_equal_bedcov(self):
        rng=random.Random(23);bed=self.root/'genes.bed'
        rows=[(0,30),(10,100),(90,100),(120,120),(500,2500),(2000,2700),(2900,3100)]
        rows += [(s,s+rng.randrange(1,400)) for s in (rng.randrange(0,2600) for _ in range(160))]
        bed.write_text(''.join(f'chr1\t{s}\t{e}\tg{i%13}\n' for i,(s,e) in enumerate(rows)))
        for bam in [self.bam,self.cram]:
            expected=subprocess.run(['samtools','bedcov','--reference',str(self.ref),str(bed),str(bam)],capture_output=True,text=True,check=True)
            actual=subprocess.run([str(native_tool('rna_cnv_sum')),str(bam),str(self.ref),'chr1',str(bed),'1'],capture_output=True,text=True,check=True)
            self.assertEqual(actual.stdout,expected.stdout)

    def test_deletion_skip_only_records_and_gene_whitespace(self):
        bam=self.root/'special.bam';bed=self.root/'special.bed';bed.write_text('chr1\t0\t200\tgene space\n')
        with pysam.AlignmentFile(str(self.bam),'rb') as template,pysam.AlignmentFile(str(bam),'wb',template=template) as out:
            for i,cigar in enumerate([[(2,20)],[(3,20)],[(4,3),(2,20)],[(1,3),(3,20)],[(2,20),(0,2)],[(0,2),(3,20)]]):
                read=pysam.AlignedSegment();read.query_name=str(i);read.reference_id=0;read.reference_start=10+i*25;read.mapping_quality=60;read.cigar=cigar
                n=sum(v for op,v in cigar if op in [0,1,4,7,8]);read.query_sequence='A'*n if n else None;out.write(read)
        pysam.index(str(bam))
        expected=subprocess.run(['samtools','bedcov',str(bed),str(bam)],capture_output=True,text=True,check=True)
        actual=subprocess.run([str(native_tool('rna_cnv_sum')),str(bam),'-','chr1',str(bed),'0'],capture_output=True,text=True,check=True)
        self.assertEqual(actual.stdout,expected.stdout)

    def test_empty_contig_and_missing_index_fail(self):
        bed=self.root/'empty.bed';bed.write_text('chr2\t0\t100\tG\n')
        tool=str(native_tool('rna_cnv_sum'));run=subprocess.run([tool,str(self.bam),str(self.ref),'chr2',str(bed),'0'],capture_output=True,text=True,check=True)
        self.assertEqual(run.stdout,'chr2\t0\t100\tG\t0\n')
        (Path(str(self.bam)+'.bai')).unlink()
        run=subprocess.run([tool,str(self.bam),str(self.ref),'chr2',str(bed),'0'],capture_output=True,text=True)
        self.assertNotEqual(run.returncode,0)

    def test_backend_fallback_only_on_unavailable_build(self):
        for mode in ['auto','native','legacy','invalid']:
            with patch.dict(os.environ,{'RNA_CNV_BACKEND':mode}),patch('integrate_wgs_rna.rna_wgs_native.native_tool',side_effect=NativeUnavailable('compiler unavailable')):
                if mode in ['auto','legacy']:self.assertIsNone(select('RNA_CNV_BACKEND','rna_cnv_sum'))
                elif mode=='native':
                    with self.assertRaises(NativeUnavailable):select('RNA_CNV_BACKEND','rna_cnv_sum')
                else:
                    with self.assertRaises(ValueError):select('RNA_CNV_BACKEND','rna_cnv_sum')

    def test_worker_budget_accounts_for_htslib_io_and_decode(self):
        for cpu in [1,2,5,6,17,68]:
            with patch.dict(os.environ,{'MFM_CPU_BUDGET':str(cpu),'MFM_SAMPLE_CPU_BUDGET':str(cpu),'MFM_CPU_POOL_STAGE_MAX':str(cpu)}):value=plan(44)
            self.assertLessEqual(value['workers_upper']*value['cpu_per_worker']+(cpu>1),cpu)
            self.assertEqual(value['cpu_per_worker'],1+value['decode_threads']+int(value['decode_threads']>0))


if __name__=='__main__':unittest.main()
