"""독립 Python 기준과 통합 계수·참조 오류·실제 spool 예약을 대조함."""
from pathlib import Path
from unittest.mock import patch
import importlib.util,json,os,random,struct,subprocess,sys,tempfile,unittest
import pysam
import test_rna_wgs_native as fixture
from integrate_wgs_rna.rna_ase_shards import _native_reader,_read_shard_python,fused_context,spool_sizes,assign_from_bams
from integrate_wgs_rna.record_stats_parallel import native_tool
from rna_allelic_common import build_exon_index,split_bucket

ROOT=Path(__file__).resolve().parents[1]
if not (ROOT/'src').exists():ROOT=fixture.ROOT
STAGE='21_allelic_counts.py' if (ROOT/'09_personal_fasta.sh').exists() else '19_allelic_counts.py'
spec=importlib.util.spec_from_file_location('fused_reference',ROOT/STAGE)
reference=importlib.util.module_from_spec(spec);sys.modules[spec.name]=reference;spec.loader.exec_module(reference)


def write_shard(path,rows,refs):
    data=bytearray(b'MFMASE1\0'+struct.pack('<I',len(refs)))
    for ref in refs:
        encoded=ref.encode();data+=struct.pack('<I',len(encoded))+encoded
    for name,flag,mq,has,score,contig,blocks in rows:
        name=name.encode();data+=struct.pack('<IHBBqiI',len(name),flag,mq,has,score,refs.index(contig),len(blocks))+name
        for start,end in blocks:data+=struct.pack('<qq',start,end)
    path.write_bytes(data)


def expected(paths,ex1,ex2,both,evidence,mapq,duplicates,scale):
    A,Am,d1,n1=_read_shard_python(paths[0],mapq,duplicates);B,Bm,d2,n2=_read_shard_python(paths[1],mapq,duplicates)
    union=n1|n2;removed=reference.drop_duplicate_fragments(union,A,Am,B,Bm) if not duplicates else 0
    corrected=sum(reference.remove_length_term(t,scale) for t in (A,Am,B,Bm))
    counts,stats=reference.assign_fragments(A,B,ex1,ex2,both,evidence);mm,mstats=reference.assign_fragments(Am,Bm,ex1,ex2,both,evidence)
    audit=dict(dup1=d1,dup2=d2,dup_union=len(union),dup_inconsistent=len(n1^n2),removed=removed,corrected=corrected,total=len(set(A)|set(B)))
    return counts,stats,mm,mstats,audit


class FusedTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.refs=['chr1','chr2','chrΔ'];raw={c:[(-10,100,'Gα'),(50,250,'G2'),(300,500,'G3'),(16370,16410,'Gbin')] for c in self.refs}
        self.ex1=build_exon_index(raw);self.ex2=build_exon_index(raw)
        self.evidence=reference.HapEvidence.__new__(reference.HapEvidence);self.evidence.idx={}
        for h in [1,2]:
            for c in self.refs:
                self.evidence.idx['good',h,c]=([30,180,16385],[70,210,16395],40)
                self.evidence.idx['bad',h,c]=([400],[420],20)
        self.paths=[self.root/'a.bin',self.root/'b.bin'];self.native=_native_reader().__self__
    def engine(self,both,mapq,duplicates,scale):
        return self.native.prepare_assign(self.ex1,self.ex2,self.evidence.idx,both,mapq,duplicates,scale,split_bucket)

    def test_무작위_모든count_audit_동점_Unicode_64bit가동일함(self):
        rng=random.Random(826);rows=[[],[]]
        for h in [0,1]:
            for i in range(650):
                name='read'+str(rng.randrange(240)) if i%9 else '샘플α'+str(i%7)
                start=rng.choice([-5,50,190,310,410,16380]);blocks=[(start,start+rng.choice([0,10,30])),(start+40,start+60)]
                rows[h].append((name,rng.choice([64,128,0])+rng.choice([0,0,1024]),rng.choice([0,1,20,254,255]),rng.choice([0,1,1,1]),rng.choice([-(1<<63),(1<<63)-1,-10,0,20]),rng.choice(self.refs),blocks))
            for mate in [64,128]:rows[h].append(('huge',mate,255,1,(1<<63)-1-h,'chr1',[(50,100)]))
            rows[h]+=[('zero',64,255,1,10+h,'chr1',[(50,50)]),('reverse',128,255,1,10+h,'chr1',[(60,50)])]
            write_shard(self.paths[h],rows[h],self.refs if h==0 else list(reversed(self.refs)))
        for mapq in [0,20,255]:
            for duplicates in [False,True]:
                for both in [False,True]:
                    for scale in [0,.25,-.25]:
                        actual=self.native.assign_shards(str(self.paths[0]),str(self.paths[1]),self.engine(both,mapq,duplicates,scale))
                        self.assertEqual(actual[:5],expected(self.paths,self.ex1,self.ex2,both,self.evidence,mapq,duplicates,scale),(mapq,duplicates,both,scale))

    def test_손상_빈block_없는파일은조용히처리하지않음(self):
        rows=[('empty',64,255,1,1,'chr1',[])]
        for p in self.paths:write_shard(p,rows,self.refs)
        for method in [lambda:expected(self.paths,self.ex1,self.ex2,False,self.evidence,255,False,.25),lambda:self.native.assign_shards(str(self.paths[0]),str(self.paths[1]),self.engine(False,255,False,.25))]:
            with self.assertRaises(ValueError):method()
        p=self.paths[0];p.write_bytes(p.read_bytes()[:-1])
        with self.assertRaises(ValueError):self.native.assign_shards(str(p),str(self.paths[1]),self.engine(False,255,False,0))
        p.unlink()
        with self.assertRaises(FileNotFoundError):self.native.assign_shards(str(p),str(self.paths[1]),self.engine(False,255,False,0))

    def test_알려진정책만허용하고큰scale는auto에서기준경로를씀(self):
        args=(reference.assign_fragments,reference.drop_duplicate_fragments,reference.remove_length_term,self.ex1,self.ex2,False,self.evidence,255,False,.25)
        with patch.dict(os.environ,{'RNA_ASE_ASSIGN_BACKEND':'native'}):
            self.assertIsNotNone(fused_context(*args,'validated_hap_v1'))
            with self.assertRaises(ValueError):fused_context(*args,None)
        with patch.dict(os.environ,{'RNA_ASE_ASSIGN_BACKEND':'auto'}):self.assertIsNone(fused_context(*args[:-1],1e20,'validated_hap_v1'))


class CnvAndMemoryTest(unittest.TestCase):
    def test_실제BAM_CLI_계수_sidecar_정리후lease가동일함(self):
        from test_allelic_outputs import write_bam,write_hap_variants
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);a=root/'a.bam';b=root/'b.bam'
            write_bam(a,[('h1',50),('h2',40),('tie',45),('only1',45)])
            write_bam(b,[('h1',40),('h2',50),('tie',45)])
            (root/'exons').write_text('G1\tchr1\t0\t1000\n');write_hap_variants(root/'variants',[{}])
            pool=root/'pool';(pool/'leases').mkdir(parents=True);(pool/'waiters').mkdir();lease=pool/'leases/own.json'
            rows=[]
            for backend in ['tables','native']:
                out=root/backend;out.mkdir();lease.write_text(json.dumps(dict(cpus=4,total_cpus=4,branch='RNA',stage='21',memory_bytes=64<<30)))
                env={**os.environ,'RNA_ASE_ASSIGN_BACKEND':backend,'RNA_ALLELIC_INPUT_BACKEND':'partitioned','RNA_ALLELIC_JOBS':'4',
                     'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_CPU_BUDGET':'4','MFM_SAMPLE_CPU_BUDGET':'4','MFM_CPU_POOL_STAGE_MAX':'4','MFM_RAM_LIMIT_BYTES':str(128<<30)}
                command=[sys.executable,str(ROOT/STAGE),'--bam1',str(a),'--bam2',str(b),'--exons1',str(root/'exons'),'--exons2',str(root/'exons'),'--hap-variants',str(root/'variants'),'--out',str(out/'counts.tsv')]
                done=subprocess.run(command,env=env,capture_output=True,text=True,timeout=60)
                self.assertEqual(done.returncode,0,done.stdout+done.stderr)
                meta=json.loads((out/'counts.meta.json').read_text());execution=meta['resource_execution']['assignment']
                self.assertEqual(execution['assignment_backend'],backend)
                self.assertTrue(execution['spool_released_after_cleanup'])
                self.assertLess(execution['spool_actual']['accounted_bytes'],meta['resource_execution']['input']['spool_estimate_bytes'])
                after=json.loads(lease.read_text());self.assertEqual(after['cpus'],1);self.assertEqual(after['memory_bytes'],execution['post_cleanup_parent_bytes'])
                rows.append((out,meta))
            for name in ['counts.tsv','counts.both_haps.tsv','counts.split_half.tsv','counts.mmap_adjusted.tsv']:
                self.assertEqual((rows[0][0]/name).read_bytes(),(rows[1][0]/name).read_bytes())
            for key in rows[0][1]:
                if key not in {'resource_execution','output_sha256'} and not key.endswith('_path'):self.assertEqual(rows[0][1][key],rows[1][1][key],key)

    def test_CNV전체해독과선택해독이BAM_CRAM에서동일하고잘못된참조를거부함(self):
        f=fixture.WgsNativeTests();f.setUp();self.addCleanup(f.temp.cleanup)
        bed=f.root/'regions.bed';bed.write_text('chr1\t0\t3000\tG\nchr1\t70\t130\tH\n')
        tool=str(native_tool('rna_cnv_sum'))
        for source in [f.bam,f.cram]:
            outputs=[]
            for mode in ['full','coverage']:
                p=subprocess.run([tool,str(source),str(f.ref),'chr1',str(bed),'1'],env={**os.environ,'RNA_CNV_CRAM_FIELDS':mode},capture_output=True,text=True,check=True);outputs.append(p.stdout)
            self.assertEqual(*outputs)
        wrong=f.root/'wrong.fa';wrong.write_text('>chr1\n'+'C'*3000+'\n>chr2\n'+'C'*100+'\n');pysam.faidx(str(wrong))
        for mode in ['full','coverage']:
            p=subprocess.run([tool,str(f.cram),str(wrong),'chr1',str(bed),'0'],env={**os.environ,'RNA_CNV_CRAM_FIELDS':mode,'REF_PATH':str(f.root/'none'),'REF_CACHE':str(f.root/'none')},capture_output=True,text=True)
            self.assertNotEqual(p.returncode,0,p.stdout+p.stderr)

    def test_spool실제논리_할당량을모두계상함(self):
        with tempfile.TemporaryDirectory() as raw:
            files=[Path(raw)/str(i) for i in range(4)]
            for i,p in enumerate(files):p.write_bytes(b'x'*(i+10))
            info=spool_sizes([((files[0],files[1]),),((files[2],files[3]),)])
            self.assertEqual(info['logical_bytes'],46);self.assertEqual(info['files'],4)
            self.assertEqual(info['accounted_bytes'],max(46,sum(p.stat().st_blocks*512 for p in files)))


if __name__=='__main__':unittest.main()
