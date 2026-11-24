"""원래 영역 occurrence·BED 합산순서 및 ASE 이진 해석의 완전 동등성을 확인함."""
from pathlib import Path
from unittest.mock import patch
import json,os,struct,subprocess,sys,tempfile,unittest
import pysam
import test_rna_wgs_native as fixture
ROOT=fixture.ROOT;VARIANT=fixture.VARIANT
from integrate_wgs_rna.rna_task_chunks import partition_rows
from integrate_wgs_rna.rna_ase_shards import read_shard,_read_shard_python,_native_reader
from integrate_wgs_rna.record_stats_parallel import NativeUnavailable
from rna_cnv_depth import collect_depth


def shard(path):
    refs=['chr1','chrΔ'];data=bytearray(b'MFMASE1\0'+struct.pack('<I',len(refs)))
    for ref in refs:
        raw=ref.encode();data+=struct.pack('<I',len(raw))+raw
    rows=[('샘플α',64,255,1,-8,1,[(1,5),(9,12)]),('샘플α',64,255,1,-8,0,[(40,60)]),
          ('샘플α',64,255,1,-7,0,[(5,10)]),('샘플α',128,255,1,-10,1,[]),
          ('zero',64,0,1,0,0,[]),('low',128,20,1,-(1<<63),0,[(1,2)]),
          ('high',64,255,1,(1<<63)-1,0,[(2,3)]),('dup',1088,0,0,0,0,[]),
          ('noAS',128,255,0,0,0,[]),('dup2',1152,255,1,5,1,[(3,8)])]
    for name,flag,mq,has,score,tid,blocks in rows:
        raw=name.encode();data+=struct.pack('<IHBBqiI',len(raw),flag,mq,has,score,tid,len(blocks))+raw
        for begin,end in blocks:data+=struct.pack('<qq',begin,end)
    path.write_bytes(data);return bytes(data)


class RnaChunkTest(unittest.TestCase):
    def test_분할은긴행_중복_원래번호를보존함(self):
        rows=[(100,105),(0,1),(8,40),(3,5),(100,105),(50,50)]
        with patch.dict(os.environ,{'RNA_WGS_TASK_LAYOUT':'chunks','RNA_WGS_CHUNK_BP':'10'}):parts=partition_rows(rows,lambda r:r)
        self.assertEqual(sorted(i for part in parts for i in part),list(range(len(rows))))
        for part in parts:
            width=max(rows[i][1] for i in part)-min(rows[i][0] for i in part)
            self.assertTrue(width<=10 or len(part)==1)
        self.assertIn([2],parts)
        with patch.dict(os.environ,{'RNA_WGS_TASK_LAYOUT':'contig'}):self.assertEqual(partition_rows(rows,lambda r:r),[list(range(len(rows)))])
        with patch.dict(os.environ,{'RNA_WGS_TASK_LAYOUT':'bad'}),self.assertRaises(ValueError):partition_rows(rows,lambda r:r)

    def test_다른작업RAM을빼고fork_worker상한을계산함(self):
        from integrate_wgs_rna.rna_wgs_native import plan
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();own=root/'leases/own.json'
            own.write_text(json.dumps(dict(cpus=17,total_cpus=68,memory_bytes=16<<30)))
            (root/'leases/other.json').write_text(json.dumps(dict(cpus=0,memory_bytes=88<<30)))
            env={'MFM_CPU_POOL_LEASE_FILE':str(own),'MFM_CPU_BUDGET':'68','MFM_SAMPLE_CPU_BUDGET':'17','MFM_CPU_POOL_STAGE_MAX':'68','MFM_RAM_LIMIT_BYTES':str(128<<30)}
            with patch.dict(os.environ,env):result=plan(300,decode_threads=0)
            self.assertEqual(result['workers_upper'],5)
            self.assertLessEqual(result['peak_ram_upper_bytes'],result['shared_ram_budget_bytes'])

    def test_CNV의BAM_CRAM_겹침행과유전자합이동일함(self):
        f=fixture.WgsNativeTests();f.setUp();self.addCleanup(f.temp.cleanup)
        rows=[('chr1',1400,2000,'G'),('chr1',0,100,'G'),('chr1',50,800,'H'),('chr1',50,800,'H'),('chr1',300,300,'zero'),('chr2',0,99,'other')]
        inputs=[dict(rows=rows,bam=bam,reference=str(f.ref)) for bam in [f.bam,f.cram]]
        result=[]
        for mode in ['contig','chunks']:
            execution={}
            with patch.dict(os.environ,{'RNA_CNV_BACKEND':'native','RNA_WGS_TASK_LAYOUT':mode,'RNA_WGS_CHUNK_BP':'300','MFM_CPU_BUDGET':'6','MFM_SAMPLE_CPU_BUDGET':'6','MFM_CPU_POOL_STAGE_MAX':'6'}):
                result.append(collect_depth(inputs,'samtools','auto',execution=execution))
            self.assertEqual(execution['layout'],mode)
        self.assertEqual(result[0],result[1])

    def test_전체변이CLI의긴DEL중복occurrence와모든열이동일함(self):
        f=fixture.WgsNativeTests();f.setUp();self.addCleanup(f.temp.cleanup)
        gtf=f.root/'genes.gtf';gtf.write_text(''.join(f'chr1\tfixture\texon\t{s+1}\t{e}\t.\t+\t.\tgene_id "G{i}";\n' for i,(s,e) in enumerate([(40,100),(200,250),(300,350),(600,700),(1400,1600)])))
        vcf=f.root/'variants.vcf';vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=3000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'+''.join(f'chr1\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0|1\n' for pos,ref,alt in [(60,'A','T'),(201,'A'*150,'A'),(650,'A','T'),(1500,'A','T')]))
        gz=Path(pysam.tabix_index(str(vcf),preset='vcf',force=True))
        for h in [1,2]:
            d=f.root/f'HAP{h}';d.mkdir();seq=list('A'*3000)
            if h==2:
                for p in [60,650,1500]:seq[p-1]='T'
            (d/'chr1.fa').write_text('>chr1\n'+''.join(seq)+'\n');pysam.faidx(str(d/'chr1.fa'))
            (d/'chr1.chain').write_text('chain 0 chr1 3000 + 0 3000 chr1 3000 + 0 3000 1\n3000\n')
        tables=[]
        for mode in ['contig','chunks']:
            out=f.root/(mode+'.tsv')
            command=[sys.executable,str(ROOT/VARIANT),'--vcf',str(gz),'--wgs-bam',str(f.cram),'--reference',str(f.ref),'--hap-dir1',str(f.root/'HAP1'),'--hap-dir2',str(f.root/'HAP2'),'--gtf',str(gtf),'--contigs','chr1','--flank','0','--threads','6','--out',str(out),'--qc-json',str(out)+'.json']
            done=subprocess.run(command,env={**os.environ,'RNA_HAP_VARIANT_BACKEND':'native','RNA_WGS_TASK_LAYOUT':mode,'RNA_WGS_CHUNK_BP':'250','MFM_CPU_BUDGET':'6','MFM_SAMPLE_CPU_BUDGET':'6','MFM_CPU_POOL_STAGE_MAX':'6'},capture_output=True,text=True,timeout=90)
            self.assertEqual(done.returncode,0,done.stdout+done.stderr);tables.append(out.read_bytes())
        self.assertEqual(*tables);self.assertEqual(sum(line.split(b'\t')[1]==b'201' for line in tables[0].splitlines()[1:]),2)


class AseReaderTest(unittest.TestCase):
    def test_CAPI해석은모든객체_Unicode_동점_mate_duplicate_MAPQ가동일함(self):
        with tempfile.TemporaryDirectory() as raw:
            path=Path(raw)/'part.bin';data=shard(path)
            for mapq in [0,20,255]:
                for duplicates in [False,True]:
                    self.assertEqual(_native_reader()(str(path),mapq,duplicates),_read_shard_python(path,mapq,duplicates))
            for cut in [0,7,11,len(data)-1,len(data)-3]:
                path.write_bytes(data[:cut])
                for reader in [_native_reader(),_read_shard_python]:
                    with self.assertRaises(ValueError):reader(str(path),255,False)
            path.write_bytes(b'MFMASE1\0'+struct.pack('<I',1000001))
            for reader in [_native_reader(),_read_shard_python]:
                with self.assertRaises(ValueError):reader(str(path),255,False)

    def test_빌드불가일때만자동Python을사용함(self):
        with tempfile.TemporaryDirectory() as raw:
            path=Path(raw)/'part.bin';shard(path)
            with patch('integrate_wgs_rna.rna_ase_shards._native_reader',side_effect=NativeUnavailable('no compiler')):
                with patch.dict(os.environ,{'RNA_ASE_SHARD_READER':'auto'}):self.assertEqual(read_shard(path,255,False),_read_shard_python(path,255,False))
                with patch.dict(os.environ,{'RNA_ASE_SHARD_READER':'native'}),self.assertRaises(NativeUnavailable):read_shard(path,255,False)


if __name__=='__main__':unittest.main()
