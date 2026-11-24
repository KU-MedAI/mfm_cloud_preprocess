"""RNA의 네이티브 경로가 독립 기준선과 모든 출력에서 같은지 검증함."""
from pathlib import Path
import gzip,hashlib,json,os,subprocess,sys,tempfile,unittest
from unittest.mock import patch
import numpy as np
import pysam,zarr
ROOT=Path(__file__).resolve().parents[1]
COUNTS='21_allelic_counts.py' if (ROOT/'09_personal_fasta.sh').exists() else '19_allelic_counts.py'
sys.path[:0]=[str(ROOT),str(ROOT/'src'),str(ROOT/'tests')]
from test_allelic_input_parallel import make_bam
from test_allelic_outputs import write_hap_variants,write_bam
from test_rna_parallel import RnaParallelTest
from integrate_wgs_rna.allelic_input import read_as
from integrate_wgs_rna.record_stats_parallel import native_tool
from integrate_wgs_rna.rna_ase_shards import read_shard
from integrate_wgs_rna.rna_fastq_native import process


def _resource_job(value):return value+1


class NativeTests(unittest.TestCase):
    def test_retained_fork_memory_is_released_after_pool_join(self):
        import integrate_wgs_rna.cpu_executor as module
        actual=module.WorkerBudget;observed=[]
        class Traced(actual):
            def capacity(self,*args,**kwargs):
                observed.append(self.base_memory)
                return super().capacity(*args,**kwargs)
        with patch.object(module,'WorkerBudget',Traced),patch.dict(os.environ,{'MFM_RAM_LIMIT_BYTES':'100000000000'}):
            values=list(module.cpu_results(_resource_job,[(i,) for i in range(4)],2,elastic=True,max_workers=2,memory_per_worker_bytes=4<<30,retain_worker_memory=True))
        self.assertTrue(all(error is None for _,_,error in values));self.assertEqual(max(observed),12<<30);self.assertEqual(observed[-1],4<<30)


    def test_qname_partition_reconstructs_every_fragment(self):
        with tempfile.TemporaryDirectory() as raw:
            folder=Path(raw);source=folder/'source.bam';bam=folder/'reads.bam';make_bam(source,index=False)
            with pysam.AlignmentFile(str(source),'rb') as handle,pysam.AlignmentFile(str(bam),'wb',template=handle) as out:
                for read in handle:out.write(read)
                for name,mq,flag,cigar in [('zero',0,64,[(0,20)]),('dup_zero',0,1088,[(0,20)]),('complex',255,576,[(4,5),(0,5),(1,2),(7,3),(2,4),(8,2),(3,20),(0,5)]),('secondary',255,320,[(0,20)]),('supplementary',255,2112,[(0,20)])]:
                    read=pysam.AlignedSegment();read.query_name=name;read.flag=flag;read.reference_id=0;read.reference_start=500;read.mapping_quality=mq;read.cigar=cigar
                    read.query_sequence='A'*sum(n for op,n in cigar if op in [0,1,4,7,8]);read.set_tag('AS',-7);out.write(read)
            for duplicates in [False,True]:
                prefix=folder/'part'
                subprocess.run([str(native_tool('rna_ase_partition')),str(bam),str(prefix),'8',str(int(duplicates)),'2','255'],check=True,capture_output=True)
                reconstructed=[{}, {},0,set()]
                for i in range(8):
                    values=read_shard(Path(str(prefix)+f'.{i}.bin'),255,duplicates)
                    for index in [0,1]:reconstructed[index].update(values[index])
                    reconstructed[2]+=values[2];reconstructed[3].update(values[3])
                self.assertEqual(tuple(reconstructed),read_as(str(bam),255,duplicates))
                subprocess.run([str(native_tool('rna_ase_partition')),str(bam),str(prefix),'8',str(int(duplicates)),'0','0'],check=True,capture_output=True)
                tables=[{}, {},0,set()]
                for i in range(8):
                    values=read_shard(Path(str(prefix)+f'.{i}.bin'),0,duplicates)
                    for index in [0,1]:tables[index].update(values[index])
                    tables[2]+=values[2];tables[3].update(values[3])
                self.assertEqual(tuple(tables),read_as(str(bam),0,duplicates))
                broken=Path(str(prefix)+'.0.bin');broken.write_bytes(b'MFMASE1\0')
                with self.assertRaisesRegex(ValueError,'잘린'):read_shard(broken,255,False)
            bad=folder/'truncated.bam';bad.write_bytes(bam.read_bytes()[:-28])
            run=subprocess.run([str(native_tool('rna_ase_partition')),str(bad),str(folder/'bad'),'8','0'],capture_output=True,text=True)
            self.assertNotEqual(run.returncode,0);self.assertIn('truncated',run.stderr)

    def test_partitioned_all_tsv_and_counts_match_legacy(self):
        with tempfile.TemporaryDirectory() as raw:
            folder=Path(raw);a=folder/'a.bam';b=folder/'b.bam'
            write_bam(a,[('h1',50),('h2',40),('tie',45),('only1',45)])
            write_bam(b,[('h1',40),('h2',50),('tie',45)])
            (folder/'exons').write_text('G1\tchr1\t0\t1000\n');write_hap_variants(folder/'hv',[{}])
            metas=[]
            for mode in ['legacy','partitioned']:
                target=folder/mode;target.mkdir()
                result=subprocess.run([sys.executable,str(ROOT/COUNTS),'--bam1',str(a),'--bam2',str(b),'--exons1',str(folder/'exons'),'--exons2',str(folder/'exons'),'--hap-variants',str(folder/'hv'),'--out',str(target/'counts.tsv')],env={**os.environ,'RNA_ALLELIC_INPUT_BACKEND':mode,'RNA_ALLELIC_JOBS':'4','MFM_SAMPLE_CPU_BUDGET':'4'},capture_output=True,text=True,timeout=60)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                meta=json.loads((target/'counts.meta.json').read_text());metas.append(meta)
            for name in ['counts.tsv','counts.split_half.tsv','counts.both_haps.tsv','counts.mmap_adjusted.tsv']:
                self.assertEqual((folder/'legacy'/name).read_bytes(),(folder/'partitioned'/name).read_bytes(),name)
            for key in metas[0]:
                if key not in {'resource_execution','output_sha256'} and not key.endswith('_path'):self.assertEqual(metas[0][key],metas[1][key],key)

    def test_fastq_copy_counts_gzip_and_corruption(self):
        with tempfile.TemporaryDirectory() as raw:
            folder=Path(raw);data=b'@read comment/1\nACNT\n+unchanged\nIJKL\n@read2/1\nAAA\n+\nFFF\n'
            source=folder/'reads.gz';source.write_bytes(gzip.compress(data));out=folder/'out.fq'
            info=process(source,out,1)
            self.assertEqual(info['backend'],'native');self.assertEqual(info['records'],2);self.assertEqual(info['max_read_length'],4)
            self.assertEqual(out.read_bytes(),data.replace(b'comment/1\n',b'comment\n').replace(b'read2/1\n',b'read2\n'))
            self.assertEqual(process(out)['records'],2)
            source.write_bytes(gzip.compress(data)[:-7])
            with self.assertRaises(ValueError):process(source)
            source.write_bytes(b'@r\nACG\n+\nII\n')
            with self.assertRaises(ValueError):process(source)

    def test_native_depth_gpu_and_buffer_contract(self):
        fixture=RnaParallelTest('test_genome_serial_parallel_and_raw_filter');fixture.setUp()
        try:
            expected=None
            modes=['legacy','native']
            if os.environ.get('RNA_TEST_GPU')=='1':modes.append('gpu')
            for mode in modes:
                output=fixture.base/(mode+'.zarr');windows=fixture.base/(mode+'_windows.zarr')
                env={'RNA_DEPTH_BACKEND':'legacy' if mode=='legacy' else 'native','MFM_COVERAGE_BACKEND':'cupy' if mode=='gpu' else 'cpu'}
                with patch.dict(os.environ,env):
                    fixture.run_script('18_coverage_zarr.py',['--bam',fixture.bam,'--fai',fixture.fai,'--out',output,'--windows-tsv',fixture.windows,'--window-out',windows,'--contigs','chr1 chr2 chr3','--sample','S1','--hap','HAP1','--jobs','2','--chunk','64'])
                store=zarr.open_group(str(output),mode='r');win=zarr.open_group(str(windows),mode='r')
                actual={c:store[c][:] for c in store.array_keys()};actual.update({'window_'+k:win[k][:] for k in win.array_keys()})
                if expected is None:expected=actual
                else:
                    self.assertEqual(expected.keys(),actual.keys())
                    for key in actual:np.testing.assert_array_equal(actual[key],expected[key],err_msg=key)
        finally:fixture.tearDown()


if __name__=='__main__':unittest.main()
