"""단계별 실제 결과 동등성과 CPU/RAM 반환·압축 병렬 예산을 검증한다."""
from pathlib import Path
from unittest.mock import patch
from dataclasses import replace
from types import SimpleNamespace
import json,os,tempfile,unittest
import pysam
from integrate_wgs_rna.memory_policy import GIB
from integrate_wgs_rna.window_index import WindowIndexConfig,run_window_index_generation
from integrate_wgs_rna.compact.resources import plan,fit_memory
from integrate_wgs_rna.realign.align_resources import alignment_plan
from integrate_wgs_rna.personal_fasta import PersonalFastaConfig,process_one_sample

class StageOptimizationTests(unittest.TestCase):
    def sample(self,root,sample='S'):
        for hap in ['HAP1','HAP2']:
            folder=root/sample/'individual_fasta'/hap;folder.mkdir(parents=True)
            for contig in ['chr1','chr2']:
                sequence='ACGT'*600
                if hap=='HAP2':sequence=sequence[:512]+'N'*32+sequence[544:]
                p=folder/(contig+'.fa');p.write_text('>'+contig+'\n'+sequence+'\n');pysam.faidx(str(p))

    def config(self,root,samples=None):
        blacklist=root/'blacklist.bed';blacklist.write_text('chr1\t32\t64\nchr2\t500\t530\n')
        return WindowIndexConfig(cohort='T',output_root=root,samples=samples or ['S'],
            chroms=['chr1','chr2'],sequence_variant_policy='snp_only',window_size=256,
            offset_unit=32,blacklist_bed=blacklist,jobs=2)

    def test_window_serial_parallel_preserve_exact_rows_offsets_and_stats(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);left=root/'serial';right=root/'parallel';self.sample(left);self.sample(right)
            env={'MFM_CPU_BUDGET':'4','MFM_RAM_LIMIT_BYTES':str(64*GIB),'MFM_CPU_UTILIZATION':'maximum'}
            with patch.dict(os.environ,dict(env,MFM_WINDOW_CONTIG_JOBS='1'),clear=True):
                a=run_window_index_generation(self.config(left))
            with patch.dict(os.environ,dict(env,MFM_WINDOW_CONTIG_JOBS='4'),clear=True):
                b=run_window_index_generation(self.config(right))
            for filename in ['windows.tsv','windows_blacklist.tsv']:
                self.assertEqual((left/'S/windows'/filename).read_bytes(),(right/'S/windows'/filename).read_bytes())
            first=json.loads((left/'S/windows/windows.json').read_text());second=json.loads((right/'S/windows/windows.json').read_text())
            self.assertEqual(first['haplotypes'],second['haplotypes'])
            self.assertEqual(a['ok_count'],1);self.assertEqual(b['ok_count'],1)
            # 유효 캐시는 parallel 실행을 다시 하지 않는다.
            with patch('integrate_wgs_rna.window_parallel.cpu_results',side_effect=AssertionError('캐시 재계산')):
                self.assertEqual(run_window_index_generation(self.config(right))['skip_count'],1)

    def test_window_ram_prediction_uses_indexed_bases_not_compressed_bytes(self):
        from integrate_wgs_rna.window_parallel import _sequence_bases
        with tempfile.TemporaryDirectory() as raw:
            p=Path(raw)/'genome.fa.gz';p.write_bytes(b'fixture')
            index=Path(str(p)+'.fai');index.write_text('1\t3000000000\t0\t60\t61\n2\t123\t0\t60\t61\n')
            self.assertEqual(_sequence_bases(dict(fasta=p,contig='chr1')),3000000000)
            index.unlink()
            with self.assertRaises(ValueError):_sequence_bases(dict(fasta=p,contig='chr1'))

    def test_multi_sample_window_does_not_create_nested_process_pool(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);self.sample(root,'S');self.sample(root,'S2')
            with patch('integrate_wgs_rna.window_parallel.cpu_results',side_effect=AssertionError('중첩pool')):
                result=run_window_index_generation(self.config(root,['S','S2']))
            self.assertEqual(result['ok_count'],2)

    def test_auto_compression_uses_reserved_cpu_and_accounts_extra_workspace(self):
        with patch.dict(os.environ,{'MFM_CPU_BUDGET':'68','MFM_BLOSC_THREADS':'1','MFM_BLOSC_THREADS_REQUESTED':'auto'},clear=True):
            execution=plan(22,22,hap_split=True)
            self.assertEqual((execution['workers'],execution['cpu_per_worker'],execution['alignment_threads'],execution['blosc_threads']),(17,4,2,2))
            forecast=dict(worker_ceiling_bytes=32*GIB,worker_reservation_bytes=8*GIB,estimated_peak_per_worker_bytes=4*GIB,memory_limit_bytes=1650000000000)
            fit_memory(22,22,forecast,hap_split=True)
            self.assertEqual(forecast['worker_reservation_bytes'],18*GIB+128*2**20)
        for cpu in [1,2,3,4,8,17,32,68]:
            with patch.dict(os.environ,{'MFM_CPU_BUDGET':str(cpu),'MFM_BLOSC_THREADS':'1','MFM_BLOSC_THREADS_REQUESTED':'auto'},clear=True):
                for split in [True,False]:
                    p=plan(22,22,hap_split=split)
                    self.assertLessEqual(p['workers']*p['cpu_per_worker'],cpu)
        with patch.dict(os.environ,{'MFM_CPU_BUDGET':'68','MFM_BLOSC_THREADS':'1','MFM_BLOSC_THREADS_REQUESTED':'1'},clear=True):
            self.assertEqual(plan(22,22,hap_split=True)['blosc_threads'],1)

    def test_explicit_codec_threads_can_use_one_hap_when_two_do_not_fit_cpu(self):
        with patch.dict(os.environ,{'MFM_CPU_BUDGET':'8','MFM_BLOSC_THREADS':'8',
                                   'MFM_BLOSC_THREADS_REQUESTED':'8'},clear=True):
            p=plan(22,22,hap_split=True)
            self.assertEqual((p['hap_jobs'],p['blosc_threads'],p['cpu_per_worker']),(1,8,8))
            with patch.dict(os.environ,MFM_COMPACT_HAP_JOBS='2'):
                with self.assertRaises(ValueError):plan(22,22,hap_split=True)

    def test_native_workspace_is_added_to_each_hap_address_limit(self):
        forecast=dict(worker_ceiling_bytes=8*GIB,worker_reservation_bytes=8*GIB,
                      estimated_peak_per_worker_bytes=4*GIB,memory_limit_bytes=1650000000000)
        with patch.dict(os.environ,{'MFM_CPU_BUDGET':'68','MFM_BLOSC_THREADS':'1',
                                   'MFM_BLOSC_THREADS_REQUESTED':'auto'},clear=True):
            execution=fit_memory(1,1,forecast,hap_split=True)
        native=execution['extra_decode_buffers_bytes']+execution['extra_codec_buffers_bytes']
        self.assertGreater(native,0)
        self.assertEqual(forecast['per_hap_address_space_bytes'],8*GIB+native//2)
        self.assertEqual(forecast['worker_address_space_bytes'],2*forecast['per_hap_address_space_bytes']+2*GIB)

    def test_accounted_native_threads_allow_cpu_growth_without_changing_gpu_gate(self):
        from integrate_wgs_rna.cpu_executor import elastic_cpu_enabled
        env={'MFM_CPU_POOL_LEASE_FILE':'/tmp/unused','MFM_COMPACT_JOBS_REQUESTED':'auto',
             'MFM_CPU_UTILIZATION':'maximum','MFM_BLOSC_THREADS':'2'}
        with patch.dict(os.environ,env,clear=True):
            self.assertFalse(elastic_cpu_enabled('MFM_COMPACT_JOBS_REQUESTED'))
            self.assertTrue(elastic_cpu_enabled('MFM_COMPACT_JOBS_REQUESTED',thread_cost_accounted=True))
            with patch.dict(os.environ,MFM_ZARR_BACKEND='cupy'):
                self.assertFalse(elastic_cpu_enabled('MFM_COMPACT_JOBS_REQUESTED',thread_cost_accounted=True))

    def test_three_hifi_balance_profiles_and_memory_limits(self):
        c=SimpleNamespace(platform='hifi',threads=30,sort_threads=4,hap_jobs=2,sort_mem='2G')
        env={'MFM_WGS_ALIGN_REQUESTED':'auto','MFM_WGS_SORT_REQUESTED':'auto','MFM_WGS_HAP_REQUESTED':'auto'}
        lengths={1:2876676610,2:2876520048}
        for mode,expected in [('legacy',(17,8,4)),('mapping',(23,4,2)),('mapping_more',(26,2,1)),('auto',(26,2,1))]:
            p=alignment_plan(c,lengths,68,1485000000000,dict(env,MFM_WGS_CPU_BALANCE=mode))
            self.assertEqual((p['threads'],p['sort_threads'],p['input_decode_threads']),expected)
            self.assertEqual(p['reserved_cpus'],68)
        p=alignment_plan(c,lengths,68,120*GIB,env)
        self.assertEqual(p['reserved_cpus'],68);self.assertLessEqual(p['reserved_memory_bytes'],120*GIB)

    def test_vcf_returns_worker_resources_before_actual_global_merge(self):
        from test_public_vcf_parallel import PublicVcfParallelTests
        from integrate_wgs_rna.public_vcf import prepare_public_vcf
        from integrate_wgs_rna import public_vcf_parallel as vcf
        with tempfile.TemporaryDirectory() as raw:
            cfg=PublicVcfParallelTests().fixture(raw);pool=Path(raw)/'pool';(pool/'leases').mkdir(parents=True)
            lease=pool/'leases/a.json';lease.write_text(json.dumps(dict(branch='WGS',stage='03',cpus=8,memory_bytes=64*GIB)))
            original=vcf.run_command;observed=[]
            def checked(command,*a,**kw):
                if command[1]=='concat':
                    row=json.loads(lease.read_text());observed.append(row)
                    self.assertEqual(row['cpus'],2);self.assertEqual(row['memory_bytes'],4*GIB)
                return original(command,*a,**kw)
            with patch.dict(os.environ,{'MFM_CPU_BUDGET':'8','MFM_CPU_POOL_GRANTED':'8','MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_PUBLIC_VCF_JOBS':'auto'},clear=True),patch.object(vcf,'run_command',side_effect=checked):
                prepare_public_vcf(cfg)
            self.assertTrue(observed)
            self.assertEqual(json.loads(lease.read_text())['cpus'],2)

    def test_fasta_serial_validation_returns_compute_cpus_before_commands(self):
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir();(root/'waiters').mkdir();lease=root/'leases/a.json'
            lease.write_text(json.dumps(dict(cpus=68,total_cpus=68,memory_bytes=128*GIB)))
            cfg=PersonalFastaConfig(cohort='T',seq_vcf_dir=root,output_root=root,reference_fasta=root/'ref.fa',chroms=['chr1'],sequence_variant_policy='snp_only',parallel_contigs=True)
            def stop(*args):
                self.assertEqual(json.loads(lease.read_text())['cpus'],1)
                raise RuntimeError('직렬 준비만 검증')
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_RAM_LIMIT_BYTES':str(512*GIB),'MFM_FASTA_CONTIG_JOBS':'auto'},clear=True),patch('integrate_wgs_rna.personal_fasta.require_vcf',side_effect=stop):
                result=process_one_sample(cfg,'S',root,root,root/'ref.fa',['chr1'],'snp_only',False)
            self.assertEqual(result.status,'FAIL')

if __name__=='__main__':unittest.main()
