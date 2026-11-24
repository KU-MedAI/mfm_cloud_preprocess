"""회수 가능한 캐시·축소된 실제CPU grant·일시 변경과HTTP403 복구를 검증함."""
from pathlib import Path
import json
import os
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import Mock, patch
from urllib.error import HTTPError
from urllib.request import Request

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'src')]
from integrate_wgs_rna.cgroup_memory import pressure_usage
from integrate_wgs_rna import alignment_record_stats as stats
from handoff import http_download as http, queue_resources
import rna_cache


class MemoryTests(unittest.TestCase):
    def test_완료입력버퍼를활성다운로드RAM으로중복예약하지않음(self):
        gib=1024**3
        with tempfile.TemporaryDirectory() as raw:
            with patch.object(queue_resources,'cgroup_limits',return_value=(None,None,None)), \
                 patch.object(queue_resources,'effective_memory',return_value=100*gib), \
                 patch.object(queue_resources,'active_ram',return_value=89*gib), \
                 patch.object(queue_resources,'live_leases',return_value={}):
                self.assertTrue(queue_resources.admit({},Path(raw),[],prefetch=[{}]*8)[0])
                with patch.object(queue_resources,'live_leases',return_value={'D':[{'branch':'DOWNLOAD','memory_bytes':3*gib}]}):
                    self.assertFalse(queue_resources.admit({},Path(raw),[],prefetch=[{}]*8)[0])

    def test_고정자원후처리도유휴대여를허용(self):
        from integrate_wgs_rna.pool_backfill import eligible_waiter
        for branch,stage in [('WGS','23'),('WGS','24'),('RNA','14')]:
            self.assertTrue(eligible_waiter({'pid':os.getpid(),'branch':branch,'stage':stage,
                                           'layout':'original','cpu_utilization':'maximum'}))
        self.assertFalse(eligible_waiter({'pid':os.getpid(),'branch':'RNA','stage':'99'}))
    def test_유휴CPU대여는캐시제외실측과예약을함께확인(self):
        from test_cpu_backfill import BackfillTests
        fixture=BackfillTests()
        with tempfile.TemporaryDirectory() as raw:
            root,observed,observer,guard=fixture.fixture(raw)
            observed.update(memory_current=9900,memory_working_bytes=100)
            self.assertTrue(fixture.lend(root,observed,observer,guard)['launched'])
    def test_비활성파일의일부만회수하고공유메모리와쓰기캐시는보호(self):
        with tempfile.TemporaryDirectory() as raw:
            folder = Path(raw)
            (folder/'memory.stat').write_text('file 700\nshmem 100\nfile_dirty 20\nfile_writeback 10\ninactive_file 500\n')
            self.assertEqual(pressure_usage(900, 1000, folder), 477)
            (folder/'memory.stat').write_text('file 700\nshmem 600\nfile_dirty 20\nfile_writeback 10\ninactive_file 500\n')
            self.assertEqual(pressure_usage(900, 1000, folder), 837)

    def test_stat누락시원래전체사용량으로되돌아감(self):
        with tempfile.TemporaryDirectory() as raw:
            self.assertEqual(pressure_usage(900, 1000, raw), 900)
            self.assertIsNone(pressure_usage(None, 1000, raw))

    def test_다운로드입장도캐시제외실측을사용하지만예약은지킴(self):
        with tempfile.TemporaryDirectory() as raw:
            config={'data_root':raw}
            with patch.object(queue_resources,'cgroup_limits',return_value=(9900,10000,None)), \
                 patch.object(queue_resources,'pressure_usage',return_value=100), \
                 patch.object(queue_resources,'effective_memory',return_value=1000), \
                 patch.object(queue_resources,'memory_plan',return_value={'reservation_bytes':100}), \
                 patch.object(queue_resources,'live_leases',return_value={}):
                self.assertTrue(queue_resources.download_admit(config,Path(raw))[0])
                with patch.object(queue_resources,'live_leases',return_value={'busy':[{'memory_bytes':850}]}):
                    self.assertFalse(queue_resources.download_admit(config,Path(raw))[0])


class GrantTests(unittest.TestCase):
    def test_CPUcoverage는부모와worker3CPU부터시작하고최대는유지(self):
        import importlib.util
        spec=importlib.util.spec_from_file_location('coverage_entry_pool',ROOT/'scripts/cpu_pool.py')
        pool=importlib.util.module_from_spec(spec);spec.loader.exec_module(pool)
        with patch.dict(os.environ,{'MFM_CPU_POOL_PROFILE':'throughput','MFM_CPU_UTILIZATION':'maximum','MFM_COVERAGE_BACKEND':'cpu','MFM_CPU_POOL_LAYOUT':'original'}):
            minimum,maximum=pool.stage_cap('RNA','18',68)
            self.assertEqual(minimum,3);self.assertGreaterEqual(maximum,3)
            with patch.dict(os.environ,{'MFM_COVERAGE_BACKEND':'cupy'}):self.assertEqual(pool.stage_cap('RNA','18',68)[0],4)

    def test_실제sort버퍼로계산하고고스레드는피크를올림(self):
        from integrate_wgs_rna.prepare_memory import transcode_plan
        self.assertEqual(transcode_plan(19,8,'2G')['reservation_bytes'],32*1024**3)
        self.assertGreater(transcode_plan(68,33,'2G')['reservation_bytes'],64*1024**3)

    def test_단일longBAM만예약을줄이고CPU는보존(self):
        from integrate_wgs_rna.prepare_memory import reserve_transcode
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'leases').mkdir()
            lease=root/'leases/self.json';lease.write_text(json.dumps({'memory_bytes':64*1024**3,'cpus':19}))
            contract=root/'input.json';contract.write_text(json.dumps({'platform':'hifi','source':{'alignment_states':['aligned'],'files':[{'format':'bam'}]}}))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_PREPARE_INPUT_CONTRACT':str(contract)}), \
                 patch('integrate_wgs_rna.prepare_memory.read_pool_records',return_value=[]), \
                 patch('integrate_wgs_rna.prepare_memory.usable_memory',return_value=128*1024**3):
                reserve_transcode(19,8,'2G')
            value=json.loads(lease.read_text());self.assertEqual(value['memory_bytes'],32*1024**3);self.assertEqual(value['cpus'],19)

    def test_과거환경16CPU대신현재예약2CPU로전수집계(self):
        with tempfile.TemporaryDirectory() as raw:
            path=Path(raw)/'read.bam';path.write_bytes(b'test')
            lease=Path(raw)/'lease.json';lease.write_text(json.dumps({'cpus':2,'memory_bytes':1<<30}))
            counts={'records':1,'mapped_records':1,'primary_records':1,'secondary_records':0,
                    'supplementary_records':0,'duplicate_records':0,'primary_duplicate_records':0,
                    'hp_tagged_records':0,'ps_tagged_records':0}
            budget=SimpleNamespace(capacity=Mock(return_value=2))
            with patch.dict(os.environ,{'MFM_CPU_POOL_LEASE_FILE':str(lease),'MFM_RECORD_STATS_SCAN_MODE':'single_pass'}), \
                 patch.object(stats,'cpu_budget',return_value=16), \
                 patch.object(stats,'WorkerBudget',return_value=budget), \
                 patch('integrate_wgs_rna.record_stats_parallel.native_tool',return_value=Path('/fake')), \
                 patch('integrate_wgs_rna.record_stats_parallel.run_scan',return_value={'counts':counts,'execution':{}}) as scan:
                self.assertEqual(stats.native_record_stats(path),counts)
                self.assertEqual(scan.call_args.args[1],2)
                self.assertEqual(budget.capacity.call_args.args[0],2)


class FingerprintTests(unittest.TestCase):
    def test_한번변경된파일은새버전으로다시읽어일관성확인(self):
        with tempfile.TemporaryDirectory() as raw:
            path=Path(raw)/'index';path.write_bytes(b'stable')
            original=rna_cache._stat(path);changed={**original,'ctime_ns':original['ctime_ns']+1}
            with patch.object(rna_cache,'_stat',side_effect=[original,changed,changed,changed]),patch.object(rna_cache.time,'sleep'):
                self.assertEqual(rna_cache.file_fingerprint(path)['ctime_ns'],changed['ctime_ns'])

    def test_계속변하는파일은3회후에도오류로유지(self):
        with tempfile.TemporaryDirectory() as raw:
            path=Path(raw)/'index';path.write_bytes(b'stable')
            original=rna_cache._stat(path);changed={**original,'ctime_ns':original['ctime_ns']+1}
            with patch.object(rna_cache,'_stat',side_effect=[original,changed]*3) as calls,patch.object(rna_cache.time,'sleep'):
                with self.assertRaises(rna_cache.FingerprintChanged):rna_cache.file_fingerprint(path)
                self.assertEqual(calls.call_count,6)


class HTTPTests(unittest.TestCase):
    def test_대체호스트403시원래공식호스트를다시확인(self):
        url='https://ftp.sra.ebi.ac.uk/vol1/run/ERR1/file'
        with patch.dict(http._EBI_HTTPS_MIRRORS,{'ftp.sra.ebi.ac.uk':'ftp.ebi.ac.uk'},clear=True), \
             patch.object(http,'urlopen',side_effect=[HTTPError(url,403,'',{},None),object()]) as opened:
            http._open(Request(url,headers={'Range':'bytes=1-2','If-Match':'v1'}),10)
            self.assertEqual(opened.call_args_list[1].args[0].full_url,url)
            self.assertNotIn('ftp.sra.ebi.ac.uk',http._EBI_HTTPS_MIRRORS)

    def test_ENA403만제한횟수재시도하고404와다른호스트는거부(self):
        url='https://ftp.sra.ebi.ac.uk/vol1/run/ERR1/file'
        self.assertTrue(http._retry_http(HTTPError(url,403,'',{},None),url))
        self.assertFalse(http._retry_http(HTTPError(url,404,'',{},None),url))
        self.assertFalse(http._retry_http(HTTPError(url,403,'',{},None),'https://other.example/file'))
        with patch.object(http,'_identify_once',side_effect=HTTPError(url,403,'',{},None)) as calls,patch.object(http.time,'sleep'):
            with self.assertRaises(HTTPError):http.identify(url,attempts=2)
            self.assertEqual(calls.call_count,2)


class ResumeTests(unittest.TestCase):
    def test_다운로드재개요청은일회만실패횟수초기화(self):
        from test_ready_queue import Harness
        with tempfile.TemporaryDirectory() as raw:
            h=Harness(raw,['S']);h.queue.RECOVERY_STAGES={'01'};s=h.scheduler()
            s.ledger['downloads']['S']={'attempts':3,'deferred':True,'next_retry_at':0,'last_error':'HTTP403'}
            h.save(h.base/'작업재개요청.json',{'samples':{'S':{'stage':'01','requested_at':1,'retry_download':True}}})
            s.step();self.assertIn('S',s.downloads)
            self.assertEqual(s.ledger['downloads']['S']['attempts'],1)
            s.downloads['S'][1].result=1;s.step()
            self.assertEqual(s.ledger['downloads']['S']['attempts'],1)
            self.assertEqual(len(s.ledger['downloads']['S']['recovery_history']),1)

    def test_수신중재개요청은프로세스를유지하고예산만갱신(self):
        from test_ready_queue import Harness
        with tempfile.TemporaryDirectory() as raw:
            h=Harness(raw,['S']);h.queue.RECOVERY_STAGES={'01'};s=h.scheduler();s.step()
            pid=s.downloads['S'][1].pid;s.ledger['downloads']['S'].update(attempts=3)
            h.save(h.base/'작업재개요청.json',{'samples':{'S':{'stage':'01','requested_at':2,'retry_download':True}}})
            s.step();self.assertEqual(s.downloads['S'][1].pid,pid)
            self.assertEqual(s.ledger['downloads']['S']['attempts'],1)

    def test_이관샘플도현재RNA시작단계를보존(self):
        from handoff.state_bindings import StateBindings
        with tempfile.TemporaryDirectory() as raw:
            root=Path(raw);(root/'references.json').write_text('{}')
            prior={'run_id':'old','state_root':str(root),'sample_pipeline_env':{'S':{'RNA_START_AT':'13'}}}
            source=root/'old.json';source.write_text(json.dumps(prior))
            binding=StateBindings.__new__(StateBindings)
            binding.items={'S':{'config':str(source),'selection':str(root/'selection.tsv'),'task_index':0}}
            binding.config={'cpu_pool_limit':68,'memory_limit_bytes':1000,'sample_cpus':68,'local_jobs':1,
                            'resource_resolution':{},'resources_resolved':True,
                            'sample_pipeline_env':{'S':{'RNA_START_AT':'15'}}}
            def save(path,value):path.write_text(json.dumps(value))
            path,_,_=binding.launch('S','worker',{'stage':'23'},SimpleNamespace(atomic_json=save),root/'launch')
            config=json.loads(path.read_text())
            self.assertEqual(config['sample_pipeline_env']['S']['RNA_START_AT'],'15')
            self.assertEqual(config['sample_pipeline_env']['S']['WGS_START_AT'],'23')
            self.assertEqual(config['run_id'],'old')

    def test_최종검증단계부터완료WGS를재사용할수있음(self):
        from handoff.queue_recovery import stages,request_recovery
        self.assertIn('23',stages(True));self.assertIn('23',stages(False))
        config={'sample_pipeline_env':{'S':{'RNA_START_AT':'15'}}}
        request_recovery(config,'S',{'stage':'23'})
        self.assertEqual(config['sample_pipeline_env']['S'],{'RNA_START_AT':'15','WGS_START_AT':'23'})


if __name__=='__main__':unittest.main()
