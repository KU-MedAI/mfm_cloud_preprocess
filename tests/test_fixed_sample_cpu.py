"""고정 샘플 CPU의 수명·하위 배분·전역 RAM 및 처리 배치의 실제 동등성을 검증함."""
import importlib.util
import json
import os
from pathlib import Path
import random
import shlex
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'src')]
from handoff.resources import resolve
from integrate_wgs_rna.sample_cpu import sample_slot, process_start, global_claim, owner_alive
from integrate_wgs_rna.memory_policy import usable_memory, GIB
from integrate_wgs_rna.cpu_executor import WorkerBudget
from integrate_wgs_rna.processing_buffers import apply, minimap_arguments, star_arguments


class FixedSampleCpuTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.env = dict(os.environ, PYTHONPATH=str(ROOT / 'src') + ':' + str(ROOT),
                        MFM_CPU_UTILIZATION='maximum', MFM_CPU_POOL_LAYOUT='original',
                        MFM_RAM_LIMIT_BYTES=str(512*GIB), MFM_CPU_POOL_PROFILE='throughput')
        for key in ['MFM_FIXED_SAMPLE_TOKEN', 'MFM_FIXED_SAMPLE_CPUS', 'MFM_CPU_POOL_LEASE_FILE', 'MFM_SAMPLE_RAM_LIMIT_BYTES']:
            self.env.pop(key, None)

    def config(self, cpus=2, total=4):
        return dict(cpu_allocation_mode='fixed_sample', sample_cpus=cpus, cpu_pool_limit=total,
                    memory_limit_bytes=512*GIB, state_root=str(self.root), run_id='test',
                    pipeline_env={'MFM_CPU_POOL_DIR': str(self.root / 'pool')})

    def wait_file(self, file, process, timeout=10):
        until = time.monotonic() + timeout
        while not file.exists() and process.poll() is None and time.monotonic() < until:
            time.sleep(.02)
        if not file.exists():
            output = process.communicate(timeout=2)[0] if process.poll() is not None else ''
            self.fail(f'준비신호 미수신: {file} {output}')

    def launch_holder(self, name):
        marker, release = self.root/(name+'.ready'), self.root/(name+'.release')
        code = ('import os,json,sys,time; from pathlib import Path; '
                'Path(sys.argv[1]).write_text(json.dumps({"cpus":os.environ["MFM_FIXED_SAMPLE_CPUS"],'
                '"ids":sorted(os.sched_getaffinity(0))})); '
                '\nwhile not Path(sys.argv[2]).exists(): time.sleep(.02)')
        proc = subprocess.Popen([sys.executable, '-m', 'integrate_wgs_rna.sample_cpu',
            '--sample', name, '--sample-cpus', '2', '--total-cpus', '4', '--pool-dir', str(self.root/'pool'),
            '--', sys.executable, '-c', code, str(marker), str(release)],
            env=self.env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        self.addCleanup(lambda: self.stop(proc, release))
        return proc, marker, release

    @staticmethod
    def stop(proc, release):
        release.touch()
        try:
            proc.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill(); proc.communicate()

    def test_two_samples_hold_fixed_cpus_third_waits_and_release(self):
        a, am, ar = self.launch_holder('A'); self.wait_file(am, a)
        b, bm, br = self.launch_holder('B'); self.wait_file(bm, b)
        self.assertTrue(set(json.loads(am.read_text())['ids']).isdisjoint(json.loads(bm.read_text())['ids']))
        c, cm, cr = self.launch_holder('C')
        time.sleep(.25)
        self.assertFalse(cm.exists())
        records = [json.loads(p.read_text()) for p in (self.root/'pool/leases').glob('*.json')]
        self.assertEqual(4, sum(global_claim(r) for r in records))
        ar.touch(); a.communicate(timeout=5)
        self.wait_file(cm, c)
        self.assertEqual('2', json.loads(bm.read_text())['cpus'])
        br.touch(); cr.touch(); b.communicate(timeout=5); c.communicate(timeout=5)
        self.assertEqual([], list((self.root/'pool/leases').glob('sample_*.json')))

    def test_stage_releases_do_not_release_sample_and_no_double_cpu_accounting(self):
        with patch.dict(os.environ, self.env, clear=True), sample_slot(self.config(), 'A') as env:
            pool = self.root/'pool'
            for stage, expected in [('11', 1), ('03', 2)]:
                output = self.root/(stage+'.json')
                code = 'import os,json,sys; from pathlib import Path; Path(sys.argv[1]).write_text(json.dumps({"grant":os.environ["MFM_CPU_POOL_GRANTED"],"pool":os.environ["MFM_CPU_BUDGET"],"affinity":len(os.sched_getaffinity(0))}))'
                result = subprocess.run([sys.executable, str(ROOT/'scripts/cpu_pool.py'), '--pool-dir', str(pool),
                    '--total-cpus', '4', '--branch', 'WGS', '--stage', stage, '--sample', 'A', '--',
                    sys.executable, '-c', code, str(output)], env=dict(self.env, **env, MFM_WINDOW_CONTIG_JOBS='1'), capture_output=True, text=True, timeout=12)
                self.assertEqual(0, result.returncode, result.stdout+result.stderr)
                data = json.loads(output.read_text())
                self.assertEqual(expected, int(data['grant']))
                self.assertEqual(('2', 2), (data['pool'], data['affinity']))
                owners = [json.loads(p.read_text()) for p in (pool/'leases').glob('*.json')]
                self.assertEqual(2, sum(global_claim(r) for r in owners))
                self.assertEqual(1, len(owners))

    def test_failure_returns_sample_slot(self):
        result = subprocess.run([sys.executable, '-m', 'integrate_wgs_rna.sample_cpu', '--sample', 'A',
            '--sample-cpus', '2', '--total-cpus', '4', '--pool-dir', str(self.root/'pool'), '--',
            sys.executable, '-c', 'raise SystemExit(7)'], env=self.env, capture_output=True, text=True, timeout=10)
        self.assertEqual(7, result.returncode, result.stdout+result.stderr)
        self.assertEqual([], list((self.root/'pool/leases').glob('sample_*.json')))

    def test_worker_uses_own_cpu_domain_and_global_ram(self):
        with patch.dict(os.environ, self.env, clear=True), sample_slot(self.config(4, 8), 'A') as env:
            pool = self.root/'pool'
            own = pool/'leases/task.json'
            own.write_text(json.dumps(dict(cpus=1, total_cpus=8, memory_bytes=4*GIB,
                fixed_sample_token=env['MFM_FIXED_SAMPLE_TOKEN'])))
            (pool/'leases/other.json').write_text(json.dumps(dict(cpus=4, total_cpus=8, memory_bytes=16*GIB,
                fixed_sample_token='B', sample_cpu_owner=True, manager_pid=os.getpid(), manager_start=process_start(os.getpid()), pool_dir=str(pool))))
            values=dict(env, MFM_CPU_POOL_LEASE_FILE=str(own))
            with patch.dict(os.environ, values):
                self.assertEqual(512*GIB*9//10, usable_memory(pool))
                waiters=pool/'waiters';waiters.mkdir(exist_ok=True)
                (waiters/'other.json').write_text(json.dumps(dict(fixed_sample_token='B',min_cpus=1,memory_bytes=0)))
                self.assertEqual(4, WorkerBudget(1, 8, True).capacity(8))
                self.assertEqual(4, json.loads(own.read_text())['cpus'])
                foreign = pool/'leases/other.json'
                row = json.loads(foreign.read_text()); row['memory_bytes'] = 440*GIB; foreign.write_text(json.dumps(row))
                self.assertEqual(1, WorkerBudget(1,8,True,memory_per_worker=8*GIB).capacity(8))
                (waiters/'other.json').unlink()
            own.unlink(); (pool/'leases/other.json').unlink()

    def test_sample_entry_honors_logical_ceiling_with_waiting_parent(self):
        from integrate_wgs_rna.sample_cpu import available
        pool=self.root/'pool';(pool/'leases').mkdir(parents=True)
        (pool/'policy.json').write_text(json.dumps({'cpu_backfill':{'enabled':True,'max_extra_cpus':16}}))
        rows=[dict(cpus=66),dict(cpus=16,cpu_backfill={'cpus':16}),dict(cpus=1,realign_phase='align:waiting')]
        for index,row in enumerate(rows):
            (pool/'leases'/f'{index}.json').write_text(json.dumps(row))
        with patch.dict(os.environ,self.env,clear=True):
            self.assertEqual(1,available(pool,68))

    def test_owner_kept_until_detached_stage_exits(self):
        pool=self.root/'pool';(pool/'leases').mkdir(parents=True)
        record=dict(sample_cpu_owner=True, fixed_sample_token='A', manager_pid=999999999,
                    manager_start=1, pool_dir=str(pool))
        child=pool/'leases/child.json'
        child.write_text(json.dumps(dict(fixed_sample_token='A', manager_pid=os.getpid())))
        self.assertTrue(owner_alive(record))
        child.unlink()
        self.assertFalse(owner_alive(record))

    def test_pool_prune_keeps_fixed_owner_of_detached_stage(self):
        from integrate_wgs_rna.pool_leases import prune_leases
        pool=self.root/'pool';(pool/'leases').mkdir(parents=True)
        owner=pool/'leases/sample_A.json'
        owner.write_text(json.dumps(dict(sample_cpu_owner=True,fixed_sample_token='A',
            manager_pid=999999999,manager_start=1,pool_dir=str(pool),cpus=2)))
        child=pool/'leases/child.json'
        child.write_text(json.dumps(dict(fixed_sample_token='A',manager_pid=os.getpid(),cpus=1)))
        self.assertEqual([],prune_leases(pool));self.assertTrue(owner.exists())
        child.unlink();self.assertEqual(['sample_A.json'],prune_leases(pool))

    def test_fixed_config_and_buffer_validation(self):
        config=dict(self.config(32,68), resource_profile='throughput', local_jobs='auto', download_jobs=1,
                    slurm={'cpus':32,'concurrency':2})
        result=resolve(config,cpu_count=68,available_memory_bytes=1024*GIB,free_storage_bytes=1<<50)
        self.assertEqual((32,2),(result['sample_cpus'],result['local_jobs']))
        with self.assertRaises(ValueError):resolve(dict(config,local_jobs=3),cpu_count=68,available_memory_bytes=1024*GIB,free_storage_bytes=1<<50)
        with self.assertRaises(ValueError):resolve(dict(config,sample_cpus='auto'),cpu_count=68,available_memory_bytes=1024*GIB,free_storage_bytes=1<<50)
        values=apply({'processing_buffers':{'sort_memory_per_thread':'1G','minimap_batch_bases':250000000,'windows_per_chunk':1,'ase_batch_fragments':25000}})
        self.assertEqual('1G', values['pipeline_env']['STAGE04_SORT_MEM'])
        self.assertEqual(' -K 250000000',minimap_arguments(values['pipeline_env']))
        self.assertEqual(['--limitIObufferSize','8000000','8000000'],star_arguments({'MFM_STAR_READ_BUFFER_BYTES':'8000000'}))
        with self.assertRaises(ValueError):apply({'processing_buffers':{'minimap_batch_bases':0}})
        with self.assertRaises(ValueError):apply({'processing_buffers':{'minimap_batch_bases':500000001}})
        with self.assertRaises(ValueError):apply({'processing_buffers':{'unknown':1}})
        with self.assertRaises(ValueError):apply({'processing_buffers':{'sort_memory_per_thread':'1G'},'pipeline_env':{'WGS_SORT_MEM':'2G'}})

    def test_cli_override_selects_fixed_mode_and_slot_limit(self):
        config=dict(self.config(16,68),cpu_allocation_mode='dynamic',resource_profile='throughput',
                    local_jobs='auto',download_jobs=1,slurm={'cpus':16,'concurrency':4})
        value=resolve(config,{'cpu_allocation_mode':'fixed_sample'},cpu_count=68,
                      available_memory_bytes=1024*GIB,free_storage_bytes=1<<50)
        self.assertEqual(('fixed_sample',16,4),(value['cpu_allocation_mode'],value['sample_cpus'],value['local_jobs']))

    def test_shell_direct_entry_uses_same_sample_slot(self):
        output=self.root/'entry.json'
        code='import os,json,sys; from pathlib import Path; Path(sys.argv[1]).write_text(json.dumps({"mode":os.environ["MFM_CPU_ALLOCATION_MODE"],"cpus":os.environ["MFM_SAMPLE_CPU_BUDGET"],"affinity":len(os.sched_getaffinity(0))}))'
        script='source ./lib/fixed_cpu_entry.sh\npipeline_fixed_cpu_entry A '+shlex.join([sys.executable,'-c',code,str(output)])+'\n'
        env=dict(self.env,MFM_CPU_ALLOCATION_MODE='fixed_sample',MFM_SAMPLE_CPU_BUDGET='2',
                 MFM_CPU_BUDGET='4',MFM_CPU_POOL_DIR=str(self.root/'pool'),PYTHON_BIN=sys.executable)
        run=subprocess.run(['bash'],input=script,cwd=ROOT,env=env,text=True,capture_output=True,timeout=10)
        self.assertEqual(0,run.returncode,run.stderr)
        self.assertEqual({'mode':'fixed_sample','cpus':'2','affinity':2},json.loads(output.read_text()))

    def test_star_large_buffers_are_in_ram_forecast(self):
        from integrate_wgs_rna.rna_memory import estimate_star_memory
        normal=estimate_star_memory('alignReads',64,index_bytes=30*GIB)
        large=estimate_star_memory('alignReads',64,index_bytes=30*GIB,io_buffer_bytes=2*GIB)
        self.assertGreater(large['memory_bytes'],normal['memory_bytes'])
        self.assertEqual(128*GIB,large['components']['thread_buffers'])

    def test_minimap_smaller_batches_preserve_real_alignments_and_tags(self):
        import pysam
        rng=random.Random(42);sequence=''.join(rng.choice('ACGT') for _ in range(12000))
        reference=self.root/'ref.fa';reference.write_text('>chr1\n'+sequence+'\n')
        fastq=self.root/'reads.fq'
        fastq.write_text(''.join(f'@r{i} HP:i:1 PS:i:17\n{sequence[i*500:i*500+400]}\n+\n'+('I'*400)+'\n' for i in range(20)))
        outputs=[]
        for batch in [1000,500000000]:
            sam=self.root/f'{batch}.sam'
            with sam.open('w') as out:
                subprocess.run(['minimap2','-ax','map-hifi','-y','-t','1','-K',str(batch),str(reference),str(fastq)],
                    stdout=out,stderr=subprocess.PIPE,check=True)
            with pysam.AlignmentFile(str(sam),'r') as handle:
                outputs.append(sorted((r.query_name,r.flag,r.reference_start,r.mapping_quality,r.cigarstring,r.query_sequence,tuple(sorted(r.tags))) for r in handle))
        self.assertEqual(outputs[0],outputs[1])
        self.assertEqual(20,len(outputs[0]))


if __name__ == '__main__':
    unittest.main()
