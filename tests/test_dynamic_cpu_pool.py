from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch

POOL = Path(__file__).resolve().parents[1] / "scripts" / "cpu_pool.py"
SPEC = importlib.util.spec_from_file_location("cpu_pool", POOL)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class DynamicCpuPoolTests(unittest.TestCase):
    def test_대기중68에서40으로바뀌면도구도최신예산을상속(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'leases').mkdir();(root/'waiters').mkdir()
            MODULE.atomic_json(root/'policy.json',{'live_limits':True,'cpu_pool_limit':68,'cpu_backfill':{'enabled':False}})
            active=root/'leases/active.json'
            MODULE.atomic_json(active,{'manager_pid':os.getpid(),'cpus':68,'memory_bytes':0})
            output=root/'environment.json'
            audit=root/'audit.jsonl'
            command='import json,os,sys; from pathlib import Path; Path(sys.argv[1]).write_text(json.dumps({k:os.environ[k] for k in ["MFM_CPU_BUDGET","MFM_CPU_POOL_GRANTED"]}))'
            env={**os.environ,'MFM_CPU_UTILIZATION':'maximum','MFM_CPU_POOL_LAYOUT':'original'}
            proc=subprocess.Popen([sys.executable,str(POOL),'--pool-dir',str(root),'--total-cpus','68',
                '--branch','WGS','--stage','03','--sample','S','--audit',str(audit),'--',sys.executable,'-c',command,str(output)],
                env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            try:
                deadline=time.monotonic()+5
                while not list((root/'waiters').glob('*.json')) and proc.poll() is None and time.monotonic()<deadline:
                    time.sleep(.01)
                self.assertTrue(list((root/'waiters').glob('*.json')))
                with (root/'pool.lock').open('a+') as lock:
                    MODULE.lock_pool(lock)
                    MODULE.atomic_json(root/'policy.json',{'live_limits':True,'cpu_pool_limit':40,'cpu_backfill':{'enabled':False}})
                    MODULE.atomic_json(active,{'manager_pid':os.getpid(),'cpus':36,'memory_bytes':0})
                logs,_=proc.communicate(timeout=8)
                self.assertEqual(proc.returncode,0,logs)
                recorded=json.loads(output.read_text())
                self.assertEqual(recorded,{'MFM_CPU_BUDGET':'40','MFM_CPU_POOL_GRANTED':'4'})
                usage=json.loads(audit.read_text())
                self.assertEqual((usage['total_cpus'],usage['maximum_cpus']),(40,40))
            finally:
                if proc.poll() is None:proc.kill();proc.communicate()

    def test_nvcomp_prepare_reserves_six_cpu_minimum(self):
        with patch.dict(os.environ, {'MFM_IO_BACKEND':'nvcomp','MFM_CPU_POOL_STAGE_MAX':'16','MFM_CPU_POOL_PROFILE':'throughput'}):
            self.assertEqual((6, 16), MODULE.stage_cap('WGS', '01', 68))
        with patch.dict(os.environ, {'MFM_IO_BACKEND':'cpu','MFM_CPU_POOL_STAGE_MAX':'16','MFM_CPU_POOL_PROFILE':'throughput'}):
            self.assertEqual((4, 16), MODULE.stage_cap('WGS', '01', 68))

    def setUp(self):
        # 기존 bounded/명시 thread 계약도 함께 유지되는지 검사한다.
        self.environment = patch.dict(os.environ, MFM_CPU_UTILIZATION='bounded')
        self.environment.start()
        self.addCleanup(self.environment.stop)

    def test_unset_profile_defaults_to_maximum_compute_caps(self):
        with patch.dict(os.environ, MFM_CPU_UTILIZATION='maximum', MFM_CPU_POOL_LAYOUT='original'):
            os.environ.pop('MFM_CPU_POOL_PROFILE', None)
            self.assertEqual((8, 72), MODULE.stage_cap('WGS', '10', 72))
            self.assertEqual((8, 72), MODULE.stage_cap('RNA', '16', 72))

    def test_illumina_realign_preparation_releases_unused_initial_cpus(self):
        with patch.dict(os.environ, MFM_CPU_UTILIZATION='maximum', MFM_CPU_POOL_LAYOUT='original',
                        MFM_REALIGN_PREPARE_MAX_CPUS='4'):
            self.assertEqual((4, 4), MODULE.stage_cap('WGS', '10', 72))
            self.assertEqual((2, 72), MODULE.stage_cap('WGS', '12', 72))
        with patch.dict(os.environ, MFM_CPU_UTILIZATION='maximum', MFM_CPU_POOL_LAYOUT='legacy',
                        MFM_REALIGN_PREPARE_MAX_CPUS='4'):
            self.assertEqual((2, 72), MODULE.stage_cap('WGS', '10', 72))

    def test_single_sample_window_contigs_use_pool_and_explicit_serial_stays_one(self):
        with patch.dict(os.environ, MFM_CPU_UTILIZATION='maximum', MFM_CPU_POOL_PROFILE='throughput', MFM_CPU_POOL_LAYOUT='original',MFM_WINDOW_CONTIG_JOBS='auto'):
            self.assertEqual((2,72),MODULE.stage_cap('WGS','11',72))
            with patch.dict(os.environ,MFM_WINDOW_CONTIG_JOBS='1'):
                self.assertEqual((1,1),MODULE.stage_cap('WGS','11',72))
        with patch.dict(os.environ, MFM_CPU_UTILIZATION='maximum', MFM_CPU_POOL_PROFILE='throughput', MFM_CPU_POOL_LAYOUT='legacy',MFM_WINDOW_CONTIG_JOBS='auto'):
            self.assertEqual((2,72),MODULE.stage_cap('WGS','09',72))

    def test_minimum_budget_starts_jobs_without_deadlock_or_oversubscription(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            audit = root / 'audit.jsonl'
            def launch(sample):
                env = {**os.environ, 'MFM_CPU_UTILIZATION': 'maximum', 'MFM_CPU_POOL_LAYOUT': 'original'}
                return subprocess.run([sys.executable, str(POOL), '--pool-dir', str(root),
                    '--total-cpus', '12', '--branch', 'WGS', '--stage', '10', '--sample', sample,
                    '--audit', str(audit), '--', sys.executable, '-c', 'import time; time.sleep(.2)'],
                    env=env, timeout=8, check=False).returncode
            with ThreadPoolExecutor(max_workers=3) as pool:
                self.assertEqual([0, 0, 0], list(pool.map(launch, ['A', 'B', 'C'])))
            rows = [json.loads(line) for line in audit.read_text().splitlines()]
            self.assertTrue(all(row['granted_cpus'] >= 8 for row in rows))
            for row in rows:
                moment = row['started_at'] + .01
                concurrent = sum(other['granted_cpus'] for other in rows
                                 if other['started_at'] <= moment < other['finished_at'])
                self.assertLessEqual(concurrent, 12)

    def test_throughput_scales_compute_stages_but_preserves_serial_caps(self):
        from unittest.mock import patch
        with patch.dict(os.environ, MFM_CPU_POOL_PROFILE="throughput", MFM_CPU_POOL_LAYOUT="original"):
            self.assertEqual(72, MODULE.stage_cap("WGS", "06", 72)[1])
            self.assertEqual(2, MODULE.stage_cap("WGS", "06", 72)[0])
            self.assertEqual(72, MODULE.stage_cap("RNA", "15", 72)[1])
            self.assertEqual(72, MODULE.stage_cap("RNA", "16", 72)[1])
            self.assertEqual(4, MODULE.stage_cap("RNA", "19", 72)[1])
        with patch.dict(os.environ, MFM_CPU_POOL_PROFILE="balanced"):
            self.assertEqual(8, MODULE.stage_cap("WGS", "06", 72)[1])

    def test_live_pool_policy_reaches_stage_environment_and_threads(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "policy.json").write_text(json.dumps({"profile": "throughput"}))
            result = subprocess.run([
                sys.executable, str(POOL), "--pool-dir", str(root), "--total-cpus", "10",
                "--branch", "WGS", "--stage", "06", "--sample", "S", "--",
                sys.executable, "-c", "import os; assert os.environ['MFM_HAPLOTAG_CPU_BUDGET']=='10'",
            ], check=False)
            self.assertEqual(0, result.returncode)

    def test_concurrent_leases_respect_total_budget(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            audit = root / "audit.jsonl"

            def launch(sample: str) -> int:
                return subprocess.run([
                    sys.executable, str(POOL), "--pool-dir", str(root / "pool"),
                    "--total-cpus", "12", "--branch", "WGS", "--stage", "10",
                    "--sample", sample, "--audit", str(audit), "--",
                    sys.executable, "-c", "import time; time.sleep(.5)",
                ], check=False).returncode

            with ThreadPoolExecutor(max_workers=3) as executor:
                self.assertEqual([0, 0, 0], list(executor.map(launch, ["S1", "S2", "S3"])))
            events = [json.loads(line) for line in audit.read_text().splitlines()]
            self.assertEqual(3, len(events))
            # 서버 부하로 순차 시작할 수 있으므로 전체 실행의 합계 대신 동시 예약을 검사함.
            self.assertTrue(all(4 <= row["granted_cpus"] <= 12 for row in events))
            for event in events:
                moment = event["started_at"]
                concurrent = sum(row["granted_cpus"] for row in events
                                 if row["started_at"] <= moment < row["finished_at"])
                self.assertLessEqual(concurrent, 12)
            self.assertTrue(all(row["exit_code"] == 0 for row in events))

    def test_rna_explicit_threads_are_capped_by_lease(self):
        previous = os.environ.get("MFM_CPU_POOL_LAYOUT")
        os.environ["MFM_CPU_POOL_LAYOUT"] = "legacy"
        try:
            command = ["python", "13_star_index.py", "--threads", "24"]
            self.assertEqual(["python", "13_star_index.py", "--threads", "6"],
                             MODULE.adjust_argv(command, "RNA", "13", 6))
            command = ["python", "11_rna_to_fastq.py", "--collate-threads", "12",
                       "--fastq-threads", "6"]
            adjusted = MODULE.adjust_argv(command, "RNA", "11", 8)
            self.assertLessEqual(int(adjusted[adjusted.index("--collate-threads") + 1]) +
                                 int(adjusted[adjusted.index("--fastq-threads") + 1]), 8)
            hap_command = ["python", "18_hap_variants.py", "--threads", "22"]
            self.assertEqual(["python", "18_hap_variants.py", "--threads", "7"],
                             MODULE.adjust_argv(hap_command, "RNA", "18", 7))
            self.assertEqual(22, MODULE.stage_cap("RNA", "18", 72)[1])
        finally:
            if previous is None:
                os.environ.pop("MFM_CPU_POOL_LAYOUT", None)
            else:
                os.environ["MFM_CPU_POOL_LAYOUT"] = previous

    def test_rna_auto_threads_grow_to_new_stage_grant(self):
        previous_layout = os.environ.get("MFM_CPU_POOL_LAYOUT")
        previous_fixed = os.environ.get("MFM_DYNAMIC_RNA_FIXED_STAGES")
        os.environ["MFM_CPU_POOL_LAYOUT"] = "legacy"
        os.environ["MFM_DYNAMIC_RNA_FIXED_STAGES"] = ""
        try:
            self.assertEqual(["python", "star.py", "--threads", "18"],
                             MODULE.adjust_argv(["python", "star.py", "--threads", "8"],
                                                "RNA", "14", 18))
            fastq = ["python", "fastq.py", "--collate-threads", "6",
                     "--fastq-threads", "2"]
            adjusted = MODULE.adjust_argv(fastq, "RNA", "11", 18)
            self.assertEqual(18,
                int(adjusted[adjusted.index("--collate-threads") + 1]) +
                int(adjusted[adjusted.index("--fastq-threads") + 1]))
            os.environ["MFM_DYNAMIC_RNA_FIXED_STAGES"] = "14"
            self.assertEqual(["python", "star.py", "--threads", "8"],
                             MODULE.adjust_argv(["python", "star.py", "--threads", "8"],
                                                "RNA", "14", 18))
        finally:
            if previous_layout is None:os.environ.pop("MFM_CPU_POOL_LAYOUT", None)
            else:os.environ["MFM_CPU_POOL_LAYOUT"] = previous_layout
            if previous_fixed is None:os.environ.pop("MFM_DYNAMIC_RNA_FIXED_STAGES", None)
            else:os.environ["MFM_DYNAMIC_RNA_FIXED_STAGES"] = previous_fixed

    def test_six_elastic_workers_share_72_cpu_budget(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            audit = root / "audit.jsonl"

            def launch(sample: str) -> int:
                env = os.environ.copy()
                env.update(MFM_LOCAL_JOBS="6", MFM_RESOURCE_PROFILE="auto",
                           MFM_ELASTIC_SAMPLE_SLOTS="true", MFM_RAM_LIMIT_BYTES="1500000000000")
                return subprocess.run([
                    sys.executable, str(POOL), "--pool-dir", str(root / "pool"),
                    "--total-cpus", "72", "--branch", "WGS", "--stage", "10",
                    "--sample", sample, "--audit", str(audit), "--",
                    sys.executable, "-c", "import time; time.sleep(.4)",
                ], env=env, check=False).returncode

            samples = [f"S{i}" for i in range(6)]
            with ThreadPoolExecutor(max_workers=6) as executor:
                self.assertEqual([0] * 6, list(executor.map(launch, samples)))
            events = [json.loads(line) for line in audit.read_text().splitlines()]
            # RAM 입장 제어로 일부 작업이 다음 차례에 시작할 수 있다.
            # 누적 grant 합 대신 같은 시점의 실제 예약 총합을 검증한다.
            for row in events:
                active = [other for other in events if other['started_at'] <= row['started_at'] + .001 < other['finished_at']]
                self.assertLessEqual(sum(other['granted_cpus'] for other in active),
                                     MODULE.pool_cpu_limit(root / 'pool', 72))
                self.assertLessEqual(sum(other['final_reserved_memory_bytes'] for other in active),
                                     MODULE.usable_memory(root / 'pool'))
            self.assertTrue(all(row["granted_cpus"] >= row["minimum_cpus"] for row in events))
            self.assertTrue(all(row["local_jobs"] == 6 for row in events))
            self.assertTrue(all(row["elastic_sample_slots"] for row in events))

    def test_failed_command_releases_lease_and_stale_waiter_is_pruned(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "pool"
            stale = root / "waiters" / "stale.json"
            stale.parent.mkdir(parents=True)
            stale.write_text(json.dumps({"pid": 2**30, "queued_at": 0,
                                         "min_cpus": 1, "max_cpus": 1,
                                         "branch": "WGS", "stage": "10", "sample": "lost"}))
            audit = Path(tmp) / "audit.jsonl"
            failed = subprocess.run([
                sys.executable, str(POOL), "--pool-dir", str(root), "--total-cpus", "6",
                "--branch", "WGS", "--stage", "10", "--sample", "FAIL",
                "--audit", str(audit), "--", sys.executable, "-c", "raise SystemExit(7)",
            ], check=False)
            self.assertEqual(7, failed.returncode)
            self.assertFalse(stale.exists())
            succeeded = subprocess.run([
                sys.executable, str(POOL), "--pool-dir", str(root), "--total-cpus", "6",
                "--branch", "WGS", "--stage", "10", "--sample", "NEXT",
                "--audit", str(audit), "--", sys.executable, "-c", "pass",
            ], check=False)
            self.assertEqual(0, succeeded.returncode)
            rows = [json.loads(line) for line in audit.read_text().splitlines()]
            self.assertEqual([7, 0], [row["exit_code"] for row in rows])

    def test_legacy_layout_maps_old_rna_threaded_stages(self):
        previous = os.environ.get("MFM_CPU_POOL_LAYOUT")
        os.environ["MFM_CPU_POOL_LAYOUT"] = "legacy"
        try:
            self.assertEqual(["python", "star.py", "--threads", "5"],
                             MODULE.adjust_argv(["python", "star.py", "--threads", "24"],
                                                "RNA", "14", 5))
            self.assertEqual(24, MODULE.stage_cap("RNA", "14", 72)[1])
            self.assertEqual(2, MODULE.stage_cap("RNA", "11", 72)[0])
        finally:
            if previous is None:
                os.environ.pop("MFM_CPU_POOL_LAYOUT", None)
            else:
                os.environ["MFM_CPU_POOL_LAYOUT"] = previous


    def test_hifi_adaptive_entry_minimum_two_in_both_layouts(self):
        plan={'adaptive_realign':True,'stage_estimates_bytes':{'personal_realign_entry':4*1024**3}}
        for layout,stage in [('original','10'),('legacy','08')]:
            with patch.dict(os.environ,{'MFM_CPU_POOL_PROFILE':'throughput','MFM_CPU_POOL_LAYOUT':layout,'MFM_LONGREAD_MEMORY_PLAN':json.dumps(plan),'MFM_REALIGN_PHASE_RESOURCES':'auto'}):
                self.assertEqual(MODULE.stage_cap('WGS',stage,68),(2,68))
                if layout=='original':
                    with patch.dict(os.environ,MFM_REALIGN_PREPARE_MAX_CPUS='6'):
                        self.assertEqual(MODULE.stage_cap('WGS',stage,68),(2,6))
                    with patch.dict(os.environ,MFM_REALIGN_PHASE_RESOURCES='0'):
                        self.assertEqual(MODULE.stage_cap('WGS',stage,68),(8,68))
        with patch.dict(os.environ,{'MFM_CPU_POOL_PROFILE':'throughput','MFM_CPU_POOL_LAYOUT':'original','MFM_LONGREAD_MEMORY_PLAN':'{}'}):
            self.assertEqual(MODULE.stage_cap('WGS','10',68),(8,68))



if __name__ == "__main__":
    unittest.main()
