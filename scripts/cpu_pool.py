#!/usr/bin/env python3
"""단계 경계에서 CPU 예산을 공정하게 임대하고 사용량을 기록한다.

실행 중인 정렬기의 thread 수를 강제로 바꾸지는 않는다. 동시에 시작한
단계들끼리 가용 코어를 나누고, 단계가 끝나면 임대를 반환한다.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import resource
import signal
import subprocess
import sys
import tempfile
import time
import uuid
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from integrate_wgs_rna.pool_leases import prune_leases, process_alive, fitting_waiters, pool_cpu_limit, regular_cpu_claim, native_allocator_ready
from integrate_wgs_rna.memory_policy import GIB, memory_limit, usable_memory, protect_pool_memory, available_memory
from integrate_wgs_rna.cpu_loan_policy import available_loan_cpus, available_regular_cpus
from integrate_wgs_rna.resource_fairness import priority_request, fair_cpu_capacity, request_key


SETTLE_SECONDS = 0.35
POLL_SECONDS = 0.10


def lock_pool(handle):
    """공유 FS의 일시 EAGAIN도 blocking lock 의미에 맞게 재시도함."""
    while True:
        try:
            fcntl.flock(handle,fcntl.LOCK_EX)
            return
        except BlockingIOError:
            time.sleep(POLL_SECONDS)


def stage_memory(branch: str, stage: str) -> int:
    """단계 시작 전 예상 RAM을 예약한다. compact 세부 예측은 입력 준비 후 갱신한다."""
    if branch.upper() == 'BENCH' and 'MFM_BENCH_ENTRY_MEMORY_BYTES' in os.environ:
        value = int(os.environ['MFM_BENCH_ENTRY_MEMORY_BYTES'])
        if value < 4 * GIB:
            raise ValueError('BENCH 시작 RAM은 최소4GiB여야 함')
        return value
    if branch.upper() == 'IO':
        return 4 * GIB  # bounded IO 배치·CUDA host 예상1GB 미만에 여유 포함함.
    legacy = os.environ.get('MFM_CPU_POOL_LAYOUT') == 'legacy'
    if branch.upper() == 'WGS' and stage == ('07' if legacy else '09'):
        declared = os.environ.get('MFM_LONGREAD_MEMORY_PLAN')
        if declared and json.loads(declared).get('schema') != 'mfm_longread_stage_memory_v1':
            raise ValueError('HiFi 단계 RAM 계획 schema가 다름')
        # FASTA 부모4GiB로 입장하며 contig 계산은 2CPU/4GiB씩 다시 예약함.
        return 4 * GIB
    if branch.upper()=='WGS' and os.environ.get('MFM_WGS_PAYLOAD')=='sequence_only':
        if stage=='03':return 20*GIB  # 22chr×768MiB+parent1GiB와10%여유를 포함함.
        if stage in {'23','24'}:
            return 4*GIB  # 부모 입장 RAM. 독립 서열 검증은 worker별 2GiB를 계산 전에 재예약함.
    if branch.upper()=='RNA' and stage==('14' if legacy else '16') and os.environ.get('MFM_RNA_BACKEND','cpu')=='parabricks':
        return (int(os.environ.get('MFM_GPU_HOST_MEMORY_BYTES',str(128*GIB)))+8*GIB)*11//10
    declared = os.environ.get('MFM_LONGREAD_MEMORY_PLAN')
    if declared:
        plan = json.loads(declared)
        if plan.get('schema') != 'mfm_longread_stage_memory_v1':
            raise ValueError('HiFi 단계 RAM 계획 schema가 다름')
        estimates = plan['stage_estimates_bytes']
        if branch.upper() == 'WGS' and stage == '03':
            from integrate_wgs_rna.public_vcf_parallel import entry_memory
            cpus = min(int(plan['cpu_upper']), int(os.environ.get('MFM_CPU_BUDGET', plan['cpu_upper'])))
            return entry_memory(cpus)
        if branch.upper() == 'WGS':
            mapping = {'01':'prepare_alignment',
                       '06':('haplotag' if 'haplotag' in estimates else 'public_vcf_and_haplotag'),
                       ('07' if legacy else '09'):'consensus',
                       ('08' if legacy else '10'):('personal_realign_entry' if 'personal_realign_entry' in estimates else 'personal_minimap_sort_qc'),
                       ('10' if legacy else '12'):'compact_entry'}
            label = mapping.get(stage)
        else:
            label = ('rna_star_index' if stage==('13' if legacy else '15') else
                     'rna_star_align' if stage==('14' if legacy else '16') else None)
            if label and label not in estimates:label='rna_star_and_sort'
        if label in estimates:
            value = estimates[label]
            if isinstance(value,bool) or not isinstance(value,int) or value < 1:
                raise ValueError('HiFi 단계 RAM 예측은 양수여야 함')
            return value
    if branch.upper() == 'WGS':
        if (stage==('08' if legacy else '10') and os.environ.get('MFM_REALIGN_PREPARE_MAX_CPUS')
                and os.environ.get('MFM_REALIGN_PHASE_RESOURCES','auto')!='0'
                and os.environ.get('WGS_REALIGN_JOBS','1')=='1'):
            # 준비 진입 뒤 Python이 입력 worker 수·버퍼 RAM을 다시 예약한다.
            return 4*GIB
        if stage == ('10' if legacy else '12'):
            return 4 * GIB
        if stage == ('09' if legacy else '11'):
            return 8 * GIB
        if stage in {'04', '08', '10'} or (legacy and stage == '07'):
            return (384 if os.environ.get('MFM_MEMORY_PROFILE') == 'high' else 192) * GIB
        return 64 * GIB
    if stage == ('15' if legacy else '17'):
        from integrate_wgs_rna.markdup_resources import markdup_memory_bytes
        return markdup_memory_bytes(int(os.environ.get('RNA_MARKDUP_HEAP_GB',os.environ.get('MFM_RNA_MARKDUP_HEAP_GB_DEFAULT','8'))))
    threaded = {'13', '14'} if legacy else {'15', '16'}
    return (128 if stage in threaded else 64) * GIB


def alive(pid: int | None) -> bool:
    return process_alive(pid)


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


def atomic_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(data, stream, ensure_ascii=False, sort_keys=True)
            stream.write("\n")
        os.replace(tmp, path)
    finally:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass


def append_jsonl(path: Path, row: dict) -> None:
    if not path:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        stream.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
        stream.flush()
        fcntl.flock(stream, fcntl.LOCK_UN)


def stage_cap(branch: str, stage: str, total: int) -> tuple[int, int]:
    """I/O/짧은 단계가 CPU를 독점하지 않게 하고, 큰 단계는 동시요청에 나눈다."""
    label = branch.upper()
    fixed = bool(os.environ.get('MFM_FIXED_SAMPLE_TOKEN')) or os.environ.get('MFM_CPU_ALLOCATION_MODE')=='fixed_sample'
    if fixed:
        total=min(total,int(os.environ.get('MFM_FIXED_SAMPLE_CPUS',os.environ.get('MFM_SAMPLE_CPU_BUDGET',total))))
    if label == 'IO':
        requested = min(total,int(os.environ.get('MFM_CPU_POOL_STAGE_MAX',total))) if fixed else (2 if stage == 'cram_crc' else (6 if stage == 'bam_roundtrip' else 4))
        return min(total,requested), min(total,requested)
    legacy = os.environ.get("MFM_CPU_POOL_LAYOUT") == "legacy"
    if label == "WGS":
        caps = {"01": 8, "02": 8, "03": 16, "04": 48, "05": 4,
                "06": 8, "07": 48, "08": total, "09": 16, "10": 36,
                "11": 16, "12": total, "23": 4, "24": 4}
        if legacy:
            caps.update({"07": 16, "08": 36, "09": 16, "10": total})
    elif label == "RNA":
        caps = {"13": 24, "14": 2, "15": 24, "16": 24, "17": 4,
                "18": 48, "19": 4, "20": 22, "21": 4, "22": 24,
                "23": 4, "24": 4, "25": 4, "26": 4}
        if legacy:
            caps = {"11": 24, "12": 2, "13": 24, "14": 24, "15": 4,
                    "16": 48, "17": 4, "18": 22, "19": 4, "20": 24,
                    "21": 4, "22": 4, "23": 4, "24": 4}
    else:
        caps = {}
    if fixed:
        # 고정 샘플의 작은 단계 상한을 없애되 동시 HAP/DAG의 합산 배분은 유지함.
        caps={number:total for number in caps}
    # 처리량 우선 정책은 실제 병렬 계산 단계가 남은 전체 CPU를 사용할 수 있게 한다.
    # 짧거나 직렬인 단계의 상한은 유지한다. 기존 balanced 정책도 선택 가능하다.
    utilization = os.environ.get("MFM_CPU_UTILIZATION", "maximum")
    profile = os.environ.get("MFM_CPU_POOL_PROFILE", "throughput" if utilization == "maximum"
                             else os.environ.get("MFM_RESOURCE_PROFILE", "balanced"))
    if profile == "throughput":
        scalable = ({"01", "02", "03", "04", "06", "07", "08", "09", "10", "11", "12"}
                    if label == "WGS" else
                    ({"11", "13", "14", "16", "18", "19", "20"} if legacy else
                     {"13", "15", "16", "18", "20", "21", "22"}) if label == "RNA" else set())
        caps.update({number: total for number in scalable})
    from integrate_wgs_rna.parabricks_bam import enabled as parabricks_bam_enabled
    gpu_markdup=label=='RNA' and stage==('15' if legacy else '17') and parabricks_bam_enabled()
    if gpu_markdup and profile=='throughput':caps[stage]=total
    maximum = max(1, min(total, int(os.environ.get("MFM_CPU_POOL_STAGE_MAX", caps.get(stage, total)))))
    # 단일 샘플도 염색체×HAP process 병렬화가 가능하다. 명시 직렬/균형 정책은 유지함.
    if label == "WGS" and stage == ("09" if legacy else "11") and (
            os.environ.get('MFM_WINDOW_CONTIG_JOBS','auto')=='1' or profile!='throughput'):
        maximum = 1
    if (not fixed and label == "WGS" and stage == "10" and not legacy
            and os.environ.get("MFM_REALIGN_PREPARE_MAX_CPUS")):
        # Illumina 준비 구간은 적게 빌리고, 정렬 직전 RealignPhases가 재산정함.
        maximum = min(maximum, max(1, int(os.environ["MFM_REALIGN_PREPARE_MAX_CPUS"])))
    if label=='RNA' and stage==('14' if legacy else '16') and os.environ.get('MFM_RNA_BACKEND','cpu')=='parabricks':
        maximum=min(maximum,max(1,int(os.environ.get('MFM_GPU_HOST_CPUS_PER_HAP','16'))))
    paired_stage = "11" if legacy else "13"
    minimum = 2 if label == "RNA" and stage == paired_stage else 1
    if profile == "throughput":
        # 장시간 병렬 계산이 한 CPU로 고정된 채 시작되지 않게 한다.
        compute_minimum = ({"01": 4, "02": 2, "03": 2, "04": 8, "06": 2,
                            "07": 2, "08": 8, "09": 2, "10": 8, "11": 2, "12": 2}
                           if label == "WGS" else
                           ({"11": 2, "13": 8, "14": 8, "16": 4, "18": 4, "20": 2}
                            if legacy else {"13": 2, "15": 8, "16": 8, "18": 4, "20": 4, "22": 2})
                           if label == "RNA" else {})
        if label == "WGS" and legacy:
            compute_minimum.update({"07": 8, "08": 2, "09": 2, "10": 2})
        if label == 'RNA' and stage == ('16' if legacy else '18') and os.environ.get('MFM_COVERAGE_BACKEND','cpu') == 'cpu':
            compute_minimum[stage] = 3  # coverage worker 최대2CPU + 부모1CPU, 최대값은 유지함.
        realign_stage = '08' if legacy else '10'
        if label == 'WGS' and stage == realign_stage:
            plan = json.loads(os.environ.get('MFM_LONGREAD_MEMORY_PLAN') or '{}')
            if (plan.get('adaptive_realign') is True
                    and plan.get('stage_estimates_bytes', {}).get('personal_realign_entry') == 4 * GIB
                    and os.environ.get('MFM_REALIGN_PHASE_RESOURCES', 'auto') != '0'):
                compute_minimum[realign_stage] = 2
        minimum = min(total, compute_minimum.get(stage, minimum))
    if label == "WGS" and stage == "01" and os.environ.get("MFM_IO_BACKEND", "cpu") == "nvcomp":
        minimum = max(minimum, min(total, 6))
    if label == "WGS" and stage == "01" and os.environ.get("MFM_IO_BACKEND", "cpu") == "nvcomp":
        minimum = max(minimum, min(total, 6))
    if label=='WGS' and stage=='03' and profile=='throughput' and os.environ.get('MFM_WGS_PAYLOAD')=='sequence_only':
        minimum=max(minimum,min(total,6))
    if gpu_markdup:minimum=max(minimum,min(2,total))
    minimum = min(minimum, maximum)
    return minimum, maximum


def prune(root: Path) -> None:
    leases = root / "leases"
    waiters = root / "waiters"
    prune_leases(root)
    for path in waiters.glob("*.json"):
        row = read_json(path)
        if not alive(row.get("pid")):
            path.unlink(missing_ok=True)


def allocate(root: Path, total: int) -> None:
    """대기열 순서대로 현재 빈 코어를 나눈다. 기존 작업의 grant는 회수하지 않는다."""
    total = pool_cpu_limit(root, total)
    prune(root)
    leases = root / "leases"
    waiters = root / "waiters"
    rows = [(p, read_json(p)) for p in sorted(waiters.glob("*.json"))]
    active_rows=[read_json(p) for p in leases.glob("*.json")]
    policy=read_json(root/'policy.json')
    free = available_regular_cpus(total,policy.get('cpu_backfill',{}),active_rows)
    free_memory = available_memory(root,active_rows)
    priority = priority_request(root,total,active_rows,waiters=rows)
    if priority is not None and priority['kind']=='internal':
        free = max(0, free-priority['additional_cpus'])
    elif priority is not None and free < priority['minimum_cpus']:
        free = 0
    if priority is not None and priority['kind']=='stage':
        rows.sort(key=lambda item:request_key(item[1])!=priority['key'])
    eligible = [(p, r) for p, r in rows
                if time.time() - float(r.get("queued_at", time.time())) >= SETTLE_SECONDS]
    eligible = fitting_waiters(eligible, free, free_memory,
                               stage_limits=policy.get('stage_concurrency'),active=active_rows,
                               prepare_read_only=policy.get('prepare_phased') is True)
    for index, (path, row) in enumerate(eligible):
        if free <= 0:
            break
        remaining = len(eligible) - index
        fair_share = max(1, free // remaining)
        # 균등 배분 값이 최소 예산보다 작아도 앞선 작업을 시작해 교착을 피한다.
        grant = min(int(row["max_cpus"]), free, max(int(row["min_cpus"]), fair_share))
        if grant < int(row["min_cpus"]):
            continue
        requested_memory = int(row.get('memory_bytes', 0))
        if requested_memory > free_memory:
            continue
        lease = {
            "ticket": path.stem,
            "manager_pid": int(row["pid"]),
            "child_pid": None,
            "cpus": grant,
            "total_cpus": total,
            "branch": row["branch"],
            "stage": row["stage"],
            "sample": row["sample"],
            "acquired_at": time.time(),
            "memory_bytes": requested_memory,
        }
        if row.get('fixed_sample_token'):
            lease['fixed_sample_token'] = row['fixed_sample_token']
        atomic_json(leases / path.name, lease)
        path.unlink(missing_ok=True)
        free -= grant
        free_memory -= requested_memory


def adjust_argv(command: list[str], branch: str, stage: str, grant: int) -> list[str]:
    """자동 RNA thread는 grant까지 늘리고, 사용자 고정값은 상한으로 보존한다."""
    args = list(command)
    legacy = os.environ.get("MFM_CPU_POOL_LAYOUT") == "legacy"
    fixed_stages = {x for x in os.environ.get("MFM_DYNAMIC_RNA_FIXED_STAGES", "").split(",") if x}
    allow_growth = stage not in fixed_stages
    fastq_stage = "11" if legacy else "13"
    threaded_stages = {"13", "14", "18"} if legacy else {"15", "16", "20"}
    if branch.upper() == "RNA" and stage == fastq_stage:
        keys = {"--collate-threads", "--fastq-threads"}
        found = []
        for i, arg in enumerate(args[:-1]):
            if arg in keys:
                found.append(i + 1)
        if len(found) == 2:
            values = [max(1, int(args[i])) for i in found]
            target = max(2, grant if allow_growth else min(grant, sum(values)))
            first = max(1, min(target - 1, round(target * values[0] / sum(values))))
            second = target - first
            args[found[0]], args[found[1]] = str(first), str(second)
    elif branch.upper() == "RNA" and stage in threaded_stages:
        for i, arg in enumerate(args[:-1]):
            if arg == "--threads":
                requested = max(1, int(args[i + 1]))
                args[i + 1] = str(grant if allow_growth else min(requested, grant))
    return args



def reconcile_launch_memory(root, lease_path, expected, minimum):
    """과거 할당기의 lease도 자식 시작 전에 최신 CPU·요청 RAM으로 재확정함."""
    initial=read_json(lease_path)
    maximum=int(initial['cpus'])
    while True:
        with (root/'pool.lock').open('a+') as lock:
            lock_pool(lock);prune(root);row=read_json(lease_path)
            if row.get('child_pid'):raise ValueError('실행중인자식의시작예약은변경할수없음')
            total=pool_cpu_limit(root,int(row['total_cpus']))
            policy=read_json(root/'policy.json')
            revoked=(row.get('cpu_backfill') and policy.get('cpu_backfill',{}).get('enabled') is not True)
            if revoked:row['revoked_cpu_backfill']=row.pop('cpu_backfill')
            others=[read_json(p) for p in (root/'leases').glob('*.json') if p!=lease_path]
            available=(available_loan_cpus(total,policy.get('cpu_backfill',{}),others,requested=maximum)
                       if row.get('cpu_backfill') else available_regular_cpus(total,policy.get('cpu_backfill',{}),others))
            ram=available_memory(root,others)
            units=min(maximum,total,max(0,available))
            if (not revoked and int(row.get('memory_bytes',0))==expected
                    and int(row['total_cpus'])==total and int(row['cpus'])==units
                    and ram>=expected and units>=minimum):return row
            row['total_cpus']=total
            row['requested_memory_bytes']=expected
            row['initial_allocator_memory_bytes']=initial.get('memory_bytes',0)
            if ram>=expected and units>=minimum:
                row.update(cpus=units,memory_bytes=expected,memory_plan_reconciled=True)
                atomic_json(lease_path,row);return row
            # 아직 도구를 시작하지 않은 대기자끼리 CPU를 붙잡는 교착을 막음.
            row.update(cpus=0,memory_bytes=4*GIB,memory_plan_reconciled=False)
            atomic_json(lease_path,row)
        time.sleep(POLL_SECONDS)


def run(args: argparse.Namespace, command: list[str]) -> int:
    from integrate_wgs_rna.stage_gpu_config import apply_stage_environment
    gpu_pool_config = apply_stage_environment(os.environ)
    if gpu_pool_config:
        print('[GPU-POOL-CONFIG] '+json.dumps(gpu_pool_config,ensure_ascii=False), file=sys.stderr, flush=True)
    from integrate_wgs_rna.gpu_compare_resources import enabled,run_unreserved
    if enabled():return run_unreserved(args,command,adjust_argv,lambda row:append_jsonl(args.audit,row))
    root = args.pool_dir
    args.total_cpus = pool_cpu_limit(root, args.total_cpus)
    (root / "leases").mkdir(parents=True, exist_ok=True)
    (root / "waiters").mkdir(parents=True, exist_ok=True)
    ticket = f"{time.time_ns()}_{os.getpid()}_{uuid.uuid4().hex[:8]}"
    waiter_path = root / "waiters" / f"{ticket}.json"
    lease_path = root / "leases" / f"{ticket}.json"
    queued = time.time()
    # 실행 중인 shell도 다음 단계 경계에서 공유 풀의 새 정책을 읽는다.
    with (root / 'pool.lock').open('a+') as lock:
        lock_pool(lock)
        policy = read_json(root / 'policy.json')
        policy['memory_limit_bytes'] = memory_limit(root)
        if os.environ.get('MFM_CPU_UTILIZATION','maximum')=='maximum':
            policy.setdefault('cpu_backfill',{'enabled':True})
        atomic_json(root / 'policy.json', policy)
    if os.environ.get("MFM_CPU_UTILIZATION", "maximum") == "maximum":
        os.environ["MFM_CPU_POOL_PROFILE"] = "throughput"
    elif policy.get("profile") in {"balanced", "throughput", "conservative"}:
        os.environ["MFM_CPU_POOL_PROFILE"] = policy["profile"]
    min_cpus, max_cpus = stage_cap(args.branch, args.stage, args.total_cpus)
    lock_path = root / "pool.lock"
    if stage_memory(args.branch, args.stage) > usable_memory(root):
        raise ValueError("단계의 RAM 예상 예약량이 전체 예산보다 큽니다. 입력/메모리 설정을 확인하세요.")
    atomic_json(waiter_path, {"pid": os.getpid(), "queued_at": queued, "ticket": ticket,
                              "branch": args.branch, "stage": args.stage, "sample": args.sample,
                              "min_cpus": min_cpus, "max_cpus": max_cpus,
                              "memory_bytes": stage_memory(args.branch, args.stage),
                              "cpu_utilization":os.environ.get('MFM_CPU_UTILIZATION','maximum'),
                              "fixed_sample_token":os.environ.get("MFM_FIXED_SAMPLE_TOKEN",""),
                              "layout":os.environ.get('MFM_CPU_POOL_LAYOUT','original')})
    try:
        while not lease_path.exists():
            if not native_allocator_ready(root):
                with lock_path.open("a+") as lock:
                    lock_pool(lock)
                    allocate(root, args.total_cpus)
                    fcntl.flock(lock, fcntl.LOCK_UN)
            if not lease_path.exists():
                time.sleep(POLL_SECONDS)
        lease = reconcile_launch_memory(root,lease_path,stage_memory(args.branch,args.stage),min_cpus)
        args.total_cpus = int(lease['total_cpus'])
        max_cpus = min(max_cpus, args.total_cpus)
        grant = int(lease["cpus"])
        effective = adjust_argv(command, args.branch, args.stage, grant)
        env = os.environ.copy()
        env.update({"MFM_CPU_BUDGET": str(args.total_cpus),
                    "MFM_SAMPLE_CPU_BUDGET": str(grant),
                    "MFM_CPU_POOL_GRANTED": str(grant),
                    "MFM_CPU_POOL_STAGE": args.stage,
                    "MFM_CPU_POOL_LEASE_FILE": str(lease_path),
                    "MFM_DYNAMIC_CPU_POOL_ACTIVE": "1",
                    "MFM_RAM_LIMIT_BYTES": str(memory_limit(root)),
                    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"})
        from integrate_wgs_rna.stage_cpu_config import apply as apply_stage_cpus
        apply_stage_cpus(env,args.branch,args.stage,grant)
        acquired = float(lease["acquired_at"])
        print(f"[CPU-POOL] branch={args.branch} stage={args.stage} sample={args.sample} "
              f"grant={grant}/{args.total_cpus} queue_wait={acquired - queued:.2f}s "
              f"command={' '.join(effective)}", file=sys.stderr, flush=True)
        usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
        tree_peak_rss = 0
        sampled_lease_seconds = 0.0
        last_sample_at = time.monotonic()
        last_grant = grant
        try:
            proc = subprocess.Popen(effective, env=env, start_new_session=True)
        except PermissionError:
            # 공유 파일시스템에 noexec가 걸려 있어도 Bash stage는 shebang과
            # 동일한 interpreter로 실행할 수 있다. Python/바이너리는 우회하지 않는다.
            if effective[0].endswith(".sh"):
                effective = ["bash", *effective]
                proc = subprocess.Popen(effective, env=env, start_new_session=True)
            else:
                raise
        lease["child_pid"] = proc.pid
        atomic_json(lease_path, lease)

        def forward(signum: int, _frame: object) -> None:
            try:
                os.killpg(proc.pid, signum)
            except ProcessLookupError:
                pass

        old_int = signal.signal(signal.SIGINT, forward)
        old_term = signal.signal(signal.SIGTERM, forward)
        started = time.time()
        try:
            while True:
                try:
                    rc = proc.wait(timeout=0.5)
                    break
                except subprocess.TimeoutExpired:
                    guard = read_json(root/'memory_guard.json')
                    if not (0 <= time.time()-float(guard.get('time',0)) < 2
                            and guard.get('limit_bytes') == memory_limit(root)):
                        with lock_path.open('a+') as lock:
                            lock_pool(lock)
                            guard = protect_pool_memory(root)
                    tree_peak_rss = max(tree_peak_rss, int(guard.get('usage_by_lease', {}).get(lease_path.name, 0)))
                    sample_at = time.monotonic()
                    sampled_lease_seconds += (sample_at - last_sample_at) * last_grant
                    last_sample_at = sample_at
                    last_grant = int(read_json(lease_path).get('cpus', grant))
                    if guard.get('action') == 'stop_largest_stage':
                        print('[RAM-GUARD] ' + json.dumps(guard, ensure_ascii=False), file=sys.stderr, flush=True)
        finally:
            signal.signal(signal.SIGINT, old_int)
            signal.signal(signal.SIGTERM, old_term)
        sampled_lease_seconds += (time.monotonic() - last_sample_at) * last_grant
        usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_seconds = max(0.0, usage_after.ru_utime + usage_after.ru_stime - usage_before.ru_utime - usage_before.ru_stime)
        finished = time.time()
        final_lease = read_json(lease_path)
        append_jsonl(args.audit, {"schema": "mfm_cpu_lease_v1", "sample": args.sample,
            "branch": args.branch, "stage": args.stage, "total_cpus": args.total_cpus,
            "granted_cpus": grant, "minimum_cpus": min_cpus, "maximum_cpus": max_cpus,
            "local_jobs": int(os.environ.get("MFM_LOCAL_JOBS", "1")),
            "resource_profile": os.environ.get("MFM_RESOURCE_PROFILE", ""),
            "elastic_sample_slots": os.environ.get("MFM_ELASTIC_SAMPLE_SLOTS", "false") == "true",
            "queued_at": queued, "acquired_at": acquired, "started_at": started,
            "finished_at": finished, "queue_wait_sec": round(acquired - queued, 3),
            "elapsed_sec": round(finished - started, 3), "exit_code": rc,
            "actual_cpu_seconds": round(cpu_seconds, 6),
            "average_used_cores": round(cpu_seconds / max(0.000001, finished - started), 4),
            "sampled_leased_cpu_seconds": round(sampled_lease_seconds, 6),
            "cpu_fraction_of_lease": round(cpu_seconds / max(0.000001, sampled_lease_seconds), 4),
            "peak_tree_rss_observed_bytes": tree_peak_rss,
            "peak_single_process_rss_bytes": int(usage_after.ru_maxrss) * 1024,
            "resource_measurement_scope": "wait한 자식·후손 CPU 누계; tree RSS는 RAM guard 시각 기반 표본이며 잠금 대기로 간격이 늘어날 수 있음; 실제 순간 피크와 구분함",

            "peak_granted_cpus": max(grant, int(final_lease.get("peak_cpus", grant))),
            "memory_limit_bytes": memory_limit(root),
            "cpu_allocation_mode": os.environ.get("MFM_CPU_ALLOCATION_MODE","dynamic"),
            "sample_cpu_budget": os.environ.get("MFM_FIXED_SAMPLE_CPUS"),
            "affinity_cpu_ids": sorted(os.sched_getaffinity(0)),
            "final_reserved_memory_bytes": int(final_lease.get('memory_bytes', 0)),
            "cpu_backfill": final_lease.get('cpu_backfill'),
            "command": command, "effective_command": effective,
            **({"gpu_pool_config": gpu_pool_config} if gpu_pool_config else {})})
        return rc
    finally:
        try:
            with lock_path.open("a+") as lock:
                lock_pool(lock)
                waiter_path.unlink(missing_ok=True)
                lease_path.unlink(missing_ok=True)
                allocate(root, args.total_cpus)
                fcntl.flock(lock, fcntl.LOCK_UN)
        except OSError:
            pass


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool-dir", type=Path, required=True)
    parser.add_argument("--total-cpus", type=int, required=True)
    parser.add_argument("--branch", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--audit", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.total_cpus < 1:
        parser.error("--total-cpus는 양의 정수여야 합니다")
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("-- 뒤에 실행할 command가 필요합니다")
    return run(args, command)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
