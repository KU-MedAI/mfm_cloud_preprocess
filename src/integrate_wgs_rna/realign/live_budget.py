"""실행 중 BWA/HiFi minimap2의 PID·인자·reference·RSS를 검증한 예약 감소만 허용함."""
from pathlib import Path
from types import SimpleNamespace
import fcntl
import hashlib
import json
import os
import shlex
import time

from .align_memory import alignment_memory
from .phases import _atomic
from ..memory_policy import GIB


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def tree_snapshot(pid):
    processes = {}
    for path in Path('/proc').iterdir():
        if not path.name.isdigit():
            continue
        try:
            raw = (path / 'stat').read_text()
            fields = raw[raw.rfind(')') + 2:].split()
            if fields[0] == 'Z':
                continue
            processes[int(path.name)] = {'parent': int(fields[1]), 'start': int(fields[19]),
                'rss': int(fields[21]) * os.sysconf('SC_PAGE_SIZE'),
                'argv': (path / 'cmdline').read_bytes().decode(errors='replace').split('\0')}
        except (OSError, ValueError, IndexError):
            continue
    if pid not in processes:
        raise ValueError('단계 manager PID가 생존하지 않음')
    found = {pid}
    while True:
        extra = {key for key, value in processes.items() if value['parent'] in found} - found
        if not extra:
            break
        found |= extra
    return {'manager_start': processes[pid]['start'],
            'rss_bytes': sum(processes[key]['rss'] for key in found),
            'bwa': [{'pid': key, **processes[key]} for key in found
                    if 'mem' in processes[key]['argv'] and
                    'bwa-mem2' in Path(processes[key]['argv'][0]).name],
            'minimap2': [{'pid': key, **processes[key]} for key in found
                         if Path(processes[key]['argv'][0]).name == 'minimap2'],
            'sort': [{'pid': key, **processes[key]} for key in found
                     if Path(processes[key]['argv'][0]).name == 'samtools'
                     and len(processes[key]['argv']) > 1 and processes[key]['argv'][1] == 'sort']}


def propose(lease_path, cohort_root):
    lease_path = Path(lease_path)
    row = json.loads(lease_path.read_text())
    if row.get('realign_phase') != 'align_restore_sort_qc' or row.get('branch') != 'WGS':
        raise ValueError('현재 실제 정렬 중인 lease만 보정 가능함')
    sample = row['sample']
    if Path(sample).name != sample or sample in {'.', '..'}:
        raise ValueError('잘못된 sample ID')
    effective = row['effective_alignment']
    haps = int(effective['hap_jobs']);threads = int(effective['threads']);sort_threads = int(effective['sort_threads'])
    tree = tree_snapshot(int(row['manager_pid']))
    hifi = bool(tree.get('minimap2'))
    decoder = int(effective.get('input_decode_threads', 0)) if hifi else 0
    overhead = 5 if hifi else 3
    if haps not in {1, 2} or int(row['cpus']) != haps * (threads + sort_threads + overhead + decoder):
        raise ValueError('실제 예약과 HAP/thread 계획 불일치')
    aligners = tree['minimap2'] if hifi else tree['bwa']
    if not aligners or len(aligners) > haps or (hifi and (len(aligners) != haps or tree['bwa'])):
        raise ValueError('현재 HAP 계획과 생존 정렬 도구 수가 일치하지 않음')
    for process in aligners:
        argv = process['argv']
        if '-t' not in argv or int(argv[argv.index('-t')+1]) != threads:
            raise ValueError('실제 정렬 thread와 lease 계획 불일치')
        if hifi and ('-ax' not in argv or argv[argv.index('-ax')+1] != 'map-hifi'):
            raise ValueError('검증된 HiFi minimap2 명령만 보정 가능함')
    root = Path(cohort_root) / sample
    lengths = {}; references = {}
    for hap in (1, 2):
        total = 0
        for number in range(1, 23):
            path = root / 'individual_fasta' / f'HAP{hap}' / f'chr{number}.fa.fai'
            text = path.read_text();parts = [line.split('\t') for line in text.splitlines()]
            if len(parts) != 1 or parts[0][0] != f'chr{number}':
                raise ValueError('상염색체 reference FAI 계약 불일치')
            total += int(parts[0][1]);references[str(path)] = hashlib.sha256(text.encode()).hexdigest()
        lengths[hap] = total
    sort_values = set()
    scripts = list((root / 'logs/personal_realign').glob(sample + '.HAP*.realign.sh'))
    for script in scripts:
        for line in script.read_text().splitlines():
            if ' sort -@ ' in line:
                args = shlex.split(line)
                if int(args[args.index('-@')+1]) != sort_threads:
                    # 완료 HAP의 예전 계획은 현재 계획과 다를 수 있으므로 현재 BWA HAP만 검사함.
                    continue
                sort_values.add(args[args.index('-m')+1])
    if len(sort_values) != 1:
        raise ValueError('현재 sort 버퍼를 단일 값으로 검증할 수 없음')
    sort_mem = sort_values.pop()
    if hifi:
        sorters = tree.get('sort', [])
        if len(sorters) != haps:
            raise ValueError('현재 HAP별 실제 samtools sort를 모두 검증할 수 없음')
        for process in sorters:
            argv = process['argv']
            if ('-@' not in argv or '-m' not in argv
                    or int(argv[argv.index('-@')+1]) != sort_threads
                    or argv[argv.index('-m')+1] != sort_mem):
                raise ValueError('실제 sort 인자와 현재 계획 불일치')
        states = row.get('hap_resources', {})
        if len(states) != haps or any(value.get('phase') != 'align_restore_sort_qc' for value in states.values()):
            raise ValueError('동시 HAP 하위 단계가 변경되어 재검토가 필요함')
    config = SimpleNamespace(sort_mem=sort_mem, platform='hifi' if hifi else 'illumina',
                             input_decode_threads=decoder)
    model = alignment_memory(config, lengths, haps, threads, sort_threads,
                             {'MFM_REALIGN_MEMORY_MODEL': 'reference_v2'})
    observed_guard = (tree['rss_bytes']*2 + GIB-1)//GIB*GIB
    new = max(model['reserved_memory_bytes'], observed_guard)
    if new >= int(row['memory_bytes']):
        raise ValueError('감소시킬 보수 예약 여유가 없음')
    return {'schema': 'mfm_live_alignment_rebalance_v1', 'lease_path': str(lease_path),
            'cohort_root': str(cohort_root), 'sample': sample,
            'manager_pid': row['manager_pid'], 'manager_start': tree['manager_start'],
            'lease_fingerprint': fingerprint(row), 'references': references,
            'bwa_processes': {str(p['pid']): p['start'] for p in tree['bwa']},
            'minimap2_processes': {str(p['pid']): p['start'] for p in tree.get('minimap2', [])},
            'sort_processes': {str(p['pid']): p['start'] for p in tree.get('sort', [])} if hifi else {},
            'rss_bytes': tree['rss_bytes'], 'rss_guard_factor': 2, 'model': model,
            'old_memory_bytes': row['memory_bytes'], 'new_memory_bytes': new,
            'cpu_unchanged': row['cpus'], 'status': 'planned'}


def apply(proposal):
    path = Path(proposal['lease_path'])
    with (path.parent.parent / 'pool.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        row = json.loads(path.read_text())
        if fingerprint(row) != proposal['lease_fingerprint']:
            raise ValueError('계획 이후 lease 변경. 재검토 필요함')
        current = propose(path, proposal['cohort_root'])
        for key in ['manager_pid', 'manager_start', 'references', 'bwa_processes', 'minimap2_processes', 'sort_processes']:
            if current[key] != proposal[key]:
                raise ValueError('PID/reference/실행 중 정렬 도구가 변경됨')
        if current['new_memory_bytes'] > proposal['new_memory_bytes']:
            raise ValueError('계획 이후 RSS 증가로 보정 여유 부족함')
        row['memory_bytes'] = proposal['new_memory_bytes']
        states = row.get('hap_resources', {})
        if states:
            per_hap, remainder = divmod(row['memory_bytes'], len(states))
            row['hap_resources'] = {hap: dict(value, memory_bytes=per_hap + (index < remainder))
                                    for index, (hap, value) in enumerate(states.items())}
        row['effective_alignment'] = dict(row['effective_alignment'],
            reserved_memory_bytes=row['memory_bytes'], memory_model=proposal['model'])
        row['memory_rebalance'] = {'time': time.time(), 'old_bytes': proposal['old_memory_bytes'],
            'new_bytes': row['memory_bytes'], 'manager_start': current['manager_start'],
            'rss_bytes': current['rss_bytes'], 'model': 'reference_v2'}
        _atomic(path, row)
    return dict(proposal, status='applied', applied_at=time.time())
