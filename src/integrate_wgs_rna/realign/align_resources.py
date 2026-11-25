"""정렬을 시작하기 직전 가용 CPU/RAM으로 HAP 동시 수와 thread를 결정함."""
import os
from ..memory_policy import GIB
from .align_memory import alignment_memory


def parse_memory(value):
    units = {'K': 1024, 'M': 1024 ** 2, 'G': GIB, 'T': 1024 ** 4}
    text = str(value).upper()
    if text[-1:] in units:
        return int(float(text[:-1]) * units[text[-1]])
    return int(text)


def alignment_cpu_overhead(config):
    # minimap2의추가I/O thread + view/fastq + sort main + unmapped awk를포함함.
    platform=str(getattr(config,'platform','')).lower()
    return 5 if platform in {'hifi','ont','pacbio_hifi','nanopore'} else 3


def _fixed_cpu_ids(env, cpus):
    """살아 있는 샘플 소유권과 실제 affinity가 같은 경우만 겹침을 허용함."""
    if (env.get('MFM_CPU_ALLOCATION_MODE') != 'fixed_sample'
            or not env.get('MFM_FIXED_SAMPLE_TOKEN')):
        return None
    from ..sample_cpu import owner, owner_alive
    try:
        _, row = owner(env['MFM_CPU_POOL_DIR'], env['MFM_FIXED_SAMPLE_TOKEN'])
        ids = set(row['cpu_ids'])
        if (int(row['cpus']) == cpus == len(ids) and owner_alive(row)
                and os.sched_getaffinity(0) == ids):
            return sorted(ids)
    except (KeyError, OSError, ValueError, TypeError, AttributeError):
        pass
    return None


def _alignment_candidate(config, lengths, available_cpus, available_memory, environment=None, *, sort_cap=None):
    env = os.environ if environment is None else environment
    cpu = max(0, int(available_cpus))
    if env.get('MFM_WGS_BACKEND','cpu')=='parabricks':
        cpu=min(cpu,max(1,int(env.get('MFM_GPU_HAP_JOBS','2')))*max(4,int(env.get('MFM_GPU_HOST_CPUS_PER_HAP',str(cpu)))))
    if not lengths or cpu < 4:
        return None
    overhead=alignment_cpu_overhead(config)
    automatic = env.get('WGS_AUTO_THREADS', '1') == '1'
    align_request = env.get('MFM_WGS_ALIGN_REQUESTED')
    sort_request = env.get('MFM_WGS_SORT_REQUESTED')
    hap_request = env.get('MFM_WGS_HAP_REQUESTED')
    auto_align = automatic and align_request == 'auto'
    auto_sort = automatic and sort_request == 'auto'
    auto_haps = hap_request == 'auto'
    balance=env.get('MFM_WGS_CPU_BALANCE','auto')
    if balance not in {'auto','legacy','mapping','mapping_more'}:
        raise ValueError('MFM_WGS_CPU_BALANCE는 auto|legacy|mapping|mapping_more이어야 함')
    platform=str(getattr(config,'platform','')).lower()
    profile='legacy'
    if auto_align and auto_sort and overhead==5 and env.get('MFM_WGS_BACKEND','cpu')=='cpu':
        if balance in {'mapping','mapping_more'}:profile=balance
        elif balance=='auto' and platform in {'hifi','pacbio_hifi'}:profile='mapping_more'
    mapping=profile!='legacy'
    default_minimum = '4' if env.get('MFM_CPU_UTILIZATION','maximum')=='maximum' else '1'
    minimum_threads = int(env.get('MFM_WGS_MIN_ALIGN_THREADS', default_minimum)) if auto_align else 1
    if minimum_threads < 1:
        raise ValueError('MFM_WGS_MIN_ALIGN_THREADS는 양수여야 함')
    fixed_threads = int(align_request) if align_request not in {None, 'auto'} else config.threads
    fixed_sort = int(sort_request) if sort_request not in {None, 'auto'} else None
    max_haps = min(2, len(lengths)) if auto_haps else min(int(hap_request) if hap_request in {'1','2'} else config.hap_jobs, len(lengths))
    if auto_haps and env.get('MFM_HAP_PARALLEL_ALLOWED', '1') == '0':
        max_haps = 1
    if env.get('MFM_WGS_BACKEND','cpu')=='parabricks':max_haps=min(max_haps,max(1,int(env.get('MFM_GPU_HAP_JOBS','2'))))
    fixed_ids = (_fixed_cpu_ids(env, cpu)
                 if cpu == 16 and auto_align and auto_sort and auto_haps
                 and profile == 'mapping_more' and balance == 'auto'
                 and env.get('MFM_MINIMAP_DECODE_THREADS', 'auto') == 'auto' else None)
    candidates = range(max_haps, 0, -1) if auto_haps else [max_haps]
    for haps in candidates:
        per_hap = cpu // haps
        decoder=0
        if overhead==5:
            request=env.get('MFM_MINIMAP_DECODE_THREADS','auto')
            if request!='auto' and (not request.isdigit()):raise ValueError('MFM_MINIMAP_DECODE_THREADS는auto또는0이상정수여야함')
            decoder=min(1 if profile=='mapping_more' else 2 if mapping else 4,max(0,(per_hap-overhead)//6)) if request=='auto' else int(request)

        if fixed_sort is not None:
            auxiliary = fixed_sort
        elif auto_sort:
            auxiliary = (max(0,(per_hap-overhead-decoder)//(12 if profile=='mapping_more' else 6 if mapping else 3)) if overhead==5 else min(4,max(0,(per_hap-overhead)//4)))
        else:
            auxiliary = min(config.sort_threads, max(0, per_hap - fixed_threads - overhead - decoder))
        if auto_sort and sort_cap is not None:
            auxiliary = min(auxiliary, sort_cap)
        if auto_align:
            if auto_sort and per_hap-auxiliary-overhead-decoder<minimum_threads:
                auxiliary=max(0,per_hap-overhead-decoder-minimum_threads)
            threads = per_hap - auxiliary - overhead - decoder
        else:
            threads = fixed_threads
            if fixed_sort is None:
                auxiliary = min(auxiliary, max(0, per_hap - threads - overhead - decoder))
        overlap = fixed_ids is not None and haps == 2
        if overlap:
            # 두 파이프라인의 유휴 구간을 같은 소유 CPU 집합에서 서로 채움.
            threads, auxiliary, decoder = 8, min(1, sort_cap) if sort_cap is not None else 1, 1
        if threads < minimum_threads or auxiliary < 0:
            continue
        logical = haps * (threads + auxiliary + overhead + decoder)
        cpus = cpu if overlap else logical
        if cpus > cpu:
            continue
        memory_model = alignment_memory(config, lengths, haps, threads, auxiliary, dict(env,MFM_MINIMAP_DECODE_THREADS_EFFECTIVE=str(decoder)))
        memory = memory_model['reserved_memory_bytes']
        if env.get('MFM_WGS_BACKEND','cpu')=='parabricks' and env.get('MFM_GPU_PHASE_RESOURCES','1')=='1':
            # 도구 실행 시에는 StageBudget이 host 피크와 GPU를 함께 재확보한다.
            memory_model=dict(memory_model,gpu_compute_peak_bytes=memory,
                              reservation_phase='input_prepare',reserved_memory_bytes=4*GIB*haps)
            memory=memory_model['reserved_memory_bytes']
        if memory > available_memory:
            continue
        return {'hap_jobs': haps, 'threads': threads, 'sort_threads': auxiliary,
                'fastq_threads': 0, 'input_decode_threads':decoder, 'cpu_balance': profile,
                'reserved_cpus': cpus, 'reserved_memory_bytes': memory, 'memory_model': memory_model,
                **({'cpu_execution_policy': 'fixed_affinity_overlap',
                    'logical_pipeline_cpus': logical, 'cpu_affinity_ids': fixed_ids} if overlap else {})}
    return None


def alignment_plan(config, lengths, available_cpus, available_memory, environment=None):
    """CPU를 유지하며 자동 sort 비중을 조정하고, 불가능할 때만 CPU 예산을 낮춤."""
    env = os.environ if environment is None else environment
    available = max(0, int(available_cpus))
    automatic = env.get('WGS_AUTO_THREADS', '1') == '1' and env.get('MFM_WGS_ALIGN_REQUESTED') == 'auto'
    adjustable = automatic and env.get('MFM_WGS_BACKEND', 'cpu') == 'cpu'
    budgets = range(available, 3, -1) if adjustable else [available]
    plan = None
    rebalanced = False
    selected = available
    for selected in budgets:
        plan = _alignment_candidate(config, lengths, selected, available_memory, env)
        if plan is None and adjustable and env.get('MFM_WGS_SORT_REQUESTED') == 'auto':
            for cap in range(selected // 3, -1, -1):
                plan = _alignment_candidate(config, lengths, selected, available_memory, env, sort_cap=cap)
                if plan is not None:
                    rebalanced = True
                    break
        if plan is not None:
            break
    if plan is None:
        return None
    reason = 'RAM에 맞춘 sort/정렬 CPU 재배분' if rebalanced else '기존 선호 구성'
    if plan.get('cpu_execution_policy') == 'fixed_affinity_overlap':
        reason = '고정 16CPU 집합에서 HAP 파이프라인 겹침'
    if selected < available:
        reason = 'RAM에 맞춘 자동 CPU 축소; ' + reason
    return dict(plan, allocation_diagnostics={
        'available_cpus': available, 'selected_cpu_budget': selected,
        'available_memory_bytes': int(available_memory),
        'memory_limited_cpu': selected < available,
        'memory_rebalanced_sort': rebalanced, 'reason': reason,
    })
