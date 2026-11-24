"""ASE 입력 순서를 보존하며 반복 reference 이름 생성을 줄임."""
import gc,json,os,time
from pathlib import Path
import pysam
from .cpu_executor import WorkerBudget
from .memory_policy import GIB,memory_limit

def _read_as_impl(bam_path: str, mapq_min: int, include_duplicates: bool = False, *, cache_reference_names=True):
    """primary 정렬을 MAPQ 로 갈라 두 벌로 읽는다.

    반환 (uniq, mmap). 각각 query_name → [mate0, mate1], 원소는 (AS, contig, blocks) 또는 None.
      uniq: MAPQ >= mapq_min        — 기존 카운트 경로
      mmap: 0 < MAPQ <  mapq_min    — multimapper 보정용(AlleleSeq2 의 mmapreads BAM 대응)
    두 집합은 서로소다. MAPQ 0(정렬 위치 불확정)은 양쪽 모두에서 제외한다.

    duplicate 는 기본적으로 제외한다 — AlleleSeq2 는 카운트 전에 rmdup 으로 제거한다.
    어느 fragment 가 duplicate 인지는 호출부가 두 HAP 의 합집합으로 다시 정한다
    (drop_duplicate_fragments). 여기서는 이 BAM 에서 DUP flag 가 붙은 이름만 모아 돌려준다.
    """
    uniq: dict = {}
    mmap: dict = {}
    n_dup_skipped = 0
    dup_names: set = set()
    with pysam.AlignmentFile(bam_path, "rb") as bf:
        references=bf.references
        for r in bf.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if r.is_duplicate:
                n_dup_skipped += 1
                dup_names.add(r.query_name)
                if not include_duplicates:
                    continue
            mq = r.mapping_quality
            if mq >= mapq_min:
                target = uniq
            elif mq > 0:
                target = mmap
            else:
                continue
            try:
                a = int(r.get_tag("AS"))
            except KeyError:
                continue
            i = 0 if r.is_read1 else 1
            name=r.query_name
            slot = target.get(name)
            if slot is None:
                slot = target[name] = [None, None]
            prev = slot[i]
            if prev is None or a > prev[0]:
                slot[i] = (a, references[r.reference_id] if cache_reference_names and 0 <= r.reference_id < len(references) else r.reference_name, r.get_blocks())
    return uniq, mmap, n_dup_skipped, dup_names






def read_as(bam_path,mapq_min,include_duplicates=False,*,cache_reference_names=True):
    # 생성하는 dict/list/tuple/scalar는 순환 참조가 없는 입력 table임.
    # 구축 중 반복 cyclic GC만 일시 중지하고 호출 전 상태로 복원함.
    enabled=gc.isenabled();gc.disable()
    try:return _read_as_impl(bam_path,mapq_min,include_duplicates,cache_reference_names=cache_reference_names)
    finally:
        if enabled:gc.enable()


def read_pair_as(bam1,bam2,mapq,duplicates=False,*,workers=1,benchmark=False,scratch=None):
    started=time.monotonic();paths=[str(bam1),str(bam2)];raw=sum(Path(p).stat().st_size for p in paths)
    estimate=max(2*GIB,raw*32)*(2 if benchmark else 1)
    if estimate+4*GIB>memory_limit()*9//10:raise ValueError('ASE 입력 예상 RAM이 예약 예산보다 큼')
    # 대형 table을 IPC/pickle로 복제하지 않음. 이 직렬 단계는 한 CPU만 유지함.
    budget=WorkerBudget(1,1,False,memory_per_worker=estimate)
    while budget.capacity(1)<1:time.sleep(.2)
    values=[read_as(path,mapq,duplicates) for path in paths]
    report={'mode':'serial_optimized','workers':1,'input_bytes':raw,'estimated_table_bytes':estimate,
            'seconds':time.monotonic()-started,'scratch_bytes':0}
    if benchmark:
        comparison=time.monotonic();legacy=[_read_as_impl(path,mapq,duplicates,cache_reference_names=False) for path in paths]
        legacy_seconds=time.monotonic()-comparison
        if values!=legacy:raise ValueError('ASE 입력 최적화 전체 record·duplicate 불일치')
        report.update(serial_equal=True,legacy_seconds=legacy_seconds)
        report['input_read_speedup']=report['legacy_seconds']/max(.000001,report['seconds'])
    print('[ASE-INPUT] '+json.dumps(report,ensure_ascii=False),flush=True)
    return values,report
