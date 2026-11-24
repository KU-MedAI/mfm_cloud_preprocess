"""출력 BAM/CRAM을 contig별로 나누어 정확한 QC와 실행 자원을 기록함."""
from __future__ import annotations
from contextlib import ExitStack,contextmanager
from concurrent.futures import ProcessPoolExecutor
from itertools import chain
from pathlib import Path
import gzip,json,multiprocessing,os,resource,tempfile,time
import pysam
from ..memory_policy import memory_limit,limit_worker_address_space

MIB=1024**2
GIB=1024**3


def output_scan_plan(config,output,references,indexed,requested,resource_limits=None):
    affinity=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else os.cpu_count() or 1
    caps=[affinity];lease={}
    if os.environ.get('MFM_CPU_POOL_LEASE_FILE'):
        lease=json.loads(Path(os.environ['MFM_CPU_POOL_LEASE_FILE']).read_text())
    for key in ['MFM_CPU_BUDGET','SLURM_CPUS_PER_TASK']:
        if os.environ.get(key):caps.append(int(os.environ[key]))
    if lease.get('cpus'):caps.append(int(lease['cpus']))
    else:
        for key in ['MFM_CPU_POOL_GRANTED','MFM_SAMPLE_CPU_BUDGET']:
            if os.environ.get(key):caps.append(int(os.environ[key]))
    active_haps=lease.get('effective_alignment',{}).get('hap_jobs',getattr(config,'hap_jobs',1))
    if resource_limits:
        caps.append(int(resource_limits['cpus']));active_haps=1
    concurrency=max(1,int(getattr(config,'jobs',1)))*max(1,int(active_haps))
    cpu=max(1,min(caps)//concurrency)
    ram=memory_limit()
    for key,factor in [('SLURM_MEM_PER_NODE',1),('SLURM_MEM_PER_CPU',min(caps))]:
        if int(os.environ.get(key,'0'))>0:ram=min(ram,int(os.environ[key])*MIB*factor)
    ram=ram*9//10
    if lease.get('memory_bytes'):ram=min(ram,int(lease['memory_bytes']))
    if resource_limits:ram=min(ram,int(resource_limits['memory_bytes']))
    ram//=concurrency
    longest=max((x[1] for x in references),default=0)
    per_worker=max(512*MIB,2*longest+256*MIB)
    parent=128*MIB
    ram_workers=max(0,(ram*10//11-parent)//per_worker)
    if ram_workers<1:raise ValueError('출력 QC의 예측 RAM이 예약 예산보다 큼')
    input_bytes=output.stat().st_size
    selected=min(cpu,len(references)+1,ram_workers)
    if requested is not None:
        if int(requested)<1:raise ValueError('QC workers는 양수여야 함')
        selected=min(selected,int(requested))
    elif input_bytes<32*MIB:selected=1
    if not indexed:selected=1
    count=max(1,selected)
    decode=max(0,min(4,int(config.sort_threads),cpu-1)) if count==1 else 0
    return {'mode':'indexed_contigs' if count>1 else 'serial_full_scan','input_bytes':input_bytes,
            'reference_contigs':len(references),'largest_reference_bases':longest,'indexed':indexed,
            'cpu_budget':cpu,'concurrent_hap_samples':concurrency,'workers_requested':requested,
            'workers_used':count,'decode_threads':decode,'ram_budget_bytes':ram,
            'estimated_peak_rss_bytes':(parent+count*per_worker)*11//10,
            'worker_address_space_bytes':min(ram,max(4*GIB,per_worker+2*GIB))}


def _measured_partition(output,ref,regions,decode):
    start=time.process_time()
    counts,flags=_scan_output_partition(output,ref,regions,decode)
    return counts,flags,{'cpu_seconds':time.process_time()-start,
                         'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024}


def _empty_output_scan() -> tuple[dict[str, int], dict[str, dict[str, int]]]:
    counts = {
        "primary": 0,
        "secondary": 0,
        "supplementary": 0,
        "paired_primary": 0,
        "proper_pair_primary": 0,
        "duplicate_primary": 0,
        "od_tag_primary": 0,
        "routed_primary": 0,
        "routed_mapped": 0,
        "routed_hp1": 0,
        "routed_hp2": 0,
        "routed_untagged": 0,
        "rescue_primary": 0,
        "rescue_mapped": 0,
        "nm_sum": 0,
        "mapped_bp": 0,
    }
    flagstat = {
        group: {
            key: 0 for key in (
                "total", "primary", "secondary", "supplementary", "duplicates",
                "primary_duplicates", "mapped", "primary_mapped", "paired",
                "read1", "read2", "properly_paired", "self_mate_mapped",
                "singletons", "different_chr", "different_chr_mapq5",
            )
        }
        for group in ("passed", "failed")
    }
    return counts, flagstat


@contextmanager
def _output_scan_reader(output: str, ref: str | None, regions: list[str] | None, decode_threads: int):
    """CRAM 미배치 컨테이너 전체를 읽도록 임시 CRAI만 정규화한다.

    일부 HTSlib의 '*' 조회는 동일 좌표의 마지막 미배치 CRAI 행에서 시작한다.
    미배치 행을 첫 컨테이너 한 행으로 묶으면 EOF까지 모든 미배치 read를 읽는다.
    원 CRAM과 CRAI는 수정하지 않으며, BAM과 순차 전체 scan은 그대로 둔다.
    """
    with ExitStack() as stack:
        kwargs = {"reference_filename": ref, "threads": decode_threads}
        if regions and "*" in regions and output.lower().endswith(".cram"):
            source = next((p for p in (Path(output + ".crai"), Path(output).with_suffix(".crai")) if p.exists()), None)
            if source is not None:
                with source.open("rb") as handle:
                    compressed = handle.read(2) == b"\x1f\x8b"
                with (gzip.open(source, "rt") if compressed else source.open()) as handle:
                    rows = handle.readlines()
                unmapped = [row for row in rows if row.startswith("-1\t")]
                if len(unmapped) > 1:
                    first = min(unmapped, key=lambda row: int(row.split("\t")[3]))
                    temp = Path(stack.enter_context(tempfile.TemporaryDirectory(prefix="mfm_cram_tail_"))) / "scan.crai"
                    with gzip.open(temp, "wt") as handle:
                        handle.writelines(row for row in rows if not row.startswith("-1\t"))
                        handle.write(first)
                    kwargs["index_filename"] = str(temp)
        yield stack.enter_context(pysam.AlignmentFile(output, "rb", **kwargs))


def _scan_output_partition(
    output: str, ref: str | None, regions: list[str] | None, decode_threads: int,
) -> tuple[dict[str, int], dict[str, dict[str, int]]]:
    counts, flagstat = _empty_output_scan()
    with _output_scan_reader(output, ref, regions, decode_threads) as bam:
        reads = (
            bam.fetch(until_eof=True) if regions is None
            else chain.from_iterable(bam.fetch(region) for region in regions)
        )
        for read in reads:
            qc_group = flagstat["failed" if read.is_qcfail else "passed"]
            qc_group["total"] += 1
            is_primary = not (read.is_secondary or read.is_supplementary)
            if is_primary:
                qc_group["primary"] += 1
            if read.is_secondary:
                qc_group["secondary"] += 1
            elif read.is_supplementary:
                qc_group["supplementary"] += 1
            if read.is_duplicate:
                qc_group["duplicates"] += 1
                if is_primary:
                    qc_group["primary_duplicates"] += 1
            if not read.is_unmapped:
                qc_group["mapped"] += 1
                if is_primary:
                    qc_group["primary_mapped"] += 1
            if is_primary and read.is_paired:
                qc_group["paired"] += 1
                qc_group["read1"] += int(read.is_read1)
                qc_group["read2"] += int(read.is_read2)
                qc_group["properly_paired"] += int(read.is_proper_pair and not read.is_unmapped)
                if not read.is_unmapped and not read.mate_is_unmapped:
                    qc_group["self_mate_mapped"] += 1
                    if read.reference_id != read.next_reference_id:
                        qc_group["different_chr"] += 1
                        if read.mapping_quality >= 5:
                            qc_group["different_chr_mapq5"] += 1
                elif not read.is_unmapped and read.mate_is_unmapped:
                    qc_group["singletons"] += 1
            if read.is_secondary:
                counts["secondary"] += 1
                continue
            if read.is_supplementary:
                counts["supplementary"] += 1
                continue
            counts["primary"] += 1
            if read.is_paired:
                counts["paired_primary"] += 1
            if read.is_proper_pair:
                counts["proper_pair_primary"] += 1
            counts["duplicate_primary"] += int(read.is_duplicate)
            counts["od_tag_primary"] += int(read.has_tag("OD"))
            xu = read.has_tag("XU") and int(read.get_tag("XU")) == 1
            if xu:
                counts["rescue_primary"] += 1
                if not read.is_unmapped:
                    counts["rescue_mapped"] += 1
            else:
                counts["routed_primary"] += 1
                if not read.is_unmapped:
                    counts["routed_mapped"] += 1
                if read.has_tag("HP"):
                    hp = int(read.get_tag("HP"))
                    if hp == 1:
                        counts["routed_hp1"] += 1
                    elif hp == 2:
                        counts["routed_hp2"] += 1
                    else:
                        counts["routed_untagged"] += 1
                else:
                    counts["routed_untagged"] += 1
            if not read.is_unmapped:
                if read.has_tag("NM"):
                    counts["nm_sum"] += int(read.get_tag("NM"))
                counts["mapped_bp"] += int(read.query_length or len(read.query_sequence or ""))
    return counts, flagstat


def scan_realigned_output(config, output, target_fa, hap, *, workers=1, execution=None, resource_limits=None):
    """직렬과 동일한 집계를 반환하고 실행 정보는 별도 dict에 기록함."""
    output=Path(output);ref=str(target_fa) if config.output_format=='cram' else None
    with pysam.AlignmentFile(str(output),'rb',reference_filename=ref) as bam:
        indexed=bam.has_index();references=list(zip(bam.references,bam.lengths))
    plan=output_scan_plan(config,output,references,indexed,workers,resource_limits=resource_limits)
    print('[OUTPUT-QC] '+json.dumps({'phase':'plan',**plan}),flush=True)
    started=time.monotonic();count=plan['workers_used']
    info=execution if execution is not None else {}
    info.update(plan,status='running')
    try:
        if count==1:
            results=[_measured_partition(str(output),ref,None,plan['decode_threads'])]
        else:
            buckets=[[] for _ in range(count)];loads=[0]*count
            for contig,length in sorted(references,key=lambda x:x[1],reverse=True):
                i=min(range(count),key=lambda i:loads[i]);buckets[i].append(contig);loads[i]+=length
            buckets[min(range(count),key=lambda i:loads[i])].append('*')
            with ProcessPoolExecutor(max_workers=count,mp_context=multiprocessing.get_context('spawn'),
                                     initializer=limit_worker_address_space,initargs=(plan['worker_address_space_bytes'],)) as pool:
                results=list(pool.map(_measured_partition,[str(output)]*count,[ref]*count,buckets,[0]*count))
        counts,flagstat=_empty_output_scan()
        for part,flags,stats in results:
            for key,value in part.items():counts[key]+=value
            for group in flagstat:
                for key,value in flags[group].items():flagstat[group][key]+=value
        counts['routed_mapped_frac']=round(counts['routed_mapped']/max(counts['routed_primary'],1),6)
        counts['rescue_mapped_frac']=round(counts['rescue_mapped']/max(counts['rescue_primary'],1),6)
        counts['paired_primary_frac']=round(counts['paired_primary']/max(counts['primary'],1),6)
        counts['nm_per_base_all_mapped']=round(counts['nm_sum']/max(counts['mapped_bp'],1),8)
        counts['flagstat_text']=format_flagstat_text(flagstat)
        info.update(status='pass',wall_seconds=round(time.monotonic()-started,6),
                    worker_cpu_seconds=sum(x[2]['cpu_seconds'] for x in results),
                    worker_peak_rss_sum_bytes=sum(x[2]['peak_rss_bytes'] for x in results))
        return counts
    except BaseException:
        info.update(status='fail',wall_seconds=round(time.monotonic()-started,6));raise
    finally:print('[OUTPUT-QC] '+json.dumps({'phase':'complete',**info}),flush=True)


def format_flagstat_text(stats: dict[str, dict[str, int]]) -> str:
    """Format samtools flagstat-compatible text from the existing QC pass."""
    passed, failed = stats["passed"], stats["failed"]

    def pct(group: dict[str, int], numerator: str, denominator: str) -> str:
        total = group[denominator]
        return "N/A" if total == 0 else f"{100.0 * group[numerator] / total:.2f}%"

    def line(key: str, label: str, percent: tuple[str, str] | None = None) -> str:
        value = f"{passed[key]} + {failed[key]} {label}"
        if percent is not None:
            value += f" ({pct(passed, *percent)} : {pct(failed, *percent)})"
        return value

    rows = [
        line("total", "in total (QC-passed reads + QC-failed reads)"),
        line("primary", "primary"),
        line("secondary", "secondary"),
        line("supplementary", "supplementary"),
        line("duplicates", "duplicates"),
        line("primary_duplicates", "primary duplicates"),
        line("mapped", "mapped", ("mapped", "total")),
        line("primary_mapped", "primary mapped", ("primary_mapped", "primary")),
        line("paired", "paired in sequencing"),
        line("read1", "read1"),
        line("read2", "read2"),
        line("properly_paired", "properly paired", ("properly_paired", "paired")),
        line("self_mate_mapped", "with itself and mate mapped"),
        line("singletons", "singletons", ("singletons", "paired")),
        line("different_chr", "with mate mapped to a different chr"),
        line("different_chr_mapq5", "with mate mapped to a different chr (mapQ>=5)"),
    ]
    return "\n".join(rows) + "\n"


