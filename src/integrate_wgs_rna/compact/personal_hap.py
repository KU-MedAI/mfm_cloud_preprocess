"""한 HAP의 독립 Zarr만 생성함. 부모는 두 작업 종료 후 기존 형식으로 결합함."""
from pathlib import Path
import os,time,json
from numcodecs import blosc
from ..memory_policy import limit_worker_address_space, process_memory_snapshot

def write_hap(job,label,single_root,*,isolated=False):
    from ..direct_compact_writer import (ExportConfig,DirectTiming,hap_paths,prepare_hap_context,
        _window_index_tsv,_personal_alignment_reference,DirectCompactChrWriter,DEFAULT_COMPACT_COMPRESSOR,DEFAULT_COMPACT_COMPRESSION_LEVEL,
        build_paired_hap_window_artifacts,_merge_timing_dict)
    os.environ['MFM_COMPACT_ALIGNMENT_THREADS']=str(job.get('alignment_threads',1))
    os.environ['MFM_SAMPLE_CPU_BUDGET']=str(job.get('hap_cpu_cost',1))
    if isolated:os.environ.pop('MFM_CPU_POOL_LEASE_FILE',None)  # 부모가 두 HAP 전체를 예약함.
    blosc.set_nthreads(int(job.get('blosc_threads',1)))
    if isolated:limit_worker_address_space(int(job.get('hap_memory_ceiling_bytes',0)))
    began_wall=time.perf_counter();began_cpu=time.process_time()
    timing=DirectTiming();config=ExportConfig(**job['config']);sample=str(job['sample']);contig=str(job['contig'])
    profile_direct=bool(job.get('profile_direct_compact'));single_root=Path(single_root)
    alignment=Path(job[f'hap{label[-1]}_alignment']);alignments={label:alignment};alignment_refs={label:_personal_alignment_reference(alignment)}
    fasta,_=hap_paths(job,label)
    with timing.step('prepare_hap_context_'+label):
        context=prepare_hap_context(sample=sample,label=label,contig=contig,individual_fasta=fasta,chain_path=None,
            config=config,window_index_path=_window_index_tsv(job.get('window_index_dir','')))
    total=int(context['summary']['coordinate_eligible_window_count_total'])
    if job.get('max_windows') is not None:total=min(total,int(job['max_windows']))
    if total<=0:raise ValueError('HAP에 사용 가능한 window 없음')
    contexts={label:context};hap_window_counts={label:total}
    writer = DirectCompactChrWriter(
        single_root,
        final_target=single_root,
        sample=sample,
        contig=contig,
        hap_window_counts={label: hap_window_counts[label]},
        window_size=config.window_size,
        source_alignment=alignments[label],
        reference_fasta=alignment_refs[label],
        sample_metadata=job.get("sample_metadata", {}),
        compact_compressor=str(job.get("compact_compressor") or DEFAULT_COMPACT_COMPRESSOR),
        compact_compression_level=int(job.get("compact_compression_level") or DEFAULT_COMPACT_COMPRESSION_LEVEL),
        attrs={
            "read_coordinate_mode": "personal_realigned",
            "readset_layout": "hap_split",
            "source_bam_path": str(alignments[label]),
            "reference_fasta": str(alignment_refs[label]) if alignment_refs[label] else "",
            "window_index_dir": str(job.get("window_index_dir", "")),
            "source_hap1_fasta": str(job["hap1_fasta"]),
            "source_hap1_chain": "",
            "source_hap2_fasta": str(job["hap2_fasta"]),
            "source_hap2_chain": "",
            "config": config.to_dict(),
            "builder": "direct_compact_personal_realigned_hap_split",
        },
    )
    chunk_size = int(job["windows_per_chunk"])
    total = hap_window_counts[label]
    with timing.step(f"process_chunks_total_{label}"):
        for chunk_index, start_index in enumerate(range(0, total, chunk_size)):
            chunk_timing: dict[str, Any] | None = None
            chunk_started = time.perf_counter()
            if profile_direct:
                chunk_timing = {
                    "chunk_index": chunk_index,
                    "label": label,
                    "window_start_index": int(start_index),
                    "requested_windows": int(min(chunk_size, total - start_index)),
                }
            artifact_timing: dict[str, Any] = chunk_timing if chunk_timing is not None else {}
            artifacts = build_paired_hap_window_artifacts(
                sample=sample,
                contig=contig,
                bam_path=alignments[label],
                hap_contexts={label: contexts[label]},
                reference_fasta=alignment_refs[label],
                max_windows=min(chunk_size, total - start_index),
                window_start_index=start_index,
                config=config,
                timing=artifact_timing,
            )
            _merge_timing_dict(timing, artifact_timing)
            label_timing: dict[str, Any] = {"label": label}
            writer.append_hap_chunk(label, artifacts[label], timing=label_timing)
            _merge_timing_dict(timing, label_timing)
            next_index = start_index + min(chunk_size, total - start_index)
            if next_index < total and os.environ.get("MFM_COMPACT_INDEX_MODE", "bounded") == "bounded":
                boundary = int(contexts[label]["eligible_rows"][next_index][2].raw_start)
                writer.read_index.retire_before(boundary)
            # 다음 배치를 만들기 전에 이전 read 사전·문자열·배열 참조를 반환한다.
            del artifacts
            if chunk_timing is not None:
                chunk_timing["elapsed_seconds"] = round(time.perf_counter() - chunk_started, 6)
                timing.add_chunk(chunk_timing)
    meta = writer.finalize(validation_level=job["validation_level"], timing=timing)
    elapsed=time.perf_counter()-began_wall;cpu_seconds=time.process_time()-began_cpu
    print('[COMPACT-WORKER-MEMORY] ' + json.dumps(dict(process_memory_snapshot(),
        sample=sample,contig=contig,haplotype=label,wall_seconds=elapsed,
        cpu_seconds=cpu_seconds,actual_mean_cpu_cores=cpu_seconds/max(elapsed,1e-9),
        allocated_cpu_per_hap=int(job.get('hap_cpu_cost',1)))), flush=True)
    del writer
    return {'label':label,'metadata':meta,'timing':timing.to_dict()}

def parallel_haps(tasks):
    from concurrent.futures import ProcessPoolExecutor,as_completed
    import multiprocessing
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
    # CPU/RAM은 외부 contig worker가 이미 두 HAP 합계로 예약함.
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn'),max_tasks_per_child=1) as pool:
        futures={pool.submit(write_hap,*args,isolated=True):i for i,args in enumerate(tasks)}
        for future in as_completed(futures):
            try:result,error=future.result(),None
            except Exception as exc:result,error=None,exc
            yield futures[future],result,error
