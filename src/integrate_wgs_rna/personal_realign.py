from __future__ import annotations

import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

import pysam
from .processing_buffers import minimap_arguments

from .prepare_alignment import bwa_mem2_version, ensure_bwa_mem2_index, reference_fingerprint
from .realign.phases import RealignPhases, reference_bases
from .realign.config import (
    DEFAULT_BWA_MEM2,
    DEFAULT_MINIMAP2,
    DEFAULT_SAMTOOLS,
    HAP_READ_MODES,
    HAP_READ_MODE_SETTINGS,
    ILLUMINA_DUPLICATE_MODES,
    ILLUMINA_HAP_SPECIFIC_INPUT_MODE,
    ILLUMINA_PAIRED_INPUT_MODE,
    ILLUMINA_ALIGNERS,
    MIN_ILLUMINA_PAIRED_PRIMARY_FRAC,
    RATIO_QC_POLICIES,
    READ_ROUTING_POLICIES,
    SCHEMA,
    UNKNOWN_READ_POLICIES,
    PersonalRealignConfig,
    PreparedPairedInput,
)
from .realign.policy import (
    alignment_contract,
    effective_hap_read_mode,
    is_hap_specific_mode,
    is_illumina_platform,
    normalize_platform,
    platform_minimap2_preset,
    resolve_hap_read_mode,
    sam_rg_platform,
    selected_aligner,
)
from .realign.qc import (
    add_alignment_input_metrics,
    carry_source_duplicate,
    classify_illumina_fragment,
    count_illumina_input_reads,
    count_input_reads,
    empty_input_counts,
    expected_routed_primary,
    htslib_decode_threads,
    illumina_prepared_path,
    iter_qname_groups,
    nm_per_base,
    prepare_illumina_hap_specific_inputs,
    prepared_qname_intersection_count,
    primary_hp,
    require_count_within_loss_tolerance,
    samtools_count,
    scan_realigned_output,
    select_read_for_haplotype,
    verify_realign_counts,
    write_realign_counts,
)


_BWA_INDEX_BUILD_LOCK = threading.Lock()


def gpu_wgs_enabled():
    backend=os.environ.get('MFM_WGS_BACKEND','cpu')
    if backend not in {'cpu','parabricks'}:raise ValueError('MFM_WGS_BACKEND는 cpu|parabricks여야 함')
    return backend=='parabricks'



def alignment_provenance(config: PersonalRealignConfig) -> dict[str, str]:
    """Metadata에 실제 정렬 입력 선택 방식을 기록한다."""
    if config.alignment is not None:
        return {
            "alignment_source": "explicit_path",
            "alignment_pattern": "",
            "alignment_path": str(config.alignment.resolve()),
        }
    return {
        "alignment_source": "pattern",
        "alignment_pattern": config.alignment_pattern,
        "alignment_path": "",
    }


def run_personal_realign_generation(config: PersonalRealignConfig) -> dict[str, Any]:
    started = time.time()
    resolved_hap_read_mode = resolve_hap_read_mode(config)
    validate_config(config)
    requested_preset = config.preset
    config.preset = platform_minimap2_preset(config.platform, config.preset)
    samples = read_samples(config.sample_list, config.samples)
    config = replace(config, jobs=min(config.jobs, max(1, len(samples))))
    contigs = expand_contigs(config.chroms)
    config.output_root.mkdir(parents=True, exist_ok=True)
    summary_tsv = config.summary_tsv or (config.output_root / f"{config.cohort}.personal_realign.summary.tsv")
    metadata_json = config.metadata_json or (config.output_root / f"{config.cohort}.personal_realign.metadata.json")

    if config.jobs <= 1:
        rows = [process_one_sample(config, sample, contigs) for sample in samples]
    else:
        rows = []
        with ProcessPoolExecutor(max_workers=config.jobs) as executor:
            futures = {executor.submit(process_one_sample, config, sample, contigs): sample for sample in samples}
            for future in as_completed(futures):
                rows.append(future.result())
        rows.sort(key=lambda row: row["sample"])

    write_summary(summary_tsv, rows)
    metadata = {
        "schema": SCHEMA,
        "cohort": config.cohort,
        "sample_count": len(samples),
        "chroms": contigs,
        "chromset": " ".join(contigs),
        "output_root": str(config.output_root),
        "sequence_variant_policy": config.sequence_variant_policy,
        "dataset_role": config.dataset_role,
        **alignment_provenance(config),
        "output_format": config.output_format,
        "platform": config.platform,
        "requested_preset": requested_preset,
        "preset": config.preset,
        "hap_read_mode": resolved_hap_read_mode,
        "read_routing_policy": config.read_routing_policy,
        "unknown_read_policy": config.unknown_read_policy,
        "include_unmapped": config.include_unmapped,
        "verify_counts": config.verify_counts,
        "count_tolerance_abs": config.count_tolerance_abs,
        "min_mapped_frac": config.min_mapped_frac,
        "samtools_bin": config.samtools_bin,
        "minimap2_bin": config.minimap2_bin,
        "bwa_mem2_bin": config.bwa_mem2_bin,
        "illumina_aligner": config.illumina_aligner,
        "bwa_mem2_index_root": str(config.bwa_mem2_index_root) if config.bwa_mem2_index_root else "",
        "illumina_duplicate_mode": config.illumina_duplicate_mode,
        "hap_jobs": config.hap_jobs,
        "route_jobs": config.route_jobs,
        "threads": {"align_per_hap": config.threads, "sort_per_hap": config.sort_threads},
        "dry_run": config.dry_run,
        "started_at_unix": started,
        "completed_at_unix": time.time(),
        "elapsed_sec": round(time.time() - started, 3),
        "ok_count": sum(1 for row in rows if row["status"] == "OK"),
        "skip_count": sum(1 for row in rows if row["status"] == "SKIP"),
        "fail_count": sum(1 for row in rows if row["status"] == "FAIL"),
        "summary_tsv": str(summary_tsv),
    }
    metadata_json.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    if metadata["fail_count"]:
        raise RuntimeError(f"personal realign failed for {metadata['fail_count']} sample(s); see {summary_tsv}")
    return metadata


def validate_config(config: PersonalRealignConfig) -> None:
    realigned_subdir = Path(config.realigned_subdir)
    if realigned_subdir.is_absolute() or ".." in realigned_subdir.parts or str(realigned_subdir) in {"", "."}:
        raise ValueError("--realigned-subdir must be a non-empty relative path without '..'")

    if config.hap_read_mode and config.hap_read_mode not in HAP_READ_MODES:
        raise ValueError(f"--hap-read-mode must be one of: {', '.join(sorted(HAP_READ_MODES))}")
    if config.sequence_source not in {"vcf_consensus", "assembly_direct"}:
        raise ValueError("--sequence-source must be vcf_consensus or assembly_direct")
    if (
        config.sequence_source == "assembly_direct"
        and effective_hap_read_mode(config) != "all_to_both"
        and not config.assembly_phase_labels_validated
    ):
        raise ValueError(
            "assembly_direct defaults to all_to_both; HP 기반 routing은 HAP1/HAP2 orientation 검증 뒤 "
            "--assembly-phase-labels-validated가 필요합니다"
        )
    if config.output_format != "cram":
        raise ValueError("--format must be cram: Stage 12 (legacy Stage70) verifies personal target.fa against CRAM @SQ M5")
    if config.read_routing_policy not in READ_ROUTING_POLICIES:
        raise ValueError(f"--read-routing-policy must be one of: {', '.join(sorted(READ_ROUTING_POLICIES))}")
    if config.unknown_read_policy not in UNKNOWN_READ_POLICIES:
        raise ValueError(f"--unknown-read-policy must be one of: {', '.join(sorted(UNKNOWN_READ_POLICIES))}")
    if config.illumina_aligner not in ILLUMINA_ALIGNERS:
        raise ValueError(f"--illumina-aligner must be one of: {', '.join(sorted(ILLUMINA_ALIGNERS))}")
    if config.illumina_duplicate_mode not in ILLUMINA_DUPLICATE_MODES:
        raise ValueError(f"--illumina-duplicate-mode must be one of: {', '.join(sorted(ILLUMINA_DUPLICATE_MODES))}")
    if config.count_tolerance_abs < 0:
        raise ValueError("--count-tolerance-abs must be >= 0")
    if config.hap_jobs not in {1, 2}:
        raise ValueError("--hap-jobs must be 1 or 2")
    if config.route_jobs < 1:
        raise ValueError("--route-jobs must be >= 1")
    if config.threads < 1 or config.sort_threads < 1:
        raise ValueError("--threads and --sort-threads must be >= 1")
    if not 0.0 <= config.min_mapped_frac <= 1.0:
        raise ValueError("--min-mapped-frac must be in [0, 1]")
    if config.ratio_qc_policy not in RATIO_QC_POLICIES:
        raise ValueError("--ratio-qc-policy must be record_only or strict")
    if config.sequence_variant_policy == "snp_only" and not config.allow_non_pretraining:
        raise ValueError("personal realign is disabled for snp_only by default; pass --allow-non-pretraining for validation-only runs")
    # post_training도 WGS를 사용하므로 dataset_role에 따른 realignment 차단은 하지 않는다.
    # post_training cleanup은 RNA 완료 전까지 80단계에서 자동으로 유예된다.
    if not shutil.which(config.samtools_bin) and not Path(config.samtools_bin).exists():
        raise FileNotFoundError(f"samtools not found: {config.samtools_bin}")
    use_bwa = is_illumina_platform(config.platform) and config.illumina_aligner == "bwa-mem2"
    if use_bwa and not shutil.which(config.bwa_mem2_bin) and not Path(config.bwa_mem2_bin).exists():
        raise FileNotFoundError(f"bwa-mem2 not found: {config.bwa_mem2_bin}")
    if not use_bwa and not shutil.which(config.minimap2_bin) and not Path(config.minimap2_bin).exists():
        raise FileNotFoundError(f"minimap2 not found: {config.minimap2_bin}")
    if config.alignment is not None and len(read_samples(config.sample_list, config.samples)) != 1:
        raise ValueError("--alignment may be used only with exactly one --sample")
    if (config.hap1_fasta_dir is None) ^ (config.hap2_fasta_dir is None):
        raise ValueError("--hap1-fasta-dir and --hap2-fasta-dir must be provided together")
    if config.hap1_fasta_dir is not None and len(read_samples(config.sample_list, config.samples)) != 1:
        raise ValueError("--hap1-fasta-dir/--hap2-fasta-dir overrides are single-sample options")
    if config.alignment is not None and config.alignment.suffix.lower() == ".cram" and config.reference_fasta is None:
        raise ValueError("CRAM input requires --reference-fasta")


def process_one_sample(config: PersonalRealignConfig, sample: str, contigs: list[str]) -> dict[str, Any]:
    started = time.time()
    timing_started = time.perf_counter()
    sample_root = config.output_root / sample
    realigned_root = sample_root / config.realigned_subdir
    log_dir = sample_root / "logs" / "personal_realign"
    scratch_root = os.environ.get("MFM_STAGE11_SCRATCH_ROOT", "").strip()
    scratch_key = hashlib.sha256(str(sample_root.resolve()).encode()).hexdigest()[:12]
    work_dir = (
        Path(scratch_root).resolve() / scratch_key / sample / "personal_realign"
        if scratch_root
        else sample_root / "work" / "personal_realign"
    )
    log_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    timing_path = sample_root / "qc" / f"{sample}.stage11_timing.json"
    timing_path.parent.mkdir(parents=True, exist_ok=True)
    timing: dict[str, Any] = {
        "schema": "mfm_stage11_timing_v1",
        "sample": sample,
        "hap_jobs": config.hap_jobs,
        "route_jobs": config.route_jobs,
        "threads": {"align_per_hap": config.threads, "sort_per_hap": config.sort_threads},
        "work_dir": str(work_dir),
        "scratch_enabled": bool(scratch_root),
        "steps_sec": {},
        "haplotypes": {},
        "status": "running",
    }
    phases = None

    try:
        hap_dirs = {
            1: config.hap1_fasta_dir or (sample_root / "individual_fasta" / "HAP1"),
            2: config.hap2_fasta_dir or (sample_root / "individual_fasta" / "HAP2"),
        }
        alignment = config.alignment or (config.output_root / config.alignment_pattern.format(sample=sample))
        if not alignment.exists():
            raise FileNotFoundError(f"alignment not found: {alignment}")
        if not has_index(alignment):
            raise FileNotFoundError(f"alignment index not found: {alignment}")

        illumina_input = is_illumina_platform(config.platform)
        needs_realign = (
            {
                hap: needs_illumina_paired_realign(config, sample, hap, realigned_root, contigs)
                for hap in (1, 2)
            }
            if illumina_input
            else {}
        )
        fuse_input_count = bool(
            illumina_input
            and config.verify_counts
            and config.read_routing_policy == "all_to_both"
            and any(needs_realign.values())
        )
        phase_enabled = os.environ.get('MFM_REALIGN_PHASE_RESOURCES','auto')
        long_phase=(not illumina_input and str(config.platform).lower() in {'hifi','pacbio_hifi'} and config.jobs==1 and not config.dry_run and phase_enabled!='0')
        count_config=config
        overlapped=None
        overlap_flag=os.environ.get('MFM_GPU_BAM_INDEX_OVERLAP','1')
        if overlap_flag not in {'0','1'}:raise ValueError('MFM_GPU_BAM_INDEX_OVERLAP는0/1이어야 함')
        if (overlap_flag=='1' and illumina_input and gpu_wgs_enabled() and config.jobs==1
                and not config.dry_run and phase_enabled!='0' and fuse_input_count
                and not is_hap_specific_mode(config)):
            from .realign.gpu_prepare_overlap import prepare as prepare_overlapped
            phases=RealignPhases();overlap_metrics={};overlap_started=time.perf_counter()
            overlapped=prepare_overlapped(config,sample,alignment,contigs,hap_dirs,realigned_root,
                work_dir,log_dir,[hap for hap in (1,2) if needs_realign[hap]],phases,overlap_metrics)
            timing['gpu_prepare_overlap']=overlap_metrics
            timing['steps_sec']['gpu_bam_index_prepare']=round(time.perf_counter()-overlap_started,3)
        if illumina_input and config.jobs==1 and not config.dry_run and phase_enabled!='0' and overlapped is None:
            phases=phases or RealignPhases()
            grant=phases.reserve('input_prepare',max_units=phases.total,min_units=min(4,phases.total),
                                 memory_per_unit=2*1024**3)
            # 실제 lease를 읽는 bounded 집계·collate가 준비 구간 전체 CPU를 사용한다.
            route_jobs=grant if os.environ.get('MFM_WGS_ROUTE_REQUESTED')=='auto' else min(config.route_jobs,grant)
            config=replace(config,route_jobs=route_jobs,fastq_threads=min(config.fastq_threads,max(0,grant-3)))
            count_config=replace(config,sort_threads=min(config.sort_threads,max(0,grant-1)))
        if long_phase:
            phases=RealignPhases()
            count_max=(1 if not config.verify_counts else min(phases.total,1+min(4,config.sort_threads))) if os.environ.get('MFM_REALIGN_INPUT_BACKEND','auto')=='serial' or not config.verify_counts else phases.total
            count_grant=phases.reserve('input_counts',max_units=count_max,memory_per_unit=2*1024**3)
            count_config=replace(config,sort_threads=min(config.sort_threads,max(0,count_grant-1)))
        step_started = time.perf_counter()
        input_counts = overlapped.input_counts if overlapped is not None else (empty_input_counts() if fuse_input_count else count_input_reads(count_config, alignment, contigs))
        timing["steps_sec"]["count_input_reads"] = round(time.perf_counter() - step_started, 3)
        timing["input_count_mode"] = "fused_with_prepare" if fuse_input_count else "separate_scan"
        prepared_by_hap: dict[int, PreparedPairedInput | None] = {1: overlapped, 2: overlapped}
        created_prepared: list[Path] = [overlapped.path] if overlapped is not None else []
        routing_qc: dict[str, Any] = {
            "schema": SCHEMA,
            "sample": sample,
            "platform": config.platform,
            "hap_read_mode": effective_hap_read_mode(config),
            "resolved_read_routing_policy": config.read_routing_policy,
            "resolved_unknown_read_policy": config.unknown_read_policy,
            "resolved_include_unmapped": config.include_unmapped,
            "input_record_counts": input_counts,
        }
        step_started = time.perf_counter()
        if illumina_input:
            routing_qc["routing_granularity"] = "QNAME_fragment" if is_hap_specific_mode(config) else "primary_read"
            if any(needs_realign.values()) and overlapped is None:
                if is_hap_specific_mode(config):
                    prepared_by_hap, strict_qc = prepare_illumina_hap_specific_inputs(
                        config, sample, alignment, contigs, work_dir
                    )
                    routing_qc["illumina_fragment_routing"] = strict_qc
                    created_prepared.extend(prepared.path for prepared in prepared_by_hap.values())
                elif config.read_routing_policy == "all_to_both":
                    prepared = prepare_illumina_paired_input(
                        config, sample, 1, alignment, contigs, work_dir,
                        None if fuse_input_count else input_counts, "all_to_both"
                    )
                    if prepared.input_counts is not None:
                        input_counts = prepared.input_counts
                    prepared_by_hap = {1: prepared, 2: prepared}
                    created_prepared.append(prepared.path)
                else:
                    for hap in (1, 2):
                        if needs_realign[hap]:
                            prepared = prepare_illumina_paired_input(
                                config,
                                sample,
                                hap,
                                alignment,
                                contigs,
                                work_dir,
                                input_counts,
                                f"hp{hap}_{config.unknown_read_policy}",
                            )
                            prepared_by_hap[hap] = prepared
                            created_prepared.append(prepared.path)
            routing_qc["input_record_counts"] = input_counts
            routing_qc["prepared_inputs"] = {
                f"HAP{hap}": prepared.stats if prepared is not None else {"status": "not_rebuilt"}
                for hap, prepared in prepared_by_hap.items()
            }
        else:
            routing_qc["routing_granularity"] = "primary_read"
        timing["steps_sec"]["prepare_routed_input"] = round(time.perf_counter() - step_started, 3)
        routing_qc_path = sample_root / "qc" / f"{sample}.hap_read_routing.json"
        routing_qc_path.parent.mkdir(parents=True, exist_ok=True)
        if created_prepared or not routing_qc_path.exists():
            routing_qc_path.write_text(json.dumps(routing_qc, indent=2, sort_keys=True) + "\n")
        prebuilt = {}
        alignment_request_config = config
        phase_enabled = os.environ.get('MFM_REALIGN_PHASE_RESOURCES', 'auto')
        if phase_enabled not in {'auto', '0', '1'}:
            raise ValueError('MFM_REALIGN_PHASE_RESOURCES는 auto|0|1이어야 함')
        if (illumina_input and not gpu_wgs_enabled() and selected_aligner(config) == 'bwa-mem2'
                and config.jobs == 1 and not config.dry_run and phase_enabled != '0'
                and any(needs_realign.values())
                and all(value is None or value.name_collated for value in prepared_by_hap.values())):
            if phases is None:
                phases = RealignPhases()
            lengths = reference_bases(hap_dirs, contigs)
            required = [hap for hap in (1, 2) if needs_realign[hap]]
            jobs = phases.index_jobs({hap: lengths[hap] for hap in required},
                                     os.environ.get('MFM_BWA_INDEX_HAP_JOBS', 'auto'))
            # HAP별 cache flock은 유지함. 서로 다른 HAP의 전역 생성 lock만 우회함.
            def build_index(hap):
                target = realigned_root / f'HAP{hap}' / 'target.fa'
                target.parent.mkdir(parents=True, exist_ok=True)
                began = time.time()
                build_target_fasta(config, hap_dirs[hap], target, contigs, ' '.join(contigs))
                prepared_at = time.time()
                resource_output = log_dir / f'{sample}.HAP{hap}.bwa_index.resource.json'
                with (log_dir / f'{sample}.HAP{hap}.bwa_index.err').open('ab') as handle:
                    index = ensure_bwa_mem2_index(target, bwa_mem2=config.bwa_mem2_bin,
                        index_root=config.bwa_mem2_index_root or (work_dir / 'bwa_mem2_index'), log=handle,
                        resource_output=resource_output, address_space_bytes=phases.index_address_space)
                execution = {'hap': hap, 'started_at': began, 'finished_at': time.time(),
                             'target_seconds': prepared_at - began, 'index_status': index['status']}
                if index['status'] == 'created':
                    measured = json.loads(resource_output.read_text())
                    execution.update(command_started_at=measured['started_at_unix'],
                        command_finished_at=measured['started_at_unix'] + measured['elapsed_sec'],
                        actual_cpu_seconds=measured['user_cpu_sec'] + measured['system_cpu_sec'],
                        max_rss_bytes=measured['max_rss_kib'] * 1024)
                return hap, index, execution
            executions = []
            with ThreadPoolExecutor(max_workers=jobs) as executor:
                for future in as_completed([executor.submit(build_index, hap) for hap in required]):
                    hap, index, execution = future.result()
                    prebuilt[hap] = index
                    executions.append(execution)
            timing['bwa_index_execution'] = sorted(executions, key=lambda row: row['hap'])
            timing['bwa_index_jobs'] = jobs
            config = phases.align_config(config, {hap: lengths[hap] for hap in required})
            timing['effective_threads'] = {'hap_jobs': config.hap_jobs, 'align_per_hap': config.threads,
                'sort_aux_per_hap': config.sort_threads, 'fastq_aux_per_hap': config.fastq_threads}
        elif long_phase or (gpu_wgs_enabled() and config.jobs==1 and not config.dry_run and phase_enabled!='0'):
            if phases is None:
                phases = RealignPhases()
            lengths=reference_bases(hap_dirs,contigs);required=[1,2];needs_realign={1:True,2:True}
            phases.reserve('target_fasta',max_units=2,memory_per_unit=2*1024**3)
            def prepare_target(hap):
                target=realigned_root/f'HAP{hap}'/'target.fa'
                build_target_fasta(config,hap_dirs[hap],target,contigs,' '.join(contigs))
            with ThreadPoolExecutor(max_workers=min(2,phases.resource_limits['cpus'])) as executor:
                list(executor.map(prepare_target,(1,2)))
            config=phases.align_config(alignment_request_config,lengths)
            timing['effective_threads']={'hap_jobs':config.hap_jobs,'align_per_hap':config.threads,'sort_aux_per_hap':config.sort_threads,'fastq_aux_per_hap':config.fastq_threads,'input_decode_per_hap':config.input_decode_threads}
        hap_resources=None
        if phases is not None and config.hap_jobs>1:
            from .realign.hap_resources import HaplotypeResources
            hap_resources=HaplotypeResources(phases,config,[1,2],local_root=log_dir/'gpu_host' if gpu_wgs_enabled() else None)
        def run_hap(hap: int, hap_config=None) -> dict[str, Any]:
            try:
                return run_one_haplotype(
                        hap_config if hap_config is not None else config,
                        sample,
                        hap,
                        hap_dirs[hap],
                        alignment,
                        contigs,
                        realigned_root,
                        work_dir,
                        log_dir,
                        input_counts,
                        prepared_by_hap[hap],
                        prebuilt.get(hap),
                        phase_controller=(hap_resources.controller(hap) if hap_resources is not None else phases),
                    )
            finally:
                if hap_resources is not None:hap_resources.release(hap)

        if config.hap_jobs == 1:
            hap_results = []
            for hap in (1, 2):
                hap_config = config
                if phases is not None and needs_realign.get(hap):
                    # 앞 HAP의 정렬·검증이 모두 끝난 안전한 경계에서만 재예약함.
                    hap_config = phases.align_config(alignment_request_config, {hap: lengths[hap]})
                    timing.setdefault('effective_threads_by_hap', {})[f'HAP{hap}'] = {
                        'hap_jobs': 1, 'align_per_hap': hap_config.threads,
                        'sort_aux_per_hap': hap_config.sort_threads}
                hap_results.append(run_hap(hap, hap_config))
        else:
            by_hap: dict[int, dict[str, Any]] = {}
            with ThreadPoolExecutor(max_workers=config.hap_jobs) as executor:
                futures = {executor.submit(run_hap, hap): hap for hap in (1, 2)}
                for future in as_completed(futures):
                    by_hap[futures[future]] = future.result()
                    if hap_resources is None and phases is not None:
                        remaining = {hap: lengths[hap] for hap in required if hap not in by_hap}
                        phases.hap_completed(config, remaining)
            hap_results = [by_hap[hap] for hap in (1, 2)]
        if hap_resources is not None:timing['hap_resource_phases']=hap_resources.events
        timing["haplotypes"] = {
            f"HAP{hap}": hap_results[hap - 1].get("timing_sec", {}) for hap in (1, 2)
        }
        if not config.dry_run and all(row["status"] in {"OK", "SKIP"} for row in hap_results):
            for path in sorted(set(created_prepared)):
                if path.exists():
                    path.unlink()
                from .realign.prepared_cache import discard
                discard(path)
            # --bwa-mem2-index-root를 명시하지 않은 production 기본은 sample scratch의
            # 고유 index다. 두 HAP 성공 뒤 즉시 제거해 수천 sample의 index 누적을 막는다.
            if config.bwa_mem2_index_root is None:
                ephemeral_index_root = work_dir / "bwa_mem2_index"
                if ephemeral_index_root.exists():
                    shutil.rmtree(ephemeral_index_root)
                timing["ephemeral_bwa_index_removed"] = True
        status = "SKIP" if all(row["status"] == "SKIP" for row in hap_results) else "OK"
        timing["status"] = status.lower()
        timing["elapsed_sec"] = round(time.perf_counter() - timing_started, 3)
        timing_path.write_text(json.dumps(timing, indent=2, sort_keys=True) + "\n")
        return {
            "sample": sample,
            "status": status,
            "hap_read_mode": effective_hap_read_mode(config),
            "read_routing_policy": config.read_routing_policy,
            "unknown_read_policy": config.unknown_read_policy,
            "routing_qc_json": str(routing_qc_path),
            "realigned_dir": str(realigned_root),
            "hap1_output": hap_results[0].get("output", ""),
            "hap2_output": hap_results[1].get("output", ""),
            "hap1_mapped_frac": hap_results[0].get("routed_mapped_frac", ""),
            "hap2_mapped_frac": hap_results[1].get("routed_mapped_frac", ""),
            "hap1_rescue_frac": hap_results[0].get("rescue_mapped_frac", ""),
            "hap2_rescue_frac": hap_results[1].get("rescue_mapped_frac", ""),
            "elapsed_sec": round(time.time() - started, 3),
            "timing_json": str(timing_path),
            "message": "ok",
        }
    except Exception as exc:  # noqa: BLE001 - row-level failure is persisted in TSV
        timing["status"] = "fail"
        timing["error"] = f"{type(exc).__name__}: {exc}"
        timing["elapsed_sec"] = round(time.perf_counter() - timing_started, 3)
        timing_path.write_text(json.dumps(timing, indent=2, sort_keys=True) + "\n")
        return {
            "sample": sample,
            "status": "FAIL",
            "hap_read_mode": effective_hap_read_mode(config),
            "read_routing_policy": config.read_routing_policy,
            "unknown_read_policy": config.unknown_read_policy,
            "routing_qc_json": "",
            "realigned_dir": str(realigned_root),
            "hap1_output": "",
            "hap2_output": "",
            "hap1_mapped_frac": "",
            "hap2_mapped_frac": "",
            "hap1_rescue_frac": "",
            "hap2_rescue_frac": "",
            "elapsed_sec": round(time.time() - started, 3),
            "timing_json": str(timing_path),
            "message": f"{type(exc).__name__}: {exc}",
        }
    finally:
        if phases is not None:
            phases.release()
            timing['resource_phases'] = phases.events
            timing_path.write_text(json.dumps(timing, indent=2, sort_keys=True) + "\n")


def run_one_haplotype(
    config: PersonalRealignConfig,
    sample: str,
    hap: int,
    hap_fasta_dir: Path,
    alignment: Path,
    contigs: list[str],
    realigned_root: Path,
    work_dir: Path,
    log_dir: Path,
    input_counts: dict[str, int],
    prepared_input: PreparedPairedInput | None,
    prebuilt_index: dict[str, Any] | None = None,
    phase_controller=None,
) -> dict[str, Any]:
    hap_started = time.perf_counter()
    hap_timing: dict[str, float] = {}
    hap_label = f"HAP{hap}"
    out_dir = realigned_root / hap_label
    out_dir.mkdir(parents=True, exist_ok=True)
    target_fa = out_dir / "target.fa"
    chromset = " ".join(contigs)
    step_started = time.perf_counter()
    build_target_fasta(config, hap_fasta_dir, target_fa, contigs, chromset)
    chain_qc = chain_provenance(config, hap_fasta_dir, contigs)
    hap_timing["target_and_chain"] = round(time.perf_counter() - step_started, 3)
    contract = alignment_contract(config)

    prefix = config.output_prefix.format(sample=sample, hap=hap, haplotype=hap_label)
    output = out_dir / f"{prefix}.{config.output_format}"
    index_candidates = output_index_candidates(output)
    realign_json = out_dir / f"{sample}.realign.json"
    if output.exists() and any(path.exists() for path in index_candidates) and realign_json.exists() and not config.overwrite:
        previous = json.loads(realign_json.read_text())
        previous_paired = previous.get("input_preparation", {}).get("mode") in {
            ILLUMINA_PAIRED_INPUT_MODE, ILLUMINA_HAP_SPECIFIC_INPUT_MODE
        }
        if previous_matches_hap_read_mode(previous, config) and previous_matches_qc_contract(previous, config) and previous.get("chromset") == chromset and (
            not is_illumina_platform(config.platform) or previous_paired
        ) and (not is_illumina_platform(config.platform) or (
            previous.get("alignment_contract") == contract and previous.get("chain_provenance") == chain_qc
        )):
            hap_timing["total"] = round(time.perf_counter() - hap_started, 3)
            return {"status": "SKIP", "output": str(output), "timing_sec": hap_timing, **previous.get("routed", {}), **previous.get("rescue_unmapped", {})}
        # 새 실행 실패 후 이전 JSON이 새 CRAM provenance로 오인되지 않게 한다.
        realign_json.unlink(missing_ok=True)

    if is_illumina_platform(config.platform) and prepared_input is None:
        raise RuntimeError(f"HAP{hap} Illumina realignment requires a prepared paired input")

    if config.overwrite:
        for path in [output, realign_json, out_dir / f"{sample}.realign_counts.tsv", out_dir / f"{sample}.flagstat.txt", *index_candidates]:
            if path.exists():
                path.unlink()

    bwa_index: dict[str, Any] | None = prebuilt_index
    if not gpu_wgs_enabled() and selected_aligner(config) == "bwa-mem2" and bwa_index is None:
        step_started = time.perf_counter()
        index_root = config.bwa_mem2_index_root or (work_dir / "bwa_mem2_index")
        index_log = log_dir / f"{sample}.HAP{hap}.bwa_index.err"
        with index_log.open("ab") as log_handle:
            # 두 개인 HAP index의 동시 생성은 peak RSS를 크게 올릴 수 있다.
            # cache 조회/생성을 직렬화하고 실제 정렬만 병렬 실행한다.
            with _BWA_INDEX_BUILD_LOCK:
                bwa_index = ensure_bwa_mem2_index(
                    target_fa, bwa_mem2=config.bwa_mem2_bin, index_root=index_root, log=log_handle
                )
        hap_timing["bwa_index"] = round(time.perf_counter() - step_started, 3)
    step_started = time.perf_counter()
    from .realign.recover_output import completed_output_proof
    recovered = completed_output_proof(config, sample, hap, alignment, target_fa, output, log_dir)
    script = Path(recovered['script']) if recovered else write_realign_script(
        config, sample, hap, alignment, target_fa, output, contigs, work_dir, log_dir,
        prepared_input, bwa_index,
    )
    hap_timing["write_script"] = round(time.perf_counter() - step_started, 3)
    if config.dry_run:
        info = {
            "schema": SCHEMA,
            "sample": sample,
            "haplotype": hap_label,
            "status": "dry_run",
            "script": str(script),
            "target": str(target_fa),
            "output": str(output),
            "chromset": chromset,
            "hap_read_mode": effective_hap_read_mode(config),
            "output_qc_execution": {"status": "not_run", "reason": "dry_run"},
            "qc_contract": {
                "min_mapped_frac": config.min_mapped_frac,
                "status": "planned",
            },
            "read_routing_policy": config.read_routing_policy,
            "unknown_read_policy": config.unknown_read_policy,
            "alignment_contract": contract,
            "bwa_mem2_index": bwa_index,
            "chain_provenance": chain_qc,
            "input_preparation": (
                {
                    "mode": (
                        ILLUMINA_HAP_SPECIFIC_INPUT_MODE
                        if prepared_input is not None and is_hap_specific_mode(config)
                        else ILLUMINA_PAIRED_INPUT_MODE
                    ),
                    "records": prepared_input.stats if prepared_input is not None else {},
                }
                if prepared_input is not None
                else {"mode": "legacy_single_read_stream"}
            ),
        }
        realign_json.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n")
        hap_timing["total"] = round(time.perf_counter() - hap_started, 3)
        return {"status": "OK", "output": str(output), "routed_mapped_frac": "", "rescue_mapped_frac": "", "timing_sec": hap_timing}

    step_started = time.perf_counter()
    resource_usage_path = log_dir / f"{sample}.HAP{hap}.resource.json"
    resource_runner = Path(__file__).resolve().parents[2] / "scripts" / "resource_runner.py"
    child_env=dict(os.environ)
    if phase_controller is not None and getattr(getattr(phase_controller,'owner',None),'local',False):
        child_env['MFM_GPU_LOCAL_LEASE_FILE']=str(phase_controller.owner.path)
    if recovered is None:
        subprocess.run(
            [sys.executable, str(resource_runner), "--output", str(resource_usage_path), "--", "bash", str(script)],
            check=True, env=child_env,
        )
    hap_timing["align_restore_sort"] = round(time.perf_counter() - step_started, 3)
    step_started = time.perf_counter()
    index_aux = min(4, max(0, config.threads + config.sort_threads + 2))
    if phase_controller is not None:
        index_cpu = phase_controller.reserve('output_index', max_units=min(5,phase_controller.total),
                                             base_memory=8*1024**3)
        index_aux = max(0,index_cpu-1)
    subprocess.run([config.samtools_bin, "index", f"-@{index_aux}", str(output)], check=True)
    hap_timing["output_index"] = round(time.perf_counter() - step_started, 3)
    step_started = time.perf_counter()
    scan_execution = {}
    if phase_controller is not None:
        phase_controller.reserve('output_qc', max_units=min(len(contigs)+1,phase_controller.total),
                                 memory_per_unit=4*1024**3)
    scan = scan_realigned_output(config, output, target_fa, hap, workers=None, execution=scan_execution, resource_limits=phase_controller.resource_limits if phase_controller is not None else None)
    write_flagstat(
        config,
        output,
        out_dir / f"{sample}.flagstat.txt",
        target_fa if config.output_format == "cram" else None,
        scan=scan,
    )
    write_realign_counts(out_dir / f"{sample}.realign_counts.tsv", scan, hap)
    verify_realign_counts(config, hap, input_counts, scan, prepared_input)
    hap_timing["scan_and_count_qc"] = round(time.perf_counter() - step_started, 3)
    hap_timing["scan_workers_used"] = scan_execution["workers_used"]
    step_started = time.perf_counter()
    axis_config=config
    if phase_controller is not None:
        axis_grant=phase_controller.reserve('axis_qc', max_units=min(2,phase_controller.total),base_memory=8*1024**3)
        axis_config=replace(config,sort_threads=min(config.sort_threads,max(0,axis_grant-1)))
    nm_ref = nm_per_base(axis_config, alignment, config.reference_fasta, contigs[0], config.inspect_reads)
    nm_personal = nm_per_base(axis_config, output, target_fa if config.output_format == "cram" else None, contigs[0], config.inspect_reads)
    header_check = validate_personal_axis(output, target_fa, config.output_format, contigs)
    hap_timing["nm_and_axis_qc"] = round(time.perf_counter() - step_started, 3)
    hap_timing["total"] = round(time.perf_counter() - hap_started, 3)

    result = {
        "schema": SCHEMA,
        "sample": sample,
        "cohort": config.cohort,
        "policy": config.sequence_variant_policy,
        "dataset_role": config.dataset_role,
        "haplotype": hap_label,
        "target": str(target_fa),
        "chromset": chromset,
        "output": str(output),
        "format": config.output_format,
        "hap_read_mode": effective_hap_read_mode(config),
        "routing_granularity": "fragment" if is_illumina_platform(config.platform) and is_hap_specific_mode(config) else "primary_read",
        "execution_backend": 'parabricks' if gpu_wgs_enabled() else 'cpu',
        "aligner": ('Parabricks4.8 GPU '+('BWA-MEM' if is_illumina_platform(config.platform) else 'minimap2')
                    if gpu_wgs_enabled() else bwa_mem2_version(config.bwa_mem2_bin)
                    if selected_aligner(config) == "bwa-mem2" else minimap2_version(config.minimap2_bin)),
        "platform": config.platform,
        "preset": config.preset,
        "aligner_args": contract["args"],
        "alignment_contract": contract,
        "target_fingerprint": reference_fingerprint(target_fa),
        "bwa_mem2_index": bwa_index,
        "chain_provenance": chain_qc,
        "input_preparation": (
            {
                "mode": (
                    ILLUMINA_HAP_SPECIFIC_INPUT_MODE
                    if is_hap_specific_mode(config)
                    else ILLUMINA_PAIRED_INPUT_MODE
                ),
                "source": str(prepared_input.path),
                "route_label": prepared_input.route_label,
                "records": prepared_input.stats,
                "command": "samtools fastq -T HP,PS,XU,OD -0 /dev/null -s /dev/null -n -o - <qname-collated BAM>",
                "duplicate_mode": config.illumina_duplicate_mode,
            }
            if prepared_input is not None
            else {"mode": "legacy_single_read_stream"}
        ),
        "output_qc_execution": scan_execution,
        "qc_contract": {
            "min_mapped_frac": config.min_mapped_frac,
            "observed_routed_mapped_frac": scan["routed_mapped_frac"],
            "routed_mapped_fraction_below_reference": scan["routed_mapped_frac"] < config.min_mapped_frac,
            "observed_paired_primary_frac": scan["paired_primary_frac"],
            "paired_primary_reference": MIN_ILLUMINA_PAIRED_PRIMARY_FRAC,
            "paired_primary_fraction_below_reference": (
                scan["paired_primary_frac"] < MIN_ILLUMINA_PAIRED_PRIMARY_FRAC
            ),
            "ratio_qc_policy": config.ratio_qc_policy,
            "qc_action": config.ratio_qc_policy,
            "status": "pass",
        },
        "read_routing_policy": config.read_routing_policy,
        "unknown_read_policy": config.unknown_read_policy,
        "include_unmapped": config.include_unmapped,
        "routing": routing_description(config, hap),
        "input_reads": {
            "hp1": input_counts["hp1"],
            "hp2": input_counts["hp2"],
            "untagged": input_counts["untagged"],
            "all_primary_mapped": input_counts["all_primary_mapped"],
            "expected_routed": expected_routed_primary(config, hap, input_counts),
            "unmapped": input_counts["unmapped"],
        },
        "records": {
            "primary": scan["primary"],
            "secondary": scan["secondary"],
            "supplementary": scan["supplementary"],
            "paired_primary": scan["paired_primary"],
            "proper_pair_primary": scan["proper_pair_primary"],
            "paired_primary_frac": scan["paired_primary_frac"],
            "duplicate_primary": scan["duplicate_primary"],
            "od_tag_primary": scan["od_tag_primary"],
        },
        "routed": {
            "primary": scan["routed_primary"],
            "mapped": scan["routed_mapped"],
            "mapped_frac": scan["routed_mapped_frac"],
            "hp_tagged": scan[f"routed_hp{hap}"],
            "hp1": scan["routed_hp1"],
            "hp2": scan["routed_hp2"],
            "untagged": scan["routed_untagged"],
        },
        "rescue_unmapped": {
            "primary": scan["rescue_primary"],
            "newly_mapped": scan["rescue_mapped"],
            "frac": scan["rescue_mapped_frac"],
            "note": "참조축 unmapped read가 개인 서열에 새로 붙은 비율",
        },
        "contig_len_vs_reference": header_check,
        "nm_per_base_all_mapped": scan["nm_per_base_all_mapped"],
        "nm_per_base_sample": {
            "contig": contigs[0],
            "reference_axis": nm_ref["nm_per_base"],
            "personal_axis": nm_personal["nm_per_base"],
            "reads_ref": nm_ref["reads"],
            "reads_personal": nm_personal["reads"],
            "verdict": "개선" if nm_personal["nm_per_base"] < nm_ref["nm_per_base"] else "비개선",
        },
        "filters_applied_here": "none; MAPQ/secondary/supplementary filtering is deferred to read-set export",
        "command_script": str(script),
        "resource_usage_file": str(resource_usage_path) if resource_usage_path.exists() else "",
        "timing_sec": hap_timing,
    }
    if recovered is not None:
        result['recovered_completed_alignment'] = recovered
    realign_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if config.output_format == "bam" and not config.keep_target:
        for path in [target_fa, Path(str(target_fa) + ".fai"), out_dir / "target.chromset"]:
            if path.exists():
                path.unlink()
    return {
        "status": "OK",
        "output": str(output),
        "routed_mapped_frac": scan["routed_mapped_frac"],
        "rescue_mapped_frac": scan["rescue_mapped_frac"],
        "timing_sec": hap_timing,
    }


def build_target_fasta(config: PersonalRealignConfig, hap_fasta_dir: Path, target_fa: Path, contigs: list[str], chromset: str) -> None:
    stamp = target_fa.parent / "target.chromset"
    if target_fa.exists() and Path(str(target_fa) + ".fai").exists() and stamp.exists() and stamp.read_text() == chromset:
        return
    target_fa.parent.mkdir(parents=True, exist_ok=True)
    part = target_fa.with_suffix(target_fa.suffix + ".part")
    for path in [target_fa, Path(str(target_fa) + ".fai"), stamp, part]:
        if path.exists():
            path.unlink()
    with part.open("w") as out:
        for contig in contigs:
            src = hap_fasta_dir / f"{contig}.fa"
            if not src.exists():
                raise FileNotFoundError(f"missing personal FASTA: {src}")
            with src.open() as fh:
                shutil.copyfileobj(fh, out, length=1024 * 1024)
    part.rename(target_fa)
    subprocess.run([config.samtools_bin, "faidx", str(target_fa)], check=True)
    validate_target_lengths(hap_fasta_dir, target_fa, contigs)
    stamp.write_text(chromset)


def write_realign_script(
    config: PersonalRealignConfig,
    sample: str,
    hap: int,
    alignment: Path,
    target_fa: Path,
    output: Path,
    contigs: list[str],
    work_dir: Path,
    log_dir: Path,
    prepared_input: PreparedPairedInput | None,
    bwa_index: dict[str, Any] | None,
) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    if gpu_wgs_enabled():
        from .gpu_wgs import write_script
        return write_script(config,sample,hap,alignment,target_fa,output,contigs,work_dir,log_dir,prepared_input)
    part = output.with_suffix(output.suffix + ".part")
    if prepared_input is not None:
        return write_illumina_paired_realign_script(
            config, sample, hap, prepared_input, target_fa, output, work_dir, log_dir, bwa_index
        )
    contig_args = " ".join(q(contig) for contig in contigs)
    view_ref = ""
    if config.reference_fasta is not None:
        view_ref = f"VIEW_REF=(-T {q(config.reference_fasta)})"
    else:
        view_ref = "VIEW_REF=()"
    fmt = ""
    if config.output_format == "cram":
        fmt = f"--output-fmt cram --reference {q(target_fa)}"
    read_group = (
        f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:{sam_rg_platform(config.platform)}"
        f"\\tDS:realigned_to_personal_HAP{hap};routing={config.read_routing_policy};unknown={config.unknown_read_policy}"
    )
    script = f"""#!/usr/bin/env bash
set -uo pipefail
SAMTOOLS={q(config.samtools_bin)}
MINIMAP2={q(config.minimap2_bin)}
INPUT={q(alignment)}
TARGET={q(target_fa)}
OUTPUT={q(output)}
PART={q(part)}
SORT_PREFIX={q(work_dir / f'sort.{sample}.HAP{hap}')}
{view_ref}
UNTAGGED_EXPR='[HP]==0 || ![HP]'
rm -f "$PART"
(
  grp() {{
    local code=$1
    shift
    "$SAMTOOLS" view -u "${{VIEW_REF[@]}}" "$@" | "$SAMTOOLS" fastq -T HP,PS -
    local s=("${{PIPESTATUS[@]}}")
    [[ "${{s[0]}}" == 0 && "${{s[1]}}" == 0 ]] || exit "$code"
  }}
{routing_fastq_commands(config, hap, contig_args)}"""
    if config.include_unmapped:
        script += """  "$SAMTOOLS" view -u "${VIEW_REF[@]}" -F 0x900 "$INPUT" '*' | "$SAMTOOLS" fastq -T HP,PS - | awk 'NR%4==1{print $0 "\\tXU:i:1"; next} {print}'
  s=("${PIPESTATUS[@]}")
  [[ "${s[0]}" == 0 && "${s[1]}" == 0 && "${s[2]}" == 0 ]] || exit 93
"""
    script += f"""  exit 0
) 2>{q(log_dir / f'{sample}.HAP{hap}.fastq.err')} \\
  | (
    "$MINIMAP2" -ax {q(config.preset)} -y -Y --MD -t {int(config.threads)}{minimap_arguments(os.environ)} \\
      -R {q(read_group)} \\
      "$TARGET" - 2>{q(log_dir / f'{sample}.HAP{hap}.mm2.err')}
    mapper_status=$?
    if [[ -n "${{MFM_CPU_POOL_LEASE_FILE:-}}" ]]; then
      {q(sys.executable)} -m integrate_wgs_rna.realign.stream_resources --hap {hap} --sort-threads {config.sort_threads} >&2 || true
    fi
    exit "$mapper_status"
  ) \\
  | "$SAMTOOLS" sort -@ {int(config.sort_threads)} -m {q(config.sort_mem)} -T "$SORT_PREFIX" {fmt} -o "$PART" - 2>{q(log_dir / f'{sample}.HAP{hap}.sort.err')}
st=("${{PIPESTATUS[@]}}")
if [[ "${{st[0]}}" != 0 || "${{st[1]}}" != 0 || "${{st[2]}}" != 0 ]]; then
  echo "samtools_fastq=${{st[0]}} minimap2=${{st[1]}} sort=${{st[2]}}" >&2
  rm -f "$PART"
  exit 1
fi
mv "$PART" "$OUTPUT"
"""
    if selected_aligner(config)=='minimap2' and config.input_decode_threads>0:
        script=script.replace('"$SAMTOOLS" view -u ',f'"$SAMTOOLS" view -@ {config.input_decode_threads} -u ')
    script_path = log_dir / f"{sample}.HAP{hap}.realign.sh"
    script_path.write_text(script)
    os.chmod(script_path, 0o755)
    return script_path


def write_illumina_paired_realign_script(
    config: PersonalRealignConfig,
    sample: str,
    hap: int,
    prepared_input: PreparedPairedInput,
    target_fa: Path,
    output: Path,
    work_dir: Path,
    log_dir: Path,
    bwa_index: dict[str, Any] | None,
) -> Path:
    """Write the Illumina paired pipeline with an explicit aligner contract.

    ``samtools fastq`` requires QNAME-collated input. The prepared BAM has
    already applied routing and provenance tags. Only complete R1/R2 fragments
    enter the interleaved aligner input; singleton and READ_OTHER are discarded.
    """
    part = output.with_suffix(output.suffix + ".part")
    prepared_alignment = prepared_input.path
    fmt = f"--output-fmt cram --reference {q(target_fa)}" if config.output_format == "cram" else ""
    read_group = (
        f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
        f"\\tDS:realigned_to_personal_HAP{hap};routing={config.read_routing_policy};unknown={config.unknown_read_policy};paired=preserved"
    )
    if config.illumina_aligner == "bwa-mem2":
        if not bwa_index or not bwa_index.get("index_prefix"):
            raise RuntimeError("BWA-MEM2 index metadata is missing")
        tool_lines = f"BWA_MEM2={q(config.bwa_mem2_bin)}"
        align_command = (
            f'"$BWA_MEM2" mem -p -Y -C -t {int(config.threads)} '
            f'-R {q(read_group)} {q(bwa_index["index_prefix"])} - '
            f'2>{q(log_dir / f"{sample}.HAP{hap}.bwa_mem2.err")}'
        )
    else:
        tool_lines = f"MINIMAP2={q(config.minimap2_bin)}"
        align_command = (
            f'"$MINIMAP2" -ax {q(config.preset)} -y -Y --MD -t {int(config.threads)} '
            f'-R {q(read_group)} "$TARGET" - '
            f'2>{q(log_dir / f"{sample}.HAP{hap}.mm2.err")}'
        )
    tool_lines += (
        f"\nPYTHON={q(sys.executable)}\n"
        f"DUP_FILTER={q(Path(__file__).with_name('restore_duplicate_flags.py'))}"
    )
    duplicate_filter = (
        f'| "$PYTHON" "$DUP_FILTER" --mode {q(config.illumina_duplicate_mode)} '
        f'2>{q(log_dir / f"{sample}.HAP{hap}.duplicate_restore.json")}'
    )
    script = f"""#!/usr/bin/env bash
set -uo pipefail
SAMTOOLS={q(config.samtools_bin)}
{tool_lines}
INPUT={q(prepared_alignment)}
TARGET={q(target_fa)}
OUTPUT={q(output)}
PART={q(part)}
COLLATE_PREFIX={q(work_dir / f'collate.{sample}.HAP{hap}')}
SORT_PREFIX={q(work_dir / f'sort.{sample}.HAP{hap}')}
rm -f "$PART"
{illumina_fastq_source_command(prepared_input, log_dir, sample, hap, config.fastq_threads)} \\
  | {align_command} \\
  {duplicate_filter} \\
  | "$SAMTOOLS" sort -@ {int(config.sort_threads)} -m {q(config.sort_mem)} -T "$SORT_PREFIX" {fmt} -o "$PART" - 2>{q(log_dir / f'{sample}.HAP{hap}.sort.err')}
st=("${{PIPESTATUS[@]}}")
for code in "${{st[@]}}"; do
  if [[ "$code" != 0 ]]; then
    echo "Illumina realignment pipeline failed: PIPESTATUS=${{st[*]}}" >&2
    rm -f "$PART"
    exit 1
  fi
done
mv "$PART" "$OUTPUT"
"""
    script_path = log_dir / f"{sample}.HAP{hap}.realign.sh"
    script_path.write_text(script)
    os.chmod(script_path, 0o755)
    return script_path


def illumina_fastq_source_command(
    prepared_input: PreparedPairedInput, log_dir: Path, sample: str, hap: int, threads: int = 2
) -> str:
    fastq = (
        f'"$SAMTOOLS" fastq -@ {int(threads)} -T HP,PS,XU,OD -0 /dev/null -s /dev/null '
        f'-n -o - "$INPUT" 2>{q(log_dir / f"{sample}.HAP{hap}.fastq.err")}'
    )
    if prepared_input.name_collated:
        return fastq
    return (
        f'"$SAMTOOLS" collate -u -O -T "$COLLATE_PREFIX" "$INPUT" '
        f'2>{q(log_dir / f"{sample}.HAP{hap}.collate.err")} | {fastq}'
    )


def routing_description(config: PersonalRealignConfig, hap: int) -> str:
    if config.read_routing_policy == "all_to_both":
        routed = "all mapped primary reads"
    elif config.unknown_read_policy == "include_both":
        routed = f"HP:i:{hap} + untagged/HP0"
    else:
        routed = f"HP:i:{hap}; untagged/HP0 dropped"
    if config.include_unmapped:
        routed += " + unmapped(XU:i:1)"
    return routed


def routing_fastq_commands(config: PersonalRealignConfig, hap: int, contig_args: str) -> str:
    if config.read_routing_policy == "all_to_both":
        return f'''  grp 91 -F 0x904 "$INPUT" {contig_args}\n'''
    commands = [f'''  grp 91 -F 0x904 -e '[HP]=={hap}' "$INPUT" {contig_args}''']
    if config.unknown_read_policy == "include_both":
        commands.append(f'''  grp 92 -F 0x904 -e "$UNTAGGED_EXPR" "$INPUT" {contig_args}''')
    return "\n".join(commands) + "\n"


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_personal_chain(
    chain: Path, contig: str, reference_length: int, consensus_length: int
) -> dict[str, Any]:
    if not chain.exists() or not chain.stat().st_size:
        raise RuntimeError(f"chain 파일이 없거나 비었습니다: {chain}")
    with chain.open("r", encoding="ascii", errors="replace") as handle:
        header = handle.readline().strip().split()
    if len(header) < 13 or header[0] != "chain":
        raise RuntimeError(f"chain header가 올바르지 않습니다: {chain}")
    qc = {
        "chain": str(chain),
        "target_name": header[2], "target_size": int(header[3]), "target_strand": header[4],
        "query_name": header[7], "query_size": int(header[8]), "query_strand": header[9],
    }
    if qc["target_size"] != reference_length:
        raise RuntimeError(f"chain tSize가 reference 길이와 다릅니다: {qc}")
    if qc["query_size"] != consensus_length:
        raise RuntimeError(f"chain qSize가 개인 FASTA 길이와 다릅니다: {qc}")
    if qc["target_strand"] != "+" or qc["query_strand"] != "+":
        raise RuntimeError(f"chain strand가 plus/plus가 아닙니다: {qc}")
    if qc["target_name"] != contig or qc["query_name"] != contig:
        raise RuntimeError(f"chain contig 이름이 예상과 다릅니다: {qc}, expected={contig}")
    qc["status"] = "pass"
    return qc


def chain_provenance(
    config: PersonalRealignConfig, hap_fasta_dir: Path, contigs: list[str]
) -> dict[str, Any]:
    if config.sequence_variant_policy != "snv_indel":
        return {"required": False, "direction": "not_applicable", "contigs": []}
    if config.reference_fasta is None:
        raise ValueError("snv_indel chain 검증에는 --reference-fasta가 필요합니다")
    reference_fai = Path(str(config.reference_fasta) + ".fai")
    if not reference_fai.exists():
        raise FileNotFoundError(f"reference FASTA index not found: {reference_fai}")
    reference_lengths = fasta_lengths(reference_fai)
    rows: list[dict[str, Any]] = []
    for contig in contigs:
        fasta = hap_fasta_dir / f"{contig}.fa"
        chain = hap_fasta_dir / f"{contig}.chain"
        personal_lengths = fasta_lengths(Path(str(fasta) + ".fai"))
        qc = validate_personal_chain(chain, contig, reference_lengths[contig], personal_lengths[contig])
        rows.append({**qc, "md5": md5_file(chain)})
    return {
        "required": True,
        "direction": "GRCh38/reference(target) -> personal_HAP(query)",
        "reference_fasta": str(config.reference_fasta),
        "contigs": rows,
    }


def needs_illumina_paired_realign(
    config: PersonalRealignConfig, sample: str, hap: int, realigned_root: Path, contigs: list[str]
) -> bool:
    if config.overwrite:
        return True
    hap_label = f"HAP{hap}"
    out_dir = realigned_root / hap_label
    prefix = config.output_prefix.format(sample=sample, hap=hap, haplotype=hap_label)
    output = out_dir / f"{prefix}.{config.output_format}"
    realign_json = out_dir / f"{sample}.realign.json"
    if not output.exists() or not any(path.exists() for path in output_index_candidates(output)) or not realign_json.exists():
        return True
    try:
        previous = json.loads(realign_json.read_text())
    except json.JSONDecodeError:
        return True
    return not (
        previous_matches_hap_read_mode(previous, config) and previous_matches_qc_contract(previous, config)
        and previous.get("alignment_contract") == alignment_contract(config)
        and previous.get("chromset") == " ".join(contigs)
        and previous.get("input_preparation", {}).get("mode") in {
            ILLUMINA_PAIRED_INPUT_MODE, ILLUMINA_HAP_SPECIFIC_INPUT_MODE
        }
    )


def previous_matches_hap_read_mode(previous: dict[str, Any], config: PersonalRealignConfig) -> bool:
    """Accept old all-to-both/hp-routed metadata, but never old strict Illumina routing."""
    if previous.get('execution_backend','cpu') != ('parabricks' if gpu_wgs_enabled() else 'cpu'):
        return False
    wanted = effective_hap_read_mode(config)
    previous_mode = previous.get("hap_read_mode")
    if previous_mode:
        return previous_mode == wanted
    if wanted == "all_to_both":
        return previous.get("read_routing_policy") == "all_to_both"
    if wanted == "hp_routed_with_unknown":
        return (
            previous.get("read_routing_policy") == "hp_routed"
            and previous.get("unknown_read_policy") == "include_both"
        )
    return False



def previous_matches_qc_contract(previous: dict[str, Any], config: PersonalRealignConfig) -> bool:
    """Do not reuse record-only outputs when legacy strict gating is requested."""
    try:
        contract = previous["qc_contract"]
        recorded = float(contract["min_mapped_frac"])
        observed = float(previous["routed"]["mapped_frac"])
    except (KeyError, TypeError, ValueError):
        return False
    recorded_policy = contract.get("ratio_qc_policy", "strict")
    if config.ratio_qc_policy == "record_only":
        return recorded_policy in RATIO_QC_POLICIES
    return (
        recorded_policy == "strict"
        and recorded >= config.min_mapped_frac
        and observed >= config.min_mapped_frac
    )

def prepare_illumina_paired_input(config, sample, hap, alignment, contigs, work_dir,
                                  input_counts, route_label):
    from .realign.prepared_cache import prepare
    if os.environ.get('MFM_WGS_BACKEND')=='parabricks' and os.environ.get('MFM_GPU_SCRATCH_ROOT'):
        from .gpu_wgs import scratch_directory
        work_dir=scratch_directory(work_dir)
    output = illumina_prepared_path(work_dir, sample, route_label)
    return prepare(config, sample, hap, alignment, contigs, work_dir, input_counts,
                   route_label, output, _prepare_illumina_paired_input_uncached)


def _prepare_illumina_paired_input_uncached(
    config: PersonalRealignConfig,
    sample: str,
    hap: int,
    alignment: Path,
    contigs: list[str],
    work_dir: Path,
    input_counts: dict[str, int] | None,
    route_label: str,
) -> PreparedPairedInput:
    """Materialize routed Illumina records before name collation.

    Reads are selected with the same flag/HP rules as the former samtools-view
    commands. Original-unmapped records receive XU:i:1 before collation, so
    the paired fastq stream retains their provenance without header rewriting.
    """
    output = illumina_prepared_path(work_dir, sample, route_label)
    part = output.with_suffix(output.suffix + ".part")
    if part.exists():
        part.unlink()
    if output.exists():
        output.unlink()

    # all_to_both는 HP 판정이 필요 없으므로 indexed contig별 독립 scan이 가능하다.
    # 각 part를 source header의 contig 순서로 결합한 뒤 기존 name-collate/QC를
    # 그대로 적용한다. 따라서 read 선택·tag·검증 계약은 유지하면서 Python의
    # 단일 iterator 병목만 제거한다. HP 기반 routing은 fragment 의미 보존을 위해
    # 기존 직렬 경로를 유지한다.
    if config.route_jobs > 1 and config.read_routing_policy == "all_to_both":
        return prepare_illumina_paired_input_parallel(
            config, sample, hap, alignment, contigs, work_dir,
            input_counts, route_label, output, part,
        )

    contig_set = set(contigs)
    stats = {
        "routed_primary": 0,
        "rescue_unmapped_primary": 0,
        "routed_hp1": 0,
        "routed_hp2": 0,
        "routed_untagged": 0,
        "source_duplicate_primary": 0,
    }
    collect_input_counts = input_counts is None
    observed_input_counts = empty_input_counts()
    open_kwargs: dict[str, Any] = {"threads": htslib_decode_threads(config)}
    if alignment.suffix.lower() == ".cram" and config.reference_fasta is not None:
        open_kwargs["reference_filename"] = str(config.reference_fasta)
    with pysam.AlignmentFile(str(alignment), "rb", **open_kwargs) as source:
        with pysam.AlignmentFile(str(part), "wb", header=source.header) as dest:
            for read in source.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary:
                    continue
                if read.is_unmapped:
                    if collect_input_counts:
                        observed_input_counts["unmapped"] += 1
                    if config.include_unmapped:
                        read.set_tag("XU", 1, value_type="i", replace=True)
                        carry_source_duplicate(config, read, stats)
                        dest.write(read)
                        stats["rescue_unmapped_primary"] += 1
                    continue
                if read.reference_name not in contig_set:
                    continue
                if collect_input_counts:
                    observed_input_counts["all_primary_mapped"] += 1
                    observed_hp = int(read.get_tag("HP")) if read.has_tag("HP") else 0
                    if observed_hp == 1:
                        observed_input_counts["hp1"] += 1
                    elif observed_hp == 2:
                        observed_input_counts["hp2"] += 1
                    else:
                        observed_input_counts["untagged"] += 1
                if not select_read_for_haplotype(config, read, hap):
                    continue
                carry_source_duplicate(config, read, stats)
                dest.write(read)
                stats["routed_primary"] += 1
                hp = int(read.get_tag("HP")) if read.has_tag("HP") else 0
                if hp == 1:
                    stats["routed_hp1"] += 1
                elif hp == 2:
                    stats["routed_hp2"] += 1
                else:
                    stats["routed_untagged"] += 1
    effective_input_counts = observed_input_counts if collect_input_counts else input_counts
    if effective_input_counts is None:
        raise RuntimeError("internal error: Illumina input counts were not resolved")
    if config.verify_counts:
        require_count_within_loss_tolerance(
            f"prepared HAP{hap} routed primary",
            stats["routed_primary"],
            expected_routed_primary(config, hap, effective_input_counts),
            0,
        )
        require_count_within_loss_tolerance(
            f"prepared HAP{hap} XU rescue primary",
            stats["rescue_unmapped_primary"],
            effective_input_counts["unmapped"] if config.include_unmapped else 0,
            0,
        )
    collate_prefix = work_dir / f"collate.prepare.{sample}.{route_label}"
    subprocess.run(
        [config.samtools_bin, "collate", "-l", "1", "-@", str(max(1, config.sort_threads)),
         "-o", str(output), "-T", str(collate_prefix), str(part)],
        check=True,
    )
    part.unlink()
    add_alignment_input_metrics(output, stats, threads=htslib_decode_threads(config),
                                    samtools_bin=config.samtools_bin, concurrency=config.jobs)
    return PreparedPairedInput(
        path=output, route_label=route_label, stats=stats,
        input_counts=effective_input_counts, name_collated=True,
    )


def _route_illumina_region(task: tuple[str, str, str, str, bool, bool, str]) -> dict[str, int]:
    """Write one indexed mapped contig (or ``*`` unmapped bin) to a BAM part."""
    alignment_s, reference_s, region, part_s, include_unmapped, collect_counts, duplicate_mode = task
    stats = {
        "routed_primary": 0,
        "rescue_unmapped_primary": 0,
        "routed_hp1": 0,
        "routed_hp2": 0,
        "routed_untagged": 0,
        "source_duplicate_primary": 0,
        "input_hp1": 0,
        "input_hp2": 0,
        "input_untagged": 0,
        "input_all_primary_mapped": 0,
        "input_unmapped": 0,
    }
    open_kwargs: dict[str, Any] = {"threads": 1}
    if reference_s:
        open_kwargs["reference_filename"] = reference_s
    with pysam.AlignmentFile(alignment_s, "rb", **open_kwargs) as source:
        with pysam.AlignmentFile(part_s, "wb", header=source.header, threads=1) as dest:
            iterator = source.fetch("*") if region == "*" else source.fetch(region)
            for read in iterator:
                if read.is_secondary or read.is_supplementary:
                    continue
                if region == "*":
                    if not read.is_unmapped:
                        continue
                    if collect_counts:
                        stats["input_unmapped"] += 1
                    if not include_unmapped:
                        continue
                    read.set_tag("XU", 1, value_type="i", replace=True)
                    stats["rescue_unmapped_primary"] += 1
                else:
                    if read.is_unmapped:
                        continue
                    hp = int(read.get_tag("HP")) if read.has_tag("HP") else 0
                    if collect_counts:
                        stats["input_all_primary_mapped"] += 1
                        if hp == 1:
                            stats["input_hp1"] += 1
                        elif hp == 2:
                            stats["input_hp2"] += 1
                        else:
                            stats["input_untagged"] += 1
                    stats["routed_primary"] += 1
                    if hp == 1:
                        stats["routed_hp1"] += 1
                    elif hp == 2:
                        stats["routed_hp2"] += 1
                    else:
                        stats["routed_untagged"] += 1
                if read.is_duplicate:
                    stats["source_duplicate_primary"] += 1
                    if duplicate_mode == "preserve_source":
                        read.set_tag("OD", 1, value_type="i", replace=True)
                dest.write(read)
    return stats


def prepare_illumina_paired_input_parallel(
    config: PersonalRealignConfig,
    sample: str,
    hap: int,
    alignment: Path,
    contigs: list[str],
    work_dir: Path,
    input_counts: dict[str, int] | None,
    route_label: str,
    output: Path,
    part: Path,
) -> PreparedPairedInput:
    """Parallel all-to-both routing with deterministic source-contig ordering."""
    collect_input_counts = input_counts is None
    reference_s = (
        str(config.reference_fasta)
        if alignment.suffix.lower() == ".cram" and config.reference_fasta is not None
        else ""
    )
    open_kwargs: dict[str, Any] = {"threads": 1}
    if reference_s:
        open_kwargs["reference_filename"] = reference_s
    with pysam.AlignmentFile(str(alignment), "rb", **open_kwargs) as source:
        wanted = set(contigs)
        ordered_regions = [name for name in source.references if name in wanted]
    missing = sorted(set(contigs) - set(ordered_regions))
    if missing:
        raise RuntimeError(f"indexed routing contig missing from alignment header: {', '.join(missing)}")
    if config.include_unmapped or collect_input_counts:
        ordered_regions.append("*")

    region_parts = [work_dir / f"{sample}.{route_label}.route.{idx:03d}.bam" for idx in range(len(ordered_regions))]
    for path in region_parts:
        if path.exists():
            path.unlink()
    tasks = [
        (
            str(alignment), reference_s, region, str(region_parts[idx]),
            config.include_unmapped, collect_input_counts, config.illumina_duplicate_mode,
        )
        for idx, region in enumerate(ordered_regions)
    ]
    try:
        with ProcessPoolExecutor(max_workers=min(config.route_jobs, len(tasks))) as executor:
            region_stats = list(executor.map(_route_illumina_region, tasks))
        keys = (
            "routed_primary", "rescue_unmapped_primary", "routed_hp1", "routed_hp2",
            "routed_untagged", "source_duplicate_primary",
        )
        stats = {key: sum(row[key] for row in region_stats) for key in keys}
        stats["route_jobs"] = min(config.route_jobs, len(tasks))
        observed_input_counts = {
            "hp1": sum(row["input_hp1"] for row in region_stats),
            "hp2": sum(row["input_hp2"] for row in region_stats),
            "untagged": sum(row["input_untagged"] for row in region_stats),
            "all_primary_mapped": sum(row["input_all_primary_mapped"] for row in region_stats),
            "unmapped": sum(row["input_unmapped"] for row in region_stats),
        }
        effective_input_counts = observed_input_counts if collect_input_counts else input_counts
        if effective_input_counts is None:
            raise RuntimeError("internal error: parallel Illumina input counts were not resolved")
        if config.verify_counts:
            require_count_within_loss_tolerance(
                f"prepared HAP{hap} routed primary", stats["routed_primary"],
                expected_routed_primary(config, hap, effective_input_counts), 0,
            )
            require_count_within_loss_tolerance(
                f"prepared HAP{hap} XU rescue primary", stats["rescue_unmapped_primary"],
                effective_input_counts["unmapped"] if config.include_unmapped else 0, 0,
            )
        # cat 결과를 중간 BAM으로 다시 쓰지 않고 collate stdin으로 직접 보낸다.
        # contig part와 결합 BAM의 동시 보관 및 36GB급 추가 write/read를 피한다.
        collate_prefix = work_dir / f"collate.prepare.{sample}.{route_label}"
        collate_threads = max(config.sort_threads, min(8, config.route_jobs))
        cat_proc = subprocess.Popen(
            [config.samtools_bin, "cat", *map(str, region_parts)],
            stdout=subprocess.PIPE,
        )
        if cat_proc.stdout is None:
            raise RuntimeError("samtools cat stdout pipe was not created")
        try:
            collate_proc = subprocess.run(
                [config.samtools_bin, "collate", "-l", "1", "-@", str(max(1, collate_threads)),
                 "-o", str(output), "-T", str(collate_prefix), "-"],
                stdin=cat_proc.stdout,
                check=False,
            )
        finally:
            cat_proc.stdout.close()
        cat_rc = cat_proc.wait()
        if cat_rc != 0 or collate_proc.returncode != 0:
            raise subprocess.CalledProcessError(
                cat_rc or collate_proc.returncode,
                f"samtools cat | samtools collate (cat={cat_rc}, collate={collate_proc.returncode})",
            )
        add_alignment_input_metrics(output, stats, threads=htslib_decode_threads(config),
                                    samtools_bin=config.samtools_bin, concurrency=config.jobs)
        return PreparedPairedInput(
            path=output, route_label=route_label, stats=stats,
            input_counts=effective_input_counts, name_collated=True,
        )
    finally:
        if part.exists():
            part.unlink()
        for path in region_parts:
            if path.exists():
                path.unlink()




def validate_personal_axis(output: Path, target_fa: Path, output_format: str, contigs: list[str]) -> dict[str, int]:
    target_lengths = fasta_lengths(Path(str(target_fa) + ".fai"))
    ref = str(target_fa) if output_format == "cram" else None
    same = 0
    diff = 0
    with pysam.AlignmentFile(
        str(output), "rb", reference_filename=ref, threads=4
    ) as bam:
        for contig in contigs:
            header_len = bam.get_reference_length(contig)
            target_len = target_lengths[contig]
            if header_len != target_len:
                raise RuntimeError(f"{contig} output header length {header_len} != target FASTA length {target_len}")
            same += 1
    return {"validated_contigs": same, "length_mismatch": diff}


def validate_target_lengths(hap_fasta_dir: Path, target_fa: Path, contigs: list[str]) -> None:
    target_lengths = fasta_lengths(Path(str(target_fa) + ".fai"))
    for contig in contigs:
        src_fai = hap_fasta_dir / f"{contig}.fa.fai"
        if not src_fai.exists():
            raise FileNotFoundError(f"missing FASTA index: {src_fai}")
        src_lengths = fasta_lengths(src_fai)
        if target_lengths.get(contig) != src_lengths.get(contig):
            raise RuntimeError(f"{contig} target length {target_lengths.get(contig)} != source length {src_lengths.get(contig)}")


def fasta_lengths(fai: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with fai.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out[parts[0]] = int(parts[1])
    return out


def write_flagstat(
    config: PersonalRealignConfig,
    output: Path,
    flagstat: Path,
    reference_fasta: Path | None,
    *,
    scan: dict[str, Any] | None = None,
) -> None:
    if scan is not None and isinstance(scan.get("flagstat_text"), str):
        flagstat.write_text(scan["flagstat_text"])
        return
    cmd = [config.samtools_bin, "flagstat", "-@4"]
    if reference_fasta is not None:
        cmd.extend(["--input-fmt-option", f"reference={reference_fasta}"])
    cmd.append(str(output))
    result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode == 0:
        flagstat.write_text(result.stdout)
    else:
        flagstat.write_text(f"[WARN] flagstat failed code={result.returncode}\n{result.stderr}\n")


def output_index_candidates(output: Path) -> list[Path]:
    if output.suffix == ".cram":
        return [Path(str(output) + ".crai"), output.with_suffix(".crai")]
    return [Path(str(output) + ".bai"), output.with_suffix(".bai")]


def has_index(alignment: Path) -> bool:
    return any(path.exists() for path in output_index_candidates(alignment))


def read_samples(sample_list: Path | None, samples: list[str]) -> list[str]:
    values: list[str] = []
    if sample_list is not None:
        with sample_list.open() as fh:
            values.extend(line.strip() for line in fh if line.strip() and not line.startswith("#"))
    values.extend(samples)
    out: list[str] = []
    seen = set()
    for sample in values:
        if sample not in seen:
            seen.add(sample)
            out.append(sample)
    if not out:
        raise ValueError("no samples provided; use --sample-list or --sample")
    return out


def expand_contigs(values: list[str]) -> list[str]:
    expanded: list[str] = []
    for value in values:
        for part in str(value).split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                left, right = part.split("-", 1)
                left_i = int(normalize_contig(left).removeprefix("chr"))
                right_i = int(normalize_contig(right).removeprefix("chr"))
                expanded.extend(f"chr{i}" for i in range(left_i, right_i + 1))
            else:
                expanded.append(normalize_contig(part))
    seen = set()
    out: list[str] = []
    for contig in expanded:
        if contig not in seen:
            seen.add(contig)
            out.append(contig)
    return out


def normalize_contig(value: str) -> str:
    value = value.strip()
    return value if value.startswith("chr") else f"chr{value}"


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "sample",
        "status",
        "hap_read_mode",
        "read_routing_policy",
        "unknown_read_policy",
        "routing_qc_json",
        "realigned_dir",
        "hap1_output",
        "hap2_output",
        "hap1_mapped_frac",
        "hap2_mapped_frac",
        "hap1_rescue_frac",
        "hap2_rescue_frac",
        "elapsed_sec",
        "message",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def minimap2_version(minimap2_bin: str) -> str:
    try:
        result = subprocess.run([minimap2_bin, "--version"], check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return f"minimap2 {result.stdout.strip()}"
    except Exception:  # noqa: BLE001
        return "minimap2 unknown"


def q(value: str | Path) -> str:
    return shlex.quote(str(value))
