"""입력 정렬과 처리 배치로 compact 계산의 최대 RAM을 사전에 예측한다."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import pysam

from ..memory_policy import GIB, memory_limit


def alignment_summary(path: Path, sample: str, contig: str) -> dict:
    qc = path.parent / f'{sample}.realign.json'
    total = None
    if qc.exists():
        data = json.loads(qc.read_text())
        total = data.get('routed', {}).get('primary')
    reference = path.parent / 'target.fa'
    with pysam.AlignmentFile(str(path), reference_filename=str(reference) if reference.exists() else None) as handle:
        lengths = dict(zip(handle.references, handle.lengths))
        read_lengths = []
        for read in handle.fetch(contig, 0, min(lengths[contig], 1_000_000)):
            if read.query_length:
                read_lengths.append(read.query_length)
            if len(read_lengths) >= 256:
                break
    source = '정렬 QC의 primary 읽기 수'
    if not total:
        total = max(1, path.stat().st_size // 40)
        source = 'QC 없음: 압축 파일 크기 기반 보수적 추정'
    return {'total_reads': int(total), 'total_reference_bases': sum(lengths.values()),
            'contig_lengths': lengths, 'read_length': max(read_lengths, default=150),
            'count_source': source, 'read_length_samples': len(read_lengths)}


def estimate_peak(summary: dict, contig: str, windows_per_chunk: int, window_size: int,
                  bounded_index: bool = True) -> dict:
    length = summary['contig_lengths'][contig]
    density = summary['total_reads'] / max(1, summary['total_reference_bases']) * 1.5
    reads = math.ceil(length * density)
    batch_reads = math.ceil(min(length, windows_per_chunk * window_size) * density)
    # read 객체·투영 사전·중복 키·패킹 임시 배열은 읽기당 4KB와 SEQ 복사 여유로 잡는다.
    batch = batch_reads * (4000 + 2 * summary['read_length'])
    # 누적 HP(1), PS(8), 길이(4), reference_end(8). 전체 키를 유지하면 Python 객체도 포함한다.
    index = reads * (21 if bounded_index else 621)
    reference_buffers = length * 4
    peak = math.ceil((GIB + batch + index + reference_buffers) * 1.3)
    # 긴 read의 국소 밀집은 전장 평균으로 예측되지 않으므로 추가 여유를 예약함.
    # 주소 공간 상한은 유지하고, 통계가 확인된 입력의 RSS 예약만 별도로 계산함.
    minimum = 32 * GIB if summary['read_length'] >= 1000 else 8 * GIB
    ceiling = max(minimum, math.ceil((peak + 2 * GIB) / (256 << 20)) * (256 << 20))
    model = os.environ.get('MFM_COMPACT_RESERVATION_MODEL', 'input_v2')
    if model not in {'input_v2', 'legacy'}:
        raise ValueError('MFM_COMPACT_RESERVATION_MODEL은 input_v2|legacy여야 함')
    trusted = (summary.get('count_source') == '정렬 QC의 primary 읽기 수'
               and summary.get('read_length_samples', 0) > 0)
    reservation = ceiling
    if model == 'input_v2' and trusted:
        # 기존 입력 피크+2GiB 여유/8GiB 기본값을 사용하되 VMS 한도는 유지한다.
        reservation = max(8 * GIB, math.ceil((peak + 2 * GIB) / (256 << 20)) * (256 << 20))
    return {'contig': contig, 'estimated_peak_bytes': peak, 'worker_ceiling_bytes': ceiling,
            'memory_reservation_bytes': reservation, 'address_space_limit_bytes': ceiling,
            'reservation_basis': '입력 통계+기존 여유' if reservation < ceiling else '기존 보수적 상한',
            'batch_reads': batch_reads, 'estimated_contig_reads': reads,
            'batch_bytes': batch, 'index_bytes': index, 'reference_buffer_bytes': reference_buffers}


def forecast_compact(jobs: list[dict], windows_per_chunk: int, window_size: int,
                     bounded_index: bool = True) -> dict:
    cached = {}
    estimates = []
    for job in jobs:
        if job.get('read_coordinate_mode') != 'personal_realigned':
            # 참조 좌표 투영은 원본 파일 수·긴 read에 따른 인덱스까지 보수적으로 예약한다.
            estimates.append({'contig': job['contig'], 'estimated_peak_bytes': 48 * GIB,
                              'worker_ceiling_bytes': 64 * GIB, 'memory_reservation_bytes': 64 * GIB,
                              'address_space_limit_bytes': 64 * GIB, 'reservation_basis': '참조 좌표 보수적 상한'})
            continue
        haps = []
        for hap in (1, 2):
            path = Path(job[f'hap{hap}_alignment'])
            if path not in cached:
                cached[path] = alignment_summary(path, job['sample'], job['contig'])
            haps.append(estimate_peak(cached[path], job['contig'], windows_per_chunk,
                                      window_size, bounded_index))
        representative = dict(max(haps, key=lambda row: row['estimated_peak_bytes']))
        for key in ['worker_ceiling_bytes', 'memory_reservation_bytes', 'address_space_limit_bytes']:
            representative[key] = max(row[key] for row in haps)
        estimates.append(representative)
    ceiling = max((row['worker_ceiling_bytes'] for row in estimates), default=8 * GIB)
    reservation = max((row['memory_reservation_bytes'] for row in estimates), default=8 * GIB)
    limit = memory_limit()
    usable = limit * 9 // 10
    if reservation + 4 * GIB > usable:
        raise ValueError('입력 기반 작업 1개의 RAM 예상량이 전체 상한을 초과합니다. 배치를 줄여야 합니다.')
    return {'schema': 'mfm_compact_memory_preflight_v2', 'memory_limit_bytes': limit,
            'windows_per_chunk': windows_per_chunk, 'window_size': window_size,
            'index_mode': 'bounded' if bounded_index else 'full',
            'reservation_model': os.environ.get('MFM_COMPACT_RESERVATION_MODEL', 'input_v2'),
            'processing_buffers': {
                'validation_record_chunk': int(os.environ.get('MFM_VALIDATION_RECORD_CHUNK', '262144')),
                'validation_cache_chunks': int(os.environ.get('MFM_VALIDATION_CACHE_CHUNKS', '4')),
                'trace_decoded_cache_mib': int(os.environ.get('MFM_TRACE_DECODED_CACHE_MIB', '64')),
                'trace_packed_batch_mib': int(os.environ.get('MFM_TRACE_PREFETCH_PACKED_MIB', '16')),
            },
            'worker_ceiling_bytes': ceiling, 'worker_reservation_bytes': reservation,
            'worker_address_space_bytes': ceiling, 'estimated_peak_per_worker_bytes': max(row['estimated_peak_bytes'] for row in estimates),
            'single_stage_memory_worker_cap': max(1, (usable - 4 * GIB) // reservation),
            'input_summaries': {str(path): value for path, value in cached.items()},
            'contigs': estimates, 'uncertainty': '평균 정렬 밀도 1.5배와 객체/버퍼 여유를 포함한 예상량; 실측 감시 병행'}
