"""reference·HAP·스레드·sort·QC 버퍼를 합산한 정렬 RAM 예약 모델임."""
import os
from ..memory_policy import GIB


def alignment_memory(config, lengths, haps, threads, sort_threads, environment=None):
    env = os.environ if environment is None else environment
    mode = env.get('MFM_REALIGN_MEMORY_MODEL', 'reference_v2')
    if mode not in {'reference_v2', 'legacy192'}:
        raise ValueError('MFM_REALIGN_MEMORY_MODEL은 reference_v2|legacy192여야 함')
    if not lengths or haps < 1 or haps > len(lengths) or threads < 1 or sort_threads < 0:
        raise ValueError('정렬 RAM 예상의 HAP/길이/스레드 값이 올바르지 않음')
    from .align_resources import parse_memory
    longest = max(int(value) for value in lengths.values())
    if longest < 1:
        raise ValueError('reference 길이는 양수여야 함')
    sort_unit = parse_memory(config.sort_mem)
    if sort_unit < 1:
        raise ValueError('sort 버퍼는 양수여야 함')
    reference = longest * 16 * haps
    minimap=str(getattr(config,'platform','')).lower() in {'hifi','ont','pacbio_hifi','nanopore'}
    shared_buffers = (28 if minimap else 24) * GIB * haps
    decoder=int(env.get('MFM_MINIMAP_DECODE_THREADS_EFFECTIVE',getattr(config,'input_decode_threads',0)))
    threads_buffer = (threads+decoder) * (GIB // 4) * haps
    pipeline_buffer = 4 * GIB * haps
    # samtools sort는 -@ 0/1에서 1개, 그 이상은 -@ 수만큼 버퍼를 할당함.
    sort_buffer = sort_unit * max(1, sort_threads) * haps
    raw = reference + shared_buffers + threads_buffer + pipeline_buffer + sort_buffer
    if mode == 'legacy192':
        estimated = max(192 * GIB, (reference + shared_buffers + sort_buffer) * 11 // 10)
    else:
        estimated = max(32 * GIB, (raw * 11 + 9) // 10)
    # GiB 단위로 올림하며 입력 크기·스레드가 증가할 때 단조 증가함.
    reserved = (estimated + GIB - 1) // GIB * GIB
    if env.get('MFM_WGS_BACKEND','cpu')=='parabricks':
        estimated=max(estimated,haps*(int(env.get('MFM_GPU_HOST_MEMORY_BYTES',str(128*GIB)))+8*GIB)*11//10)
        reserved=(estimated+GIB-1)//GIB*GIB
    return {'schema': 'mfm_alignment_memory_v2', 'mode': mode,
            'reference_bases_per_hap': longest, 'active_haps': haps,
            'reference_index_bytes': reference, 'hap_qc_buffers_bytes': shared_buffers,
            'thread_buffers_bytes': threads_buffer, 'input_decode_threads':decoder, 'pipeline_buffers_bytes': pipeline_buffer,
            'sort_buffers_bytes': sort_buffer, 'margin_percent': 10,
            'raw_estimated_bytes': raw, 'model_margin_bytes': max(0, reserved-raw),
            'estimation_basis': 'reference/thread/sort 입력 기반 보수적 모델; 실측 보정 전',
            'reserved_memory_bytes': reserved}
