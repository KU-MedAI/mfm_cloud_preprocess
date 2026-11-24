"""read별 홀수 nibble padding을 유지하는 bounded 4bit 일괄 변환."""
from __future__ import annotations
import numpy as np
import os

BATCH_BYTES = 4 << 20
BATCH_RECORDS = 65536


def pack_sequences(sequences, allowed, global_start_index, invalid_policy, lookup, symbols, fallback):
    lengths = np.fromiter((len(s) for s in sequences), dtype=np.int32, count=len(sequences))
    counts = (lengths.astype(np.int64) + 1) // 2
    offsets = np.empty(len(sequences) + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])
    packed = np.empty(int(offsets[-1]), dtype=np.uint8)
    valid = np.zeros(256, dtype=bool)
    for base in allowed:
        valid[int(lookup[ord(base)])] = True
    valid[255] = False
    invalid_counts = {}
    first = 0
    while first < len(sequences):
        end, size = first, 0
        while end < len(sequences) and end-first < BATCH_RECORDS:
            n = int(lengths[end])
            if end > first and size+n+(n & 1) > BATCH_BYTES:
                break
            size += n+(n & 1)
            end += 1
        try:
            data = b''.join(s.encode('ascii') + (b'A' if len(s) & 1 else b'') for s in sequences[first:end])
        except UnicodeEncodeError:
            # 이전 구현과 같은 최초 오류(record/position)를 유지함.
            fallback(sequences[first:end], allowed, global_start_index+first, invalid_policy)
            raise
        raw = np.frombuffer(data, dtype=np.uint8)
        if os.environ.get('MFM_ZARR_BACKEND','cpu')=='cupy':
            from ..gpu_arrays import pack_codes
            local_lengths=lengths[first:end].astype(np.int64)
            ends=np.cumsum(local_lengths+(local_lengths&1));odd=np.flatnonzero(local_lengths&1)
            encoded,bad,hist=pack_codes(raw,lookup,valid,ends[odd]-1,invalid_policy)
            if bad>=0:
                record=int(np.searchsorted(ends,bad,side='right'));start=int(ends[record-1]) if record else 0
                raise ValueError(f"invalid base for 4bit_acgtnp: base={chr(int(raw[bad]))!r} "
                                 f"record_index={global_start_index+first+record} position={bad-start} allowed={symbols}")
            for key,value in hist.items():invalid_counts[key]=invalid_counts.get(key,0)+value
            packed[int(offsets[first]):int(offsets[end])]=encoded
            first=end
            continue
        codes = lookup[raw]
        invalid = ~valid[codes]
        local_lengths = lengths[first:end].astype(np.int64)
        ends = np.cumsum(local_lengths + (local_lengths & 1))
        odd = np.flatnonzero(local_lengths & 1)
        if odd.size:
            invalid[ends[odd]-1] = False
        if invalid.any():
            if invalid_policy == 'fail':
                pos = int(np.flatnonzero(invalid)[0])
                record = int(np.searchsorted(ends, pos, side='right'))
                start = int(ends[record-1]) if record else 0
                raise ValueError(f"invalid base for 4bit_acgtnp: base={chr(int(raw[pos]))!r} "
                                 f"record_index={global_start_index+first+record} position={pos-start} allowed={symbols}")
            values, frequencies = np.unique(raw[invalid], return_counts=True)
            for value, frequency in zip(values, frequencies):
                key = chr(int(value)).upper()
                invalid_counts[key] = invalid_counts.get(key, 0) + int(frequency)
            codes[invalid] = int(lookup[ord('N')])
        packed[int(offsets[first]):int(offsets[end])] = (codes[0::2] << 4) | codes[1::2]
        first = end
    return packed, offsets, lengths, int(lengths.sum(dtype=np.int64)), invalid_counts
