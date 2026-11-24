from __future__ import annotations

from collections import OrderedDict
from typing import Any

import numpy as np

from .runtime import VALIDATION_CACHE_CHUNKS, VALIDATION_RECORD_CHUNK


class ReadValidationChunkCache:
    """연속 metadata/CIGAR Zarr slice를 제한된 메모리에 보관한다.

    full validation은 같은 read를 CIGAR 무결성 검사와 window depth 재계산에서
    다시 참조한다. record 단위 scalar 접근 대신 연속 chunk를 읽어 I/O 호출 수를
    줄이되 LRU chunk 수를 제한해 대형 contig에서도 메모리를 고정한다.
    """

    _META_NAMES = (
        "cigar_offset",
        "cigar_len",
        "ref_start",
        "ref_end",
        "lengths",
        "segment_type",
    )

    def __init__(
        self,
        read_arrays: dict[str, Any],
        cigar_ops: Any,
        record_count: int,
        *,
        chunk_size: int = VALIDATION_RECORD_CHUNK,
        max_chunks: int = VALIDATION_CACHE_CHUNKS,
    ) -> None:
        self.read_arrays = read_arrays
        self.cigar_ops = cigar_ops
        self.record_count = int(record_count)
        self.chunk_size = max(1, int(chunk_size))
        self.max_chunks = int(max_chunks)
        if self.max_chunks < 0:
            raise ValueError("캐시 chunk 수는 음수일 수 없습니다.")
        self._cache: OrderedDict[int, dict[str, Any]] = OrderedDict()

    def _load(self, chunk_index: int) -> dict[str, Any]:
        lo = chunk_index * self.chunk_size
        hi = min(self.record_count, lo + self.chunk_size)
        metadata = {
            name: np.asarray(self.read_arrays[name][lo:hi])
            for name in self._META_NAMES
        }
        offsets = metadata["cigar_offset"].astype(np.int64, copy=False)
        lengths = metadata["cigar_len"].astype(np.int64, copy=False)
        if offsets.size:
            cigar_lo = int(offsets.min())
            cigar_hi = int(np.max(offsets + lengths))
        else:
            cigar_lo = cigar_hi = 0
        metadata.update(
            {
                "record_lo": lo,
                "cigar_lo": cigar_lo,
                "cigar_ops": np.asarray(
                    self.cigar_ops[cigar_lo:cigar_hi], dtype=np.uint32
                ),
            }
        )
        return metadata

    def record(self, read_id: int) -> tuple[dict[str, Any], int]:
        if read_id < 0 or read_id >= self.record_count:
            raise ValueError(f"full validation invalid read_id={read_id}")
        chunk_index = read_id // self.chunk_size
        if not self.max_chunks:
            chunk = self._load(chunk_index)
            return chunk, read_id - int(chunk["record_lo"])
        chunk = self._cache.pop(chunk_index, None)
        if chunk is None:
            chunk = self._load(chunk_index)
        self._cache[chunk_index] = chunk
        while len(self._cache) > self.max_chunks:
            self._cache.popitem(last=False)
        return chunk, read_id - int(chunk["record_lo"])

    def decode(
        self, read_id: int
    ) -> tuple[list[tuple[int, int]], int, int, int, int]:
        chunk, local_index = self.record(read_id)
        cigar_offset = int(chunk["cigar_offset"][local_index])
        cigar_len = int(chunk["cigar_len"][local_index])
        if cigar_len <= 0:
            raise ValueError(f"full validation empty CIGAR for read_id={read_id}")
        op_lo = cigar_offset - int(chunk["cigar_lo"])
        ops = chunk["cigar_ops"][op_lo : op_lo + cigar_len]
        if ops.shape[0] != cigar_len:
            raise ValueError(f"full validation truncated CIGAR for read_id={read_id}")
        ref_start = int(chunk["ref_start"][local_index])
        ref_end = int(chunk["ref_end"][local_index])
        read_length = int(chunk["lengths"][local_index])
        segment_type = int(chunk["segment_type"][local_index])
        ref_pos = ref_start
        reference_consumed = 0
        query_consumed = 0
        intervals: list[tuple[int, int]] = []
        for packed_value in ops:
            packed = int(packed_value)
            op = packed & 0x0F
            length = packed >> 4
            if op not in {0, 1, 2, 3, 6, 7, 8} or length <= 0:
                raise ValueError(
                    f"full validation invalid stored CIGAR op={op} "
                    f"length={length} read_id={read_id}"
                )
            if op in {0, 1, 7, 8}:
                query_consumed += length
            if op in {0, 2, 3, 7, 8}:
                reference_consumed += length
            if op in {0, 7, 8}:
                intervals.append((ref_pos, ref_pos + length))
                ref_pos += length
            elif op in {2, 3}:
                ref_pos += length
        if query_consumed != read_length:
            raise ValueError(
                f"full validation CIGAR/query mismatch read_id={read_id}: "
                f"cigar={query_consumed} stored={read_length}"
            )
        expected_reference = ref_end - ref_start
        if reference_consumed != expected_reference:
            raise ValueError(
                f"full validation CIGAR/reference mismatch read_id={read_id}: "
                f"cigar={reference_consumed} stored={expected_reference}"
            )
        return (
            intervals,
            query_consumed,
            reference_consumed,
            read_length,
            segment_type,
        )
