# 검토용 원문 발췌임. 단독 실행 모듈이 아니며 import/호출 문맥은 전체 전달본을 사용함.

# 원본: vendor/MFM-DNARNAFM/pre_training/preprocess_raw_readset/shared_hap_writer.py:146-156
def _clipping_removed_query_sequence(read: pysam.AlignedSegment) -> tuple[str, int, int]:
    """Return mapped query sequence and its original query-coordinate bounds.

    pysam's query_alignment_start/end exclude soft clips. Hard clips are absent
    from query_sequence. The returned sequence preserves M/=/X and insertion
    bases while removing clipping bases.
    """
    sequence = (read.query_sequence or "").upper()
    start = int(read.query_alignment_start or 0)
    end = int(read.query_alignment_end or len(sequence))
    return sequence[start:end], start, end

# 원본: vendor/MFM-DNARNAFM/pre_training/preprocess_raw_readset/shared_hap_writer.py:159-180
def _pack_cigartuples(cigartuples: list[tuple[int, int]] | tuple[tuple[int, int], ...] | None) -> list[int]:
    """Pack CIGAR as uint32 values matching the stored clipping-removed query.

    HTSlib encodes one CIGAR op as (length << 4) | op_code. The stored read
    sequence in this writer has soft/hard clipping removed, so clipping ops are
    omitted here as well. This keeps query_start/query_end and CIGAR operations
    in the same query coordinate system used by the zarr read_sequence.
    """
    if not cigartuples:
        return []
    packed: list[int] = []
    for op_code, length in cigartuples:
        op_code = int(op_code)
        length = int(length)
        if op_code in (4, 5):  # soft/hard clipping, not present in stored query sequence
            continue
        if length < 0:
            raise ValueError(f"negative CIGAR length: {length}")
        if length >= (1 << 28):
            raise ValueError(f"CIGAR length too large for uint32 packing: {length}")
        packed.append((length << 4) | op_code)
    return packed

# 원본: vendor/MFM-DNARNAFM/pre_training/preprocess_raw_readset/shared_hap_writer.py:701-733
def _record_from_projected(
    read: pysam.AlignedSegment,
    projected,
    window,
    read_sequence: str,
    aln_query_start: int,
) -> dict[str, Any] | None:
    query_start = int(projected.query_start) - int(aln_query_start)
    query_end = int(projected.query_end) - int(aln_query_start)
    if query_start < 0 or query_end <= query_start or query_end > len(read_sequence):
        return None
    read_start = int(projected.window_start) - int(window.raw_start)
    read_end = int(projected.window_end) - int(window.raw_start)
    if read_start < 0 or read_end <= read_start or read_end > int(window.actual_length):
        return None
    return {
        "read_key": _read_key(read),
        "shared_read_key": _alignment_segment_key(read),
        "read_sequence": read_sequence,
        "read_haplotype": _read_haplotype(read),
        "read_phase_set": _read_phase_set(read),
        "read_index": -1,
        "query_start": query_start,
        "query_end": query_end,
        "read_start": read_start,
        "read_end": read_end,
        "read_strand": -1 if read.is_reverse else 1,
        "reference_start": int(read.reference_start),
        "reference_end": int(read.reference_end or read.reference_start),
        "cigar": read.cigarstring or "",
        "cigar_ops_packed": _pack_cigartuples(read.cigartuples),
        "segment_type": 1 if read.is_supplementary else 0,
    }

# 원본: vendor/MFM-DNARNAFM/pre_training/preprocess_raw_readset/shared_hap_writer.py:736-820
def _project_identity_read_to_window_segment(
    read: pysam.AlignedSegment,
    *,
    window_start: int,
    window_end: int,
) -> ProjectedWindowSegment | None:
    """Fast path for SNP-only/identity HAP coordinates.

    In identity coordinate mode, individual raw coordinates are equal to BAM
    reference coordinates. We can therefore avoid dense coord.ref_to_raw lookup
    and project CIGAR operations directly against the target window.
    """
    if not read.query_sequence or not read.cigartuples:
        return None

    query_position = 0
    reference_position = int(read.reference_start)
    raw_min: int | None = None
    raw_max: int | None = None
    previous_raw: int | None = None
    selected_raw_min: int | None = None
    selected_raw_max: int | None = None
    selected_query_min: int | None = None
    selected_query_max: int | None = None

    def _select_query(q0: int, q1: int, raw0: int, raw1: int) -> None:
        nonlocal selected_query_min, selected_query_max
        nonlocal selected_raw_min, selected_raw_max
        if q1 <= q0:
            return
        selected_query_min = q0 if selected_query_min is None else min(selected_query_min, q0)
        selected_query_max = q1 - 1 if selected_query_max is None else max(selected_query_max, q1 - 1)
        selected_raw_min = raw0 if selected_raw_min is None else min(selected_raw_min, raw0)
        selected_raw_max = raw1 if selected_raw_max is None else max(selected_raw_max, raw1)

    for operation, length_value in read.cigartuples:
        length = int(length_value)
        if operation in (0, 7, 8):  # M, =, X
            block_ref_start = reference_position
            block_ref_end = reference_position + length
            if length > 0:
                raw_min = block_ref_start if raw_min is None else min(raw_min, block_ref_start)
                raw_max = block_ref_end - 1 if raw_max is None else max(raw_max, block_ref_end - 1)
                previous_raw = block_ref_end - 1
            overlap_start = max(block_ref_start, int(window_start))
            overlap_end = min(block_ref_end, int(window_end))
            if overlap_start < overlap_end:
                q0 = query_position + (overlap_start - block_ref_start)
                q1 = query_position + (overlap_end - block_ref_start)
                _select_query(q0, q1, overlap_start, overlap_end - 1)
            query_position += length
            reference_position += length
        elif operation == 1:  # insertion relative to reference
            anchor = previous_raw + 1 if previous_raw is not None else reference_position
            if window_start <= anchor < window_end:
                _select_query(query_position, query_position + length, anchor, anchor)
            query_position += length
        elif operation in (2, 3):  # deletion or skipped reference
            reference_position += length
        elif operation == 4:  # soft clip
            query_position += length
        elif operation in (5, 6):  # hard clip or padding
            continue
        else:
            raise ValueError(f"Unsupported CIGAR operation: {operation}")

    if (
        raw_min is None
        or raw_max is None
        or selected_query_min is None
        or selected_query_max is None
        or selected_raw_min is None
        or selected_raw_max is None
    ):
        return None

    return ProjectedWindowSegment(
        raw_start=raw_min,
        raw_end=raw_max + 1,
        query_start=selected_query_min,
        query_end=selected_query_max + 1,
        strand="-" if read.is_reverse else "+",
        window_start=selected_raw_min,
        window_end=selected_raw_max + 1,
    )
