# 검토용 원문 발췌임. 단독 실행 모듈이 아니며 import/호출 문맥은 전체 전달본을 사용함.

# 원본: src/integrate_wgs_rna/direct_compact_writer.py:1549-1669
def _register_reads_bulk(self, records: Iterable[dict[str, Any]], timing: dict[str, Any] | None = None) -> None:
    step_started = time.perf_counter()
    pending: dict[tuple[Any, ...], dict[str, Any]] = {}
    pending_order: list[tuple[Any, ...]] = []
    for record in records:
        key = tuple(record.get("shared_read_key", record["read_key"]))
        hap = int(record["read_haplotype"])
        phase_set = int(record.get("read_phase_set", 0))
        existing = self.read_key_to_index.get(key)
        if existing is not None:
            old_hap = self.read_haplotype_by_index[existing]
            if old_hap != hap and hap != 0:
                self._replace_read_hp(existing, 0)
                self.haplotype_conflict_count += 1
            old_phase_set = self.read_phase_set_by_index[existing]
            if old_phase_set != phase_set and phase_set != 0:
                self._replace_read_ps(existing, 0)
                self.phase_set_conflict_count += 1
            continue
        if key in pending:
            old_hap = int(pending[key]["read_haplotype"])
            if old_hap != hap and hap != 0:
                pending[key]["read_haplotype"] = 0
            old_phase_set = int(pending[key].get("read_phase_set", 0))
            if old_phase_set != phase_set and phase_set != 0:
                pending[key]["read_phase_set"] = 0
            continue
        pending[key] = record
        pending_order.append(key)

    _timing_add(timing, "register_reads_dedup_conflict", time.perf_counter() - step_started)
    _timing_set_count(timing, "new_shared_reads", len(pending_order))
    if not pending_order:
        return

    old = len(self.read_lengths_by_index)
    sequences = [str(pending[key]["read_sequence"]).upper() for key in pending_order]
    step_started = time.perf_counter()
    packed, local_offsets, lengths, total_bases, invalid_counts = pack_4bit_acgtnp_chunk(
        sequences,
        READ_ALLOWED,
        old,
        "fail",
    )
    _timing_add(timing, "pack_read_sequences_4bit", time.perf_counter() - step_started)
    _timing_set_count(timing, "read_bases_packed", int(total_bases))
    if invalid_counts:
        raise ValueError(f"unexpected read invalid base counts: {invalid_counts}")
    if np.any(lengths <= 0):
        raise ValueError("empty shared read sequence")
    step_started = time.perf_counter()
    haplotypes = np.asarray([int(pending[key]["read_haplotype"]) for key in pending_order], dtype=np.uint8)
    phase_sets = np.asarray([int(pending[key].get("read_phase_set", 0)) for key in pending_order], dtype=np.int64)
    reference_starts = np.asarray([int(pending[key].get("reference_start", -1)) for key in pending_order], dtype=np.int64)
    reference_ends = np.asarray([int(pending[key].get("reference_end", -1)) for key in pending_order], dtype=np.int64)
    segment_types = np.asarray([int(pending[key].get("segment_type", 0)) for key in pending_order], dtype=np.int8)
    cigar_lists = [np.asarray(pending[key].get("cigar_ops_packed", []), dtype=np.uint32) for key in pending_order]
    cigar_lens = np.asarray([int(values.shape[0]) for values in cigar_lists], dtype=np.int32)
    cigar_offsets = self.cigar_op_count + np.concatenate(([0], np.cumsum(cigar_lens[:-1], dtype=np.int64))).astype(np.int64)
    if np.any(~np.isin(haplotypes, np.asarray([0, 1, 2], dtype=np.uint8))):
        raise ValueError("invalid HP tag value outside 0/1/2")
    if np.any(phase_sets < 0):
        raise ValueError("negative phase-set id")
    if np.any(reference_ends < reference_starts):
        raise ValueError("invalid reference interval")
    if np.any(segment_types < 0):
        raise ValueError("negative segment type")
    _timing_add(timing, "build_read_metadata_arrays", time.perf_counter() - step_started)

    step_started = time.perf_counter()
    old_base_offset = self.read_base_offset
    new_read_offsets = local_offsets[1:] + old_base_offset
    append_1d(self.read_arrays["bases_packed"], packed)
    append_1d(self.read_arrays["lengths"], lengths.astype(np.int32, copy=False))
    append_1d(self.read_arrays["offsets"], new_read_offsets)
    self.read_base_offset = old_base_offset + int(packed.shape[0])
    _timing_add_count(timing, "read_sequence_append_calls", 3)
    _timing_add(timing, "write_read_sequence_arrays", time.perf_counter() - step_started)

    step_started = time.perf_counter()
    total_new_cigar_ops = int(cigar_lens.sum())
    if total_new_cigar_ops:
        old_cigar = self.cigar_op_count
        new_cigar = old_cigar + total_new_cigar_ops
        self.cigar_ops.resize((new_cigar,))
        self.cigar_ops[old_cigar:new_cigar] = np.concatenate(cigar_lists).astype(np.uint32, copy=False)
        self.cigar_op_count = new_cigar
        _timing_add_count(timing, "cigar_array_writes", 1)
    _timing_add(timing, "write_cigar_arrays", time.perf_counter() - step_started)

    step_started = time.perf_counter()
    read_metadata_items = (
        ("hp", haplotypes),
        ("ps", phase_sets),
        ("ref_start", reference_starts),
        ("ref_end", reference_ends),
        ("segment_type", segment_types),
        ("cigar_offset", cigar_offsets),
        ("cigar_len", cigar_lens),
    )
    for name, values in read_metadata_items:
        append_1d(self.read_arrays[name], values)
    _timing_add_count(timing, "read_metadata_append_calls", len(read_metadata_items))
    _timing_add(timing, "write_read_metadata_arrays", time.perf_counter() - step_started)

    step_started = time.perf_counter()
    for offset, key in enumerate(pending_order, start=old):
        local_index = offset - old
        self.read_key_to_index[key] = offset
        self.read_haplotype_by_index.append(int(haplotypes[local_index]))
        self.read_phase_set_by_index.append(int(phase_sets[local_index]))
        self.read_lengths_by_index.append(int(lengths[local_index]))
        self.read_index.reference_ends.append(int(reference_ends[local_index]))
    self.read_hp_counts.update(int(value) for value in haplotypes)
    self.read_ps_counts.update(int(value) for value in phase_sets)
    batch_min = int(lengths.min())
    batch_max = int(lengths.max())
    self.read_length_min = batch_min if self.read_length_min is None else min(self.read_length_min, batch_min)
    self.read_length_max = batch_max if self.read_length_max is None else max(self.read_length_max, batch_max)
    self.read_length_sum += int(total_bases)
    _timing_add(timing, "update_read_indices", time.perf_counter() - step_started)

# 원본: src/integrate_wgs_rna/direct_compact_writer.py:446-479
def _cigar_ref_match_intervals(record: dict[str, Any]) -> list[tuple[int, int]]:
    """Return absolute reference/HAP intervals covered by read bases.

    CIGAR semantics for depth:
    - M/= /X consume query and reference, so these bases contribute depth.
    - D/N consume reference only, so the read has no base at those positions.
    - I/S/H/P consume no reference depth position here. S/H are normally removed
      before storage, but skipping them keeps the helper robust.
    """
    ref_pos = int(record.get("reference_start", -1))
    ref_end = int(record.get("reference_end", ref_pos))
    if ref_pos < 0 or ref_end < ref_pos:
        return []
    ops = record.get("cigar_ops_packed", [])
    if ops is None or len(ops) == 0:
        return [(ref_pos, ref_end)] if ref_end > ref_pos else []
    intervals: list[tuple[int, int]] = []
    for packed_value in ops:
        packed = int(packed_value)
        op = packed & 0xF
        length = packed >> 4
        if length <= 0:
            continue
        if op in (0, 7, 8):  # M, =, X
            intervals.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op in (2, 3):  # D, N
            ref_pos += length
        elif op in (1, 4, 5, 6):  # I, S, H, P
            continue
        else:
            # Unknown/obsolete CIGAR op: do not invent reference coverage.
            continue
    return intervals

# 원본: src/integrate_wgs_rna/direct_compact_writer.py:482-558
def _depth_arrays_for_windows(
    windows: list[Any],
    records_by_window: list[list[dict[str, Any]]],
    *,
    method: str,
    segment_types: tuple[int, ...] = (0, 1),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    """Build flat uint32 depth and per-window offsets/lengths.

    For personal_realigned data, ``method`` is CIGAR-aware and uses absolute
    HAP coordinates. For reference_projected data, exact CIGAR-to-HAP projection
    is not representable from the stored original CIGAR alone after indels in
    the HAP chain, so this writer falls back to segment span depth. v2
    production uses personal_realigned mode for depth intended for modelling.
    """
    flat_parts: list[np.ndarray] = []
    allowed_segment_types = {int(value) for value in segment_types}
    offsets = [0]
    lengths: list[int] = []
    depth_sum = 0
    depth_max = 0
    for window, records in zip(windows, records_by_window):
        # SequenceWindow.sequence is padded with "P" to the model window size,
        # while actual_length is the number of real FASTA bases before padding.
        # Depth is a model-input track and must share the same axis as the
        # stored sequence. Padding bases are stored with depth 0, while
        # CIGAR-derived coverage is accumulated only over actual_length.
        storage_length = len(str(window.sequence))
        actual_length = int(window.actual_length)
        lengths.append(storage_length)
        if storage_length <= 0:
            flat_parts.append(np.asarray([], dtype=np.uint32))
            offsets.append(offsets[-1])
            continue
        gpu_depth=os.environ.get('MFM_COVERAGE_BACKEND','cpu')=='cupy'
        starts=[];ends=[]
        diff = None if gpu_depth else np.zeros(storage_length + 1, dtype=np.int64)
        window_start = int(window.raw_start)
        actual_window_end = window_start + actual_length
        for record in records:
            if int(record.get("segment_type", 0)) not in allowed_segment_types:
                continue
            if method == "cigar_ref_match":
                intervals = _cigar_ref_match_intervals(record)
                for ref_start, ref_end in intervals:
                    start = max(int(ref_start), window_start)
                    end = min(int(ref_end), actual_window_end)
                    if end <= start:
                        continue
                    if gpu_depth:starts.append(start-window_start);ends.append(end-window_start)
                    else:
                        diff[start - window_start] += 1
                        diff[end - window_start] -= 1
            else:
                # Reference-projected fallback: uses already projected window
                # segment interval. This is kept for backwards-compatible
                # reference_projected outputs, not as the v2 preferred depth.
                start = max(int(record.get("read_start", 0)), 0)
                end = min(int(record.get("read_end", start)), actual_length)
                if end <= start:
                    continue
                if gpu_depth:starts.append(start);ends.append(end)
                else:
                    diff[start] += 1
                    diff[end] -= 1
        if gpu_depth:
            from .gpu_arrays import depth_from_intervals
            depth=depth_from_intervals(starts,ends,storage_length)
        else:
            depth = np.cumsum(diff[:-1], dtype=np.int64).astype(np.uint32, copy=False)
        if depth.size:
            depth_sum += int(depth.sum(dtype=np.uint64))
            depth_max = max(depth_max, int(depth.max()))
        flat_parts.append(depth)
        offsets.append(offsets[-1] + storage_length)
    flat = np.concatenate(flat_parts) if flat_parts else np.asarray([], dtype=np.uint32)
    return flat, np.asarray(offsets, dtype=np.int64), np.asarray(lengths, dtype=np.int32), depth_sum, depth_max
