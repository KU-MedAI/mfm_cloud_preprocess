# 검토용 원문 발췌임. 단독 실행 모듈이 아니며 import/호출 문맥은 전체 전달본을 사용함.

# 원본: src/integrate_wgs_rna/direct_compact_writer.py:1403-1449
self.readset_root.attrs.update(
    {
        "schema": "mfm_readset_4bit_acgtnp_direct_v0",
        "encoding": "4bit_acgtnp",
        "symbols": format_symbols(READ_ALLOWED),
        "sample": sample,
        "contig": contig,
        "window_size": self.window_size,
        "read_sequence_storage": "shared_clipping_removed_mapped_query_sequence",
        "read_haplotype_storage": "BAM_HP_tag_uint8_zero_for_missing_or_conflict",
        "read_phase_set_storage": "BAM_PS_tag_int64_zero_for_missing_or_conflict",
        "read_cigar_storage": "uint32_packed_htslib_cigar_ops_clipping_removed_length_shift_4_or_op",
        "read_reference_storage": "bam_reference_0based_halfopen_per_shared_alignment_segment",
        "created_at": utc_now(),
        "compressor": self.compact_compressor,
        "compression_level": self.compact_compression_level,
        **(attrs or {}),
    }
)
self.readset_root.attrs["model_segment_types"] = list(self.model_segment_types)
self.readset_root.attrs["model_segment_policy"] = "primary_only" if self.model_segment_types == (0,) else "explicit_segment_types"

reads = self.readset_root.create_group("reads")
reads.attrs.update(
    {
        "encoding": "4bit_acgtnp",
        "symbols": format_symbols(READ_ALLOWED),
        "kind": "read_sequence",
        "invalid_base_policy": "fail",
        "bases_per_byte": 2,
    }
)
self.read_arrays = {
    "bases_packed": reads.empty("bases_packed", shape=(0,), dtype=np.uint8, chunks=(PACKED_BASE_CHUNK,), compressor=comp),
    "offsets": reads.empty("offsets", shape=(1,), dtype=np.int64, chunks=(READ_META_CHUNK,), compressor=comp),
    "lengths": reads.empty("lengths", shape=(0,), dtype=np.int32, chunks=(READ_META_CHUNK,), compressor=comp),
    "hp": reads.empty("hp", shape=(0,), dtype=np.uint8, chunks=(READ_META_CHUNK,), compressor=comp),
    "ps": reads.empty("ps", shape=(0,), dtype=np.int64, chunks=(READ_META_CHUNK,), compressor=comp),
    "ref_start": reads.empty("ref_start", shape=(0,), dtype=np.int64, chunks=(READ_META_CHUNK,), compressor=comp),
    "ref_end": reads.empty("ref_end", shape=(0,), dtype=np.int64, chunks=(READ_META_CHUNK,), compressor=comp),
    "segment_type": reads.empty("segment_type", shape=(0,), dtype=np.int8, chunks=(READ_META_CHUNK,), compressor=comp),
    "cigar_offset": reads.empty("cigar_offset", shape=(0,), dtype=np.int64, chunks=(READ_META_CHUNK,), compressor=comp),
    "cigar_len": reads.empty("cigar_len", shape=(0,), dtype=np.int32, chunks=(READ_META_CHUNK,), compressor=comp),
}
self.read_arrays["offsets"][0] = 0
cigar = self.readset_root.create_group("cigar")
self.cigar_ops = cigar.empty("ops_packed", shape=(0,), dtype=np.uint32, chunks=(SEGMENT_CHUNK,), compressor=comp)
