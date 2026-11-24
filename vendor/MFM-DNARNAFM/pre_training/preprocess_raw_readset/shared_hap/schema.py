"""Shared-HAP output schema constants.

Keep schema names and storage conventions in one module so the writer, CLI,
wrapper, validation, and documentation do not drift apart.
"""

from __future__ import annotations

SHARED_SCHEMA = "shared_hap_model_input_by_chr"
SCHEMA_VERSION = "hp_ps_cigar_integrated_2026_08_10"
MANIFEST_SCHEMA = "shared_hap_model_input_manifest"
OUTPUT_LAYOUT = "caller_defined_no_v1_v2_path_split"

# Existing production data may still carry a historical schema attr after the
# PS field was migrated into the original v1 folders. Validation accepts these
# attrs as long as required HP/PS arrays are present and structurally valid.
COMPATIBLE_SHARED_SCHEMAS = {
    SHARED_SCHEMA,
    "shared_hap_model_input_by_chr_v1",
    "shared_hap_model_input_by_chr_v2",
}

HAP_LABELS = ("HAP1", "HAP2")
READ_CHUNK_SIZE = 4096
SEGMENT_CHUNK_SIZE = 8192
VALIDATION_LEVELS = ("full", "standard", "fast")

READ_HAPLOTYPE_STORAGE = "BAM_HP_tag_int8_zero_for_missing_or_nonstandard"
READ_PHASE_SET_STORAGE = "BAM_PS_tag_int64_zero_for_missing_or_conflict"
READ_SEQUENCE_STORAGE = "shared_clipping_removed_mapped_query_sequence"
READ_CIGAR_STORAGE = "uint32_packed_htslib_cigar_ops_clipping_removed_length_shift_4_or_op"
READ_REFERENCE_STORAGE = "bam_reference_0based_halfopen_per_shared_alignment_segment"
