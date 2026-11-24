"""Stage 12 compact Zarr storage constants.

These values define storage layout only. Changing them can alter artifact bytes,
chunk counts, and performance, so defaults are versioned and regression-tested.
"""

READ_META_CHUNK = 524_288
SEGMENT_CHUNK = 1_048_576
PACKED_BASE_CHUNK = 32 << 20
WINDOW_META_CHUNK = 1024
DEPTH_CHUNK = 1 << 20

FULL_DIRECT_CONTIG_SCHEMA = "mfm_integrated_direct_shared_hap_4bit_acgtnp_depth_contig_v1"
FULL_DIRECT_PERSONAL_CONTIG_SCHEMA = "mfm_integrated_direct_personal_realigned_hap_split_4bit_acgtnp_depth_contig_v1"
PERSONAL_INPUT_PROVENANCE_SCHEMA = "stage70_personal_input_provenance_v1"
PERSONAL_FULL_VALIDATION_SCHEMA = "stage70_personal_full_validation_v1"

DEFAULT_COMPACT_COMPRESSOR = "zstd"
DEFAULT_COMPACT_COMPRESSION_LEVEL = 5
SUPPORTED_COMPACT_COMPRESSORS = {"zstd", "lz4", "blosclz"}
