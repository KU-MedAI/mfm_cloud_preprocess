"""개인 HAP 재정렬의 config와 platform/routing policy."""

from .config import PersonalRealignConfig, PreparedPairedInput
from .policy import (
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

__all__ = [
    "PersonalRealignConfig",
    "PreparedPairedInput",
    "alignment_contract",
    "effective_hap_read_mode",
    "is_hap_specific_mode",
    "is_illumina_platform",
    "normalize_platform",
    "platform_minimap2_preset",
    "resolve_hap_read_mode",
    "sam_rg_platform",
    "selected_aligner",
]
