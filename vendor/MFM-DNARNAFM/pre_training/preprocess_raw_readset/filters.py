from __future__ import annotations


def rejection_reason(read, config) -> str | None:
    """Return the deterministic record-level rejection reason.

    Keep this function conservative: do not reject reads by reference mismatch,
    NM/aligned length, minor fraction, or custom discordance scores here. Those
    signals can overlap true individual variation and require cohort/window
    context. Percentile-based clipping/read-quality outlier filters should be
    handled by an explicit higher-level preset.
    """

    if config.exclude_unmapped and read.is_unmapped:
        return "unmapped"
    if config.exclude_secondary and read.is_secondary:
        return "secondary"
    if config.exclude_supplementary and read.is_supplementary:
        return "supplementary"
    if config.exclude_qcfail and read.is_qcfail:
        return "qcfail"
    if config.exclude_duplicate and read.is_duplicate:
        return "duplicate"
    if config.require_proper_pair and read.is_paired and not read.is_proper_pair:
        return "improper_pair"
    if int(read.mapping_quality) < config.min_mapping_quality:
        return "low_mapq"
    # query_length는 SEQ 문자열을 디코딩·캐시하지 않고 부재를 판정함.
    if not read.query_length:
        return "missing_sequence"
    return None
