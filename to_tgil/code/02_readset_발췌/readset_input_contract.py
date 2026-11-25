# 검토용 원문 발췌임. 단독 실행 모듈이 아니며 import/호출 문맥은 전체 전달본을 사용함.

# 원본: src/integrate_wgs_rna/full_preprocess.py:155-183
def validate_read_coordinate_mode(config: FullPreprocessConfig) -> None:
    """Validate how read coordinates will be represented in the final compact output.

    reference_projected keeps the historical layout: one shared readset.zarr stores
    reads fetched from the original reference-coordinate alignment and segment
    coordinates are projected to each HAP through chain files.

    personal_realigned stores one readset per HAP.  Reads must already be aligned
    against the corresponding HAP FASTA, so read starts/ends/CIGAR are interpreted
    directly in that HAP coordinate system.  This mode is intentionally restricted
    to direct compact output because the legacy two-step raw layout has only one
    shared readset path and cannot represent HAP-specific alignments cleanly.
    """
    mode = config.read_coordinate_mode or "reference_projected"
    if mode not in {"reference_projected", "personal_realigned"}:
        raise ValueError(f"unsupported read_coordinate_mode: {mode}")
    config.read_coordinate_mode = mode
    if mode != "personal_realigned":
        return
    if config.compact_mode != "direct":
        raise ValueError("--read-coordinate-mode personal_realigned requires --compact-mode direct")
    missing = []
    for label, path in (("--hap1-alignment", config.hap1_alignment), ("--hap2-alignment", config.hap2_alignment)):
        if path is None:
            missing.append(label)
        elif not path.exists():
            raise FileNotFoundError(f"missing {label} personal realigned alignment: {path}")
    if missing:
        raise ValueError("--read-coordinate-mode personal_realigned requires " + ", ".join(missing))
