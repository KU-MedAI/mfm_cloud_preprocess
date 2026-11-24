from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class ExportConfig:
    """Static gates for deterministic raw-read preprocessing.

    This config is intentionally limited to checks available from a single
    alignment record or a single sequence window. Sample-level percentile
    filters, such as clipping/read-quality outlier removal, need distributional
    context and should be applied by an explicit preprocessing preset or wrapper
    instead of being hidden in the low-level record filter.
    """

    window_size: int = 1 << 20
    offset_unit: int | None = None
    curriculum_unit: int = 1 << 20
    seed: int = 42
    min_mapping_quality: int = 30
    exclude_unmapped: bool = True
    exclude_secondary: bool = True
    # Keep supplementary alignments by default. For long reads, split
    # alignments can carry useful complex/repetitive-region context.
    exclude_supplementary: bool = False
    exclude_qcfail: bool = True
    # Historical F0 output preserved duplicates. The fixed training preset
    # `high_conf_read_v1` should set this to True.
    exclude_duplicate: bool = False
    # Stage 70 stores the resolved policy instead of inferring it from platform.
    filter_preset: str = "legacy_all_retained_v0"
    # segment_type: primary=0, supplementary=1. Explicit tracks prevent split
    # read evidence from silently inflating model coverage.
    depth_segment_types: tuple[int, ...] = (0, 1)
    supplementary_depth_segment_types: tuple[int, ...] = ()
    model_segment_types: tuple[int, ...] = (0, 1)
    require_proper_pair: bool = False
    max_n_ratio: float = 1.0
    window_storage_policy: str = "annotated_acgtn_v1"
    # Coordinate-coverage gate REMOVED (jhleee, 2026-08-17). 0.0 keeps every
    # sequence-valid window; the coverage value itself is still measured and
    # recorded in QC (`coordinate_coverage_min/max`, `..._threshold`).
    #
    # Why: requiring a window to map back to the reference is a *target
    # alignment* requirement. It only means something when the objective has a
    # reference-coordinate target (post-training's 1bp RNA track). Pre-training
    # learns a representation from sequence (MLM) and predicts no such track, so
    # reference correspondence is unrelated to its objective. Worse, the windows
    # this gate drops are exactly the insertion-rich ones — i.e. where the
    # individual genome differs most from the reference — so it discards the most
    # individual-specific sequence. The individual sequence IS the point.
    #
    # Measured 2026-08-17 (data/mfm_integrated/final/*/long_read/*/metadata/
    # preprocess.json): 0 windows had ever been dropped by this gate —
    # pre_training GIAB 264 hap x chrom min coverage 0.98407, post_training
    # geuvadis min 0.99692. It was non-binding for germline samples, so this
    # change reproduces existing outputs; it prevents the gate from firing on
    # cancer / SV-rich / MHC-segdup data later.
    #
    # NOTE for the read channel: windows kept with low coverage have little or no
    # reference->individual projection, so the projection-based read path leaves
    # them with sparse/empty pileup. The fix is to realign reads onto the
    # individual genome (post-training repo: process_code/06_realign_personal.sh)
    # rather than to reinstate this gate.
    min_coordinate_coverage: float = 0.0

    def effective_offset_unit(self) -> int:
        """Return the coordinate span used for deterministic random offsets.

        The default samples seed-based random starts across the full model
        window span. With the standard 1,048,576 bp window, this means
        random_offset is in [0, 1,048,576). Set `offset_unit=1` when an
        experiment intentionally needs fixed-start windows.
        """
        return self.window_size if self.offset_unit is None else int(self.offset_unit)

    def validate(self) -> None:
        if self.window_size <= 0 or self.curriculum_unit <= 0:
            raise ValueError("window sizes must be positive")
        offset_unit = self.effective_offset_unit()
        if offset_unit <= 0:
            raise ValueError("offset_unit must be positive")
        if offset_unit > self.window_size:
            raise ValueError("offset_unit cannot exceed window_size")
        if not 0.0 <= self.max_n_ratio <= 1.0:
            raise ValueError("max_n_ratio must be within [0, 1]")
        if not 0.0 <= self.min_coordinate_coverage <= 1.0:
            raise ValueError("min_coordinate_coverage must be within [0, 1]")
        for name, values, allow_empty in (
            ("depth_segment_types", self.depth_segment_types, False),
            ("supplementary_depth_segment_types", self.supplementary_depth_segment_types, True),
            ("model_segment_types", self.model_segment_types, False),
        ):
            normalized = tuple(int(value) for value in values)
            if (not allow_empty and not normalized) or any(value not in {0, 1} for value in normalized):
                raise ValueError(f"{name} must contain only primary(0)/supplementary(1) segment types")
        if set(self.depth_segment_types) & set(self.supplementary_depth_segment_types):
            raise ValueError("primary and supplementary depth segment types must not overlap")

    def to_dict(self) -> dict:
        return asdict(self)
