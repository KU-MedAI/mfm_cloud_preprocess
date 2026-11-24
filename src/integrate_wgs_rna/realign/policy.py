from __future__ import annotations

from typing import Any
import os

from .config import (
    HAP_READ_MODES,
    HAP_READ_MODE_SETTINGS,
    PersonalRealignConfig,
)


def normalize_platform(value: str) -> str:
    return (value or "").strip().lower().replace("-", "_").replace(" ", "_")


def platform_minimap2_preset(
    platform: str, requested_preset: str
) -> str:
    """sequencing platform에 맞는 minimap2 preset을 결정한다."""
    requested = (requested_preset or "auto").strip()
    if requested and requested != "auto":
        return requested
    normalized = normalize_platform(platform)
    if (
        normalized
        in {
            "ont",
            "nanopore",
            "oxfordnanopore",
            "oxford_nanopore",
            "1kg_ont_vienna",
        }
        or normalized.startswith("ont")
        or "nanopore" in normalized
    ):
        return "map-ont"
    if (
        normalized
        in {
            "hifi",
            "pacbio_hifi",
            "pacbiohifi",
            "pb_hifi",
            "ccs",
            "pacbio_ccs",
            "pacbioccs",
        }
        or "hifi" in normalized
        or "ccs" in normalized
    ):
        return "map-hifi"
    if (
        normalized in {"clr", "pacbio_clr", "pacbioclr", "pb_clr"}
        or "_clr" in normalized
    ):
        return "map-pb"
    if (
        normalized in {"illumina", "short_read", "shortread", "ilmn", "sr"}
        or "illumina" in normalized
        or normalized.startswith("short")
    ):
        return "sr"
    raise ValueError(
        "--preset auto requires --platform ont|hifi|clr|illumina, "
        "or pass --preset explicitly"
    )


def sam_rg_platform(platform: str) -> str:
    normalized = normalize_platform(platform)
    if (
        normalized
        in {
            "ont",
            "nanopore",
            "oxfordnanopore",
            "oxford_nanopore",
            "1kg_ont_vienna",
        }
        or normalized.startswith("ont")
        or "nanopore" in normalized
    ):
        return "ONT"
    if (
        normalized
        in {
            "hifi",
            "pacbio_hifi",
            "pacbiohifi",
            "pb_hifi",
            "ccs",
            "pacbio_ccs",
            "pacbioccs",
            "clr",
            "pacbio_clr",
            "pacbioclr",
            "pb_clr",
        }
        or "hifi" in normalized
        or "ccs" in normalized
        or normalized.startswith("pacbio")
    ):
        return "PACBIO"
    if (
        normalized in {"illumina", "short_read", "shortread", "ilmn", "sr"}
        or "illumina" in normalized
        or normalized.startswith("short")
    ):
        return "ILLUMINA"
    return "UNKNOWN"


def resolve_hap_read_mode(config: PersonalRealignConfig) -> str:
    """상위 HAP mode를 기존 routing primitive에 연결한다."""
    mode = (config.hap_read_mode or "").strip()
    if not mode:
        return "legacy"
    if mode not in HAP_READ_MODES:
        raise ValueError(
            f"--hap-read-mode must be one of: "
            f"{', '.join(sorted(HAP_READ_MODES))}"
        )
    policy, unknown, include_unmapped = HAP_READ_MODE_SETTINGS[mode]
    config.read_routing_policy = policy
    config.unknown_read_policy = unknown
    config.include_unmapped = include_unmapped
    return mode


def effective_hap_read_mode(config: PersonalRealignConfig) -> str:
    return (config.hap_read_mode or "").strip() or "legacy"


def is_hap_specific_mode(config: PersonalRealignConfig) -> bool:
    return effective_hap_read_mode(config) == "hap_specific"


def is_illumina_platform(platform: str) -> bool:
    return sam_rg_platform(platform) == "ILLUMINA"


def selected_aligner(config: PersonalRealignConfig) -> str:
    return (
        config.illumina_aligner
        if is_illumina_platform(config.platform)
        else "minimap2"
    )


def alignment_contract(config: PersonalRealignConfig) -> dict[str, Any]:
    aligner = selected_aligner(config)
    illumina = is_illumina_platform(config.platform)
    gpu=os.environ.get('MFM_WGS_BACKEND','cpu')=='parabricks'
    return {
        **({'execution_backend':'parabricks','cpu_baseline_aligner':aligner} if gpu else {}),
        "aligner": ('parabricks-bwa-mem' if illumina else 'parabricks-minimap2') if gpu else aligner,
        "args": (
            ('fq2bam --bwa-options=-Y -C --no-markdups' if illumina else f'minimap2 --preset {config.preset} --md --copy-comment --seeding-chaining-on-gpu on')
            if gpu else "mem -p -Y -C"
            if aligner == "bwa-mem2"
            else f"-ax {config.preset} -y -Y --MD"
        ),
        "sequence_source": config.sequence_source,
        "assembly_phase_labels_validated": bool(config.assembly_phase_labels_validated),
        "illumina_duplicate_mode": (
            config.illumina_duplicate_mode if illumina else "not_applicable"
        ),
        "paired_fastq_policy": (
            "complete_R1_R2_only; singleton_and_read_other_discarded"
            if illumina
            else "not_applicable"
        ),
    }
