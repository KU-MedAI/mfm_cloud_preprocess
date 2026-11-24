from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

DEFAULT_SAMTOOLS = "samtools"
DEFAULT_MINIMAP2 = "minimap2"
DEFAULT_BWA_MEM2 = "bwa-mem2"
SCHEMA = "mfm_integrated_personal_realign_v2"
READ_ROUTING_POLICIES = {"hp_routed", "all_to_both"}
UNKNOWN_READ_POLICIES = {"include_both", "drop_unknown"}
HAP_READ_MODES = {"all_to_both", "hp_routed_with_unknown", "hap_specific"}
HAP_READ_MODE_SETTINGS = {
    "all_to_both": ("all_to_both", "include_both", True),
    "hp_routed_with_unknown": ("hp_routed", "include_both", True),
    "hap_specific": ("hp_routed", "drop_unknown", False),
}
ILLUMINA_PAIRED_INPUT_MODE = "illumina_name_collated_paired"
ILLUMINA_HAP_SPECIFIC_INPUT_MODE = (
    "illumina_name_collated_hap_specific_fragment"
)
MIN_ILLUMINA_PAIRED_PRIMARY_FRAC = 0.90
ILLUMINA_ALIGNERS = {"bwa-mem2", "minimap2"}
ILLUMINA_DUPLICATE_MODES = {"preserve_source", "none"}
RATIO_QC_POLICIES = {"record_only", "strict"}


@dataclass
class PersonalRealignConfig:
    cohort: str
    output_root: Path
    sample_list: Path | None
    samples: list[str]
    chroms: list[str]
    sequence_variant_policy: str
    dataset_role: str
    sequence_source: str = "vcf_consensus"
    assembly_phase_labels_validated: bool = False
    alignment: Path | None = None
    alignment_pattern: str = (
        "{sample}/haplotagged/{sample}.ont.hg38.haplotagged.bam"
    )
    hap1_fasta_dir: Path | None = None
    hap2_fasta_dir: Path | None = None
    output_prefix: str = "{sample}.personal.HAP{hap}"
    realigned_subdir: str = "realigned"
    reference_fasta: Path | None = None
    output_format: str = "cram"
    platform: str = ""
    preset: str = "auto"
    threads: int = 24
    sort_threads: int = 8
    input_decode_threads: int = 0
    fastq_threads: int = 2
    sort_mem: str = "2G"
    jobs: int = 1
    hap_jobs: int = 1
    route_jobs: int = 1
    min_mapped_frac: float = 0.90
    ratio_qc_policy: str = "record_only"
    hap_read_mode: str = ""
    read_routing_policy: str = "hp_routed"
    unknown_read_policy: str = "include_both"
    include_unmapped: bool = True
    verify_counts: bool = True
    count_tolerance_abs: int = 10
    keep_target: bool = False
    allow_non_pretraining: bool = False
    samtools_bin: str = DEFAULT_SAMTOOLS
    minimap2_bin: str = DEFAULT_MINIMAP2
    bwa_mem2_bin: str = DEFAULT_BWA_MEM2
    illumina_aligner: str = "bwa-mem2"
    bwa_mem2_index_root: Path | None = None
    illumina_duplicate_mode: str = "preserve_source"
    inspect_reads: int = 5000
    overwrite: bool = False
    dry_run: bool = False
    summary_tsv: Path | None = None
    metadata_json: Path | None = None


@dataclass(frozen=True)
class PreparedPairedInput:
    """v2 WGS Stage 10용 Illumina 입력: 라우팅 후 FASTQ 변환 전 자료."""

    path: Path
    route_label: str
    stats: dict[str, Any]
    input_counts: dict[str, int] | None = None
    name_collated: bool = False
