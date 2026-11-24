from __future__ import annotations

import csv
import json
import os
import hashlib
import re
import shlex
import shutil
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pysam

from pre_training.preprocess_raw_readset.config import ExportConfig
from pre_training.preprocess_raw_readset.shared_hap_cli import _job_schedule_key, _run_contig
from pre_training.preprocess_raw_readset.shared_hap.manifest import (
    init_manifest,
    mark_manifest_complete,
    mark_manifest_failed,
    update_manifest_record,
)
from .direct_compact_writer import (
    DEFAULT_COMPACT_COMPRESSOR,
    DEFAULT_COMPACT_COMPRESSION_LEVEL,
    run_direct_compact_jobs,
    summarize_direct_compact_outputs,
)
from .legacy_shared_hap_converter import LegacyConvertConfig, directory_size, human_bytes, run_legacy_convert
from .pipeline.manifest_progress import ManifestProgressJournal, manifest_progress_signature
from .pipeline.catalog import (
    FULL_LAYOUT_NAME,
    catalog_row_from_result,
    detect_rna_personal,
    read_catalog,
    read_optional_json,
    upsert_sample_catalog,
    write_tsv,
)
from .timing import StepTimer
from .personal_fasta import reference_fingerprint
from .window_policy import persist_window_annotations
from .compact.stages import compact_step, annotation_step


RECORD_COLUMNS = ("sample", "contig", "hap1_fasta", "hap1_chain", "hap2_fasta", "hap2_chain", "bam", "window_index_dir", "read_coordinate_mode", "hap1_alignment", "hap2_alignment")
DEFAULT_BCFTOOLS = "bcftools"
DEFAULT_SAMPLE_VCF_FILTER_EXPR = 'GT="het"'
DEFAULT_CONSENSUS_FILTER_EXPR = '(FILTER="PASS" || FILTER=".")'
DEFAULT_HAPLOTAG_VCF_TYPES = "snps,indels"
DEFAULT_SEQUENCE_VARIANT_POLICY = "snv_indel"
PROJECT_SEQUENCE_POLICY = {
    "igsr_ont_wgs_only": "snv_indel",
    "igsr_ont_geuvadis": "snv_indel",
}
PROJECT_ROLE_POLICY = {
    "igsr_ont_wgs_only": "pre_training",
    "igsr_ont_geuvadis": "post_training",
}
PROJECT_READ_TYPE_POLICY = {
    "igsr_ont_wgs_only": "long_read",
    "igsr_ont_geuvadis": "long_read",
}
LONG_READ_PLATFORMS = {"ont", "nanopore", "oxfordnanopore", "oxford_nanopore", "hifi", "pacbio_hifi", "pacbiohifi", "ccs", "pacbio_ccs", "pacbioccs"}

COMMON_VCF_PATTERN_CANDIDATES = (
    # process_code policy: default common panel is the 2022 1kGP 3202-sample
    # phased SNV/INDEL/SV VCF. Older 2020 CCDG panels must be requested via
    # --common-vcf-pattern explicitly so we do not silently mix VCF releases.
    "1kGP_high_coverage_Illumina.{contig}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
)
CONSENSUS_ALPHABET = frozenset("ACGTN")


@dataclass
class FullPreprocessConfig:
    sample_id: str
    dataset_role: str
    read_type: str
    project_id: str
    hap1_fasta_dir: Path | None
    hap1_chain_dir: Path | None
    hap2_fasta_dir: Path | None
    hap2_chain_dir: Path | None
    alignment: Path
    contigs: list[str]
    output_sample_root: Path | None = None
    dataset_root: Path | None = None
    work_root: Path | None = None
    window_index_dir: Path | None = None
    read_coordinate_mode: str = "personal_realigned"
    hap1_alignment: Path | None = None
    hap2_alignment: Path | None = None
    layout_name: str = FULL_LAYOUT_NAME
    hap_read_mode: str = "all_to_both"
    source_sample_id: str = ""
    project_label: str = ""
    platform: str = ""
    reference_build: str = "hg38"
    sample_status: str = "healthy"
    disease: str = "NaN"
    storage_policy: str = "model_input_only"
    jobs: int = 1
    windows_per_chunk: int = 32
    compact_mode: str = "direct"
    validation_level: str = "full"
    compact_compressor: str = DEFAULT_COMPACT_COMPRESSOR
    compact_compression_level: int = DEFAULT_COMPACT_COMPRESSION_LEVEL
    profile_direct_compact: bool = False
    schedule: str = "longest-first"
    window_size: int = 1 << 20
    offset_unit: int | None = None
    seed: int = 42
    min_mapq: int = 30
    stage70_filter_preset: str = "auto"
    max_windows: int | None = None
    overwrite: bool = False
    resume: bool = True
    allow_untagged: bool = False
    haplotag_if_missing: bool = False
    allow_longread_vcf_haplotag: bool = False
    allow_unphased_longread: bool = False
    auto_prepare_hap_fasta: bool = False
    phased_vcf: Path | None = None
    common_vcf_dir: Path | None = None
    common_vcf_pattern: str = ""
    common_vcf_filter_expr: str = DEFAULT_SAMPLE_VCF_FILTER_EXPR
    consensus_filter_expr: str = DEFAULT_CONSENSUS_FILTER_EXPR
    haplotag_vcf_types: str = DEFAULT_HAPLOTAG_VCF_TYPES
    sequence_variant_policy: str = DEFAULT_SEQUENCE_VARIANT_POLICY
    reference_fasta: Path | None = None
    sequence_reference_fasta: Path | None = None
    personal_fasta_qc: Path | None = None
    phase_inspect_max_reads: int = 50_000
    synthesize_missing_ps: bool = False
    whatshap_bin: str = "whatshap"
    samtools_bin: str = "samtools"
    bcftools_bin: str = DEFAULT_BCFTOOLS
    source_wgs: list[Path] = field(default_factory=list)
    source_vcfs: list[Path] = field(default_factory=list)
    rna_bigwigs: list[Path] = field(default_factory=list)


@dataclass
class FullManifestConfig:
    manifest: Path
    status_tsv: Path
    dataset_root: Path | None = None
    sample_jobs: int = 1
    jobs: int = 1
    overwrite: bool = False
    resume: bool = True
    storage_policy: str = "model_input_only"
    sample_catalog: Path | None = None


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


# These presets encode the Stage 70 representation contract. They do not
# create long-read duplicate marks or turn assembly/graph policies into BAM flags.
STAGE70_FILTER_PRESET_CHOICES = (
    "auto",
    "short_read_base_v1",
    "long_read_primary_supp_evidence_v1",
)


def resolve_stage70_filter_policy(read_type: str, requested_preset: str) -> dict[str, Any]:
    preset = (requested_preset or "auto").strip()
    if preset == "auto":
        preset = "short_read_base_v1" if read_type == "short_read" else "long_read_primary_supp_evidence_v1"
    if preset == "short_read_base_v1":
        policy = {
            "read_type": "short_read", "filter_preset": preset,
            "exclude_supplementary": True, "exclude_duplicate": True,
            "depth_segment_types": (0,), "supplementary_depth_segment_types": (), "model_segment_types": (0,),
        }
    elif preset == "long_read_primary_supp_evidence_v1":
        policy = {
            "read_type": "long_read", "filter_preset": preset,
            "exclude_supplementary": False, "exclude_duplicate": False,
            "depth_segment_types": (0,), "supplementary_depth_segment_types": (1,), "model_segment_types": (0,),
        }
    else:
        raise ValueError(f"unknown Stage 70 filter preset: {requested_preset}")
    if read_type != policy["read_type"]:
        raise ValueError(f"Stage 70 filter preset {preset} is incompatible with read_type={read_type}")
    return policy


def stage70_export_config(config: FullPreprocessConfig) -> ExportConfig:
    policy = resolve_stage70_filter_policy(config.read_type, config.stage70_filter_preset)
    return ExportConfig(
        window_size=config.window_size,
        offset_unit=config.offset_unit,
        seed=config.seed,
        min_mapping_quality=config.min_mapq,
        exclude_supplementary=bool(policy["exclude_supplementary"]),
        exclude_duplicate=bool(policy["exclude_duplicate"]),
        filter_preset=str(policy["filter_preset"]),
        depth_segment_types=tuple(policy["depth_segment_types"]),
        supplementary_depth_segment_types=tuple(policy["supplementary_depth_segment_types"]),
        model_segment_types=tuple(policy["model_segment_types"]),
    )


def _replay_cli_commands(
    argv: list[str] | None = None,
    *,
    python_executable: str | None = None,
    module_root: Path | None = None,
) -> list[str]:
    """Build a standalone replay command for metadata/commands.sh."""
    args = list(sys.argv if argv is None else argv)
    if not args:
        raise ValueError("argv must contain the CLI program name")
    root = module_root or Path(__file__).resolve().parents[1]
    executable = python_executable or sys.executable
    return [
        f'export PYTHONPATH={shlex.quote(str(root))}:"${{PYTHONPATH:-}}"',
        shell_join([str(executable), "-m", "integrate_wgs_rna.cli", *args[1:]]),
    ]



def validate_personal_fasta_reference_contract(
    config: FullPreprocessConfig, contigs: list[str]
) -> dict[str, Any] | None:
    """Fail closed before Stage 70 when personal FASTA and sequence reference diverge."""
    if config.read_coordinate_mode != "personal_realigned":
        return None
    reference = effective_sequence_reference(config)
    if reference is None or not reference.is_file() or not Path(str(reference) + ".fai").is_file():
        raise ValueError("personal_realigned requires --sequence-reference-fasta and its .fai")
    qc_path = config.personal_fasta_qc
    if qc_path is None and config.hap1_fasta_dir is not None:
        qc_path = config.hap1_fasta_dir.parent / "personal_fasta_qc.json"
    if qc_path is None or not qc_path.is_file():
        raise FileNotFoundError(f"personal FASTA QC가 없습니다: {qc_path}")
    qc = json.loads(qc_path.read_text())
    if qc.get("schema") != "mfm_personal_fasta_sample_qc_v1" or qc.get("status") != "pass":
        raise ValueError(f"personal FASTA QC schema/status 오류: {qc_path}")
    if qc.get("sample") != config.sample_id:
        raise ValueError(f"personal FASTA QC sample 불일치: {qc.get('sample')} != {config.sample_id}")
    expected_fingerprint = reference_fingerprint(reference)
    if qc.get("reference_fingerprint") != expected_fingerprint:
        raise ValueError(
            "personal FASTA reference fingerprint 불일치: "
            f"qc={qc.get('reference_fingerprint')} current={expected_fingerprint}"
        )
    if Path(str(qc.get("reference_fasta", ""))).resolve() != reference.resolve():
        raise ValueError("personal FASTA QC reference path와 --sequence-reference-fasta가 다릅니다")
    hap_results = qc.get("hap_results", {})
    require_chain = config.sequence_variant_policy == "snv_indel"
    for hap in ("HAP1", "HAP2"):
        rows = {str(row.get("contig")): row for row in hap_results.get(hap, [])}
        missing = [contig for contig in contigs if contig not in rows]
        if missing:
            raise ValueError(f"personal FASTA QC {hap} contig 누락: {missing}")
        if require_chain:
            bad = [contig for contig in contigs if rows[contig].get("chain_qc", {}).get("status") != "pass"]
            if bad:
                raise ValueError(f"personal FASTA QC {hap} chain 검증 실패: {bad}")
    return {
        "status": "pass",
        "qc_path": str(qc_path.resolve()),
        "reference_fingerprint": expected_fingerprint,
        "validated_contigs": contigs,
        "chain_required": require_chain,
    }

def run_full_preprocess(config: FullPreprocessConfig) -> dict[str, Any]:
    timer = StepTimer()
    layout_path = Path(config.layout_name)
    if layout_path.is_absolute() or ".." in layout_path.parts or len(layout_path.parts) != 1:
        raise ValueError("layout_name must be a single safe directory name")
    if config.hap_read_mode not in {"all_to_both", "hp_routed_with_unknown", "hap_specific"}:
        raise ValueError("invalid hap_read_mode")
    enforce_project_policy(config)
    sample_root = resolve_sample_root(config)
    work_root = config.work_root or (sample_root / "work")
    raw_output_root = work_root / "full_shared_hap_raw"
    compact_root = sample_root / "compact" / config.layout_name
    metadata_dir = sample_root / "metadata"
    logs_dir = sample_root / "logs"
    qc_dir = sample_root / "qc"
    commands: list[str] = _replay_cli_commands()
    hap_fasta_summary: dict[str, Any] | None = None
    sample_vcf_for_haplotagging: Path | None = None
    cleanup_removed: list[str] = []
    personal_fasta_contract: dict[str, Any] | None = None

    with timer.step("prepare_sample_root"):
        prepare_sample_root(sample_root)
        if config.compact_mode == "two_step":
            (raw_output_root / config.sample_id).mkdir(parents=True, exist_ok=True)

    with timer.step("resolve_inputs"):
        contigs = expand_contigs(config.contigs)
        if config.alignment.suffix.lower() == ".cram" and config.reference_fasta is None:
            raise ValueError("CRAM input requires --reference-fasta to avoid remote reference lookup")
        validate_read_coordinate_mode(config)
        personal_fasta_contract = validate_personal_fasta_reference_contract(config, contigs)

    if config.auto_prepare_hap_fasta:
        with timer.step("prepare_hap_fasta_from_common_vcf"):
            hap_fasta_summary = prepare_common_vcf_hap_fastas(config, contigs, work_root, commands)

    with timer.step("build_input_records"):
        ensure_hap_input_dirs(config)
        records = build_records(config, contigs)
        record_tsv = metadata_dir / "full_shared_hap.records.tsv"
        write_record_tsv(record_tsv, records)

    with timer.step("inspect_alignment_phase_tags"):
        phase_summary = inspect_alignment_phase_tags(config.alignment, contigs, config.reference_fasta, config.phase_inspect_max_reads)

    alignment_for_preprocess = config.alignment
    if phase_summary["hp_present_count"] == 0 and is_long_read_config(config) and not config.allow_unphased_longread:
        raise ValueError(
            "long-read alignment has no HP tags in inspected reads. v2 production requires a validated public "
            "or explicit Clair3 phased VCF plus WhatsHap HP/PS haplotagging before compact export. Use --allow-unphased-longread only for unphased ablation."
        )

    unphased_allowed = config.allow_untagged or (is_long_read_config(config) and config.allow_unphased_longread)
    if phase_summary["hp_present_count"] == 0 and not unphased_allowed:
        if not config.haplotag_if_missing:
            raise ValueError(
                "alignment has no HP tags in inspected reads. Use --haplotag-if-missing with "
                "--phased-vcf/--reference-fasta, or pass --allow-untagged intentionally."
            )
        if config.phased_vcf is None and config.common_vcf_dir is not None:
            with timer.step("prepare_sample_vcf_from_common_vcf"):
                sample_vcf_for_haplotagging = prepare_sample_phased_vcf(config, work_root, commands)
        with timer.step("whatshap_haplotag"):
            alignment_for_preprocess = run_haplotag(config, work_root, commands, sample_vcf_for_haplotagging)
        with timer.step("inspect_haplotagged_alignment"):
            phase_summary = inspect_alignment_phase_tags(alignment_for_preprocess, contigs, config.reference_fasta, config.phase_inspect_max_reads)
            if phase_summary["hp_present_count"] == 0:
                raise RuntimeError("whatshap haplotag completed but HP tags are still absent in inspected reads")

    if alignment_for_preprocess != config.alignment:
        records = [dict(row, bam=str(alignment_for_preprocess)) for row in records]
        write_record_tsv(record_tsv, records)

    export_config = stage70_export_config(config)
    if config.compact_mode == "direct":
        with timer.step("run_direct_compact_4bit"):
            shared_manifest = compact_step(sample_root, compact_root, records,
                {"export": export_config.to_dict(), "max_windows": config.max_windows,
                 "validation_level": config.validation_level, "compressor": config.compact_compressor,
                 "compression_level": config.compact_compression_level, "force_nonce": time.time_ns() if config.overwrite or not config.resume else None},
                lambda: run_direct_compact_jobs(
                records=records,
                record_tsv=record_tsv,
                output_root=compact_root,
                windows_per_chunk=config.windows_per_chunk,
                jobs=config.jobs,
                schedule=config.schedule,
                max_windows=config.max_windows,
                validation_level=config.validation_level,
                overwrite=config.overwrite,
                resume=config.resume,
                reference_fasta=config.reference_fasta,
                window_index_dir=str(config.window_index_dir or ""),
                read_coordinate_mode=config.read_coordinate_mode,
                hap1_alignment=str(config.hap1_alignment or ""),
                hap2_alignment=str(config.hap2_alignment or ""),
                export_config=export_config,
                sample_metadata={
                    "sample_id": config.sample_id,
                    "source_sample_id": config.source_sample_id or config.sample_id,
                    "dataset_role": config.dataset_role,
                    "read_type": config.read_type,
                    "project_id": config.project_id,
                    "layout_name": config.layout_name,
                    "hap_read_mode": config.hap_read_mode,
                },
                compact_compressor=config.compact_compressor,
                compact_compression_level=config.compact_compression_level,
                profile_direct_compact=config.profile_direct_compact,
            ))
        with timer.step("summarize_direct_compact_4bit"):
            compact_summary = summarize_direct_compact_outputs(
                sample_id=config.sample_id,
                source_sample_id=config.source_sample_id or config.sample_id,
                dataset_role=config.dataset_role,
                read_type=config.read_type,
                project_id=config.project_id,
                output_sample_root=sample_root,
                compact_root=compact_root,
                contigs=contigs,
            )
    elif config.compact_mode == "two_step":
        with timer.step("run_full_shared_hap_zarr"):
            shared_manifest = run_shared_hap_jobs(config, records, record_tsv, raw_output_root)

        with timer.step("convert_raw_to_compact_4bit"):
            compact_summary = run_legacy_convert(
                LegacyConvertConfig(
                    sample_id=config.sample_id,
                    source_sample_id=config.source_sample_id or config.sample_id,
                    legacy_root=raw_output_root / config.sample_id,
                    output_sample_root=sample_root,
                    dataset_role=config.dataset_role,
                    read_type=config.read_type,
                    project_id=config.project_id,
                    contigs=contigs,
                    overwrite=config.overwrite,
                    rna_bigwigs=config.rna_bigwigs,
                    source_wgs=config.source_wgs or [config.alignment],
                    source_vcfs=effective_source_vcfs(config),
                    layout_name=config.layout_name,
                    hap_invalid_policy="to_n",
                )
            )
    else:
        raise ValueError("compact_mode must be one of: direct, two_step")

    with timer.step("persist_window_annotations"):
        window_annotations = annotation_step(sample_root, compact_root, config.window_index_dir, records, personal_fasta_contract)

    with timer.step("write_full_metadata"):
        summary = {
            "schema": "mfm_integrated_full_shared_hap_preprocess_v0",
            "sample_id": config.sample_id,
            "source_sample_id": config.source_sample_id or config.sample_id,
            "dataset_role": config.dataset_role,
            "read_type": config.read_type,
            "project_id": config.project_id,
            "project_label": config.project_label,
            "platform": config.platform,
            "reference_build": config.reference_build,
            "sample_status": config.sample_status,
            "disease": config.disease,
            "storage_policy": config.storage_policy,
            "window_annotations": window_annotations,
            "compact_mode": config.compact_mode,
            "compact_compressor": config.compact_compressor,
            "compact_compression_level": int(config.compact_compression_level),
            "profile_direct_compact": bool(config.profile_direct_compact),
            "allow_longread_vcf_haplotag": bool(config.allow_longread_vcf_haplotag),
            "allow_unphased_longread": bool(config.allow_unphased_longread),
            "layout_name": config.layout_name,
            "hap_read_mode": config.hap_read_mode,
            "output_sample_root": str(sample_root),
            "raw_output_root": str(raw_output_root),
            "compact_root": str(compact_root),
            "alignment_input": str(config.alignment),
            "alignment_used": str(alignment_for_preprocess),
            "reference_roles": reference_roles(config),
            "phase_summary": phase_summary,
            "hap_fasta_summary": hap_fasta_summary,
            "personal_fasta_reference_contract": personal_fasta_contract,
            "common_vcf_dir": str(config.common_vcf_dir) if config.common_vcf_dir else "",
            "consensus_filter_expr": config.consensus_filter_expr,
            "haplotag_vcf_types": config.haplotag_vcf_types,
            "haplotag_vcf_role": haplotag_vcf_role(config),
            "longread_haplotag_policy": "Validated public or explicit Clair3 phased VCF/HP-PS BAM is required for long-read production; Illumina panel VCF is not used for long-read FASTA or haplotagging",
            "sequence_variant_policy": config.sequence_variant_policy,
            "read_coordinate_mode": config.read_coordinate_mode,
            "hap1_alignment": str(config.hap1_alignment or ""),
            "hap2_alignment": str(config.hap2_alignment or ""),
            "project_policy_guard": {
                "expected_sequence_variant_policy": PROJECT_SEQUENCE_POLICY.get(config.project_id, ""),
                "expected_dataset_role": PROJECT_ROLE_POLICY.get(config.project_id, ""),
                "expected_read_type": PROJECT_READ_TYPE_POLICY.get(config.project_id, ""),
            },
            "sequence_vcf_role": sequence_vcf_role(config),
            "sequence_vcf_types": sequence_consensus_vcf_types(config),
            "sequence_vcf_nonref_min_ac": 1,
            "sequence_reference_fasta": str(effective_sequence_reference(config)) if effective_sequence_reference(config) else "",
            "synthesize_missing_ps": config.synthesize_missing_ps,
            "sample_vcf_for_haplotagging": str(sample_vcf_for_haplotagging) if sample_vcf_for_haplotagging else "",
            "haplotag_vcf_for_whatshap": str(sample_vcf_for_haplotagging or config.phased_vcf or ""),
            "cleanup_removed": cleanup_removed,
            "shared_hap_manifest": shared_manifest,
            "compact_summary": compact_summary,
            "timing": timer.records,
        }
        summary["size_bytes"] = directory_size(sample_root)
        summary["size_human"] = human_bytes(int(summary["size_bytes"]))
        write_full_metadata(sample_root, config, summary, records, commands)

    if config.storage_policy == "model_input_only":
        with timer.step("cleanup_raw_intermediate"):
            cleanup_removed.extend(cleanup_model_input_work(config, raw_output_root, work_root))
            summary["cleanup_removed"] = cleanup_removed
            summary["size_bytes_after_cleanup"] = directory_size(sample_root)
            summary["size_human_after_cleanup"] = human_bytes(int(summary["size_bytes_after_cleanup"]))
            summary["timing"] = timer.records
            write_json(metadata_dir / "preprocess.json", summary)
            write_json(logs_dir / "timing.json", {"timing": timer.records})

    sample_qc = {
        "sample_id": config.sample_id,
        "status": "complete",
        "contig_count": len(contigs),
        "compact_root": str(compact_root),
        "size_bytes": directory_size(sample_root),
        "size_human": human_bytes(directory_size(sample_root)),
        "timing": timer.records,
    }
    write_json(qc_dir / "sample_qc.json", sample_qc)
    return sample_qc | {"output_sample_root": str(sample_root)}


def run_full_manifest(config: FullManifestConfig) -> list[dict[str, Any]]:
    rows = read_manifest_rows(config.manifest)
    config.status_tsv.parent.mkdir(parents=True, exist_ok=True)
    journal_path = Path(f"{config.status_tsv}.progress.sqlite3")
    # 산출물에 영향을 주는 manifest-level 옵션만 signature에 넣는다. status TSV
    # 또는 catalog 경로 변경만으로 이미 완료된 sample을 불필요하게 재실행하지 않는다.
    signature = manifest_progress_signature(rows, {
        "dataset_root": config.dataset_root,
        "jobs": config.jobs,
        "overwrite": config.overwrite,
        "storage_policy": config.storage_policy,
    })
    if not config.resume:
        # --no-resume은 이 실행의 새 journal을 요구한다. sqlite sidecar만 지우며
        # final output·원본·기존 status TSV는 건드리지 않는다.
        for path in (journal_path, Path(f"{journal_path}-wal"), Path(f"{journal_path}-shm")):
            try:
                path.unlink()
            except FileNotFoundError:
                pass
    journal = ManifestProgressJournal(journal_path, signature)
    try:
        completed = journal.completed_sample_ids() if config.resume else set()
        pending_rows = [row for row in rows if row["sample_id"] not in completed]
        if completed:
            print(f"[RESUME] completed={len(completed)} pending={len(pending_rows)}")
        if config.sample_jobs <= 1:
            for row in pending_rows:
                result = run_full_manifest_row(row, config)
                print_manifest_result(result)
                journal.record(result)
        else:
            with ProcessPoolExecutor(max_workers=config.sample_jobs) as executor:
                future_to_row = {
                    executor.submit(run_full_manifest_row, row, config): row for row in pending_rows
                }
                for future in as_completed(future_to_row):
                    result = future.result()
                    print_manifest_result(result)
                    journal.record(result)
        results = journal.results()
        results.sort(key=lambda item: int(item.get("manifest_order", 0)))
        write_status_tsv(config.status_tsv, results)
        if config.sample_catalog is not None:
            complete_results = [item for item in results if item.get("status") == "complete"]
            if complete_results:
                upsert_sample_catalog(config.sample_catalog, rows, complete_results)
        journal.cleanup()
        return results
    finally:
        # 예외/interrupt에서는 journal을 남겨 다음 실행에서 완료 sample만 재사용한다.
        journal.close()


def persist_manifest_progress(
    config: FullManifestConfig,
    manifest_rows: list[dict[str, str]],
    results: list[dict[str, Any]],
) -> None:
    ordered = sorted(results, key=lambda item: int(item.get("manifest_order", 0)))
    write_status_tsv(config.status_tsv, ordered)
    if config.sample_catalog is not None:
        complete_results = [result for result in ordered if result.get("status") == "complete"]
        if complete_results:
            upsert_sample_catalog(config.sample_catalog, manifest_rows, complete_results)


def run_full_manifest_row(row: dict[str, str], manifest_config: FullManifestConfig) -> dict[str, Any]:
    start = time.monotonic()
    sample_id = row["sample_id"]
    try:
        config = config_from_row(row, manifest_config)
        result = run_full_preprocess(config)
        elapsed = time.monotonic() - start
        return {
            "manifest_order": int(row.get("_manifest_order", 0)),
            "sample_id": sample_id,
            "status": "complete",
            "elapsed_seconds": f"{elapsed:.3f}",
            "elapsed_hms": format_hms(elapsed),
            "output_sample_root": result.get("output_sample_root", ""),
            "size_bytes": result.get("size_bytes", 0),
            "size_human": result.get("size_human", ""),
            "error": "",
        }
    except Exception as exc:  # noqa: BLE001
        elapsed = time.monotonic() - start
        return {
            "manifest_order": int(row.get("_manifest_order", 0)),
            "sample_id": sample_id,
            "status": "failed",
            "elapsed_seconds": f"{elapsed:.3f}",
            "elapsed_hms": format_hms(elapsed),
            "output_sample_root": row.get("output_sample_root", ""),
            "size_bytes": 0,
            "size_human": "0.00 B",
            "error": f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}",
        }


def run_shared_hap_jobs(
    config: FullPreprocessConfig,
    records: list[dict[str, str]],
    record_tsv: Path,
    output_root: Path,
) -> dict[str, Any]:
    export_config = stage70_export_config(config)
    jobs = [
        {
            "sample": row["sample"],
            "contig": row["contig"],
            "hap1_fasta": row["hap1_fasta"],
            "hap1_chain": row["hap1_chain"],
            "hap2_fasta": row["hap2_fasta"],
            "hap2_chain": row["hap2_chain"],
            "bam": row["bam"],
            "output_root": str(output_root),
            "windows_per_chunk": config.windows_per_chunk,
            "resume": config.resume,
            "overwrite": config.overwrite,
            "max_windows": config.max_windows,
            "validation_level": config.validation_level,
            "reference_fasta": str(config.reference_fasta) if config.reference_fasta else "",
            "config": export_config.to_dict(),
        }
        for row in records
    ]
    if config.schedule == "longest-first":
        jobs = sorted(jobs, key=_job_schedule_key)
    manifest = init_manifest(
        record_tsv=record_tsv,
        output_root=output_root,
        windows_per_chunk=config.windows_per_chunk,
        jobs=config.jobs,
        schedule=config.schedule,
        max_windows=config.max_windows,
        validation_level=config.validation_level,
        config=export_config.to_dict(),
    )
    manifest_path = output_root / "manifest.json"
    write_json(manifest_path, manifest)
    try:
        if config.jobs <= 1:
            for job in jobs:
                result = _run_contig(job)
                update_manifest_record(manifest, result)
                write_json(manifest_path, manifest)
                print(json.dumps(result, ensure_ascii=False), flush=True)
        else:
            with ProcessPoolExecutor(max_workers=config.jobs) as pool:
                futures = [pool.submit(_run_contig, job) for job in jobs]
                for future in as_completed(futures):
                    result = future.result()
                    update_manifest_record(manifest, result)
                    write_json(manifest_path, manifest)
                    print(json.dumps(result, ensure_ascii=False), flush=True)
    except BaseException as exc:
        mark_manifest_failed(manifest, exc)
        write_json(manifest_path, manifest)
        raise
    mark_manifest_complete(manifest)
    write_json(manifest_path, manifest)
    return manifest


def enforce_project_policy(config: FullPreprocessConfig) -> dict[str, str]:
    """Apply process_code cohort policy guards to integrated preprocessing.

    Known IGSR projects are not allowed to fall back to a user-supplied policy.
    WGS-only and GEUVADIS both use SNV+indel HAP FASTA in v2. GEUVADIS RNA
    targets must be realigned/projected to the personal HAP coordinate system
    instead of forcing SNP-only identity coordinates.
    """
    policy = normalize_sequence_variant_policy(config.sequence_variant_policy)
    config.sequence_variant_policy = policy
    expected_policy = PROJECT_SEQUENCE_POLICY.get(config.project_id)
    if expected_policy is not None and policy != expected_policy:
        raise ValueError(
            f"project_id={config.project_id} requires sequence_variant_policy={expected_policy}, "
            f"but got {policy}"
        )
    expected_role = PROJECT_ROLE_POLICY.get(config.project_id)
    if expected_role is not None and config.dataset_role != expected_role:
        raise ValueError(
            f"project_id={config.project_id} requires dataset_role={expected_role}, "
            f"but got {config.dataset_role}"
        )
    expected_read_type = PROJECT_READ_TYPE_POLICY.get(config.project_id)
    if expected_read_type is not None and config.read_type != expected_read_type:
        raise ValueError(
            f"project_id={config.project_id} requires read_type={expected_read_type}, "
            f"but got {config.read_type}"
        )
    if config.project_id == "igsr_ont_wgs_only" and config.rna_bigwigs:
        raise ValueError("igsr_ont_wgs_only is WGS-only; RNA inputs must use a GEUVADIS/RNA project_id")
    return {
        "project_id": config.project_id,
        "sequence_variant_policy": policy,
        "dataset_role": config.dataset_role,
        "read_type": config.read_type,
    }


def normalize_platform(value: str) -> str:
    return (value or "").strip().lower().replace("-", "_").replace(" ", "_")


def is_long_read_config(config: FullPreprocessConfig) -> bool:
    if config.read_type == "long_read":
        return True
    normalized = normalize_platform(config.platform)
    return (
        normalized in LONG_READ_PLATFORMS
        or normalized.startswith("ont")
        or "nanopore" in normalized
        or "hifi" in normalized
        or "ccs" in normalized
        or normalized.startswith("pacbio")
    )


def normalize_sequence_variant_policy(value: str) -> str:
    value = (value or DEFAULT_SEQUENCE_VARIANT_POLICY).strip().lower().replace("-", "_")
    aliases = {
        "snp": "snp_only",
        "snps": "snp_only",
        "snp_only": "snp_only",
        "snv_indel": "snv_indel",
        "snps_indels": "snv_indel",
        "snps,indels": "snv_indel",
    }
    if value not in aliases:
        raise ValueError(f"unknown sequence variant policy: {value}")
    return aliases[value]


def sequence_policy_needs_chain(config: FullPreprocessConfig) -> bool:
    return normalize_sequence_variant_policy(config.sequence_variant_policy) != "snp_only"


def sequence_policy_type_expr(config: FullPreprocessConfig) -> str:
    policy = normalize_sequence_variant_policy(config.sequence_variant_policy)
    if policy == "snp_only":
        return 'TYPE="snp"'
    return '(TYPE="snp" || TYPE="indel")'


def consensus_include_expression(config: FullPreprocessConfig) -> str:
    filter_expr = (config.consensus_filter_expr or DEFAULT_CONSENSUS_FILTER_EXPR).strip()
    return f"({filter_expr}) && ({sequence_policy_type_expr(config)})"


def sequence_consensus_vcf_types(config: FullPreprocessConfig) -> str:
    policy = normalize_sequence_variant_policy(config.sequence_variant_policy)
    if policy == "snp_only":
        return "snps"
    return "snps,indels"


def vcf_type_label(value: str) -> str:
    parts = [part.strip().lower().replace("-", "_") for part in str(value).split(",") if part.strip()]
    return "_".join(parts) if parts else "variants"


def sequence_vcf_role(config: FullPreprocessConfig) -> str:
    policy = normalize_sequence_variant_policy(config.sequence_variant_policy)
    if policy == "snp_only":
        return "sequence_snp_nonref_for_consensus"
    return "sequence_snv_indel_nonref_for_consensus"


def haplotag_vcf_role(config: FullPreprocessConfig) -> str:
    label = vcf_type_label(config.haplotag_vcf_types or DEFAULT_HAPLOTAG_VCF_TYPES)
    return f"haplotag_het_{label}_for_whatshap"


def reference_roles(config: FullPreprocessConfig) -> dict[str, str]:
    return {
        "cram_decode_and_haplotag": str(config.reference_fasta or ""),
        "personal_sequence_consensus": str(effective_sequence_reference(config) or ""),
    }


def effective_sequence_reference(config: FullPreprocessConfig) -> Path | None:
    """Return the reference used to build personal HAP FASTA.

    CRAM decoding and WhatsHap haplotagging use ``reference_fasta`` (Vienna for
    IGSR ONT). Personal FASTA consensus must use the separately agreed UCSC
    hg38 sequence reference. We intentionally do not fall back to
    ``reference_fasta`` here, because that would silently build HAP FASTA on the
    CRAM/Vienna reference when ``--sequence-reference-fasta`` is omitted.
    """
    return config.sequence_reference_fasta


def parse_applied_variant_count(stderr: str) -> int:
    matches = re.findall(r"Applied\s+(\d+)", stderr or "")
    return sum(int(value) for value in matches) if matches else 0


def is_symbolic_alt(alt: str | None) -> bool:
    if alt is None:
        return True
    return alt.startswith("<") or alt in {"*", "."}


def count_bad_bases(text_value: str | None) -> int:
    if text_value is None or is_symbolic_alt(text_value):
        return 0
    return sum(1 for base in text_value.upper() if base not in CONSENSUS_ALPHABET)


def genotype_tuple(record: pysam.VariantRecord, sample_id: str) -> tuple[int | None, ...]:
    if sample_id not in record.samples:
        raise ValueError(f"VCF record does not contain sample {sample_id}: {record.contig}:{record.pos}")
    gt = record.samples[sample_id].get("GT")
    return tuple(gt or ())


def genotype_is_hom_ref(gt: tuple[int | None, ...]) -> bool:
    called = [allele for allele in gt if allele is not None]
    return bool(called) and all(allele == 0 for allele in called)


def genotype_is_hom(gt: tuple[int | None, ...]) -> bool:
    called = [allele for allele in gt if allele is not None]
    return bool(called) and len(set(called)) == 1


def summarize_sample_vcf(vcf: Path, sample_id: str) -> dict[str, int]:
    summary = {
        "variant_count": 0,
        "indel_count": 0,
        "hom_ref_count": 0,
        "hom_count": 0,
        "multiallelic_count": 0,
        "symbolic_alt_count": 0,
        "bad_base_count": 0,
    }
    with pysam.VariantFile(str(vcf)) as handle:
        if sample_id not in handle.header.samples:
            raise ValueError(f"VCF에 sample이 없습니다: sample={sample_id}, vcf={vcf}")
        for record in handle.fetch():
            summary["variant_count"] += 1
            alts = tuple(record.alts or ())
            if len(alts) > 1:
                summary["multiallelic_count"] += 1
            if any(is_symbolic_alt(alt) for alt in alts):
                summary["symbolic_alt_count"] += 1
            summary["bad_base_count"] += count_bad_bases(record.ref)
            for alt in alts:
                summary["bad_base_count"] += count_bad_bases(alt)
            if any(alt and not is_symbolic_alt(alt) and len(record.ref or "") != len(alt) for alt in alts):
                summary["indel_count"] += 1
            gt = genotype_tuple(record, sample_id)
            if genotype_is_hom_ref(gt):
                summary["hom_ref_count"] += 1
            if genotype_is_hom(gt):
                summary["hom_count"] += 1
    return summary


def validate_sequence_vcf_summary(vcf: Path, sample_id: str) -> dict[str, int]:
    summary = summarize_sample_vcf(vcf, sample_id)
    if summary["variant_count"] <= 0:
        raise RuntimeError(f"sequence VCF 변이가 0개입니다: {vcf}")
    if summary["hom_ref_count"]:
        raise RuntimeError(f"sequence VCF에 hom-ref가 남아 있습니다: {vcf} hom_ref={summary['hom_ref_count']}")
    if summary["symbolic_alt_count"]:
        raise RuntimeError(f"sequence VCF에 symbolic ALT/SV가 남아 있습니다: {vcf} symbolic={summary['symbolic_alt_count']}")
    if summary["bad_base_count"]:
        raise RuntimeError(f"sequence VCF REF/ALT에 ATGCN 외 문자가 있습니다: {vcf} bad_base={summary['bad_base_count']}")
    return summary


def validate_haplotag_vcf_summary(vcf: Path, sample_id: str) -> dict[str, int]:
    summary = summarize_sample_vcf(vcf, sample_id)
    if summary["variant_count"] <= 0:
        raise RuntimeError(f"haplotag VCF 변이가 0개입니다: {vcf}")
    if summary["hom_count"]:
        raise RuntimeError(f"haplotag VCF에 homozygous GT가 남아 있습니다: {vcf} hom={summary['hom_count']}")
    if summary["multiallelic_count"]:
        raise RuntimeError(f"haplotag VCF에 multiallelic record가 남아 있습니다: {vcf} multiallelic={summary['multiallelic_count']}")
    if summary["symbolic_alt_count"]:
        raise RuntimeError(f"haplotag VCF에 symbolic ALT/SV가 남아 있습니다: {vcf} symbolic={summary['symbolic_alt_count']}")
    if summary["bad_base_count"]:
        raise RuntimeError(f"haplotag VCF REF/ALT에 ATGCN 외 문자가 있습니다: {vcf} bad_base={summary['bad_base_count']}")
    return summary

def run_haplotag(
    config: FullPreprocessConfig,
    work_root: Path,
    commands: list[str],
    phased_vcf_override: Path | None = None,
) -> Path:
    if is_long_read_config(config):
        raise ValueError(
            "v2 policy forbids panel/VCF-based WhatsHap haplotagging for long-read data. "
            "Long-read input must already carry HP/PS from the manifest-approved public or explicit Clair3 VCF path before this compact stage."
        )
    if config.reference_fasta is None:
        raise ValueError("--reference-fasta is required for --haplotag-if-missing")
    phased_vcf = phased_vcf_override or config.phased_vcf
    if phased_vcf is None:
        phased_vcf = prepare_sample_phased_vcf(config, work_root, commands)
    out_dir = work_root / "haplotagged"
    out_dir.mkdir(parents=True, exist_ok=True)
    tagged_bam = out_dir / f"{config.sample_id}.haplotagged.bam"
    if tagged_bam.exists() and config.resume and not config.overwrite:
        if config.synthesize_missing_ps:
            normalize_zero_ps_for_hp_reads(tagged_bam, commands)
        ensure_alignment_index(config.samtools_bin, tagged_bam, commands)
        return tagged_bam
    cmd = [
        config.whatshap_bin,
        "haplotag",
        "--sample",
        config.sample_id,
        "--reference",
        str(config.reference_fasta),
        "--skip-missing-contigs",
        "--regions",
        " ".join(expand_contigs(config.contigs)),
        "--output-threads",
        str(max(1, min(int(config.jobs or 1), 8))),
        "--output",
        str(tagged_bam),
        str(phased_vcf),
        str(config.alignment),
    ]
    commands.append(shell_join(cmd))
    subprocess.run(cmd, check=True)
    if config.synthesize_missing_ps:
        normalize_zero_ps_for_hp_reads(tagged_bam, commands)
    else:
        commands.append("# preserved source PS tags; missing or zero PS remains unknown")
    ensure_alignment_index(config.samtools_bin, tagged_bam, commands)
    return tagged_bam


def normalize_zero_ps_for_hp_reads(tagged_bam: Path, commands: list[str]) -> dict[str, int]:
    summary_path = tagged_bam.with_suffix(".ps_normalize.json")
    if summary_path.exists():
        return json.loads(summary_path.read_text(encoding="utf-8"))
    normalized = tagged_bam.with_suffix(".ps_normalized.bam")
    counts = {
        "total_reads": 0,
        "hp_tagged_reads": 0,
        "ps_preserved_reads": 0,
        "ps_synthetic_reads": 0,
        "ps_missing_or_invalid_reads": 0,
    }
    with pysam.AlignmentFile(str(tagged_bam), "rb") as src, pysam.AlignmentFile(str(normalized), "wb", template=src) as dst:
        for read in src.fetch(until_eof=True):
            counts["total_reads"] += 1
            hp = read_tag_int(read, "HP")
            if hp in (1, 2):
                counts["hp_tagged_reads"] += 1
                ps = read_tag_int(read, "PS")
                if ps is not None and ps > 0:
                    counts["ps_preserved_reads"] += 1
                else:
                    if ps is None:
                        counts["ps_missing_or_invalid_reads"] += 1
                    read.set_tag("PS", 1, value_type="i")
                    counts["ps_synthetic_reads"] += 1
            dst.write(read)
    os.replace(normalized, tagged_bam)
    Path(str(tagged_bam) + ".bai").unlink(missing_ok=True)
    Path(str(tagged_bam) + ".csi").unlink(missing_ok=True)
    write_json(summary_path, {
        "schema": "mfm_haplotagged_bam_ps_normalize_v1",
        "policy": "HP-tagged reads with missing, invalid, or zero PS get synthetic PS=1. For 1000G SHAPEIT common VCF this represents a chromosome-wide phased set because the VCF has phased GT but no explicit PS blocks.",
        **counts,
    })
    commands.append("# normalized HP-tagged reads with PS<=0 to synthetic PS=1 in " + shlex.quote(str(tagged_bam)))
    return counts


def read_tag_int(read: pysam.AlignedSegment, tag: str) -> int | None:
    if not read.has_tag(tag):
        return None
    try:
        return int(read.get_tag(tag))
    except (TypeError, ValueError):
        return None


def validate_consensus_length_policy(config: FullPreprocessConfig, results: list[dict[str, Any]]) -> dict[str, Any]:
    policy = normalize_sequence_variant_policy(config.sequence_variant_policy)
    by_hap: dict[str, dict[str, int]] = {}
    for result in results:
        hap = str(result.get("haplotype", ""))
        if hap not in by_hap:
            by_hap[hap] = {"contigs": 0, "same_length": 0, "different_length": 0}
        by_hap[hap]["contigs"] += 1
        if int(result.get("reference_length", -1)) == int(result.get("consensus_length", -2)):
            by_hap[hap]["same_length"] += 1
        else:
            by_hap[hap]["different_length"] += 1
    if policy == "snp_only":
        bad = {hap: row for hap, row in by_hap.items() if row["different_length"]}
        if bad:
            raise RuntimeError(f"snp_only 정책인데 개인 FASTA 길이가 reference와 다른 hap이 있습니다: {bad}")
        expectation = "all_contigs_same_as_reference"
    else:
        bad = {hap: row for hap, row in by_hap.items() if row["different_length"] == 0}
        if bad:
            raise RuntimeError(f"snv_indel 정책인데 모든 contig 길이가 reference와 같습니다. indel 미적용 의심: {bad}")
        expectation = "at_least_one_contig_differs_per_haplotype"
    return {
        "schema": "mfm_consensus_length_policy_qc_v1",
        "sequence_variant_policy": policy,
        "expectation": expectation,
        "haplotypes": by_hap,
    }


def chain_header_qc(chain_path: Path, contig: str, reference_length: int, consensus_length: int) -> dict[str, Any]:
    if not chain_path.exists() or not chain_path.stat().st_size:
        raise RuntimeError(f"chain 파일이 없거나 비었습니다: {chain_path}")
    with chain_path.open("r", encoding="ascii", errors="replace") as handle:
        header = handle.readline().strip().split()
    if len(header) < 13 or header[0] != "chain":
        raise RuntimeError(f"chain header가 올바르지 않습니다: {chain_path}")
    qc = {
        "chain": str(chain_path),
        "target_name": header[2],
        "target_size": int(header[3]),
        "target_strand": header[4],
        "query_name": header[7],
        "query_size": int(header[8]),
        "query_strand": header[9],
    }
    if qc["target_name"] != contig or qc["query_name"] != contig:
        raise RuntimeError(f"chain contig 불일치: {chain_path} qc={qc}")
    if qc["target_size"] != reference_length:
        raise RuntimeError(f"chain tSize가 reference 길이와 다릅니다: {chain_path} qc={qc} reference={reference_length}")
    if qc["query_size"] != consensus_length:
        raise RuntimeError(f"chain qSize가 개인 FASTA 길이와 다릅니다: {chain_path} qc={qc} consensus={consensus_length}")
    if qc["target_strand"] != "+" or qc["query_strand"] != "+":
        raise RuntimeError(f"plus/plus가 아닌 chain은 지원하지 않습니다: {chain_path} qc={qc}")
    qc["status"] = "pass"
    return qc


def prepare_common_vcf_hap_fastas(
    config: FullPreprocessConfig,
    contigs: list[str],
    work_root: Path,
    commands: list[str],
) -> dict[str, Any]:
    if config.common_vcf_dir is None:
        raise ValueError("--common-vcf-dir is required for --auto-prepare-hap-fasta")
    if effective_sequence_reference(config) is None:
        raise ValueError("--sequence-reference-fasta is required for --auto-prepare-hap-fasta; --reference-fasta is reserved for CRAM decoding/WhatsHap")
    output_root = work_root / "individual_fasta"
    output_root.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    for contig in contigs:
        vcf = resolve_common_vcf_path(config.common_vcf_dir, contig, config.common_vcf_pattern)
        ensure_vcf_index(vcf)
        ensure_vcf_has_sample(vcf, config.sample_id)
        for mode in ("HAP1", "HAP2"):
            completed = load_completed_common_consensus(output_root, contig, mode, config)
            if completed is not None:
                results.append(completed | {"status": "skipped_existing"})
                continue
            result = create_common_vcf_consensus_contig(config, vcf, contig, mode, output_root, commands)
            results.append(result)
    config.hap1_fasta_dir = output_root / "HAP1"
    config.hap2_fasta_dir = output_root / "HAP2"
    if sequence_policy_needs_chain(config):
        config.hap1_chain_dir = output_root / "HAP1"
        config.hap2_chain_dir = output_root / "HAP2"
    else:
        config.hap1_chain_dir = None
        config.hap2_chain_dir = None
    length_policy_qc = validate_consensus_length_policy(config, results)
    summary = {
        "schema": "mfm_common_vcf_consensus_autoprepare_v1",
        "sample_id": config.sample_id,
        "common_vcf_dir": str(config.common_vcf_dir),
        "common_vcf_pattern": config.common_vcf_pattern,
        "sequence_variant_policy": normalize_sequence_variant_policy(config.sequence_variant_policy),
        "sequence_vcf_role": sequence_vcf_role(config),
        "sequence_vcf_types": sequence_consensus_vcf_types(config),
        "sequence_vcf_nonref_min_ac": 1,
        "consensus_filter_expr": config.consensus_filter_expr,
        "consensus_include_expression": consensus_include_expression(config),
        "reference_roles": reference_roles(config),
        "sequence_reference_fasta": str(effective_sequence_reference(config)),
        "chain_required": sequence_policy_needs_chain(config),
        "length_policy_qc": length_policy_qc,
        "haplotag_vcf_types": config.haplotag_vcf_types,
        "haplotag_vcf_role": haplotag_vcf_role(config),
            "longread_haplotag_policy": "Validated public or explicit Clair3 phased VCF/HP-PS BAM is required for long-read production; Illumina panel VCF is not used for long-read FASTA or haplotagging",
        "synthesize_missing_ps": config.synthesize_missing_ps,
        "output_root": str(output_root),
        "contigs": contigs,
        "result_count": len(results),
        "results": results,
    }
    write_json(output_root / "manifest.common_vcf.json", summary)
    return summary


def create_common_vcf_consensus_contig(
    config: FullPreprocessConfig,
    vcf: Path,
    contig: str,
    mode: str,
    output_root: Path,
    commands: list[str],
) -> dict[str, Any]:
    haplotype = "1" if mode == "HAP1" else "2"
    directory = output_root / mode
    directory.mkdir(parents=True, exist_ok=True)
    fasta_path = directory / f"{contig}.fa"
    needs_chain = sequence_policy_needs_chain(config)
    chain_path = directory / f"{contig}.chain" if needs_chain else None
    qc_path = directory / f"{contig}.qc.json"
    existing_paths = [fasta_path, qc_path] + ([chain_path] if chain_path is not None else [])
    if not config.overwrite and any(path.exists() for path in existing_paths):
        raise FileExistsError(f"기존 출력이 있습니다: {directory}/{contig}.*")
    fasta_partial = fasta_path.with_suffix(".fa.partial")
    chain_partial = chain_path.with_suffix(".chain.partial") if chain_path is not None else None
    stderr_partial = qc_path.with_suffix(".bcftools.stderr.partial")
    for path in (fasta_partial, chain_partial, stderr_partial):
        if path is not None:
            path.unlink(missing_ok=True)
    include_expression = consensus_include_expression(config)
    sequence_vcf = prepare_sequence_consensus_vcf_contig(config, vcf, contig, output_root, commands)
    cmd = [
        config.bcftools_bin,
        "consensus",
        "--samples",
        config.sample_id,
        "--haplotype",
        haplotype,
        "--include",
        include_expression,
    ]
    if chain_partial is not None:
        cmd.extend(["--chain", str(chain_partial)])
    cmd.extend(["--output", str(fasta_partial), str(sequence_vcf)])
    commands.append(
        "# sequence_vcf_for_consensus="
        + shlex.quote(str(sequence_vcf))
        + "; common_source_vcf="
        + shlex.quote(str(vcf))
    )
    commands.append(
        shell_join([config.samtools_bin, "faidx", str(effective_sequence_reference(config)), contig])
        + " | "
        + shell_join(cmd)
    )
    stderr_handle = stderr_partial.open("wb")
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=stderr_handle)
    try:
        write_reference_stream(process, require_path(effective_sequence_reference(config), "sequence_reference_fasta/reference_fasta"), contig)
        return_code = process.wait()
        stderr_handle.close()
        stderr = stderr_partial.read_text(encoding="utf-8", errors="replace") if stderr_partial.exists() else ""
        if return_code != 0:
            raise RuntimeError(f"bcftools consensus 실패({return_code}): {stderr}")
        if not fasta_partial.exists() or not fasta_partial.stat().st_size:
            raise RuntimeError("bcftools consensus가 빈 FASTA를 생성했습니다")
        if chain_partial is not None:
            if not chain_partial.exists() or not chain_partial.stat().st_size:
                raise RuntimeError("bcftools consensus가 빈 chain을 생성했습니다")
        os.replace(fasta_partial, fasta_path)
        if chain_partial is not None and chain_path is not None:
            os.replace(chain_partial, chain_path)
    except BaseException:
        if process.poll() is None:
            process.kill()
            process.wait()
        if not stderr_handle.closed:
            stderr_handle.close()
        for path in (fasta_partial, chain_partial, stderr_partial):
            if path is not None:
                path.unlink(missing_ok=True)
        raise
    finally:
        if not stderr_handle.closed:
            stderr_handle.close()
    sanitized_invalid_bases = sanitize_fasta_to_acgtn(fasta_path)
    pysam.faidx(str(fasta_path))
    qc = consensus_qc(config, sequence_vcf, contig, mode, fasta_path, chain_path, stderr)
    qc["common_source_vcf"] = str(vcf)
    qc["sanitized_invalid_bases"] = sanitized_invalid_bases
    write_json(qc_path, qc)
    stderr_partial.unlink(missing_ok=True)
    if qc["invalid_bases"]:
        raise RuntimeError(f"생성 FASTA에 허용하지 않는 base가 있습니다: {qc_path}")
    if qc.get("applied_variant_count", 0) <= 0:
        raise RuntimeError(f"consensus 적용 변이가 0개입니다. FILTER/type 정책 또는 VCF sample을 확인하세요: {qc_path}")
    return qc


def load_completed_common_consensus(output_root: Path, contig: str, mode: str, config: FullPreprocessConfig) -> dict[str, Any] | None:
    if config.overwrite or not config.resume:
        return None
    directory = output_root / mode
    fasta_path = directory / f"{contig}.fa"
    chain_path = directory / f"{contig}.chain" if sequence_policy_needs_chain(config) else None
    qc_path = directory / f"{contig}.qc.json"
    index_path = Path(str(fasta_path) + ".fai")
    paths = tuple(path for path in (fasta_path, chain_path, qc_path, index_path) if path is not None)
    if not any(path.exists() for path in paths):
        return None
    if not all(path.exists() and path.stat().st_size for path in paths):
        raise ValueError(f"불완전한 기존 consensus 출력이 있습니다: {paths}")
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    if qc.get("status") != "pass":
        raise ValueError(f"기존 consensus QC가 pass가 아닙니다: {qc_path}")
    return qc


def sanitize_fasta_to_acgtn(fasta_path: Path) -> dict[str, int]:
    invalid_counts: dict[str, int] = {}
    partial = fasta_path.with_suffix(fasta_path.suffix + ".sanitize.partial")
    allowed = CONSENSUS_ALPHABET
    with fasta_path.open("r", encoding="ascii", errors="replace") as src, partial.open("w", encoding="ascii") as dst:
        for line in src:
            if line.startswith(">"):
                dst.write(line)
                continue
            newline = "\n" if line.endswith("\n") else ""
            seq = line.rstrip("\r\n").upper()
            out = []
            for base in seq:
                if base in allowed:
                    out.append(base)
                else:
                    invalid_counts[base] = invalid_counts.get(base, 0) + 1
                    out.append("N")
            dst.write("".join(out) + newline)
    os.replace(partial, fasta_path)
    Path(str(fasta_path) + ".fai").unlink(missing_ok=True)
    return invalid_counts


def consensus_qc(
    config: FullPreprocessConfig,
    vcf: Path,
    contig: str,
    mode: str,
    fasta_path: Path,
    chain_path: Path | None,
    stderr: str,
) -> dict[str, Any]:
    reference = require_path(effective_sequence_reference(config), "sequence_reference_fasta/reference_fasta")
    with pysam.FastaFile(str(reference)) as ref:
        reference_length = ref.get_reference_length(contig)
    invalid_bases: dict[str, int] = {}
    n_count = 0
    with pysam.FastaFile(str(fasta_path)) as consensus:
        if tuple(consensus.references) != (contig,):
            raise ValueError(f"생성 FASTA contig 불일치: {consensus.references}")
        consensus_length = consensus.get_reference_length(contig)
        for start in range(0, consensus_length, 1 << 20):
            seq = consensus.fetch(contig, start, min(consensus_length, start + (1 << 20))).upper()
            n_count += seq.count("N")
            for base in seq:
                if base not in CONSENSUS_ALPHABET:
                    invalid_bases[base] = invalid_bases.get(base, 0) + 1
    overlap_skip_count = stderr.count("overlaps with another variant, skipping")
    chain_qc = chain_header_qc(chain_path, contig, reference_length, consensus_length) if chain_path else None
    return {
        "schema": "mfm_common_vcf_consensus_qc_v1",
        "status": "pass" if not invalid_bases else "fail",
        "sample": config.sample_id,
        "haplotype": mode,
        "contig": contig,
        "reference_fasta": str(reference),
        "reference_role": "personal_sequence_consensus",
        "source_vcf": str(vcf),
        "sequence_vcf_for_consensus": str(vcf),
        "sequence_vcf_role": sequence_vcf_role(config),
        "sequence_vcf_types": sequence_consensus_vcf_types(config),
        "sequence_vcf_nonref_min_ac": 1,
        "sequence_variant_policy": normalize_sequence_variant_policy(config.sequence_variant_policy),
        "consensus_filter_expr": config.consensus_filter_expr,
        "consensus_include_expression": consensus_include_expression(config),
        "reference_length": reference_length,
        "consensus_length": consensus_length,
        "n_count": n_count,
        "invalid_bases": invalid_bases,
        "overlap_skip_count": overlap_skip_count,
        "applied_variant_count": parse_applied_variant_count(stderr),
        "fasta": str(fasta_path),
        "chain": str(chain_path) if chain_path else "",
        "chain_qc": chain_qc,
        "fasta_md5": file_md5(fasta_path),
        "chain_md5": file_md5(chain_path) if chain_path else "",
    }

def prepare_sequence_consensus_vcf_contig(
    config: FullPreprocessConfig,
    source_vcf: Path,
    contig: str,
    output_root: Path,
    commands: list[str],
) -> Path:
    out_dir = output_root / "sequence_vcf"
    out_dir.mkdir(parents=True, exist_ok=True)
    role = sequence_vcf_role(config)
    split_vcf = out_dir / f"{config.sample_id}.{contig}.{role}.split.vcf.gz"
    output_vcf = out_dir / f"{config.sample_id}.{contig}.{role}.refok.vcf.gz"
    metadata_path = out_dir / f"{config.sample_id}.{contig}.{role}.json"
    if output_vcf.exists() and has_vcf_index(output_vcf) and metadata_path.exists() and config.resume and not config.overwrite:
        metadata = read_optional_json(metadata_path) or {}
        if metadata.get("schema") == "mfm_sequence_vcf_from_common_v2" and metadata.get("ref_check_policy") == "bcftools_norm_check_ref_x":
            return output_vcf
    for path in (
        split_vcf,
        Path(str(split_vcf) + ".tbi"),
        Path(str(split_vcf) + ".csi"),
        output_vcf,
        Path(str(output_vcf) + ".tbi"),
        Path(str(output_vcf) + ".csi"),
    ):
        path.unlink(missing_ok=True)
    filter_expr = (config.consensus_filter_expr or DEFAULT_CONSENSUS_FILTER_EXPR).strip()
    cmd1 = [config.bcftools_bin, "view", "-s", config.sample_id, "-Ou", str(source_vcf)]
    cmd2 = [config.bcftools_bin, "norm", "-m", "-any", "-Ou", "-"]
    cmd3 = [
        config.bcftools_bin,
        "view",
        "-v",
        sequence_consensus_vcf_types(config),
        "-c",
        "1",
        "-i",
        filter_expr,
        "-Oz",
        "-o",
        str(split_vcf),
        "-",
    ]
    run_bcftools_pipeline([cmd1, cmd2, cmd3], commands, "bcftools sequence consensus VCF split/filter 실패")
    index_vcf(config, split_vcf, commands)
    split_summary = validate_sequence_vcf_summary(split_vcf, config.sample_id)
    reference = require_path(effective_sequence_reference(config), "sequence_reference_fasta/reference_fasta")
    refcheck_cmd = [
        config.bcftools_bin,
        "norm",
        "--check-ref",
        "x",
        "-f",
        str(reference),
        "-Oz",
        "-o",
        str(output_vcf),
        str(split_vcf),
    ]
    commands.append(shell_join(refcheck_cmd))
    subprocess.run(refcheck_cmd, check=True)
    index_vcf(config, output_vcf, commands)
    refok_summary = validate_sequence_vcf_summary(output_vcf, config.sample_id)
    dropped_ref_mismatch = split_summary["variant_count"] - refok_summary["variant_count"]
    write_json(
        metadata_path,
        {
            "schema": "mfm_sequence_vcf_from_common_v2",
            "role": role,
            "sample_id": config.sample_id,
            "contig": contig,
            "common_source_vcf": str(source_vcf),
            "filter_order": "sample_subset_then_norm_split_then_type_nonref_then_check_ref",
            "multiallelic_policy": "bcftools_norm_m_any_preserve_as_split_records",
            "ref_check_policy": "bcftools_norm_check_ref_x",
            "filter_expr": filter_expr,
            "vcf_types": sequence_consensus_vcf_types(config),
            "nonref_min_ac": 1,
            "sequence_variant_policy": normalize_sequence_variant_policy(config.sequence_variant_policy),
            "reference_role": "personal_sequence_consensus",
            "sequence_reference_fasta": str(reference),
            "split_vcf": str(split_vcf),
            "split_summary": split_summary,
            "vcf": str(output_vcf),
            "refok_summary": refok_summary,
            "dropped_ref_mismatch_count": dropped_ref_mismatch,
            "size_bytes": output_vcf.stat().st_size,
        },
    )
    return output_vcf

def prepare_sample_phased_vcf(config: FullPreprocessConfig, work_root: Path, commands: list[str]) -> Path:
    if config.common_vcf_dir is None:
        raise ValueError("--phased-vcf or --common-vcf-dir is required for --haplotag-if-missing")
    out_dir = work_root / "common_vcf"
    out_dir.mkdir(parents=True, exist_ok=True)
    contigs = expand_contigs(config.contigs)
    role = haplotag_vcf_role(config)
    combined = out_dir / f"{config.sample_id}.{role}.vcf.gz"
    metadata_path = out_dir / f"{config.sample_id}.{role}.json"
    if combined.exists() and has_vcf_index(combined) and metadata_path.exists() and config.resume and not config.overwrite:
        metadata = read_optional_json(metadata_path) or {}
        if metadata.get("schema") == "mfm_haplotag_vcf_from_common_v2" and metadata.get("biallelic_policy") == "bcftools_view_m2_M2":
            validate_haplotag_vcf_summary(combined, config.sample_id)
            return combined
    combined.unlink(missing_ok=True)
    Path(str(combined) + ".tbi").unlink(missing_ok=True)
    Path(str(combined) + ".csi").unlink(missing_ok=True)
    part_paths: list[Path] = []
    for contig in contigs:
        vcf = resolve_common_vcf_path(config.common_vcf_dir, contig, config.common_vcf_pattern)
        ensure_vcf_index(vcf)
        ensure_vcf_has_sample(vcf, config.sample_id)
        part = out_dir / f"{config.sample_id}.{contig}.{role}.vcf.gz"
        part_paths.append(part)
        if part.exists() and has_vcf_index(part) and config.resume and not config.overwrite:
            try:
                validate_haplotag_vcf_summary(part, config.sample_id)
                continue
            except RuntimeError:
                pass
        for path in (part, Path(str(part) + ".tbi"), Path(str(part) + ".csi")):
            path.unlink(missing_ok=True)
        run_sample_vcf_extract(config, vcf, part, commands)
    concat_cmd = [config.bcftools_bin, "concat", "-a", "-Oz", "-o", str(combined), *map(str, part_paths)]
    commands.append(shell_join(concat_cmd))
    subprocess.run(concat_cmd, check=True)
    index_vcf(config, combined, commands)
    combined_summary = validate_haplotag_vcf_summary(combined, config.sample_id)
    for part in part_paths:
        for path in (part, Path(str(part) + ".tbi"), Path(str(part) + ".csi")):
            path.unlink(missing_ok=True)
    write_json(
        metadata_path,
        {
            "schema": "mfm_haplotag_vcf_from_common_v2",
            "role": role,
            "sample_id": config.sample_id,
            "common_vcf_dir": str(config.common_vcf_dir),
            "common_vcf_pattern": config.common_vcf_pattern,
            "filter_order": "sample_subset_then_biallelic_type_het_filter",
            "filter_expr": config.common_vcf_filter_expr,
            "vcf_types": config.haplotag_vcf_types or DEFAULT_HAPLOTAG_VCF_TYPES,
            "biallelic_policy": "bcftools_view_m2_M2",
            "reference_role": "cram_decode_and_haplotag",
            "reference_fasta": str(config.reference_fasta or ""),
            "contigs": contigs,
            "vcf": str(combined),
            "vcf_summary": combined_summary,
            "size_bytes": combined.stat().st_size,
        },
    )
    return combined


def run_sample_vcf_extract(config: FullPreprocessConfig, source_vcf: Path, output_vcf: Path, commands: list[str]) -> None:
    cmd1 = [config.bcftools_bin, "view", "-s", config.sample_id, "-Ou", str(source_vcf)]
    cmd2 = [
        config.bcftools_bin,
        "view",
        "-m2",
        "-M2",
        "-v",
        config.haplotag_vcf_types or DEFAULT_HAPLOTAG_VCF_TYPES,
        "-i",
        config.common_vcf_filter_expr,
        "-Oz",
        "-o",
        str(output_vcf),
        "-",
    ]
    run_bcftools_view_pipe(cmd1, cmd2, commands, "bcftools haplotag VCF 추출 실패")
    index_vcf(config, output_vcf, commands)
    validate_haplotag_vcf_summary(output_vcf, config.sample_id)


def run_bcftools_pipeline(pipeline: list[list[str]], commands: list[str], error_label: str) -> None:
    if len(pipeline) < 2:
        raise ValueError("bcftools pipeline requires at least two commands")
    commands.append(" | ".join(shell_join(cmd) for cmd in pipeline))
    processes: list[subprocess.Popen[bytes]] = []
    previous_stdout = None
    for index, cmd in enumerate(pipeline):
        is_last = index == len(pipeline) - 1
        process = subprocess.Popen(
            cmd,
            stdin=previous_stdout,
            stdout=subprocess.DEVNULL if is_last else subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if previous_stdout is not None:
            previous_stdout.close()
        previous_stdout = process.stdout
        processes.append(process)
    stderr_by_step = [""] * len(processes)
    for index in range(len(processes) - 1, -1, -1):
        stderr = processes[index].communicate()[1]
        stderr_by_step[index] = stderr.decode(errors="replace")
    return_codes = [process.returncode for process in processes]
    if any(code != 0 for code in return_codes):
        details = " ".join(
            f"p{idx + 1}={code} stderr={stderr_by_step[idx]}"
            for idx, code in enumerate(return_codes)
        )
        raise RuntimeError(f"{error_label}: {details}")


def run_bcftools_view_pipe(cmd1: list[str], cmd2: list[str], commands: list[str], error_label: str) -> None:
    run_bcftools_pipeline([cmd1, cmd2], commands, error_label)

def index_vcf(config: FullPreprocessConfig, vcf: Path, commands: list[str]) -> None:
    index_cmd = [config.bcftools_bin, "index", "--tbi", "--force", str(vcf)]
    commands.append(shell_join(index_cmd))
    subprocess.run(index_cmd, check=True)

def write_reference_stream(process: subprocess.Popen, reference: Path, contig: str) -> None:
    if process.stdin is None:
        raise RuntimeError("bcftools consensus stdin is not available")
    with pysam.FastaFile(str(reference)) as fasta:
        length = fasta.get_reference_length(contig)
        process.stdin.write(f">{contig}\n".encode())
        for start in range(0, length, 1 << 20):
            seq = fasta.fetch(contig, start, min(length, start + (1 << 20))).upper()
            for offset in range(0, len(seq), 60):
                process.stdin.write(seq[offset:offset + 60].encode() + b"\n")
    process.stdin.close()


def resolve_common_vcf_path(common_vcf_dir: Path, contig: str, pattern: str = "") -> Path:
    candidates: list[Path] = []
    if pattern:
        candidates.append(common_vcf_dir / pattern.format(contig=contig, chrom=contig, chromosome=contig))
    candidates.extend(common_vcf_dir / item.format(contig=contig) for item in COMMON_VCF_PATTERN_CANDIDATES)
    seen: set[Path] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"common VCF를 찾을 수 없습니다: dir={common_vcf_dir}, contig={contig}")


def ensure_vcf_has_sample(vcf: Path, sample_id: str) -> None:
    with pysam.VariantFile(str(vcf)) as handle:
        if sample_id not in handle.header.samples:
            raise ValueError(f"common VCF에 sample이 없습니다: sample={sample_id}, vcf={vcf}")


def ensure_vcf_index(vcf: Path) -> None:
    if not has_vcf_index(vcf):
        raise FileNotFoundError(f"VCF index가 없습니다: {vcf}.tbi 또는 .csi")


def has_vcf_index(vcf: Path) -> bool:
    return Path(str(vcf) + ".tbi").exists() or Path(str(vcf) + ".csi").exists()


def ensure_alignment_index(samtools_bin: str, alignment: Path, commands: list[str]) -> None:
    if Path(str(alignment) + ".bai").exists() or Path(str(alignment) + ".csi").exists() or Path(str(alignment) + ".crai").exists():
        return
    index_cmd = [samtools_bin, "index", str(alignment)]
    commands.append(shell_join(index_cmd))
    subprocess.run(index_cmd, check=True)


def ensure_hap_input_dirs(config: FullPreprocessConfig) -> None:
    for name in ("hap1_fasta_dir", "hap2_fasta_dir"):
        require_path(getattr(config, name), name)
    if sequence_policy_needs_chain(config):
        for name in ("hap1_chain_dir", "hap2_chain_dir"):
            require_path(getattr(config, name), name)


def require_path(path: Path | None, name: str) -> Path:
    if path is None:
        raise ValueError(f"{name} is required unless --auto-prepare-hap-fasta supplies it")
    return path


def cleanup_model_input_work(config: FullPreprocessConfig, raw_output_root: Path, work_root: Path) -> list[str]:
    targets = [raw_output_root, work_root / "haplotagged", work_root / "common_vcf"]
    if config.auto_prepare_hap_fasta:
        targets.append(work_root / "individual_fasta")
    removed: list[str] = []
    for target in targets:
        if target.exists():
            shutil.rmtree(target)
            removed.append(str(target))
    return removed


def effective_source_vcfs(config: FullPreprocessConfig) -> list[Path]:
    if config.source_vcfs:
        return config.source_vcfs
    if config.phased_vcf is not None:
        return [config.phased_vcf]
    return []


def file_md5(path: Path, block_size: int = 8 << 20) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect_alignment_phase_tags(alignment: Path, contigs: list[str], reference_fasta: Path | None = None, max_reads: int = 50_000) -> dict[str, Any]:
    inspected = 0
    hp_present = 0
    ps_present = 0
    ps_nonzero = 0
    ps_zero = 0
    ps_invalid = 0
    kwargs = {"reference_filename": str(reference_fasta)} if reference_fasta is not None else {}
    with pysam.AlignmentFile(str(alignment), "rb", **kwargs) as bam:
        for contig in contigs:
            bam_contig = contig if contig in bam.references else contig.removeprefix("chr")
            if bam_contig not in bam.references:
                continue
            for read in bam.fetch(bam_contig):
                if read.is_unmapped or read.is_secondary:
                    continue
                inspected += 1
                if read.has_tag("HP"):
                    hp_present += 1
                if read.has_tag("PS"):
                    ps_present += 1
                    ps_value = read_tag_int(read, "PS")
                    if ps_value is None:
                        ps_invalid += 1
                    elif ps_value > 0:
                        ps_nonzero += 1
                    else:
                        ps_zero += 1
                if inspected >= max_reads:
                    break
            if inspected >= max_reads:
                break
    return {
        "alignment": str(alignment),
        "inspected_read_count": inspected,
        "hp_present_count": hp_present,
        "ps_present_count": ps_present,
        "ps_nonzero_count": ps_nonzero,
        "ps_zero_count": ps_zero,
        "ps_invalid_count": ps_invalid,
        "hp_present_fraction": (hp_present / inspected) if inspected else 0.0,
        "ps_present_fraction": (ps_present / inspected) if inspected else 0.0,
        "ps_nonzero_fraction": (ps_nonzero / inspected) if inspected else 0.0,
    }


def build_records(config: FullPreprocessConfig, contigs: list[str]) -> list[dict[str, str]]:
    rows = []
    chain_required = sequence_policy_needs_chain(config)
    for contig in contigs:
        hap1_chain = ""
        hap2_chain = ""
        if chain_required:
            hap1_chain = str(first_existing(require_path(config.hap1_chain_dir, "hap1_chain_dir"), contig, [".chain"]))
            hap2_chain = str(first_existing(require_path(config.hap2_chain_dir, "hap2_chain_dir"), contig, [".chain"]))
        rows.append(
            {
                "sample": config.sample_id,
                "contig": contig,
                "hap1_fasta": str(first_existing(require_path(config.hap1_fasta_dir, "hap1_fasta_dir"), contig, [".fa", ".fasta", ".fa.gz", ".fasta.gz"])),
                "hap1_chain": hap1_chain,
                "hap2_fasta": str(first_existing(require_path(config.hap2_fasta_dir, "hap2_fasta_dir"), contig, [".fa", ".fasta", ".fa.gz", ".fasta.gz"])),
                "hap2_chain": hap2_chain,
                "bam": str(config.alignment),
                "window_index_dir": str(config.window_index_dir or ""),
                "read_coordinate_mode": config.read_coordinate_mode,
                "hap1_alignment": str(config.hap1_alignment or ""),
                "hap2_alignment": str(config.hap2_alignment or ""),
            }
        )
    return rows


def first_existing(directory: Path, contig: str, suffixes: list[str]) -> Path:
    candidates = [directory / f"{contig}{suffix}" for suffix in suffixes]
    candidates.extend(directory / contig / f"{contig}{suffix}" for suffix in suffixes)
    for path in candidates:
        if path.exists():
            return path
    shown = ", ".join(str(path) for path in candidates[:4])
    raise FileNotFoundError(f"missing {contig} file under {directory}; first candidates: {shown}")


def expand_contigs(values: list[str]) -> list[str]:
    expanded: list[str] = []
    for value in values:
        for part in str(value).split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                left, right = part.split("-", 1)
                left_i = int(normalize_contig(left).removeprefix("chr"))
                right_i = int(normalize_contig(right).removeprefix("chr"))
                expanded.extend(f"chr{i}" for i in range(left_i, right_i + 1))
            else:
                expanded.append(normalize_contig(part))
    seen = set()
    out = []
    for contig in expanded:
        if contig not in seen:
            seen.add(contig)
            out.append(contig)
    return out


def normalize_contig(value: str) -> str:
    value = value.strip()
    return value if value.startswith("chr") else f"chr{value}"


def resolve_sample_root(config: FullPreprocessConfig) -> Path:
    if config.output_sample_root is not None:
        return config.output_sample_root
    if config.dataset_root is None:
        raise ValueError("either output_sample_root or dataset_root is required")
    return config.dataset_root / "final" / config.dataset_role / config.read_type / config.project_id / config.sample_id


def prepare_sample_root(sample_root: Path) -> None:
    # `tracks` 는 여기서 만들지 않는다 — bigwig 이 있을 때만 write_rna_track_manifest 가
    # 스스로 mkdir 한다. STAR 개인좌표 재정렬(snv_indel) 경로에는 bigwig 이 없어서,
    # 미리 만들면 10명 전부 **빈 tracks/** 가 남고 "RNA track 이 여기 있나" 하고 찾게 된다.
    # 그 경로의 실제 위치는 rna_personal/<phase>/HAP{1,2}/coverage.zarr 다(sample.json 참조).
    for name in ("metadata", "compact", "qc", "logs", "work"):
        (sample_root / name).mkdir(parents=True, exist_ok=True)


def write_full_metadata(
    sample_root: Path,
    config: FullPreprocessConfig,
    summary: dict[str, Any],
    records: list[dict[str, str]],
    commands: list[str],
) -> None:
    metadata_dir = sample_root / "metadata"
    logs_dir = sample_root / "logs"
    rna_personal = detect_rna_personal(sample_root)
    write_json(
        metadata_dir / "sample.json",
        {
            "schema": "mfm_integrated_sample_v1",
            "sample_id": config.sample_id,
            "source_sample_id": config.source_sample_id or config.sample_id,
            "dataset_role": config.dataset_role,
            "read_type": config.read_type,
            "project_id": config.project_id,
            "project_label": config.project_label,
            "platform": config.platform,
            "reference_build": config.reference_build,
            "sample_status": config.sample_status,
            "disease": config.disease,
            # bigwig(구 경로) 또는 rna_personal(개인좌표 STAR 경로) 중 하나라도 있으면 RNA 가 있다.
            # 경로는 rna_personal 키에 그대로 실어 보내 로더가 추측하지 않게 한다(키는 추가만 — 기존 소비자 무영향).
            "has_rna": bool(config.rna_bigwigs) or bool(rna_personal),
            "rna_personal": rna_personal,
            "layout_name": config.layout_name,
            "hap_read_mode": config.hap_read_mode,
            "compact_mode": config.compact_mode,
            "reference_roles": reference_roles(config),
            "sequence_variant_policy": config.sequence_variant_policy,
            "project_policy_guard": {
                "expected_sequence_variant_policy": PROJECT_SEQUENCE_POLICY.get(config.project_id, ""),
                "expected_dataset_role": PROJECT_ROLE_POLICY.get(config.project_id, ""),
                "expected_read_type": PROJECT_READ_TYPE_POLICY.get(config.project_id, ""),
            },
        },
    )
    write_json(metadata_dir / "preprocess.json", summary)
    write_source_manifest(metadata_dir / "source_manifest.tsv", config)
    write_output_manifest(metadata_dir / "output_manifest.tsv", summary)
    write_record_tsv(metadata_dir / "full_shared_hap.records.tsv", records)
    write_commands(metadata_dir / "commands.sh", commands)
    write_json(logs_dir / "timing.json", {"timing": summary["timing"]})
    if config.rna_bigwigs:
        write_rna_track_manifest(sample_root / "tracks" / "rna_bigwig", config.sample_id, config.rna_bigwigs)


def write_source_manifest(path: Path, config: FullPreprocessConfig) -> None:
    rows = [
        {"kind": "alignment_input", "path": str(config.alignment), "source_url": "", "md5": ""},
        *[{"kind": "source_wgs", "path": str(item), "source_url": "", "md5": ""} for item in config.source_wgs],
        *[{"kind": "source_vcf", "path": str(item), "source_url": "", "md5": ""} for item in config.source_vcfs],
        *[{"kind": "rna_bigwig", "path": str(item), "source_url": "", "md5": ""} for item in config.rna_bigwigs],
    ]
    if config.phased_vcf is not None:
        rows.append({"kind": "phased_vcf_for_haplotagging", "path": str(config.phased_vcf), "source_url": "", "md5": ""})
    if config.common_vcf_dir is not None:
        rows.append({"kind": "common_phased_vcf_dir", "path": str(config.common_vcf_dir), "source_url": "", "md5": ""})
    if config.reference_fasta is not None:
        rows.append({"kind": "reference_fasta_for_cram_haplotag", "path": str(config.reference_fasta), "source_url": "", "md5": ""})
    if effective_sequence_reference(config) is not None:
        rows.append({"kind": "reference_fasta_for_personal_sequence", "path": str(effective_sequence_reference(config)), "source_url": "", "md5": ""})
    if config.window_index_dir is not None:
        rows.append({"kind": "window_index_dir", "path": str(config.window_index_dir), "source_url": "", "md5": ""})
    if config.hap1_alignment is not None:
        rows.append({"kind": "hap1_personal_alignment", "path": str(config.hap1_alignment), "source_url": "", "md5": ""})
    if config.hap2_alignment is not None:
        rows.append({"kind": "hap2_personal_alignment", "path": str(config.hap2_alignment), "source_url": "", "md5": ""})
    write_tsv(path, rows, ["kind", "path", "source_url", "md5"])


def write_output_manifest(path: Path, summary: dict[str, Any]) -> None:
    rows = []
    compact_mode = summary.get("compact_mode", "two_step")
    shared_kind = "direct_compact_4bit_contig" if compact_mode == "direct" else "raw_shared_hap_zarr"
    records = summary.get("shared_hap_manifest", {}).get("records", {})
    for key, result in sorted(records.items()):
        rows.append({"kind": shared_kind, "id": key, "path": result.get("zarr", ""), "status": result.get("status", "")})
    for contig_summary in summary.get("compact_summary", {}).get("contigs", []):
        rows.append(
            {
                "kind": "compact_full_shared_hap_4bit",
                "id": contig_summary.get("contig", ""),
                "path": contig_summary.get("target", ""),
                "status": contig_summary.get("status", ""),
            }
        )
    write_tsv(path, rows, ["kind", "id", "path", "status"])


def write_rna_track_manifest(tracks_dir: Path, sample_id: str, bigwigs: list[Path]) -> None:
    tracks_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for source in bigwigs:
        track_id = source.stem.replace(".", "_")
        rows.append(
            {
                "sample_id": sample_id,
                "track_id": track_id,
                "name": f"RNA-seq:{source.stem}",
                "assay_title": "RNA-seq",
                "strand": "unstranded",
                "ontology_curie": "NaN",
                "biosample_life_stage": "NaN",
                "tissue": "NaN",
                "cell_line": "NaN",
                "cell_type": "NaN",
                "endedness": "NaN",
                "genetically_modified": "False",
                "nonzero_mean": "NaN",
                "source_path": str(source),
            }
        )
    write_tsv(
        tracks_dir / "track_manifest.tsv",
        rows,
        [
            "sample_id",
            "track_id",
            "name",
            "assay_title",
            "strand",
            "ontology_curie",
            "biosample_life_stage",
            "tissue",
            "cell_line",
            "cell_type",
            "endedness",
            "genetically_modified",
            "nonzero_mean",
            "source_path",
        ],
    )


def write_record_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    write_tsv(path, rows, list(RECORD_COLUMNS))


def write_commands(path: Path, commands: list[str]) -> None:
    path.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(commands) + "\n", encoding="utf-8")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    try:
        partial.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)


def read_manifest_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    for index, row in enumerate(rows):
        row["_manifest_order"] = str(index)
    return rows


def config_from_row(row: dict[str, str], manifest_config: FullManifestConfig) -> FullPreprocessConfig:
    dataset_root = Path(row["dataset_root"]) if row.get("dataset_root") else manifest_config.dataset_root
    output_sample_root = Path(row["output_sample_root"]) if row.get("output_sample_root") else None
    return FullPreprocessConfig(
        sample_id=row["sample_id"],
        source_sample_id=row.get("source_sample_id", ""),
        dataset_role=row["dataset_role"],
        read_type=row["read_type"],
        project_id=row["project_id"],
        project_label=row.get("project_label", ""),
        platform=row.get("platform", ""),
        reference_build=row.get("reference_build", "hg38") or "hg38",
        sample_status=row.get("sample_status", "healthy") or "healthy",
        disease=row.get("disease", "NaN") or "NaN",
        hap1_fasta_dir=parse_optional_path(row.get("hap1_fasta_dir", "")),
        hap1_chain_dir=parse_optional_path(row.get("hap1_chain_dir", "")),
        hap2_fasta_dir=parse_optional_path(row.get("hap2_fasta_dir", "")),
        hap2_chain_dir=parse_optional_path(row.get("hap2_chain_dir", "")),
        alignment=Path(row["alignment_path"]),
        contigs=parse_contig_list(row["contigs"]),
        output_sample_root=output_sample_root,
        dataset_root=dataset_root,
        work_root=Path(row["work_root"]) if row.get("work_root") else None,
        window_index_dir=parse_optional_path(row.get("window_index_dir", "")),
        read_coordinate_mode=row.get("read_coordinate_mode", "") or "personal_realigned",
        hap1_alignment=parse_optional_path(row.get("hap1_alignment", "")),
        hap2_alignment=parse_optional_path(row.get("hap2_alignment", "")),
        storage_policy=row.get("storage_policy") or manifest_config.storage_policy,
        jobs=int(row.get("jobs") or manifest_config.jobs),
        windows_per_chunk=int(row.get("windows_per_chunk") or 32),
        compact_mode=row.get("compact_mode") or "direct",
        validation_level=row.get("validation_level") or "full",
        compact_compressor=row.get("compact_compressor") or DEFAULT_COMPACT_COMPRESSOR,
        compact_compression_level=int(row.get("compact_compression_level") or DEFAULT_COMPACT_COMPRESSION_LEVEL),
        schedule=row.get("schedule") or "longest-first",
        window_size=int(row.get("window_size") or (1 << 20)),
        offset_unit=parse_optional_int(row.get("offset_unit", "")),
        seed=int(row.get("seed") or 42),
        min_mapq=int(row.get("min_mapq") or 30),
        stage70_filter_preset=row.get("stage70_filter_preset", "") or "auto",
        max_windows=parse_optional_int(row.get("max_windows", "")),
        overwrite=manifest_config.overwrite or parse_bool(row.get("overwrite", "")),
        resume=manifest_config.resume and not parse_bool(row.get("no_resume", "")),
        allow_untagged=parse_bool(row.get("allow_untagged", "")),
        haplotag_if_missing=parse_bool(row.get("haplotag_if_missing", "")),
        allow_longread_vcf_haplotag=False,
        allow_unphased_longread=parse_bool(row.get("allow_unphased_longread", "")),
        auto_prepare_hap_fasta=parse_bool(row.get("auto_prepare_hap_fasta", "")),
        phased_vcf=parse_optional_path(row.get("phased_vcf", "")),
        common_vcf_dir=parse_optional_path(row.get("common_vcf_dir", "")),
        common_vcf_pattern=row.get("common_vcf_pattern", "") or "",
        common_vcf_filter_expr=row.get("common_vcf_filter_expr", "") or DEFAULT_SAMPLE_VCF_FILTER_EXPR,
        consensus_filter_expr=row.get("consensus_filter_expr", "") or DEFAULT_CONSENSUS_FILTER_EXPR,
        haplotag_vcf_types=row.get("haplotag_vcf_types", "") or DEFAULT_HAPLOTAG_VCF_TYPES,
        sequence_variant_policy=row.get("sequence_variant_policy", "") or DEFAULT_SEQUENCE_VARIANT_POLICY,
        reference_fasta=parse_optional_path(row.get("reference_fasta", "")),
        sequence_reference_fasta=parse_optional_path(row.get("sequence_reference_fasta", "")),
        phase_inspect_max_reads=int(row.get("phase_inspect_max_reads") or 50000),
        synthesize_missing_ps=parse_bool(row.get("synthesize_missing_ps", "")),
        whatshap_bin=row.get("whatshap_bin", "whatshap") or "whatshap",
        samtools_bin=row.get("samtools_bin", "samtools") or "samtools",
        bcftools_bin=row.get("bcftools_bin", DEFAULT_BCFTOOLS) or DEFAULT_BCFTOOLS,
        source_wgs=parse_paths(row.get("source_wgs", "")),
        source_vcfs=parse_paths(row.get("source_vcfs", "")),
        rna_bigwigs=parse_paths(row.get("rna_bigwigs", "")),
    )


def parse_contig_list(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_paths(value: str) -> list[Path]:
    return [Path(item) for item in value.split(";") if item]


def parse_optional_path(value: str) -> Path | None:
    return Path(value) if value else None


def parse_optional_int(value: str) -> int | None:
    return int(value) if value else None


def parse_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def shell_join(values: list[str]) -> str:
    return " ".join(shlex.quote(str(value)) for value in values)


def print_manifest_result(result: dict[str, Any]) -> None:
    print(
        "[FULL] "
        f"sample={result['sample_id']} status={result['status']} "
        f"elapsed={result['elapsed_seconds']}s size={result['size_human']}"
    )


def write_status_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    write_tsv(
        path,
        rows,
        ["manifest_order", "sample_id", "status", "elapsed_seconds", "elapsed_hms", "output_sample_root", "size_bytes", "size_human", "error"],
    )


def format_hms(seconds: float) -> str:
    total = int(round(seconds))
    return f"{total // 3600:02d}:{(total % 3600) // 60:02d}:{total % 60:02d}"
