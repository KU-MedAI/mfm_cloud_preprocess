from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


DEFAULT_WHATSHAP = "whatshap"
DEFAULT_SAMTOOLS = "samtools"
DEFAULT_BCFTOOLS = "bcftools"
DEFAULT_SAMPLE_VCF_FILTER = 'GT="het"'

MANIFEST_COLUMNS = [
    "sample_id",
    "source_sample_id",
    "dataset_role",
    "read_type",
    "project_id",
    "project_label",
    "platform",
    "reference_build",
    "sample_status",
    "disease",
    "hap1_fasta_dir",
    "hap1_chain_dir",
    "hap2_fasta_dir",
    "hap2_chain_dir",
    "alignment_path",
    "read_coordinate_mode",
    "hap1_alignment",
    "hap2_alignment",
    "contigs",
    "dataset_root",
    "output_sample_root",
    "work_root",
    "window_index_dir",
    "storage_policy",
    "jobs",
    "windows_per_chunk",
    "compact_mode",
    "validation_level",
    "compact_compressor",
    "compact_compression_level",
    "schedule",
    "window_size",
    "offset_unit",
    "seed",
    "min_mapq",
    "max_windows",
    "overwrite",
    "no_resume",
    "allow_untagged",
    "haplotag_if_missing",
    "allow_longread_vcf_haplotag",
    "allow_unphased_longread",
    "auto_prepare_hap_fasta",
    "phased_vcf",
    "common_vcf_dir",
    "common_vcf_pattern",
    "common_vcf_filter_expr",
    "consensus_filter_expr",
    "haplotag_vcf_types",
    "sequence_variant_policy",
    "reference_fasta",
    "sequence_reference_fasta",
    "phase_inspect_max_reads",
    "synthesize_missing_ps",
    "whatshap_bin",
    "samtools_bin",
    "bcftools_bin",
    "source_wgs",
    "source_vcfs",
    "rna_bigwigs",
]

SUMMARY_COLUMNS = [
    "sample_id",
    "alignment_path",
    "alignment_index_path",
    "cram_size_bytes",
    "cram_size_human",
    "selection_rank",
    "selection_reason",
]


@dataclass
class IgsrOntManifestConfig:
    raw_root: Path
    output_manifest: Path
    dataset_root: Path
    common_vcf_dir: Path
    reference_fasta: Path
    contigs: str = "chr1-chr22"
    dataset_role: str = "pre_training"
    read_type: str = "long_read"
    project_id: str = "igsr_ont_wgs_only"
    project_label: str = "IGSR ONT WGS-only"
    platform: str = "OxfordNanopore"
    reference_build: str = "hg38"
    sample_status: str = "healthy"
    disease: str = "NaN"
    storage_policy: str = "model_input_only"
    jobs: int = 32
    windows_per_chunk: int = 32
    compact_mode: str = "direct"
    read_coordinate_mode: str = "personal_realigned"
    hap1_alignment_pattern: str = "{sample}/realigned/HAP1/{sample}.personal.HAP1.cram"
    hap2_alignment_pattern: str = "{sample}/realigned/HAP2/{sample}.personal.HAP2.cram"
    haplotagged_alignment_pattern: str = "{sample}/haplotagged/{sample}.clair3_longread.haplotagged.bam"
    window_index_pattern: str = "{sample}/windows"
    validation_level: str = "full"
    compact_compressor: str = "zstd"
    compact_compression_level: int = 5
    schedule: str = "longest-first"
    window_size: int = 1 << 20
    offset_unit: int | None = None
    seed: int = 42
    min_mapq: int = 30
    max_windows: int | None = None
    phase_inspect_max_reads: int = 50_000
    whatshap_bin: str = DEFAULT_WHATSHAP
    samtools_bin: str = DEFAULT_SAMTOOLS
    bcftools_bin: str = DEFAULT_BCFTOOLS
    common_vcf_pattern: str = ""
    common_vcf_filter_expr: str = DEFAULT_SAMPLE_VCF_FILTER
    consensus_filter_expr: str = '(FILTER="PASS" || FILTER=".")'
    haplotag_vcf_types: str = "snps,indels"
    sequence_variant_policy: str = "snv_indel"
    sequence_reference_fasta: Path | None = None
    synthesize_missing_ps: bool = False
    allow_longread_vcf_haplotag: bool = False
    allow_unphased_longread: bool = False
    samples: list[str] = field(default_factory=list)
    sample_file: Path | None = None
    pilot_size: int | None = None
    limit: int | None = None
    summary_tsv: Path | None = None
    overwrite_manifest: bool = False


@dataclass(frozen=True)
class IgsrOntSample:
    sample_id: str
    alignment_path: Path
    alignment_index_path: Path | None
    cram_size_bytes: int


def build_igsr_ont_manifest(config: IgsrOntManifestConfig) -> dict[str, Any]:
    if config.sequence_reference_fasta is None:
        raise ValueError(
            "--sequence-reference-fasta is required for IGSR ONT manifests; "
            "use UCSC hg38 for personal HAP FASTA consensus. "
            "--reference-fasta is reserved for CRAM decoding/WhatsHap, e.g. Vienna."
        )
    if config.read_coordinate_mode not in {"reference_projected", "personal_realigned"}:
        raise ValueError("--read-coordinate-mode must be reference_projected or personal_realigned")
    if config.read_coordinate_mode == "personal_realigned" and config.compact_mode != "direct":
        raise ValueError("personal_realigned manifest requires compact_mode=direct")
    samples = scan_igsr_ont_samples(config.raw_root)
    selected = select_samples(samples, config)
    if config.output_manifest.exists() and not config.overwrite_manifest:
        raise FileExistsError(f"기존 manifest가 있습니다: {config.output_manifest}. --overwrite-manifest를 지정하세요")
    config.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    rows = [manifest_row(sample, config) for sample in selected]
    write_tsv(config.output_manifest, rows, MANIFEST_COLUMNS)
    summary_tsv = config.summary_tsv or config.output_manifest.with_suffix(config.output_manifest.suffix + ".samples.tsv")
    write_tsv(summary_tsv, summary_rows(selected), SUMMARY_COLUMNS)
    return {
        "manifest": str(config.output_manifest),
        "summary_tsv": str(summary_tsv),
        "raw_root": str(config.raw_root),
        "discovered_sample_count": len(samples),
        "selected_sample_count": len(selected),
        "selected_samples": [sample.sample_id for sample in selected],
        "total_selected_cram_size_bytes": sum(sample.cram_size_bytes for sample in selected),
        "total_selected_cram_size_human": human_bytes(sum(sample.cram_size_bytes for sample in selected)),
    }


def scan_igsr_ont_samples(raw_root: Path) -> list[IgsrOntSample]:
    if not raw_root.exists():
        raise FileNotFoundError(f"raw root가 없습니다: {raw_root}")
    out: list[IgsrOntSample] = []
    for sample_dir in sorted(path for path in raw_root.iterdir() if path.is_dir()):
        sample_id = sample_dir.name
        wgs_dir = sample_dir / "wgs"
        if not wgs_dir.exists():
            continue
        alignment = first_existing_alignment(wgs_dir, sample_id)
        if alignment is None:
            continue
        index = first_existing_index(alignment)
        out.append(IgsrOntSample(sample_id, alignment, index, alignment.stat().st_size))
    if not out:
        raise FileNotFoundError(f"ONT WGS CRAM을 찾지 못했습니다: {raw_root}")
    return out


def first_existing_alignment(wgs_dir: Path, sample_id: str) -> Path | None:
    preferred = [
        wgs_dir / f"{sample_id}.hg38.cram",
        wgs_dir / f"{sample_id}.ont.hg38.cram",
        wgs_dir / f"{sample_id}.cram",
    ]
    for path in preferred:
        if path.exists():
            return path
    candidates = sorted(wgs_dir.glob("*.cram")) + sorted(wgs_dir.glob("*.bam"))
    return candidates[0] if candidates else None


def first_existing_index(alignment: Path) -> Path | None:
    for suffix in (".crai", ".bai", ".csi"):
        candidate = Path(str(alignment) + suffix)
        if candidate.exists():
            return candidate
    sibling = alignment.with_suffix(alignment.suffix + ".crai")
    if sibling.exists():
        return sibling
    return None


def select_samples(samples: list[IgsrOntSample], config: IgsrOntManifestConfig) -> list[IgsrOntSample]:
    by_id = {sample.sample_id: sample for sample in samples}
    requested = list(config.samples)
    if config.sample_file is not None:
        requested.extend(read_sample_file(config.sample_file))
    if requested:
        missing = [sample_id for sample_id in requested if sample_id not in by_id]
        if missing:
            raise ValueError(f"raw root에서 찾지 못한 sample: {missing[:10]}")
        selected = [by_id[sample_id] for sample_id in requested]
    else:
        selected = list(samples)
        if config.pilot_size:
            selected = select_quantile_pilot(selected, config.pilot_size)
        if config.limit:
            selected = selected[: config.limit]
    return selected


def read_sample_file(path: Path) -> list[str]:
    values: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        values.append(line.split()[0])
    return values


def select_quantile_pilot(samples: list[IgsrOntSample], pilot_size: int) -> list[IgsrOntSample]:
    if pilot_size <= 0:
        raise ValueError("pilot_size는 양수여야 합니다")
    ordered = sorted(samples, key=lambda item: (item.cram_size_bytes, item.sample_id))
    if pilot_size >= len(ordered):
        return ordered
    if pilot_size == 1:
        return [ordered[len(ordered) // 2]]
    indices = [round(index * (len(ordered) - 1) / (pilot_size - 1)) for index in range(pilot_size)]
    seen: set[str] = set()
    out: list[IgsrOntSample] = []
    for index in indices:
        sample = ordered[index]
        if sample.sample_id not in seen:
            seen.add(sample.sample_id)
            out.append(sample)
    cursor = 0
    while len(out) < pilot_size and cursor < len(ordered):
        sample = ordered[cursor]
        if sample.sample_id not in seen:
            seen.add(sample.sample_id)
            out.append(sample)
        cursor += 1
    return out


def manifest_row(sample: IgsrOntSample, config: IgsrOntManifestConfig) -> dict[str, str]:
    personal_mode = config.read_coordinate_mode == "personal_realigned"
    final_root = final_sample_root(sample.sample_id, config)
    hap1_dir = final_root / "individual_fasta" / "HAP1"
    hap2_dir = final_root / "individual_fasta" / "HAP2"
    chain_required = config.sequence_variant_policy != "snp_only"
    alignment_path = haplotagged_alignment_path(sample.sample_id, config) if personal_mode else str(sample.alignment_path)
    return {
        "sample_id": sample.sample_id,
        "source_sample_id": sample.sample_id,
        "dataset_role": config.dataset_role,
        "read_type": config.read_type,
        "project_id": config.project_id,
        "project_label": config.project_label,
        "platform": config.platform,
        "reference_build": config.reference_build,
        "sample_status": config.sample_status,
        "disease": config.disease,
        "hap1_fasta_dir": str(hap1_dir) if personal_mode else "",
        "hap1_chain_dir": str(hap1_dir) if personal_mode and chain_required else "",
        "hap2_fasta_dir": str(hap2_dir) if personal_mode else "",
        "hap2_chain_dir": str(hap2_dir) if personal_mode and chain_required else "",
        "alignment_path": alignment_path,
        "read_coordinate_mode": config.read_coordinate_mode,
        "hap1_alignment": personal_alignment_path(sample.sample_id, config, 1),
        "hap2_alignment": personal_alignment_path(sample.sample_id, config, 2),
        "contigs": config.contigs,
        "dataset_root": str(config.dataset_root),
        "output_sample_root": "",
        "work_root": "",
        "window_index_dir": window_index_path(sample.sample_id, config) if personal_mode else "",
        "storage_policy": config.storage_policy,
        "jobs": str(config.jobs),
        "windows_per_chunk": str(config.windows_per_chunk),
        "compact_mode": config.compact_mode,
        "validation_level": config.validation_level,
        "compact_compressor": config.compact_compressor,
        "compact_compression_level": str(config.compact_compression_level),
        "schedule": config.schedule,
        "window_size": str(config.window_size),
        "offset_unit": "" if config.offset_unit is None else str(config.offset_unit),
        "seed": str(config.seed),
        "min_mapq": str(config.min_mapq),
        "max_windows": "" if config.max_windows is None else str(config.max_windows),
        "overwrite": "false",
        "no_resume": "false",
        "allow_untagged": "false",
        "haplotag_if_missing": "false" if personal_mode else "true",
        "allow_longread_vcf_haplotag": "false",
        "allow_unphased_longread": "true" if config.allow_unphased_longread else "false",
        "auto_prepare_hap_fasta": "false" if personal_mode else "true",
        "phased_vcf": "",
        "common_vcf_dir": str(config.common_vcf_dir),
        "common_vcf_pattern": config.common_vcf_pattern,
        "common_vcf_filter_expr": config.common_vcf_filter_expr,
        "consensus_filter_expr": config.consensus_filter_expr,
        "haplotag_vcf_types": config.haplotag_vcf_types,
        "sequence_variant_policy": config.sequence_variant_policy,
        "reference_fasta": str(config.reference_fasta),
        "sequence_reference_fasta": str(config.sequence_reference_fasta or ""),
        "phase_inspect_max_reads": str(config.phase_inspect_max_reads),
        "synthesize_missing_ps": "true" if config.synthesize_missing_ps else "false",
        "whatshap_bin": config.whatshap_bin,
        "samtools_bin": config.samtools_bin,
        "bcftools_bin": config.bcftools_bin,
        "source_wgs": str(sample.alignment_path),
        "source_vcfs": "",
        "rna_bigwigs": "",
    }


def final_sample_root(sample_id: str, config: IgsrOntManifestConfig) -> Path:
    return config.dataset_root / "final" / config.dataset_role / config.read_type / config.project_id / sample_id


def final_project_path(sample_id: str, config: IgsrOntManifestConfig, pattern: str) -> str:
    rendered = pattern.format(sample=sample_id, hap="", haplotype="")
    return str(config.dataset_root / "final" / config.dataset_role / config.read_type / config.project_id / rendered)


def haplotagged_alignment_path(sample_id: str, config: IgsrOntManifestConfig) -> str:
    return final_project_path(sample_id, config, config.haplotagged_alignment_pattern)


def window_index_path(sample_id: str, config: IgsrOntManifestConfig) -> str:
    return final_project_path(sample_id, config, config.window_index_pattern)


def personal_alignment_path(sample_id: str, config: IgsrOntManifestConfig, hap: int) -> str:
    if config.read_coordinate_mode != "personal_realigned":
        return ""
    pattern = config.hap1_alignment_pattern if hap == 1 else config.hap2_alignment_pattern
    rendered = pattern.format(sample=sample_id, hap=hap, haplotype=f"HAP{hap}")
    return str(config.dataset_root / "final" / config.dataset_role / config.read_type / config.project_id / rendered)


def summary_rows(samples: list[IgsrOntSample]) -> list[dict[str, str]]:
    ordered = sorted(samples, key=lambda item: item.cram_size_bytes)
    rank_by_id = {sample.sample_id: index for index, sample in enumerate(ordered)}
    out: list[dict[str, str]] = []
    for sample in samples:
        rank = rank_by_id[sample.sample_id]
        reason = "selected"
        if len(ordered) > 1:
            fraction = rank / (len(ordered) - 1)
            reason = f"size_quantile_{fraction:.3f}"
        out.append(
            {
                "sample_id": sample.sample_id,
                "alignment_path": str(sample.alignment_path),
                "alignment_index_path": str(sample.alignment_index_path or ""),
                "cram_size_bytes": str(sample.cram_size_bytes),
                "cram_size_human": human_bytes(sample.cram_size_bytes),
                "selection_rank": str(rank),
                "selection_reason": reason,
            }
        )
    return out


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def human_bytes(value: int) -> str:
    size = float(value)
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} TiB"
