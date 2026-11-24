from __future__ import annotations

import json
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pysam

from .cohort_vcf import DEFAULT_BCFTOOLS, expand_chroms, has_vcf_index, run_command
from .sample_vcf import normalize_sequence_variant_policy, policy_vtypes, run_pipeline


@dataclass
class PublicVcfConfig:
    sample: str
    source_sample: str
    source_vcfs: list[Path]
    reference_fasta: Path
    output_dir: Path
    callset: str
    phase_source: str
    chroms: list[str]
    sequence_variant_policy: str = "snv_indel"
    bcftools_bin: str = DEFAULT_BCFTOOLS
    overwrite: bool = False


@dataclass
class VcfStats:
    records: int = 0
    nonref_records: int = 0
    het_records: int = 0
    phased_het_records: int = 0
    ps_present_records: int = 0
    hom_ref_records: int = 0
    multiallelic_records: int = 0
    symbolic_records: int = 0

    @property
    def phased_het_fraction(self) -> float:
        return self.phased_het_records / self.het_records if self.het_records else 1.0


def prepare_public_vcf(config: PublicVcfConfig) -> dict[str, Any]:
    from .public_vcf_parallel import reserve_validation_memory
    started = time.time()
    source_vcfs = [path.resolve() for path in config.source_vcfs]
    reference = config.reference_fasta.resolve()
    output_dir = config.output_dir.resolve()
    chroms = expand_chroms(config.chroms)
    policy = normalize_sequence_variant_policy(config.sequence_variant_policy)
    phase_source = normalize_phase_source(config.phase_source)
    validate_inputs(config, source_vcfs, reference, chroms)
    output_dir.mkdir(parents=True, exist_ok=True)

    seq_vcf = output_dir / f"{config.sample}.seq.vcf.gz"
    tag_vcf = output_dir / f"{config.sample}.tag.vcf.gz"
    sidecar = output_dir / f"{config.sample}.seq.vcf.source.json"
    source_identity = {
        "source_vcfs": [str(path) for path in source_vcfs],
        "source_vcf_sizes": [path.stat().st_size for path in source_vcfs],
        "source_sample": config.source_sample,
        "sample": config.sample,
        "reference_fasta": str(reference),
        "callset": config.callset,
        "phase_source": phase_source,
        "chroms": chroms,
        "sequence_variant_policy": policy,
    }

    if not config.overwrite and outputs_match(sidecar, source_identity, seq_vcf, tag_vcf):
        reserve_validation_memory()
        seq_stats = inspect_vcf(seq_vcf, config.sample)
        tag_stats = inspect_vcf(tag_vcf, config.sample)
        validate_outputs(seq_stats, tag_stats, phase_source)
        return build_payload(source_identity, seq_vcf, tag_vcf, seq_stats, tag_stats, "skip", time.time() - started)

    cleanup_bundle(seq_vcf)
    cleanup_bundle(tag_vcf)
    sidecar.unlink(missing_ok=True)
    try:
        build_sequence_vcf(config, source_vcfs, reference, chroms, policy, seq_vcf)
        reserve_validation_memory()
        build_tag_vcf(config.bcftools_bin, seq_vcf, tag_vcf)
        seq_stats = inspect_vcf(seq_vcf, config.sample)
        tag_stats = inspect_vcf(tag_vcf, config.sample)
        validate_outputs(seq_stats, tag_stats, phase_source)
    except Exception:
        cleanup_bundle(seq_vcf)
        cleanup_bundle(tag_vcf)
        sidecar.unlink(missing_ok=True)
        raise
    payload = build_payload(source_identity, seq_vcf, tag_vcf, seq_stats, tag_stats, "created", time.time() - started)
    write_json_atomic(sidecar, payload)
    return payload


def normalize_phase_source(value: str) -> str:
    normalized = (value or "provided").strip().lower().replace("-", "_").replace(" ", "_")
    if normalized != "provided":
        raise ValueError(f"public VCF phase_source는 현재 provided만 지원합니다: {value}")
    return normalized


def validate_inputs(config: PublicVcfConfig, source_vcfs: list[Path], reference: Path, chroms: list[str]) -> None:
    if not source_vcfs:
        raise ValueError("public VCF 목록이 비어 있습니다")
    if not reference.exists() or not reference.stat().st_size:
        raise FileNotFoundError(f"VCF reference FASTA가 없습니다: {reference}")
    fai = Path(str(reference) + ".fai")
    if not fai.exists() or not fai.stat().st_size:
        raise FileNotFoundError(f"VCF reference index가 없습니다: {fai}")
    if not Path(config.bcftools_bin).exists():
        raise FileNotFoundError(f"bcftools가 없습니다: {config.bcftools_bin}")

    ref_lengths = load_fai_lengths(fai)
    available_contigs: set[str] = set()
    length_mismatch: list[str] = []
    for source_vcf in source_vcfs:
        if not source_vcf.exists() or not source_vcf.stat().st_size:
            raise FileNotFoundError(f"public VCF가 없습니다: {source_vcf}")
        if not has_vcf_index(source_vcf):
            raise FileNotFoundError(f"public VCF index가 없습니다: {source_vcf}.tbi 또는 .csi")
        with pysam.VariantFile(str(source_vcf)) as vcf:
            samples = list(vcf.header.samples)
            if config.source_sample not in samples:
                raise ValueError(f"public VCF에 sample이 없습니다: {source_vcf} sample={config.source_sample}")
            for chrom in chroms:
                if chrom not in vcf.header.contigs:
                    continue
                available_contigs.add(chrom)
                vcf_len = vcf.header.contigs[chrom].length
                if chrom not in ref_lengths:
                    continue
                if vcf_len is not None and int(vcf_len) != ref_lengths[chrom]:
                    length_mismatch.append(f"{source_vcf.name}:{chrom}:vcf={vcf_len},ref={ref_lengths[chrom]}")
    missing = [chrom for chrom in chroms if chrom not in available_contigs or chrom not in ref_lengths]
    if missing:
        raise ValueError(f"VCF/reference contig 누락: {','.join(missing)}")
    if length_mismatch:
        raise ValueError("VCF/reference contig 길이 불일치: " + ";".join(length_mismatch[:5]))


def load_fai_lengths(fai: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with fai.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def build_sequence_vcf(
    config: PublicVcfConfig, source_vcfs: list[Path], reference: Path,
    chroms: list[str], policy: str, output_vcf: Path,
) -> None:
    from .public_vcf_parallel import build_parallel
    build_parallel(config, source_vcfs, reference, chroms, policy, output_vcf, _build_sequence_vcf_serial)


def _build_sequence_vcf_serial(
    config: PublicVcfConfig,
    source_vcfs: list[Path],
    reference: Path,
    chroms: list[str],
    policy: str,
    output_vcf: Path,
) -> None:
    selected = output_vcf.with_name(output_vcf.name + ".selected.part")
    final_part = output_vcf.with_name(output_vcf.name + ".part")
    sorted_part = output_vcf.with_name(output_vcf.name + ".sorted.part")
    selected.unlink(missing_ok=True)
    final_part.unlink(missing_ok=True)
    sorted_part.unlink(missing_ok=True)
    seq_types = policy_vtypes(policy)
    temporary: list[Path] = []
    try:
        source_parts: list[Path] = []
        for index, source_vcf in enumerate(source_vcfs):
            with pysam.VariantFile(str(source_vcf)) as vcf:
                source_chroms = [chrom for chrom in chroms if chrom in vcf.header.contigs]
            if not source_chroms:
                continue
            source_view = output_vcf.with_name(f".{output_vcf.name}.source{index}.view.bcf")
            source_part = output_vcf.with_name(f".{output_vcf.name}.source{index}.bcf")
            source_view.unlink(missing_ok=True)
            source_part.unlink(missing_ok=True)
            run_command([
                config.bcftools_bin,
                "view",
                "-s",
                config.source_sample,
                "-r",
                ",".join(source_chroms),
                "-v",
                seq_types,
                "-c",
                "1",
                "-Ob",
                "-o",
                str(source_view),
                str(source_vcf),
            ])
            # 일부 public VCF는 index가 존재해도 record 좌표 순서를 보장하지 않는다.
            # concat 및 후속 norm이 요구하는 index는 정렬한 임시 BCF에만 만든다.
            run_command([config.bcftools_bin, "sort", "-Ob", "-o", str(source_part), str(source_view)])
            source_view.unlink(missing_ok=True)
            run_command([config.bcftools_bin, "index", "-f", str(source_part)])
            source_parts.append(source_part)
            temporary.extend([source_view, source_part, Path(str(source_part) + ".csi")])
        if not source_parts:
            raise RuntimeError("요청 chromosome을 포함하는 public VCF가 없습니다")

        if len(source_parts) == 1:
            normalized_input = source_parts[0]
        else:
            combined = output_vcf.with_name(f".{output_vcf.name}.combined.bcf")
            sorted_bcf = output_vcf.with_name(f".{output_vcf.name}.sorted.bcf")
            temporary.extend([combined, sorted_bcf])
            run_command([config.bcftools_bin, "concat", "-a", "-Ob", "-o", str(combined), *map(str, source_parts)])
            run_command([config.bcftools_bin, "sort", "-Ob", "-o", str(sorted_bcf), str(combined)])
            normalized_input = sorted_bcf

        pipeline = [
            [config.bcftools_bin, "norm", "-c", "e", "-f", str(reference), "-m", "-any", "-Ou", str(normalized_input)],
            [config.bcftools_bin, "view", "-v", seq_types, "-c", "1", "-Oz", "-o", str(selected), "-"],
        ]
        run_pipeline(pipeline, "public VCF sequence extraction failed")
        if config.source_sample == config.sample:
            os.replace(selected, final_part)
        else:
            rename_file = output_vcf.with_name(f".{config.sample}.reheader.tsv")
            rename_file.write_text(f"{config.source_sample}\t{config.sample}\n", encoding="utf-8")
            try:
                run_command([config.bcftools_bin, "reheader", "-s", str(rename_file), "-o", str(final_part), str(selected)])
            finally:
                rename_file.unlink(missing_ok=True)
                selected.unlink(missing_ok=True)
        # norm -m -any는 단일 입력 VCF에서도 multi-allelic record를 분리하면서 좌표 순서를
        # 보장하지 않을 수 있다. index 전 최종 산출물을 항상 정렬해, sample alias reheader 여부와
        # 무관하게 coordinate-sorted VCF 계약을 유지한다.
        run_command([config.bcftools_bin, "sort", "-Oz", "-o", str(sorted_part), str(final_part)])
        os.replace(sorted_part, output_vcf)
        run_command([config.bcftools_bin, "index", "-t", str(output_vcf)])
    except Exception:
        selected.unlink(missing_ok=True)
        final_part.unlink(missing_ok=True)
        sorted_part.unlink(missing_ok=True)
        cleanup_bundle(output_vcf)
        raise
    finally:
        for path in temporary:
            path.unlink(missing_ok=True)


def build_tag_vcf(bcftools_bin: str, seq_vcf: Path, output_vcf: Path) -> None:
    part = output_vcf.with_name(output_vcf.name + ".part")
    part.unlink(missing_ok=True)
    cmd = [
        bcftools_bin,
        "view",
        "-m2",
        "-M2",
        "-v",
        "snps,indels",
        "-i",
        'GT="het"',
        "-Oz",
        "-o",
        str(part),
        str(seq_vcf),
    ]
    try:
        run_command(cmd)
        os.replace(part, output_vcf)
        run_command([bcftools_bin, "index", "-t", str(output_vcf)])
    except Exception:
        part.unlink(missing_ok=True)
        cleanup_bundle(output_vcf)
        raise


def inspect_vcf(path: Path, sample: str) -> VcfStats:
    stats = VcfStats()
    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
        if samples != [sample]:
            raise ValueError(f"파생 VCF sample header 오류: expected={[sample]} actual={samples}")
        for record in vcf:
            stats.records += 1
            if len(record.alts or ()) != 1:
                stats.multiallelic_records += 1
            if any(alt.startswith("<") or alt == "*" or "[" in alt or "]" in alt for alt in (record.alts or ())):
                stats.symbolic_records += 1
            call = record.samples[sample]
            gt = call.get("GT")
            if gt is None or any(allele is None for allele in gt):
                continue
            if any(int(allele) > 0 for allele in gt):
                stats.nonref_records += 1
            else:
                stats.hom_ref_records += 1
            if len(set(gt)) > 1:
                stats.het_records += 1
                if call.phased:
                    stats.phased_het_records += 1
            if "PS" in call and call.get("PS") is not None:
                stats.ps_present_records += 1
    return stats


def validate_outputs(seq_stats: VcfStats, tag_stats: VcfStats, phase_source: str) -> None:
    if seq_stats.records <= 0 or tag_stats.records <= 0:
        raise RuntimeError(f"public VCF 파생 결과가 비어 있음: seq={seq_stats.records} tag={tag_stats.records}")
    if seq_stats.hom_ref_records or seq_stats.multiallelic_records or seq_stats.symbolic_records:
        raise RuntimeError(
            "sequence VCF 계약 위반: "
            f"hom_ref={seq_stats.hom_ref_records} multiallelic={seq_stats.multiallelic_records} symbolic={seq_stats.symbolic_records}"
        )
    if tag_stats.hom_ref_records or tag_stats.multiallelic_records or tag_stats.symbolic_records:
        raise RuntimeError(
            "tag VCF 계약 위반: "
            f"hom_ref={tag_stats.hom_ref_records} multiallelic={tag_stats.multiallelic_records} symbolic={tag_stats.symbolic_records}"
        )
    if tag_stats.records != tag_stats.het_records:
        raise RuntimeError(f"tag VCF에 non-het record가 있음: records={tag_stats.records} het={tag_stats.het_records}")
    if phase_source == "provided" and tag_stats.phased_het_records != tag_stats.het_records:
        raise RuntimeError(
            "public VCF의 heterozygous GT가 모두 phased가 아님: "
            f"phased={tag_stats.phased_het_records} het={tag_stats.het_records}"
        )


def outputs_match(sidecar: Path, identity: dict[str, Any], seq_vcf: Path, tag_vcf: Path) -> bool:
    if not sidecar.exists() or not seq_vcf.exists() or not tag_vcf.exists():
        return False
    if not has_vcf_index(seq_vcf) or not has_vcf_index(tag_vcf):
        return False
    try:
        previous = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return all(previous.get(key) == value for key, value in identity.items())


def build_payload(
    identity: dict[str, Any],
    seq_vcf: Path,
    tag_vcf: Path,
    seq_stats: VcfStats,
    tag_stats: VcfStats,
    action: str,
    elapsed: float,
) -> dict[str, Any]:
    return {
        "schema": "mfm_public_vcf_prepare_v1",
        "source": "public",
        **identity,
        "action": action,
        "seq_vcf": str(seq_vcf),
        "tag_vcf": str(tag_vcf),
        "seq_stats": {**asdict(seq_stats), "phased_het_fraction": round(seq_stats.phased_het_fraction, 6)},
        "tag_stats": {**asdict(tag_stats), "phased_het_fraction": round(tag_stats.phased_het_fraction, 6)},
        "missing_ps_policy": "WhatsHap uses phase block id 0 for phased heterozygous GT without PS; allowed for provided globally oriented haplotypes",
        "fallback": "disabled",
        "elapsed_sec": round(elapsed, 3),
    }


def cleanup_bundle(vcf: Path) -> None:
    for path in (vcf, Path(str(vcf) + ".tbi"), Path(str(vcf) + ".csi"), vcf.with_name(vcf.name + ".part"), vcf.with_name(vcf.name + ".selected.part")):
        path.unlink(missing_ok=True)


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(part, path)
