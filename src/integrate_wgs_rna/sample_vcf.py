from __future__ import annotations

import csv
import fcntl
import json
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from .cohort_vcf import (
    DEFAULT_BCFTOOLS,
    chrom_sort_key,
    expand_chroms,
    has_vcf_index,
    load_sample_ids,
    run_capture,
    run_command,
)


@dataclass
class PerSampleVcfConfig:
    """Configuration for process_code-style per-sample VCF generation.

    This is the implementation used by numbered Stage 03 (``03_prepare_vcf.sh``).

    The command creates two VCFs per sample:
    - ``*.seq.vcf.gz``: variants used to build personal HAP FASTA.
    - ``*.tag.vcf.gz``: heterozygous biallelic variants used by WhatsHap.
    """

    cohort: str
    common_vcf_dir: Path
    output_dir: Path
    chroms: list[str]
    sequence_variant_policy: str
    sample_list: Path | None = None
    samples: list[str] | None = None
    common_vcf_pattern: str = ""
    concat_dir: Path | None = None
    bcftools_bin: str = DEFAULT_BCFTOOLS
    jobs: int = 8
    concat_threads: int = 0
    overwrite: bool = False
    fast_validate: bool = False
    summary_tsv: Path | None = None
    metadata_json: Path | None = None
    extraction_mode: str = "chromosome"


@dataclass
class SampleVcfSummary:
    cohort: str
    sample: str
    sequence_variant_policy: str
    seq_vcf: str
    tag_vcf: str
    status: str
    seq_variant_count: int
    seq_indel_count: int
    tag_variant_count: int
    tag_indel_count: int
    seq_hom_ref_count: int
    tag_hom_count: int
    tag_complete_phased_het_count: int
    tag_unphased_het_count: int
    tag_missing_gt_count: int
    tag_other_gt_count: int
    seq_symbolic_alt_count: int
    seq_bad_base_char_count: int
    seq_size_bytes: int
    tag_size_bytes: int
    elapsed_sec: float
    message: str = ""


def run_per_sample_vcf_extract(config: PerSampleVcfConfig) -> dict[str, object]:
    """Create process_code 02-compatible per-sample VCFs."""

    common_vcf_dir = config.common_vcf_dir.resolve()
    output_dir = config.output_dir.resolve()
    concat_dir = (config.concat_dir or common_vcf_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    concat_dir.mkdir(parents=True, exist_ok=True)
    chroms = expand_chroms(config.chroms)
    samples = load_sample_ids(config.sample_list, config.samples)
    if not samples:
        raise ValueError("sample list is empty")
    policy = normalize_sequence_variant_policy(config.sequence_variant_policy)
    seq_vtypes = policy_vtypes(policy)

    started_at = time.time()
    if config.extraction_mode not in {"chromosome", "concat"}:
        raise ValueError("VCF extraction_mode must be chromosome|concat")
    concat_vcf = (prepare_concat_vcf(config=config, common_vcf_dir=common_vcf_dir,
                                   concat_dir=concat_dir, chroms=chroms)
                  if config.extraction_mode == "concat" else None)

    summaries: list[SampleVcfSummary] = []
    from .vcf_execution import sample_parallelism, SAMPLE_CONTEXT
    import threading
    parallel_samples=sample_parallelism(config.jobs,len(samples))
    sample_slots=iter(range(parallel_samples));slot_lock=threading.Lock()
    def assign_sample_slot():
        with slot_lock:SAMPLE_CONTEXT.slot=next(sample_slots)
    with ThreadPoolExecutor(max_workers=parallel_samples,initializer=assign_sample_slot) as executor:
        futures = {
            executor.submit(
                extract_one_sample,
                config,
                sample,
                concat_vcf,
                output_dir,
                policy,
                seq_vtypes,
            ): sample
            for sample in samples
        }
        for future in as_completed(futures):
            summaries.append(future.result())

    summaries.sort(key=lambda row: row.sample)
    summary_tsv = (config.summary_tsv or (output_dir / f"{config.cohort}.persample_vcf.summary.tsv")).resolve()
    metadata_json = (config.metadata_json or (output_dir / f"{config.cohort}.persample_vcf.metadata.json")).resolve()
    write_summary_tsv(summary_tsv, summaries)
    failed = [row for row in summaries if row.status == "FAIL"]
    if policy == "snp_only":
        bad = [row for row in summaries if row.status != "FAIL" and row.seq_indel_count != 0]
        if bad:
            failed_labels = ", ".join(f"{row.sample}:seq_indel={row.seq_indel_count}" for row in bad[:5])
            raise RuntimeError(f"snp_only policy violation: {failed_labels}")

    completed_at = time.time()
    metadata = {
        "schema": "mfm_persample_vcf_v1",
        "cohort": config.cohort,
        "sample_count": len(samples),
        "chroms": chroms,
        "common_vcf_dir": str(common_vcf_dir),
        "common_vcf_pattern": config.common_vcf_pattern,
        "concat_vcf": str(concat_vcf) if concat_vcf else "",
        "extraction_mode": config.extraction_mode,
        "output_dir": str(output_dir),
        "sequence_variant_policy": policy,
        "sequence_vcf_types": seq_vtypes,
        "tag_vcf_types": "snps,indels",
        "tag_filter": 'GT="het"',
        "tag_genotype_contract": "all_records_exactly_0|1_or_1|0",
        "tag_genotype_contract_applies_in_fast_validate": True,
        "tag_biallelic_policy": "bcftools_view_m2_M2",
        "sequence_multiallelic_policy": "bcftools_norm_m_any_before_type_filter",
        "bcftools_bin": config.bcftools_bin,
        "jobs": config.jobs,
        "concat_threads": config.concat_threads,
        "overwrite": config.overwrite,
        "fast_validate": config.fast_validate,
        "summary_tsv": str(summary_tsv),
        "started_at_unix": started_at,
        "completed_at_unix": completed_at,
        "elapsed_sec": round(completed_at - started_at, 3),
        "ok_count": sum(1 for row in summaries if row.status == "OK"),
        "skip_count": sum(1 for row in summaries if row.status == "SKIP"),
        "fail_count": len(failed),
    }
    write_json(metadata_json, metadata)
    if failed:
        failed_labels = ", ".join(f"{row.sample}:{row.message}" for row in failed[:5])
        raise RuntimeError(f"per-sample VCF extraction failed: {failed_labels}")
    return metadata


def prepare_concat_vcf(*, config, common_vcf_dir, concat_dir, chroms):
    concat_dir.mkdir(parents=True, exist_ok=True)
    with (concat_dir / f".{config.cohort}.concat.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return _prepare_concat_vcf_locked(config=config, common_vcf_dir=common_vcf_dir,
                                         concat_dir=concat_dir, chroms=chroms)


def _prepare_concat_vcf_locked(
    *,
    config: PerSampleVcfConfig,
    common_vcf_dir: Path,
    concat_dir: Path,
    chroms: list[str],
) -> Path:
    """Concatenate chromosome-level common VCFs once per cohort.

    process_code writes this concat VCF into common_vcf/.  ``concat_dir`` keeps
    that behavior by default but allows an explicit scratch/output directory if
    the common VCF directory is read-only.
    """

    concat_vcf = concat_dir / f"{config.cohort}.autosomes.vcf.gz"
    stamp = concat_dir / f"{config.cohort}.autosomes.vcf.gz.chromset"
    wanted = " ".join(sorted(chroms, key=chrom_sort_key)) + " "
    from .pipeline.checkpoints import file_identity, atomic_json
    inputs = [resolve_common_vcf_path(common_vcf_dir, config.cohort, c, config.common_vcf_pattern)
              for c in sorted(chroms, key=chrom_sort_key)]
    def identities():
        files=[]
        for p in inputs:
            files.append(file_identity(p))
            index=next((Path(str(p)+x) for x in (".tbi", ".csi") if Path(str(p)+x).exists()),None)
            if index is None: raise FileNotFoundError(f"VCF index 없음: {p}")
            files.append(file_identity(index))
        return files
    before_inputs=identities()
    proof=Path(str(concat_vcf)+".inputs.json")
    try:previous=json.loads(proof.read_text())
    except (OSError,ValueError):previous={}
    if (
        concat_vcf.exists()
        and has_vcf_index(concat_vcf)
        and stamp.exists()
        and stamp.read_text(encoding="utf-8").strip() == wanted.strip()
        and previous.get("inputs") == before_inputs
        and not config.overwrite
    ):
        return concat_vcf

    inputs = [
        resolve_common_vcf_path(common_vcf_dir, config.cohort, contig, config.common_vcf_pattern)
        for contig in sorted(chroms, key=chrom_sort_key)
    ]
    for path in inputs:
        if not path.exists() or not path.stat().st_size:
            raise FileNotFoundError(f"common VCF missing: {path}")
        if not has_vcf_index(path):
            raise FileNotFoundError(f"common VCF index missing: {path}.tbi or .csi")

    part = concat_vcf.with_suffix(concat_vcf.suffix + ".part")
    part.unlink(missing_ok=True)
    from .vcf_execution import compression_threads
    cmd = [
        config.bcftools_bin,
        "concat",
        "-Oz",
        "--threads",
        str(compression_threads(config)),
        "-o",
        str(part),
        *map(str, inputs),
    ]
    run_command(cmd)
    part.replace(concat_vcf)
    run_command([config.bcftools_bin, "index", "-t", str(concat_vcf)])
    stamp.write_text(wanted + "\n", encoding="utf-8")
    if before_inputs != identities():raise RuntimeError("concat 중 패널 입력 변경")
    atomic_json(proof,{"inputs":before_inputs})
    return concat_vcf


def extract_one_sample(
    config: PerSampleVcfConfig,
    sample: str,
    concat_vcf: Path,
    output_dir: Path,
    policy: str,
    seq_vtypes: str,
) -> SampleVcfSummary:
    t0 = time.time()
    seq_vcf = output_dir / f"{sample}.seq.vcf.gz"
    tag_vcf = output_dir / f"{sample}.tag.vcf.gz"
    try:
        from .vcf_execution import source_signature, validated_metrics, extract_chromosomes
        source = source_signature(config, concat_vcf)
        seq_exists = seq_vcf.exists() and has_vcf_index(seq_vcf)
        tag_exists = tag_vcf.exists() and has_vcf_index(tag_vcf)
        if seq_exists and tag_exists and not config.overwrite:
            seq_metrics, tag_metrics = validated_metrics(config, sample, seq_vcf, tag_vcf, source)
            return build_summary(
                config=config,
                sample=sample,
                policy=policy,
                seq_vcf=seq_vcf,
                tag_vcf=tag_vcf,
                status="SKIP",
                seq_metrics=seq_metrics,
                tag_metrics=tag_metrics,
                elapsed_sec=time.time() - t0,
                message="existing valid VCFs",
            )

        if (seq_vcf.exists() or tag_vcf.exists()) and not config.overwrite:
            raise RuntimeError("partial existing output; rerun with --overwrite or remove stale files")

        unlink_vcf_bundle(seq_vcf)
        unlink_vcf_bundle(tag_vcf)
        if concat_vcf is None:
            extract_chromosomes(config, sample, seq_vcf, tag_vcf, seq_vtypes)
        else:
            extract_sequence_vcf(config, sample, concat_vcf, seq_vcf, seq_vtypes)
            extract_tag_vcf(config, sample, concat_vcf, tag_vcf)
        if source != source_signature(config, concat_vcf):
            raise RuntimeError("추출 중 VCF 패널 입력 변경")
        seq_metrics, tag_metrics = validated_metrics(config, sample, seq_vcf, tag_vcf, source, fresh=True)
        return build_summary(
            config=config,
            sample=sample,
            policy=policy,
            seq_vcf=seq_vcf,
            tag_vcf=tag_vcf,
            status="OK",
            seq_metrics=seq_metrics,
            tag_metrics=tag_metrics,
            elapsed_sec=time.time() - t0,
            message="created",
        )
    except Exception as exc:  # noqa: BLE001 - per-sample errors are recorded in the summary TSV.
        return SampleVcfSummary(
            cohort=config.cohort,
            sample=sample,
            sequence_variant_policy=policy,
            seq_vcf=str(seq_vcf),
            tag_vcf=str(tag_vcf),
            status="FAIL",
            seq_variant_count=0,
            seq_indel_count=0,
            tag_variant_count=0,
            tag_indel_count=0,
            seq_hom_ref_count=0,
            tag_hom_count=0,
            tag_complete_phased_het_count=0,
            tag_unphased_het_count=0,
            tag_missing_gt_count=0,
            tag_other_gt_count=0,
            seq_symbolic_alt_count=0,
            seq_bad_base_char_count=0,
            seq_size_bytes=seq_vcf.stat().st_size if seq_vcf.exists() else 0,
            tag_size_bytes=tag_vcf.stat().st_size if tag_vcf.exists() else 0,
            elapsed_sec=round(time.time() - t0, 3),
            message=str(exc),
        )


def extract_sequence_vcf(
    config: PerSampleVcfConfig,
    sample: str,
    concat_vcf: Path,
    output_vcf: Path,
    seq_vtypes: str,
) -> None:
    part = output_vcf.with_suffix(output_vcf.suffix + ".part")
    part.unlink(missing_ok=True)
    pipeline = [
        [config.bcftools_bin, "view", "-s", sample, "-Ou", str(concat_vcf)],
        [config.bcftools_bin, "norm", "-m", "-any", "-Ou", "-"],
        [config.bcftools_bin, "view", "-v", seq_vtypes, "-c", "1", "-Oz", "-o", str(part), "-"],
    ]
    run_pipeline(pipeline, f"sequence VCF extraction failed for {sample}")
    part.replace(output_vcf)
    run_command([config.bcftools_bin, "index", "-t", str(output_vcf)])


def extract_tag_vcf(config: PerSampleVcfConfig, sample: str, concat_vcf: Path, output_vcf: Path) -> None:
    part = output_vcf.with_suffix(output_vcf.suffix + ".part")
    part.unlink(missing_ok=True)
    pipeline = [
        [config.bcftools_bin, "view", "-s", sample, "-Ou", str(concat_vcf)],
        [
            config.bcftools_bin,
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
            "-",
        ],
    ]
    run_pipeline(pipeline, f"tag VCF extraction failed for {sample}")
    part.replace(output_vcf)
    run_command([config.bcftools_bin, "index", "-t", str(output_vcf)])


def validate_seq_vcf(config: PerSampleVcfConfig, vcf: Path, policy: str) -> dict[str, int]:
    if not vcf.exists() or not vcf.stat().st_size:
        raise FileNotFoundError(f"sequence VCF missing: {vcf}")
    if not has_vcf_index(vcf):
        raise FileNotFoundError(f"sequence VCF index missing: {vcf}.tbi or .csi")
    variant_count = vcf_variant_count(config.bcftools_bin, vcf)
    if config.fast_validate:
        if variant_count <= 0:
            raise RuntimeError(f"sequence VCF has zero variants: {vcf}")
        return {
            "variant_count": variant_count,
            "indel_count": -1,
            "hom_ref_count": -1,
            "symbolic_alt_count": -1,
            "bad_base_char_count": -1,
        }
    metrics = summarize_sequence_vcf(config.bcftools_bin, vcf)
    indel_count = metrics["indel_count"]
    hom_ref_count = metrics["hom_ref_count"]
    symbolic_alt_count = metrics["symbolic_alt_count"]
    bad_base_char_count = metrics["bad_base_char_count"]
    if variant_count <= 0:
        raise RuntimeError(f"sequence VCF has zero variants: {vcf}")
    if hom_ref_count != 0:
        raise RuntimeError(f"hom-ref records remain in sequence VCF: {hom_ref_count}")
    if symbolic_alt_count != 0:
        raise RuntimeError(f"symbolic ALT records remain in sequence VCF: {symbolic_alt_count}")
    if bad_base_char_count != 0:
        raise RuntimeError(f"REF/ALT contains non-ACGTN characters: {bad_base_char_count}")
    if policy == "snp_only" and indel_count != 0:
        raise RuntimeError(f"snp_only sequence VCF contains indels: {indel_count}")
    return {
        "variant_count": variant_count,
        "indel_count": indel_count,
        "hom_ref_count": hom_ref_count,
        "symbolic_alt_count": symbolic_alt_count,
        "bad_base_char_count": bad_base_char_count,
    }


def validate_tag_vcf(config: PerSampleVcfConfig, vcf: Path) -> dict[str, int]:
    if not vcf.exists() or not vcf.stat().st_size:
        raise FileNotFoundError(f"tag VCF missing: {vcf}")
    if not has_vcf_index(vcf):
        raise FileNotFoundError(f"tag VCF index missing: {vcf}.tbi or .csi")
    variant_count = vcf_variant_count(config.bcftools_bin, vcf)
    if variant_count <= 0:
        raise RuntimeError(f"tag VCF has zero variants: {vcf}")
    # FAST 검증에서도 phase 계약은 생략하지 않는다. 이 VCF의 GT[0]/GT[1]은
    # Stage 06 haplotag와 Stage 09 HAP1/HAP2 서열의 방향을 결정하므로,
    # 일부 record만 검사하거나 단순히 het 여부만 보면 안 된다.
    gt = summarize_tag_genotypes(config.bcftools_bin, vcf, include_indels=not config.fast_validate)
    if gt["records"] != variant_count:
        raise RuntimeError(
            f"tag VCF GT/record count mismatch: GT={gt['records']} variants={variant_count}"
        )
    invalid = gt["unphased_het_count"] + gt["missing_gt_count"] + gt["other_gt_count"]
    if invalid:
        raise RuntimeError(
            "tag VCF requires complete phased diploid heterozygous GT (0|1 or 1|0): "
            f"unphased_het={gt['unphased_het_count']} missing={gt['missing_gt_count']} "
            f"other={gt['other_gt_count']}"
        )
    indel_count = int(gt.pop("indel_count", -1))
    hom_count = gt["hom_count"]
    if hom_count != 0:
        raise RuntimeError(f"hom records remain in tag VCF: {hom_count}")
    return {"variant_count": variant_count, "indel_count": indel_count, **gt}


def build_summary(
    *,
    config: PerSampleVcfConfig,
    sample: str,
    policy: str,
    seq_vcf: Path,
    tag_vcf: Path,
    status: str,
    seq_metrics: dict[str, int],
    tag_metrics: dict[str, int],
    elapsed_sec: float,
    message: str,
) -> SampleVcfSummary:
    return SampleVcfSummary(
        cohort=config.cohort,
        sample=sample,
        sequence_variant_policy=policy,
        seq_vcf=str(seq_vcf),
        tag_vcf=str(tag_vcf),
        status=status,
        seq_variant_count=seq_metrics["variant_count"],
        seq_indel_count=seq_metrics["indel_count"],
        tag_variant_count=tag_metrics["variant_count"],
        tag_indel_count=tag_metrics["indel_count"],
        seq_hom_ref_count=seq_metrics["hom_ref_count"],
        tag_hom_count=tag_metrics["hom_count"],
        tag_complete_phased_het_count=tag_metrics["complete_phased_het_count"],
        tag_unphased_het_count=tag_metrics["unphased_het_count"],
        tag_missing_gt_count=tag_metrics["missing_gt_count"],
        tag_other_gt_count=tag_metrics["other_gt_count"],
        seq_symbolic_alt_count=seq_metrics["symbolic_alt_count"],
        seq_bad_base_char_count=seq_metrics["bad_base_char_count"],
        seq_size_bytes=seq_vcf.stat().st_size,
        tag_size_bytes=tag_vcf.stat().st_size,
        elapsed_sec=round(elapsed_sec, 3),
        message=message,
    )


def resolve_common_vcf_path(common_vcf_dir: Path, cohort: str, contig: str, pattern: str = "") -> Path:
    chrom = contig.removeprefix("chr")
    candidates: list[Path] = []
    if pattern:
        if "{" in pattern:
            candidates.append(common_vcf_dir / pattern.format(cohort=cohort, contig=contig, chrom=chrom, chromosome=contig))
        elif "%s" in pattern:
            candidates.append(common_vcf_dir / (pattern % contig))
            candidates.append(common_vcf_dir / (pattern % chrom))
        else:
            candidates.append(common_vcf_dir / pattern)
    candidates.extend(
        [
            common_vcf_dir / f"{cohort}.{contig}.vcf.gz",
            common_vcf_dir / f"{cohort}.chr{chrom}.vcf.gz",
            common_vcf_dir / f"{contig}.vcf.gz",
            common_vcf_dir / f"chr{chrom}.vcf.gz",
            common_vcf_dir / f"1kGP_high_coverage_Illumina.{contig}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
        ]
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"common VCF not found for {contig}: dir={common_vcf_dir}, pattern={pattern}")


def normalize_sequence_variant_policy(value: str) -> str:
    normalized = (value or "snv_indel").strip().lower().replace("-", "_")
    aliases = {
        "snp": "snp_only",
        "snps": "snp_only",
        "snp_only": "snp_only",
        "snv_indel": "snv_indel",
        "snps_indels": "snv_indel",
        "snps,indels": "snv_indel",
    }
    if normalized not in aliases:
        raise ValueError(f"unknown sequence variant policy: {value}")
    return aliases[normalized]


def policy_vtypes(policy: str) -> str:
    if policy == "snp_only":
        return "snps"
    if policy == "snv_indel":
        return "snps,indels"
    raise ValueError(f"unknown sequence variant policy: {policy}")


def unlink_vcf_bundle(vcf: Path) -> None:
    for path in (vcf, Path(str(vcf) + ".tbi"), Path(str(vcf) + ".csi"), vcf.with_suffix(vcf.suffix + ".part")):
        path.unlink(missing_ok=True)


def vcf_variant_count(bcftools_bin: str, vcf: Path) -> int:
    return int(run_capture([bcftools_bin, "index", "-n", str(vcf)]).strip())


def summarize_sequence_vcf(bcftools_bin: str, vcf: Path) -> dict[str, int]:
    """기존 sequence-VCF QC 항목을 bcftools 한 번의 순회로 누적한다."""
    command = [bcftools_bin, "query", "-f", "%TYPE\t%REF\t%ALT[\t%GT]\n", str(vcf)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    counts = {
        "indel_count": 0,
        "hom_ref_count": 0,
        "symbolic_alt_count": 0,
        "bad_base_char_count": 0,
    }
    allowed = set("ACGTN\t\n")
    for raw in process.stdout:
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 3:
            process.kill()
            process.communicate()
            raise RuntimeError(f"bcftools sequence VCF query returned malformed row: {raw!r}")
        variant_type, ref, alt = fields[:3]
        gt = fields[3].strip() if len(fields) >= 4 else ""
        if "INDEL" in variant_type.split(","):
            counts["indel_count"] += 1
        if gt in {"0|0", "0/0"}:
            counts["hom_ref_count"] += 1
        if re.search(r"<.*>", alt):
            counts["symbolic_alt_count"] += 1
        counts["bad_base_char_count"] += sum(
            character not in allowed for character in f"{ref}\t{alt}\n"
        )
    _, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(f"bcftools sequence VCF query failed: {stderr.strip()}")
    return counts



def summarize_tag_genotypes(
    bcftools_bin: str, vcf: Path, *, include_indels: bool = True
) -> dict[str, int]:
    """한 sample tag VCF의 GT/phase와 필요 시 indel 수를 한 번 순회한다."""
    fmt = "%TYPE\t[%GT]\n" if include_indels else "[%GT]\n"
    command = [bcftools_bin, "query", "-f", fmt, str(vcf)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    counts = {
        "records": 0,
        "complete_phased_het_count": 0,
        "unphased_het_count": 0,
        "missing_gt_count": 0,
        "other_gt_count": 0,
        "hom_count": 0,
    }
    if include_indels:
        counts["indel_count"] = 0
    for raw in process.stdout:
        if include_indels:
            fields = raw.rstrip("\n").split("\t", 1)
            if len(fields) != 2:
                process.kill()
                process.communicate()
                raise RuntimeError(f"bcftools tag VCF query returned malformed row: {raw!r}")
            variant_type, value = fields
            if "INDEL" in variant_type.split(","):
                counts["indel_count"] += 1
        else:
            value = raw
        value = value.strip()
        counts["records"] += 1
        if value in {"0|1", "1|0"}:
            counts["complete_phased_het_count"] += 1
        elif "." in value or not value:
            counts["missing_gt_count"] += 1
        elif value in {"0/1", "1/0"}:
            counts["unphased_het_count"] += 1
        else:
            counts["other_gt_count"] += 1
            if value in {"0|0", "0/0", "1|1", "1/1"}:
                counts["hom_count"] += 1
    _, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(f"bcftools tag GT query failed: {stderr.strip()}")
    return counts




def run_pipeline(pipeline: list[list[str]], label: str) -> None:
    if len(pipeline) < 2:
        raise ValueError("pipeline requires at least two commands")
    processes: list[subprocess.Popen] = []
    previous_stdout = None
    stdout_data = b""
    errors: list[str] = []
    try:
        for cmd in pipeline:
            process = subprocess.Popen(cmd, stdin=previous_stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if previous_stdout is not None:
                previous_stdout.close()
            previous_stdout = process.stdout
            processes.append(process)
        # 상류 stderr도 즉시 배출해야 큰 REF 오류가 파이프를 가득 채워 멈추지 않음.
        with ThreadPoolExecutor(max_workers=len(processes)) as readers:
            try:
                stderr_tasks = [readers.submit(process.stderr.read) for process in processes]
                stdout_data = processes[-1].stdout.read()
                for process in processes:
                    process.wait()
                errors = [task.result().decode(errors="replace") for task in stderr_tasks]
            except BaseException:
                for process in processes:
                    if process.poll() is None:
                        process.kill()
                raise
        failed = [process.returncode for process in processes if process.returncode != 0]
        if failed:
            commands = " | ".join(" ".join(cmd) for cmd in pipeline)
            stdout_text = stdout_data.decode(errors="replace")
            raise RuntimeError(f"{label}: returncodes={failed}\ncmd={commands}\nstdout={stdout_text}\nstderr={''.join(errors)}")
    finally:
        for process in processes:
            if process.poll() is None:
                process.kill()
            process.wait()
            if process.stdout is not None:
                process.stdout.close()
            if process.stderr is not None:
                process.stderr.close()


def write_summary_tsv(path: Path, rows: list[SampleVcfSummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(SampleVcfSummary.__dataclass_fields__.keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def count_query_gt(bcftools_bin: str, vcf: Path, targets: set[str]) -> int:
    output = run_capture([bcftools_bin, "query", "-f", "[%GT]\n", str(vcf)])
    return sum(1 for line in output.splitlines() if line.strip() in targets)



def count_bad_ref_alt_chars(bcftools_bin: str, vcf: Path) -> int:
    output = run_capture([bcftools_bin, "query", "-f", "%REF\t%ALT\n", str(vcf)])
    allowed = set("ACGTN\t\n")
    return sum(1 for char in output if char not in allowed)



def count_vcf_lines(bcftools_bin: str, args: list[str]) -> int:
    process = subprocess.Popen([bcftools_bin, *args], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    count = sum(1 for _ in process.stdout)
    _, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(f"bcftools {' '.join(args)} failed: {stderr.strip()}")
    return count



def summarize_tag_vcf(bcftools_bin: str, vcf: Path) -> dict[str, int]:
    """Tag VCF의 기존 indel/동형접합 카운트를 한 번의 bcftools 순회로 계산한다."""
    command = [bcftools_bin, "query", "-f", "%TYPE\t[%GT]\n", str(vcf)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    counts = {"indel_count": 0, "hom_count": 0}
    hom_gt = {"0|0", "0/0", "1|1", "1/1"}
    for raw in process.stdout:
        fields = raw.rstrip("\n").split("\t", 1)
        if len(fields) != 2:
            process.kill()
            process.communicate()
            raise RuntimeError(f"bcftools tag VCF query returned malformed row: {raw!r}")
        variant_type, gt = fields
        if "INDEL" in variant_type.split(","):
            counts["indel_count"] += 1
        if gt.strip() in hom_gt:
            counts["hom_count"] += 1
    _, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(f"bcftools tag VCF query failed: {stderr.strip()}")
    return counts

