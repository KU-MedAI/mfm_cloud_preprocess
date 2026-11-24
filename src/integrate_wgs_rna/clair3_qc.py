from __future__ import annotations

import csv
import gzip
import json
import math
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .cohort_vcf import load_sample_ids


@dataclass
class Clair3QcConfig:
    cohort: str
    output_root: Path
    platform: str
    model_path: Path
    reference_fasta: Path
    sample_list: Path | None = None
    samples: list[str] | None = None
    samtools_bin: str = "samtools"
    inspect_reads: int = 50_000
    min_mq: int = 5
    min_coverage: int = 2
    qual: float = 2.0
    min_hp_fraction: float = 0.50
    output_prefix: str = "{sample}.clair3_longread"
    summary_tsv: Path | None = None
    metadata_json: Path | None = None


@dataclass
class Clair3QcRow:
    cohort: str
    sample: str
    status: str
    platform: str
    model_path: str
    model_platform_ok: int
    selected_vcf: str
    selected_vcf_kind: str
    vcf_total_records: int
    vcf_pass_records: int
    vcf_snv_records: int
    vcf_indel_records: int
    vcf_other_records: int
    vcf_multiallelic_records: int
    vcf_gt_records: int
    vcf_nonref_gt_records: int
    vcf_het_records: int
    vcf_phased_gt_records: int
    vcf_phased_het_records: int
    vcf_phased_het_fraction: float
    # phase block 파편화 — phased_het_fraction 이 "phase 됐는가"만 재는 것을 보완한다. 아래 주석 참조.
    vcf_phase_set_count: int
    vcf_phase_block_n50_bp: int
    vcf_phase_block_mean_bp: str
    vcf_phase_block_max_bp: int
    vcf_qual_missing_records: int
    vcf_qual_below_threshold_records: int
    vcf_qual_min: str
    vcf_qual_mean: str
    haplotagged_bam: str
    inspected_reads: int
    hp_present_reads: int
    ps_present_reads: int
    hp1_reads: int
    hp2_reads: int
    hp0_or_other_reads: int
    hp_fraction: float
    ps_fraction: float
    min_mq: int
    min_coverage: int
    qual_threshold: float
    parameter_notes: str
    elapsed_sec: float
    warnings: str


def run_clair3_qc(config: Clair3QcConfig) -> dict[str, Any]:
    samples = load_sample_ids(config.sample_list, config.samples)
    if not samples:
        raise ValueError("sample list is empty")
    output_root = config.output_root.resolve()
    rows = [inspect_one_sample(config, sample, output_root) for sample in samples]
    rows.sort(key=lambda row: row.sample)
    summary_tsv = (config.summary_tsv or output_root / f"{config.cohort}.clair3_output_qc.summary.tsv").resolve()
    metadata_json = (config.metadata_json or output_root / f"{config.cohort}.clair3_output_qc.metadata.json").resolve()
    write_summary_tsv(summary_tsv, rows)
    metadata = {
        "schema": "mfm_clair3_output_qc_v1",
        "cohort": config.cohort,
        "sample_count": len(samples),
        "output_root": str(output_root),
        "platform": normalize_platform(config.platform),
        "model_path": str(config.model_path.resolve()),
        "reference_fasta": str(config.reference_fasta.resolve()),
        "samtools_bin": config.samtools_bin,
        "inspect_reads": config.inspect_reads,
        "min_mq": config.min_mq,
        "min_coverage": config.min_coverage,
        "qual": config.qual,
        "min_hp_fraction": config.min_hp_fraction,
        "summary_tsv": str(summary_tsv),
        "ok_count": sum(1 for row in rows if row.status == "OK"),
        "warn_count": sum(1 for row in rows if row.status == "WARN"),
        "missing_count": sum(1 for row in rows if row.status == "MISSING"),
        "fail_count": sum(1 for row in rows if row.status == "FAIL"),
    }
    write_json(metadata_json, metadata)
    return metadata


def inspect_one_sample(config: Clair3QcConfig, sample: str, output_root: Path) -> Clair3QcRow:
    t0 = time.time()
    sample_dir = output_root / sample
    calling_dir = sample_dir / "longread_calling" / "clair3"
    hap_bam = sample_dir / "haplotagged" / f"{config.output_prefix.format(sample=sample)}.haplotagged.bam"
    qc_json = sample_dir / "qc" / "clair3_output_qc.json"
    warnings: list[str] = []
    parameter_notes: list[str] = []
    platform = normalize_platform(config.platform)
    model_ok = model_matches_platform(config.model_path, platform)
    if not model_ok:
        warnings.append("model_platform_mismatch")
    if config.min_mq <= 5:
        parameter_notes.append("permissive_min_mq")
    if config.min_coverage <= 2:
        parameter_notes.append("permissive_min_coverage")
    if config.qual <= 2:
        parameter_notes.append("permissive_qual")

    selected_vcf, vcf_kind = select_vcf(calling_dir)
    vcf_stats = empty_vcf_stats()
    hp_stats = empty_hp_stats()
    status = "OK"
    if selected_vcf is None:
        warnings.append("missing_clair3_vcf")
        status = "MISSING"
    else:
        try:
            vcf_stats = inspect_vcf(selected_vcf, config.qual)
            if vcf_stats["total_records"] <= 0:
                warnings.append("empty_vcf")
                status = "FAIL"
            if vcf_stats["pass_records"] <= 0:
                warnings.append("no_pass_records")
            if vcf_stats["het_records"] > 0 and vcf_stats["phased_het_records"] <= 0:
                warnings.append("no_phased_het")
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"vcf_parse_error:{exc}")
            status = "FAIL"

    if not hap_bam.exists() or not hap_bam.stat().st_size:
        warnings.append("missing_haplotagged_bam")
        if status == "OK":
            status = "MISSING"
    else:
        try:
            hp_stats = inspect_hp_ps_detailed(config.samtools_bin, hap_bam, config.reference_fasta, config.inspect_reads)
            hp_fraction = fraction(hp_stats["hp_present"], hp_stats["inspected"])
            if hp_stats["inspected"] <= 0:
                warnings.append("no_reads_inspected")
                status = "FAIL"
            elif hp_stats["hp_present"] <= 0:
                warnings.append("no_hp_tags")
                status = "FAIL"
            elif hp_fraction < config.min_hp_fraction:
                warnings.append("low_hp_fraction")
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"haplotag_inspect_error:{exc}")
            status = "FAIL"

    if status == "OK" and warnings:
        status = "WARN"
    row = Clair3QcRow(
        cohort=config.cohort,
        sample=sample,
        status=status,
        platform=platform,
        model_path=str(config.model_path.resolve()),
        model_platform_ok=int(model_ok),
        selected_vcf=str(selected_vcf) if selected_vcf else "",
        selected_vcf_kind=vcf_kind,
        vcf_total_records=vcf_stats["total_records"],
        vcf_pass_records=vcf_stats["pass_records"],
        vcf_snv_records=vcf_stats["snv_records"],
        vcf_indel_records=vcf_stats["indel_records"],
        vcf_other_records=vcf_stats["other_records"],
        vcf_multiallelic_records=vcf_stats["multiallelic_records"],
        vcf_gt_records=vcf_stats["gt_records"],
        vcf_nonref_gt_records=vcf_stats["nonref_gt_records"],
        vcf_het_records=vcf_stats["het_records"],
        vcf_phased_gt_records=vcf_stats["phased_gt_records"],
        vcf_phased_het_records=vcf_stats["phased_het_records"],
        vcf_phased_het_fraction=round(fraction(vcf_stats["phased_het_records"], vcf_stats["het_records"]), 6),
        vcf_phase_set_count=vcf_stats["phase_set_count"],
        vcf_phase_block_n50_bp=vcf_stats["phase_block_n50_bp"],
        vcf_phase_block_mean_bp=format_float(vcf_stats["phase_block_mean_bp"]),
        vcf_phase_block_max_bp=vcf_stats["phase_block_max_bp"],
        vcf_qual_missing_records=vcf_stats["qual_missing_records"],
        vcf_qual_below_threshold_records=vcf_stats["qual_below_threshold_records"],
        vcf_qual_min=format_float(vcf_stats["qual_min"]),
        vcf_qual_mean=format_float(vcf_stats["qual_sum"] / vcf_stats["qual_count"] if vcf_stats["qual_count"] else math.nan),
        haplotagged_bam=str(hap_bam),
        inspected_reads=hp_stats["inspected"],
        hp_present_reads=hp_stats["hp_present"],
        ps_present_reads=hp_stats["ps_present"],
        hp1_reads=hp_stats["hp1"],
        hp2_reads=hp_stats["hp2"],
        hp0_or_other_reads=hp_stats["hp0_or_other"],
        hp_fraction=round(fraction(hp_stats["hp_present"], hp_stats["inspected"]), 6),
        ps_fraction=round(fraction(hp_stats["ps_present"], hp_stats["inspected"]), 6),
        min_mq=config.min_mq,
        min_coverage=config.min_coverage,
        qual_threshold=float(config.qual),
        parameter_notes=",".join(parameter_notes),
        elapsed_sec=round(time.time() - t0, 3),
        warnings=",".join(warnings),
    )
    payload = asdict(row)
    payload["schema"] = "mfm_clair3_sample_output_qc_v1"
    write_json(qc_json, payload)
    return row


def select_vcf(calling_dir: Path) -> tuple[Path | None, str]:
    candidates = [
        (calling_dir / "phased_merge_output.vcf.gz", "phased_merge_output"),
        (calling_dir / "merge_output.vcf.gz", "merge_output"),
        (calling_dir / "pileup.vcf.gz", "pileup_gz"),
        (calling_dir / "pileup.vcf", "pileup_plain"),
    ]
    for path, kind in candidates:
        if path.exists() and path.stat().st_size:
            return path, kind
    return None, ""


def empty_vcf_stats() -> dict[str, Any]:
    return {
        "total_records": 0,
        "pass_records": 0,
        "snv_records": 0,
        "indel_records": 0,
        "other_records": 0,
        "multiallelic_records": 0,
        "gt_records": 0,
        "nonref_gt_records": 0,
        "het_records": 0,
        "phased_gt_records": 0,
        "phased_het_records": 0,
        "qual_missing_records": 0,
        "qual_below_threshold_records": 0,
        "qual_count": 0,
        "qual_sum": 0.0,
        "qual_min": math.nan,
        "phase_set_count": 0,
        "phase_block_total_bp": 0,
        "phase_block_mean_bp": math.nan,
        "phase_block_n50_bp": 0,
        "phase_block_max_bp": 0,
    }


def inspect_vcf(vcf: Path, qual_threshold: float) -> dict[str, Any]:
    stats = empty_vcf_stats()
    phase_spans: dict[tuple[str, str], list[int]] = {}
    opener = gzip.open if str(vcf).endswith(".gz") else open
    with opener(vcf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            stats["total_records"] += 1
            ref = parts[3]
            alts = [] if parts[4] in {".", ""} else parts[4].split(",")
            if len(alts) > 1:
                stats["multiallelic_records"] += 1
            filter_value = parts[6]
            if filter_value in {"PASS", "."}:
                stats["pass_records"] += 1
            classify_variant(stats, ref, alts)
            update_qual(stats, parts[5], qual_threshold)
            if len(parts) >= 10:
                update_gt(stats, parts[8], parts[9], parts[0], parts[1], phase_spans)
    summarize_phase_blocks(stats, phase_spans)
    return stats


def summarize_phase_blocks(stats: dict[str, Any], spans: dict[tuple[str, str], list[int]]) -> None:
    """phase block 개수와 길이 분포를 stats 에 채운다.

    ★ 왜 필요한가 (2026-09-04)
      `vcf_phased_het_fraction` 은 "phase 됐는가"만 재고 **"몇 조각으로 쪼개졌는가"는 안 잰다.**
      신규 10명 실측에서 이 값은 0.942–0.964 로 전부 양호했는데, 평균 block 길이는
      657 kb(HG00139) ~ **15 kb(NA20757)** 로 45배 차이가 났다. 즉 이 지표만으로는
      극단적으로 파편화된 샘플을 걸러낼 수 없다.

      개인 FASTA 는 `bcftools consensus -H 1` 로 만들고 `1|0`/`0|1` 방향은 **phase set 마다 독립**이라,
      block 이 잘게 쪼개질수록 HAP1 이라는 이름이 실제 haplotype 에서 멀어진다(모자이크).
      유전자 단위 측정에는 영향이 없지만(유전자의 99.2% 가 단일 block 안에 있다) **긴 문맥 학습에서는
      한 window 안에서 HAP1 라벨이 뒤집힌다.**
      근거: `post-training/journal/2609/ase_truth_spread_260903.md` §5·§9
    """
    lengths = sorted(hi - lo + 1 for lo, hi in spans.values())
    total = sum(lengths)
    stats["phase_set_count"] = len(lengths)
    stats["phase_block_total_bp"] = total
    stats["phase_block_mean_bp"] = (total / len(lengths)) if lengths else math.nan
    stats["phase_block_max_bp"] = lengths[-1] if lengths else 0
    # N50 — 긴 block 부터 더해 총합의 절반을 넘기는 지점의 block 길이.
    # 평균은 짧은 block 이 많으면 끌려 내려가므로 둘 다 본다.
    acc = 0
    stats["phase_block_n50_bp"] = 0
    for length in reversed(lengths):
        acc += length
        if acc * 2 >= total:
            stats["phase_block_n50_bp"] = length
            break


def classify_variant(stats: dict[str, Any], ref: str, alts: list[str]) -> None:
    if not alts:
        stats["other_records"] += 1
        return
    if all(len(ref) == 1 and len(alt) == 1 and is_plain_base(ref) and is_plain_base(alt) for alt in alts):
        stats["snv_records"] += 1
    elif any(len(ref) != len(alt) for alt in alts) and all(not is_symbolic_alt(alt) for alt in alts):
        stats["indel_records"] += 1
    else:
        stats["other_records"] += 1


def is_plain_base(value: str) -> bool:
    return value.upper() in {"A", "C", "G", "T", "N"}


def is_symbolic_alt(alt: str) -> bool:
    return alt.startswith("<") or alt.endswith(">") or "[" in alt or "]" in alt or alt == "*"


def update_qual(stats: dict[str, Any], raw_qual: str, threshold: float) -> None:
    if raw_qual in {".", ""}:
        stats["qual_missing_records"] += 1
        return
    try:
        qual = float(raw_qual)
    except ValueError:
        stats["qual_missing_records"] += 1
        return
    stats["qual_count"] += 1
    stats["qual_sum"] += qual
    if math.isnan(stats["qual_min"]) or qual < stats["qual_min"]:
        stats["qual_min"] = qual
    if qual < threshold:
        stats["qual_below_threshold_records"] += 1


def update_gt(
    stats: dict[str, Any],
    fmt: str,
    sample_field: str,
    contig: str = "",
    pos: str = "",
    phase_spans: dict[tuple[str, str], list[int]] | None = None,
) -> None:
    keys = fmt.split(":")
    values = sample_field.split(":")
    try:
        gt_index = keys.index("GT")
    except ValueError:
        gt_index = 0
    if gt_index >= len(values):
        return
    gt = values[gt_index]
    if gt in {"", ".", "./.", ".|."}:
        return
    sep = "|" if "|" in gt else "/" if "/" in gt else ""
    allele_tokens = gt.replace("|", "/").split("/") if sep else [gt]
    called = [token for token in allele_tokens if token != "."]
    if not called:
        return
    stats["gt_records"] += 1
    if any(token != "0" for token in called):
        stats["nonref_gt_records"] += 1
    is_het = len(set(called)) > 1
    is_phased = "|" in gt
    if is_het:
        stats["het_records"] += 1
    if is_phased:
        stats["phased_gt_records"] += 1
    if is_het and is_phased:
        stats["phased_het_records"] += 1
    if is_phased and phase_spans is not None:
        record_phase_set(phase_spans, keys, values, contig, pos)


def record_phase_set(
    spans: dict[tuple[str, str], list[int]],
    keys: list[str],
    values: list[str],
    contig: str,
    pos: str,
) -> None:
    """PS 태그별로 [최소 POS, 최대 POS] 를 갱신한다 — block 길이의 근거."""
    try:
        ps_index = keys.index("PS")
    except ValueError:
        return
    if ps_index >= len(values):
        return
    ps = values[ps_index]
    if ps in {"", "."}:
        return
    try:
        position = int(pos)
    except ValueError:
        return
    # PS 값은 block 첫 변이의 좌표라 contig 가 다르면 충돌할 수 있다 → (contig, PS) 로 키를 잡는다.
    key = (contig, ps)
    span = spans.get(key)
    if span is None:
        spans[key] = [position, position]
        return
    if position < span[0]:
        span[0] = position
    elif position > span[1]:
        span[1] = position


def inspect_hp_ps_detailed(samtools_bin: str, alignment: Path, reference: Path, max_reads: int) -> dict[str, int]:
    cmd = [samtools_bin, "view", "-T", str(reference), str(alignment)]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    inspected = hp_present = ps_present = hp1 = hp2 = hp0_or_other = 0
    for line in process.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12:
            continue
        try:
            flag = int(parts[1])
        except ValueError:
            flag = 0
        if flag & 0x904:  # unmapped, secondary, supplementary
            continue
        inspected += 1
        hp_value = None
        has_ps = False
        for tag in parts[11:]:
            if tag.startswith("HP:i:"):
                hp_value = tag.rsplit(":", 1)[-1]
            elif tag.startswith("PS:i:"):
                has_ps = True
        if hp_value is not None:
            hp_present += 1
            if hp_value == "1":
                hp1 += 1
            elif hp_value == "2":
                hp2 += 1
            else:
                hp0_or_other += 1
        if has_ps:
            ps_present += 1
        if inspected >= max_reads:
            process.terminate()
            break
    _, stderr = process.communicate()
    if process.returncode not in (0, -15, 143):
        raise RuntimeError(f"samtools view failed while inspecting HP/PS: {stderr.strip()}")
    return {
        "inspected": inspected,
        "hp_present": hp_present,
        "ps_present": ps_present,
        "hp1": hp1,
        "hp2": hp2,
        "hp0_or_other": hp0_or_other,
    }


def empty_hp_stats() -> dict[str, int]:
    return {"inspected": 0, "hp_present": 0, "ps_present": 0, "hp1": 0, "hp2": 0, "hp0_or_other": 0}


def normalize_platform(value: str) -> str:
    p = (value or "").strip().lower().replace("-", "_").replace(" ", "")
    if p in {"ont", "nanopore", "oxfordnanopore", "oxford_nanopore"}:
        return "ont"
    if p in {"hifi", "pacbio_hifi", "pacbiohifi", "ccs", "pacbio_ccs", "pacbioccs"}:
        return "hifi"
    return p


def model_matches_platform(model_path: Path, platform: str) -> bool:
    name = model_path.name.lower()
    if platform == "hifi":
        return "hifi" in name or "revio" in name or "sequel" in name
    if platform == "ont":
        return "ont" in name or "guppy" in name or name.startswith("r9") or name.startswith("r10")
    return False


def fraction(numerator: int, denominator: int) -> float:
    return 0.0 if denominator <= 0 else numerator / denominator


def format_float(value: float) -> str:
    if isinstance(value, float) and math.isnan(value):
        return "NaN"
    return f"{value:.6g}"


def write_summary_tsv(path: Path, rows: list[Clair3QcRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(Clair3QcRow.__dataclass_fields__.keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
