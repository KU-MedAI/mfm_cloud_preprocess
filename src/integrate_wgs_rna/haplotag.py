from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

import pysam

from .cohort_vcf import load_sample_ids, run_command
from .haplotag_contract import validate_haplotag_contract
from .haplotag_resources import reserve_tail


DEFAULT_WHATSHAP = "whatshap"
DEFAULT_SAMTOOLS = "samtools"
LONG_READ_PLATFORMS = {"ont", "nanopore", "oxfordnanopore", "oxford_nanopore", "hifi", "pacbio_hifi", "pacbiohifi", "ccs", "pacbio_ccs", "pacbioccs"}


def normalize_platform(value: str) -> str:
    return (value or "").strip().lower().replace("-", "_").replace(" ", "_")


def is_long_read_platform(value: str) -> bool:
    normalized = normalize_platform(value)
    return (
        normalized in LONG_READ_PLATFORMS
        or normalized.startswith("ont")
        or "nanopore" in normalized
        or "hifi" in normalized
        or "ccs" in normalized
        or normalized.startswith("pacbio")
    )


def validate_v2_haplotag_policy(config: "HaplotagConfig") -> None:
    if not normalize_platform(config.platform):
        raise ValueError("--platform is required in v2 haplotag command")
    if is_long_read_platform(config.platform) and not config.allow_longread_vcf_haplotag:
        raise ValueError(
            "long-read VCF haplotagging requires explicit allow_longread_vcf_haplotag opt-in; "
            "the caller must validate that the VCF is donor-matched, phased, and reference-compatible"
        )


@dataclass
class HaplotagConfig:
    """Configuration for process_code-style WhatsHap haplotagging.

    This corresponds to process_code/04_haplotag.sh, but follows the selected
    policy for already HP-tagged inputs: do not copy or rewrite the original
    alignment; record it as the downstream alignment in metadata.
    """

    cohort: str
    raw_root: Path
    tag_vcf_dir: Path
    output_root: Path
    reference_fasta: Path
    platform: str
    alignment: Path | None = None
    allow_longread_vcf_haplotag: bool = False
    sample_list: Path | None = None
    samples: list[str] | None = None
    alignment_pattern: str = "{sample}/wgs/{sample}.hg38.cram"
    output_prefix: str = "{sample}.ont.hg38.haplotagged"
    whatshap_bin: str = DEFAULT_WHATSHAP
    samtools_bin: str = DEFAULT_SAMTOOLS
    jobs: int = 6
    output_threads: int = 4
    inspect_reads: int = 50_000
    tag_supplementary: bool = False
    overwrite: bool = False
    summary_tsv: Path | None = None
    metadata_json: Path | None = None
    source_sample_id: str | None = None


@dataclass
class HaplotagSummary:
    cohort: str
    sample: str
    status: str
    action: str
    input_alignment: str
    tag_vcf: str
    downstream_alignment: str
    haplotagged_bam: str
    taglist: str
    hp_present_before: int
    ps_present_before: int
    inspected_before: int
    taglist_total_reads: int
    taglist_h1_reads: int
    taglist_h2_reads: int
    hp_present_after: int
    ps_present_after: int
    inspected_after: int
    elapsed_sec: float
    message: str = ""


def run_haplotag_generation(config: HaplotagConfig) -> dict[str, object]:
    validate_v2_haplotag_policy(config)
    raw_root = config.raw_root.resolve()
    tag_vcf_dir = config.tag_vcf_dir.resolve()
    output_root = config.output_root.resolve()
    reference = config.reference_fasta.resolve()
    if not reference.exists() or not reference.stat().st_size:
        raise FileNotFoundError(f"reference FASTA가 없습니다: {reference}")
    if not Path(str(reference) + ".fai").exists():
        raise FileNotFoundError(f"reference FASTA index가 없습니다: {reference}.fai")
    samples = load_sample_ids(config.sample_list, config.samples)
    if not samples:
        raise ValueError("sample list is empty")
    if config.alignment is not None and len(samples) != 1:
        raise ValueError("--alignment can be used only with exactly one --sample")
    if config.source_sample_id is not None and len(samples) != 1:
        raise ValueError("--source-sample-id can be used only with exactly one --sample")
    output_root.mkdir(parents=True, exist_ok=True)

    started_at = time.time()
    summaries: list[HaplotagSummary] = []
    with ThreadPoolExecutor(max_workers=max(1, config.jobs)) as executor:
        futures = {
            executor.submit(process_one_sample, config, sample, raw_root, tag_vcf_dir, output_root, reference): sample
            for sample in samples
        }
        for future in as_completed(futures):
            summaries.append(future.result())

    summaries.sort(key=lambda row: row.sample)
    summary_tsv = (config.summary_tsv or (output_root / f"{config.cohort}.haplotag.summary.tsv")).resolve()
    metadata_json = (config.metadata_json or (output_root / f"{config.cohort}.haplotag.metadata.json")).resolve()
    write_summary_tsv(summary_tsv, summaries)
    failed = [row for row in summaries if row.status == "FAIL"]
    completed_at = time.time()
    metadata = {
        "schema": "mfm_haplotag_v1",
        "cohort": config.cohort,
        "sample_count": len(samples),
        "raw_root": str(raw_root),
        "explicit_alignment": str(config.alignment.resolve()) if config.alignment is not None else "",
        "tag_vcf_dir": str(tag_vcf_dir),
        "output_root": str(output_root),
        "reference_fasta": str(reference),
        "source_sample_id": config.source_sample_id,
        "platform": config.platform,
        "allow_longread_vcf_haplotag": bool(config.allow_longread_vcf_haplotag),
        "alignment_pattern": config.alignment_pattern,
        "output_prefix": config.output_prefix,
        "already_hp_policy": (
            "retag_with_validated_public_vcf"
            if is_long_read_platform(config.platform) and config.allow_longread_vcf_haplotag
            else "skip_and_record_original_as_downstream_alignment"
        ),
        "longread_production_policy": (
            "validated public phased VCF opt-in"
            if is_long_read_platform(config.platform) and config.allow_longread_vcf_haplotag
            else "Clair3 or short-read panel policy"
        ),
        "tag_supplementary": config.tag_supplementary,
        "whatshap_bin": config.whatshap_bin,
        "samtools_bin": config.samtools_bin,
        "jobs": config.jobs,
        "output_threads": config.output_threads,
        "inspect_reads": config.inspect_reads,
        "overwrite": config.overwrite,
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
        raise RuntimeError(f"haplotag generation failed: {failed_labels}")
    return metadata


def process_one_sample(
    config: HaplotagConfig,
    sample: str,
    raw_root: Path,
    tag_vcf_dir: Path,
    output_root: Path,
    reference: Path,
) -> HaplotagSummary:
    t0 = time.time()
    alignment = config.alignment.resolve() if config.alignment is not None else resolve_alignment(raw_root, config.alignment_pattern, sample)
    tag_vcf = tag_vcf_dir / f"{sample}.tag.vcf.gz"
    hap_dir = output_root / sample / "haplotagged"
    prefix = config.output_prefix.format(sample=sample)
    bam = hap_dir / f"{prefix}.bam"
    taglist = hap_dir / f"{sample}.taglist.tsv.gz"
    qc_json = hap_dir / f"{sample}.haplotag_qc.json"
    try:
        require_alignment(alignment)
        require_vcf(tag_vcf)
        contract = validate_haplotag_contract(
            alignment, tag_vcf, reference, sample, config.source_sample_id,
        )
        hap_dir.mkdir(parents=True, exist_ok=True)
        before = inspect_hp_ps(config.samtools_bin, alignment, reference, config.inspect_reads)
        if before["hp_present"] > 0 and not config.allow_longread_vcf_haplotag:
            qc = {
                "schema": "mfm_haplotag_sample_qc_v1",
                "cohort": config.cohort,
                "sample": sample,
                "action": "skip_existing_hp",
                "policy": "input alignment already has HP tags; original alignment is recorded as downstream_alignment without copying",
                "input_alignment": str(alignment),
                "downstream_alignment": str(alignment),
                "tag_vcf": str(tag_vcf),
                "haplotagged_bam": "",
                "taglist": "",
                "platform": config.platform,
                "allow_longread_vcf_haplotag": bool(config.allow_longread_vcf_haplotag),
                "before": before,
                "after": before,
                "input_contract": contract,
                "status": "pass",
            }
            write_json(qc_json, qc)
            return HaplotagSummary(
                cohort=config.cohort,
                sample=sample,
                status="SKIP",
                action="skip_existing_hp",
                input_alignment=str(alignment),
                tag_vcf=str(tag_vcf),
                downstream_alignment=str(alignment),
                haplotagged_bam="",
                taglist="",
                hp_present_before=before["hp_present"],
                ps_present_before=before["ps_present"],
                inspected_before=before["inspected"],
                taglist_total_reads=0,
                taglist_h1_reads=0,
                taglist_h2_reads=0,
                hp_present_after=before["hp_present"],
                ps_present_after=before["ps_present"],
                inspected_after=before["inspected"],
                elapsed_sec=round(time.time() - t0, 3),
                message="input already has HP tags; downstream_alignment is original alignment",
            )

        if bam.exists() and has_alignment_index(bam) and taglist.exists() and qc_json.exists() and not config.overwrite:
            reserve_tail(config, sample, "existing_qc", 2)
            after = inspect_hp_ps(config.samtools_bin, bam, reference, config.inspect_reads)
            reserve_tail(config, sample, "tag_counts")
            tag_counts = summarize_taglist(taglist)
            if os.environ.get("MFM_HAPLOTAG_KEEP_SHARDS", "0") != "1":
                from .haplotag_cleanup import cleanup_parallel_shards
                cleanup_parallel_shards(qc_json, config.samtools_bin)
            if after["hp_present"] <= 0:
                raise RuntimeError(f"existing haplotagged BAM has no HP tags in inspected reads: {bam}")
            return build_created_summary(
                config=config,
                sample=sample,
                alignment=alignment,
                tag_vcf=tag_vcf,
                bam=bam,
                taglist=taglist,
                before=before,
                after=after,
                tag_counts=tag_counts,
                elapsed_sec=time.time() - t0,
                status="SKIP",
                message="existing valid haplotagged BAM",
            )

        if (bam.exists() or taglist.exists()) and not config.overwrite:
            raise RuntimeError("partial existing haplotag output; rerun with --overwrite or remove stale files")
        # 병렬 실행은 기존 부분 결과를 보존하고 완성 후보만 게시한다.
        parallel_cpus = int(os.environ.get("MFM_HAPLOTAG_CPU_BUDGET", "1"))
        if parallel_cpus <= 1:
            unlink_alignment_bundle(bam)
            taglist.unlink(missing_ok=True)
        run_whatshap(config, sample, alignment, tag_vcf, reference, bam, taglist)
        index_grant = reserve_tail(config, sample, "output_index", max(1, config.output_threads + 1), 64 << 20)
        index_threads = max(0, config.output_threads)
        if index_grant is not None:index_threads = min(index_threads, index_grant - 1)
        lease_file = os.environ.get("MFM_CPU_POOL_LEASE_FILE")
        if lease_file:
            index_threads = min(index_threads, max(0, int(json.loads(Path(lease_file).read_text())["cpus"]) - 1))
        run_command([config.samtools_bin, "index", f"-@{index_threads}", str(bam)])
        reserve_tail(config, sample, "inspect_tags", 2)
        after = inspect_hp_ps(config.samtools_bin, bam, reference, config.inspect_reads)
        if after["hp_present"] <= 0:
            raise RuntimeError("whatshap completed but HP tags are absent in inspected reads")
        reserve_tail(config, sample, "tag_counts")
        tag_counts = summarize_taglist(taglist)
        if tag_counts["total"] <= 0:
            raise RuntimeError("taglist is empty")
        qc = {
            "schema": "mfm_haplotag_sample_qc_v1",
            "cohort": config.cohort,
            "sample": sample,
            "action": "created_haplotagged_bam",
            "input_alignment": str(alignment),
            "downstream_alignment": str(bam),
            "tag_vcf": str(tag_vcf),
            "haplotagged_bam": str(bam),
            "taglist": str(taglist),
            "reference_fasta": str(reference),
            "platform": config.platform,
            "allow_longread_vcf_haplotag": bool(config.allow_longread_vcf_haplotag),
            "tag_supplementary": config.tag_supplementary,
            "input_contract": contract,
            "before": before,
            "after": after,
            "taglist_counts": tag_counts,
            "status": "pass",
        }
        write_json(qc_json, qc)
        if os.environ.get("MFM_HAPLOTAG_KEEP_SHARDS", "0") != "1":
            from .haplotag_cleanup import cleanup_parallel_shards
            cleanup_parallel_shards(qc_json, config.samtools_bin)
        return build_created_summary(
            config=config,
            sample=sample,
            alignment=alignment,
            tag_vcf=tag_vcf,
            bam=bam,
            taglist=taglist,
            before=before,
            after=after,
            tag_counts=tag_counts,
            elapsed_sec=time.time() - t0,
            status="OK",
            message="created",
        )
    except Exception as exc:  # noqa: BLE001 - per-sample failures are reported in summary TSV.
        return HaplotagSummary(
            cohort=config.cohort,
            sample=sample,
            status="FAIL",
            action="failed",
            input_alignment=str(alignment),
            tag_vcf=str(tag_vcf),
            downstream_alignment="",
            haplotagged_bam=str(bam),
            taglist=str(taglist),
            hp_present_before=0,
            ps_present_before=0,
            inspected_before=0,
            taglist_total_reads=0,
            taglist_h1_reads=0,
            taglist_h2_reads=0,
            hp_present_after=0,
            ps_present_after=0,
            inspected_after=0,
            elapsed_sec=round(time.time() - t0, 3),
            message=str(exc),
        )


def build_created_summary(
    *,
    config: HaplotagConfig,
    sample: str,
    alignment: Path,
    tag_vcf: Path,
    bam: Path,
    taglist: Path,
    before: dict[str, int],
    after: dict[str, int],
    tag_counts: dict[str, int],
    elapsed_sec: float,
    status: str,
    message: str,
) -> HaplotagSummary:
    return HaplotagSummary(
        cohort=config.cohort,
        sample=sample,
        status=status,
        action="created_haplotagged_bam",
        input_alignment=str(alignment),
        tag_vcf=str(tag_vcf),
        downstream_alignment=str(bam),
        haplotagged_bam=str(bam),
        taglist=str(taglist),
        hp_present_before=before["hp_present"],
        ps_present_before=before["ps_present"],
        inspected_before=before["inspected"],
        taglist_total_reads=tag_counts["total"],
        taglist_h1_reads=tag_counts["h1"],
        taglist_h2_reads=tag_counts["h2"],
        hp_present_after=after["hp_present"],
        ps_present_after=after["ps_present"],
        inspected_after=after["inspected"],
        elapsed_sec=round(elapsed_sec, 3),
        message=message,
    )


def run_whatshap(
    config: HaplotagConfig,
    sample: str,
    alignment: Path,
    tag_vcf: Path,
    reference: Path,
    output_bam: Path,
    taglist: Path,
) -> None:
    contract = validate_haplotag_contract(
        alignment, tag_vcf, reference, sample, config.source_sample_id,
    )
    parallel_cpus = int(os.environ.get("MFM_HAPLOTAG_CPU_BUDGET", "1"))
    if parallel_cpus > 1:
        from .haplotag_parallel import run_parallel
        run_parallel(config, sample, alignment, tag_vcf, reference,
                     output_bam, taglist, contract, parallel_cpus)
        return
    part = output_bam.with_suffix(output_bam.suffix + ".part")
    part.unlink(missing_ok=True)
    cmd = [
        config.whatshap_bin,
        "haplotag",
        "--sample",
        sample,
        "--reference",
        str(reference),
        "--output-threads",
        str(max(1, config.output_threads)),
        "--skip-missing-contigs",
        "--output-haplotag-list",
        str(taglist),
        "-o",
        str(part),
    ]
    if config.tag_supplementary:
        # Current WhatsHap treats --tag-supplementary as an optional-value argument.
        # Use the = form so the following positional VCF is never consumed as the option value.
        cmd.append("--tag-supplementary=copy-primary")
    if contract["ignore_read_groups"]:
        # Manifest의 단일 donor alias와 전체 @RG/VCF를 대조한 뒤에만 허용한다.
        # BAM header의 원본 SM을 바꾸지 않고 canonical donor VCF를 사용한다.
        cmd.append("--ignore-read-groups")
    cmd.extend([str(tag_vcf), str(alignment)])
    run_command(cmd)
    part.replace(output_bam)


def resolve_alignment(raw_root: Path, pattern: str, sample: str) -> Path:
    path = raw_root / pattern.format(sample=sample)
    if path.exists():
        return path
    candidates = [
        raw_root / sample / "wgs" / f"{sample}.final.cram",
        raw_root / sample / "wgs" / f"{sample}.hg38.cram",
        raw_root / sample / "wgs" / f"{sample}.ont.hg38.cram",
        raw_root / sample / "wgs" / f"{sample}.hg38.bam",
        raw_root / sample / "wgs" / f"{sample}.ont.hg38.bam",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return path


def require_alignment(alignment: Path) -> None:
    if not alignment.exists() or not alignment.stat().st_size:
        raise FileNotFoundError(f"alignment file is missing: {alignment}")


def require_vcf(vcf: Path) -> None:
    if not vcf.exists() or not vcf.stat().st_size:
        raise FileNotFoundError(f"tag VCF is missing: {vcf}")
    if not Path(str(vcf) + ".tbi").exists() and not Path(str(vcf) + ".csi").exists():
        raise FileNotFoundError(f"tag VCF index is missing: {vcf}.tbi or .csi")


def inspect_hp_ps(samtools_bin: str, alignment: Path, reference: Path, max_reads: int) -> dict[str, int]:
    cmd = [samtools_bin, "view", "-T", str(reference), str(alignment)]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    inspected = 0
    hp_present = 0
    ps_present = 0
    for line in process.stdout:
        inspected += 1
        if "\tHP:i:" in line:
            hp_present += 1
        if "\tPS:i:" in line:
            ps_present += 1
        if inspected >= max_reads:
            process.terminate()
            break
    _, stderr = process.communicate()
    if process.returncode not in (0, -15, 143):
        raise RuntimeError(f"samtools view failed while inspecting HP/PS: {stderr.strip()}")
    return {"inspected": inspected, "hp_present": hp_present, "ps_present": ps_present}


def summarize_taglist(taglist: Path) -> dict[str, int]:
    if not taglist.exists() or not taglist.stat().st_size:
        raise FileNotFoundError(f"taglist is missing: {taglist}")
    opener = gzip.open if taglist.suffix == ".gz" else open
    total = h1 = h2 = 0
    with opener(taglist, "rt", encoding="utf-8", errors="replace") as handle:
        for lineno, line in enumerate(handle, start=1):
            if lineno == 1 and line.lower().startswith("read"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            total += 1
            if parts[1] == "H1":
                h1 += 1
            elif parts[1] == "H2":
                h2 += 1
    return {"total": total, "h1": h1, "h2": h2}


def has_alignment_index(alignment: Path) -> bool:
    return Path(str(alignment) + ".bai").exists() or Path(str(alignment) + ".csi").exists() or alignment.with_suffix(".bai").exists()


def unlink_alignment_bundle(alignment: Path) -> None:
    for path in (
        alignment,
        Path(str(alignment) + ".bai"),
        Path(str(alignment) + ".csi"),
        alignment.with_suffix(".bai"),
        alignment.with_suffix(alignment.suffix + ".part"),
    ):
        path.unlink(missing_ok=True)


def write_summary_tsv(path: Path, rows: list[HaplotagSummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(HaplotagSummary.__dataclass_fields__.keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
