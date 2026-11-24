from __future__ import annotations

import csv
import gzip
import json
import os
import shlex
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .cohort_vcf import load_sample_ids
from .haplotag import inspect_hp_ps


DEFAULT_CLAIR3_BIN = "run_clair3.py"
DEFAULT_WHATSHAP = "whatshap"
DEFAULT_SAMTOOLS = "samtools"
DEFAULT_PYTHON = "python"
DEFAULT_PYPY = "pypy3"
DEFAULT_PARALLEL = "parallel"
DEFAULT_TABIX = "tabix"


PLATFORM_ALIASES = {
    "ont": "ont",
    "nanopore": "ont",
    "oxfordnanopore": "ont",
    "oxford_nanopore": "ont",
    "oxford-nanopore": "ont",
    "hifi": "hifi",
    "pacbio_hifi": "hifi",
    "pacbio-hifi": "hifi",
    "pacbiohifi": "hifi",
    "ccs": "hifi",
    "pacbio_ccs": "hifi",
    "pacbio-ccs": "hifi",
    "pacbioccs": "hifi",
}


@dataclass
class Clair3LongReadConfig:
    """Run long-read self-calling and haplotagging for v2 preprocessing.

    v2 정책에서 long-read BAM/CRAM은 Illumina statistical VCF로 직접
    haplotagging하지 않는다. 대신 같은 long-read alignment에서 Clair3로
    variant를 call하고, Clair3가 phase한 VCF를 WhatsHap haplotag 입력으로
    사용한다. 이렇게 생성된 HP/PS-tagged BAM이 이후 HAP별 personal
    realignment와 compact zarr 생성 단계의 입력이 된다.
    """

    cohort: str
    output_root: Path
    reference_fasta: Path
    model_path: Path
    platform: str
    raw_root: Path | None = None
    alignment: Path | None = None
    sample_list: Path | None = None
    samples: list[str] | None = None
    alignment_pattern: str = "{sample}/wgs/{sample}.hg38.cram"
    clair3_bin: str = DEFAULT_CLAIR3_BIN
    whatshap_bin: str = DEFAULT_WHATSHAP
    samtools_bin: str = DEFAULT_SAMTOOLS
    python_bin: str = DEFAULT_PYTHON
    pypy_bin: str = DEFAULT_PYPY
    parallel_bin: str = DEFAULT_PARALLEL
    tabix_bin: str = DEFAULT_TABIX
    jobs: int = 1
    threads: int = 32
    output_threads: int = 8
    inspect_reads: int = 50_000
    contigs: str = ""
    bed_fn: Path | None = None
    chunk_size: int = 5_000_000
    min_mq: int = 5
    min_coverage: int = 2
    qual: int = 2
    enable_long_indel: bool = False
    remove_intermediate_dir: bool = True
    tag_supplementary: bool = True
    ignore_read_groups: bool = False
    haplotag_mode: str = "external_whatshap"
    output_prefix: str = "{sample}.clair3_longread"
    export_sequence_vcf: bool = False
    sequence_vcf_dir: Path | None = None
    sequence_vcf_export_mode: str = "copy"
    allow_unphased_sequence_vcf: bool = False
    allow_empty_haplotag: bool = False
    allow_pileup_vcf_fallback: bool = False
    overwrite: bool = False
    dry_run: bool = False
    summary_tsv: Path | None = None
    metadata_json: Path | None = None
    extra_clair3_args: list[str] | None = None


@dataclass
class Clair3LongReadSummary:
    cohort: str
    sample: str
    status: str
    action: str
    input_alignment: str
    clair3_dir: str
    merged_vcf: str
    sequence_vcf: str
    sequence_vcf_source: str
    sequence_vcf_phased_het_fraction: float
    haplotagged_bam: str
    taglist: str
    hp_present_after: int
    ps_present_after: int
    inspected_after: int
    elapsed_sec: float
    message: str = ""


def run_clair3_longread_generation(config: Clair3LongReadConfig) -> dict[str, Any]:
    """Run Clair3 + WhatsHap for one or more long-read samples.

    The command is intentionally a wrapper, not a reimplementation of variant
    calling. All biological decisions are delegated to Clair3/WhatsHap, while
    this function standardizes paths, metadata, and validation for the v2
    integrated preprocessing pipeline.
    """

    validate_config(config)
    samples = resolve_samples(config)
    if not samples:
        raise ValueError("sample list is empty")
    output_root = config.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    started_at = time.time()
    summaries: list[Clair3LongReadSummary] = []
    with ThreadPoolExecutor(max_workers=max(1, config.jobs)) as executor:
        futures = {
            executor.submit(process_one_sample, config, sample, output_root): sample
            for sample in samples
        }
        for future in as_completed(futures):
            summaries.append(future.result())
    summaries.sort(key=lambda row: row.sample)
    summary_tsv = (config.summary_tsv or (output_root / f"{config.cohort}.clair3_longread.summary.tsv")).resolve()
    metadata_json = (config.metadata_json or (output_root / f"{config.cohort}.clair3_longread.metadata.json")).resolve()
    write_summary_tsv(summary_tsv, summaries)
    completed_at = time.time()
    metadata = {
        "schema": "mfm_longread_clair3_haplotag_v1",
        "cohort": config.cohort,
        "policy": "long-read BAM/CRAM is self-called with Clair3; Illumina VCF is not used for long-read haplotagging or default FASTA generation",
        "sequence_vcf_policy": "export Clair3 long-read-derived VCF to <sample>.seq.vcf.gz when export_sequence_vcf=true",
        "export_sequence_vcf": bool(config.export_sequence_vcf),
        "sequence_vcf_dir": str(config.sequence_vcf_dir.resolve()) if config.sequence_vcf_dir else "",
        "sequence_vcf_export_mode": config.sequence_vcf_export_mode,
        "allow_unphased_sequence_vcf": bool(config.allow_unphased_sequence_vcf),
        "sample_count": len(samples),
        "output_root": str(output_root),
        "reference_fasta": str(config.reference_fasta.resolve()),
        "model_path": str(config.model_path.resolve()),
        "platform_input": config.platform,
        "platform_clair3": normalize_clair3_platform(config.platform),
        "haplotag_mode": config.haplotag_mode,
        "clair3_bin": config.clair3_bin,
        "whatshap_bin": config.whatshap_bin,
        "samtools_bin": config.samtools_bin,
        "python_bin": config.python_bin,
        "pypy_bin": config.pypy_bin,
        "parallel_bin": config.parallel_bin,
        "tabix_bin": config.tabix_bin,
        "jobs": config.jobs,
        "threads_per_sample": config.threads,
        "output_threads": config.output_threads,
        "contigs": config.contigs,
        "bed_fn": str(config.bed_fn.resolve()) if config.bed_fn else "",
        "remove_intermediate_dir": bool(config.remove_intermediate_dir),
        "tag_supplementary": bool(config.tag_supplementary),
        "ignore_read_groups": bool(config.ignore_read_groups),
        "allow_empty_haplotag": bool(config.allow_empty_haplotag),
        "overwrite": bool(config.overwrite),
        "dry_run": bool(config.dry_run),
        "started_at_unix": started_at,
        "completed_at_unix": completed_at,
        "elapsed_sec": round(completed_at - started_at, 3),
        "ok_count": sum(1 for row in summaries if row.status == "OK"),
        "skip_count": sum(1 for row in summaries if row.status == "SKIP"),
        "planned_count": sum(1 for row in summaries if row.status == "PLAN"),
        "fail_count": sum(1 for row in summaries if row.status == "FAIL"),
        "summary_tsv": str(summary_tsv),
    }
    write_json(metadata_json, metadata)
    failed = [row for row in summaries if row.status == "FAIL"]
    if failed:
        labels = ", ".join(f"{row.sample}:{row.message}" for row in failed[:5])
        raise RuntimeError(f"Clair3 long-read generation failed: {labels}")
    return metadata


def validate_config(config: Clair3LongReadConfig) -> None:
    platform = normalize_clair3_platform(config.platform)
    if platform not in {"ont", "hifi"}:
        raise ValueError(f"long-read Clair3 platform must resolve to ont or hifi: {config.platform}")
    if config.haplotag_mode not in {"external_whatshap", "clair3_internal", "none"}:
        raise ValueError("--haplotag-mode must be one of: external_whatshap, clair3_internal, none")
    if config.sequence_vcf_export_mode not in {"copy", "symlink"}:
        raise ValueError("--sequence-vcf-export-mode must be copy or symlink")
    if config.export_sequence_vcf and config.sequence_vcf_dir is None:
        raise ValueError("--export-sequence-vcf requires --sequence-vcf-dir")
    if config.alignment is not None and (config.sample_list is not None or (config.samples and len(config.samples) != 1)):
        raise ValueError("--alignment is a single-sample override; provide exactly one --sample and no --sample-list")
    if config.alignment is None and config.raw_root is None:
        raise ValueError("provide --alignment for one sample or --raw-root with --alignment-pattern")
    if config.contigs and config.bed_fn:
        # Clair3 allows both, but using both makes smoke/QC interpretation ambiguous.
        raise ValueError("use either --contigs/--ctg_name or --bed-fn, not both")
    if not config.contigs and config.bed_fn is None:
        raise ValueError("provide --contigs or --bed-fn for Clair3 calling scope")
    if not config.dry_run:
        require_file(config.reference_fasta, "reference FASTA")
        require_file(Path(str(config.reference_fasta) + ".fai"), "reference FASTA index")
        require_path_exists(config.model_path, "Clair3 model path")
        require_executable(config.clair3_bin, "Clair3 run_clair3.sh")
        if config.haplotag_mode == "external_whatshap":
            require_executable(config.whatshap_bin, "WhatsHap")
            require_executable(config.samtools_bin, "samtools")


def resolve_samples(config: Clair3LongReadConfig) -> list[str]:
    if config.alignment is not None:
        samples = load_sample_ids(config.sample_list, config.samples)
        if len(samples) != 1:
            raise ValueError("--alignment mode requires exactly one sample via --sample")
        return samples
    return load_sample_ids(config.sample_list, config.samples)


def process_one_sample(config: Clair3LongReadConfig, sample: str, output_root: Path) -> Clair3LongReadSummary:
    t0 = time.time()
    alignment = resolve_alignment(config, sample)
    prefix = config.output_prefix.format(sample=sample)
    sample_dir = output_root / sample
    calling_dir = sample_dir / "longread_calling" / "clair3"
    hap_dir = sample_dir / "haplotagged"
    log_dir = sample_dir / "logs"
    hap_bam = hap_dir / f"{prefix}.haplotagged.bam"
    taglist = hap_dir / f"{prefix}.taglist.tsv.gz"
    qc_json = sample_dir / "qc" / "clair3_longread_haplotag.json"
    merged_vcf = calling_dir / "merge_output.vcf.gz"
    try:
        if not config.dry_run:
            require_file(alignment, "input alignment")
            require_alignment_index(alignment)
        calling_dir.mkdir(parents=True, exist_ok=True)
        hap_dir.mkdir(parents=True, exist_ok=True)
        qc_json.parent.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)
        commands = [build_clair3_command(config, alignment, calling_dir, sample)]
        if config.haplotag_mode == "external_whatshap":
            commands.append(build_whatshap_command(config, sample, alignment, merged_vcf, hap_bam, taglist))
            commands.append([config.samtools_bin, "index", f"-@{max(1, config.output_threads)}", str(hap_bam)])
        if config.dry_run:
            qc = build_qc_payload(config, sample, alignment, calling_dir, merged_vcf, hap_bam, taglist, commands, {}, "dry_run")
            write_json(qc_json, qc)
            return summary_row(config, sample, "PLAN", "dry_run", alignment, calling_dir, merged_vcf, hap_bam, taglist, {}, time.time() - t0, "planned only")
        if existing_outputs_valid(config, merged_vcf, hap_bam, taglist):
            haplotag_vcf = select_clair3_haplotag_vcf(config, calling_dir, merged_vcf, log_dir)
            sequence_export = export_sequence_vcf_if_requested(config, sample, haplotag_vcf, log_dir)
            after = inspect_optional_hp_ps(config, hap_bam)
            return summary_row(config, sample, "SKIP", "existing_valid", alignment, calling_dir, merged_vcf, hap_bam, taglist, after, time.time() - t0, "existing outputs", sequence_export=sequence_export)
        if not config.overwrite and any(path.exists() for path in [merged_vcf, hap_bam, taglist, qc_json]):
            raise RuntimeError("partial existing Clair3/haplotag output; rerun with --overwrite or remove stale files")
        if config.overwrite:
            cleanup_known_outputs(calling_dir, hap_bam, taglist, qc_json)
            calling_dir.mkdir(parents=True, exist_ok=True)
            hap_dir.mkdir(parents=True, exist_ok=True)
            qc_json.parent.mkdir(parents=True, exist_ok=True)
        run_logged_command(commands[0], log_dir / "clair3.log")
        haplotag_vcf = select_clair3_haplotag_vcf(config, calling_dir, merged_vcf, log_dir)
        if config.haplotag_mode == "external_whatshap":
            ensure_vcf_index(haplotag_vcf, config, log_dir)
            commands[1] = build_whatshap_command(config, sample, alignment, haplotag_vcf, hap_bam, taglist)
            run_logged_command(commands[1], log_dir / "whatshap_haplotag.log")
            finalize_part_output(hap_bam)
            run_logged_command(commands[2], log_dir / "samtools_index_haplotagged.log")
            after = inspect_optional_hp_ps(config, hap_bam)
            if after.get("hp_present", 0) <= 0 and not config.allow_empty_haplotag:
                raise RuntimeError("WhatsHap haplotag completed but inspected reads have no HP tags")
        else:
            after = {}
        sequence_export = export_sequence_vcf_if_requested(config, sample, haplotag_vcf, log_dir)
        qc = build_qc_payload(config, sample, alignment, calling_dir, merged_vcf, hap_bam, taglist, commands, after, "complete", sequence_export=sequence_export)
        write_json(qc_json, qc)
        return summary_row(config, sample, "OK", "created", alignment, calling_dir, merged_vcf, hap_bam, taglist, after, time.time() - t0, "created", sequence_export=sequence_export)
    except Exception as exc:  # noqa: BLE001 - collect per-sample failure in summary TSV.
        return summary_row(config, sample, "FAIL", "failed", alignment, calling_dir, merged_vcf, hap_bam, taglist, {}, time.time() - t0, str(exc))


def build_clair3_command(config: Clair3LongReadConfig, alignment: Path, output_dir: Path, sample: str) -> list[str]:
    entrypoint = [str(config.clair3_bin)]
    if str(config.clair3_bin).endswith(".py"):
        entrypoint = [str(config.python_bin), str(config.clair3_bin)]
    cmd = [
        *entrypoint,
        f"--bam_fn={alignment}",
        f"--ref_fn={config.reference_fasta.resolve()}",
        f"--threads={max(1, config.threads)}",
        f"--platform={normalize_clair3_platform(config.platform)}",
        f"--model_path={config.model_path.resolve()}",
        f"--output={output_dir}",
        f"--sample_name={sample}",
        f"--chunk_size={max(1, config.chunk_size)}",
        f"--min_mq={max(0, config.min_mq)}",
        f"--min_coverage={max(1, config.min_coverage)}",
        f"--qual={max(0, config.qual)}",
        f"--samtools={config.samtools_bin}",
        f"--python={config.python_bin}",
        f"--pypy={config.pypy_bin}",
        f"--parallel={config.parallel_bin}",
        f"--whatshap={config.whatshap_bin}",
        "--use_whatshap_for_final_output_phasing=True",
    ]
    if config.haplotag_mode == "clair3_internal":
        cmd.append("--use_whatshap_for_final_output_haplotagging=True")
    if config.contigs:
        cmd.append(f"--ctg_name={normalize_contig_argument(config.contigs)}")
    if config.bed_fn:
        cmd.append(f"--bed_fn={config.bed_fn.resolve()}")
    if config.enable_long_indel:
        cmd.append("--enable_long_indel=True")
    if config.remove_intermediate_dir:
        cmd.append("--remove_intermediate_dir=True")
    for extra in config.extra_clair3_args or []:
        cmd.append(extra)
    return cmd


def build_whatshap_command(
    config: Clair3LongReadConfig,
    sample: str,
    alignment: Path,
    phased_vcf: Path,
    output_bam: Path,
    taglist: Path,
) -> list[str]:
    part = output_bam.with_suffix(output_bam.suffix + ".part")
    cmd = [
        config.whatshap_bin,
        "haplotag",
        "--sample",
        sample,
        "--reference",
        str(config.reference_fasta.resolve()),
        "--output-threads",
        str(max(1, config.output_threads)),
        "--skip-missing-contigs",
        "--output-haplotag-list",
        str(taglist),
        "-o",
        str(part),
    ]
    if config.ignore_read_groups:
        cmd.append("--ignore-read-groups")
    if config.tag_supplementary:
        cmd.append("--tag-supplementary=copy-primary")
    cmd.extend([str(phased_vcf), str(alignment)])
    return cmd


def existing_outputs_valid(config: Clair3LongReadConfig, merged_vcf: Path, hap_bam: Path, taglist: Path) -> bool:
    if config.overwrite:
        return False
    if not merged_vcf.exists() or not merged_vcf.stat().st_size:
        return False
    if config.haplotag_mode == "external_whatshap":
        if not hap_bam.exists() or not hap_bam.stat().st_size:
            return False
        if not Path(str(hap_bam) + ".bai").exists() and not hap_bam.with_suffix(".bai").exists():
            return False
        if not taglist.exists() or not taglist.stat().st_size:
            return False
    return True


def select_clair3_haplotag_vcf(config: Clair3LongReadConfig, calling_dir: Path, merged_vcf: Path, log_dir: Path) -> Path:
    """Return the VCF that should be used for downstream haplotagging.

    Production order:
    1. phased_merge_output.vcf.gz, when Clair3 created a phased final VCF.
    2. merge_output.vcf.gz, the standard Clair3 final VCF.

    Smoke-only fallback:
    Some restricted BED/sliced-BAM runs can stop after pileup calling and leave
    only pileup.vcf. That file is useful for pipeline plumbing smoke but is not
    an adequate production long-read phasing result. Therefore fallback is only
    enabled when allow_pileup_vcf_fallback is explicitly true.
    """

    phased = calling_dir / "phased_merge_output.vcf.gz"
    if phased.exists() and phased.stat().st_size:
        return phased
    if merged_vcf.exists() and merged_vcf.stat().st_size:
        return merged_vcf
    if not config.allow_pileup_vcf_fallback:
        raise FileNotFoundError(f"Clair3 final VCF not found: {merged_vcf}")

    pileup_gz = calling_dir / "pileup.vcf.gz"
    pileup_plain = calling_dir / "pileup.vcf"
    if pileup_gz.exists() and pileup_gz.stat().st_size:
        shutil.copy2(pileup_gz, merged_vcf)
        for suffix in (".tbi", ".csi"):
            idx = Path(str(pileup_gz) + suffix)
            if idx.exists():
                shutil.copy2(idx, Path(str(merged_vcf) + suffix))
        ensure_vcf_index(merged_vcf, config, log_dir)
        return merged_vcf
    if pileup_plain.exists() and pileup_plain.stat().st_size:
        bgzip_plain_vcf(pileup_plain, merged_vcf, config, log_dir)
        ensure_vcf_index(merged_vcf, config, log_dir)
        return merged_vcf
    raise FileNotFoundError(f"Clair3 final VCF and pileup fallback VCF not found: {merged_vcf}")


def bgzip_plain_vcf(src_vcf: Path, dst_vcf_gz: Path, config: Clair3LongReadConfig, log_dir: Path) -> None:
    bgzip = Path(config.tabix_bin).with_name("bgzip") if "/" in config.tabix_bin else None
    bgzip_bin = str(bgzip) if bgzip and bgzip.exists() else shutil.which("bgzip")
    if not bgzip_bin:
        raise FileNotFoundError("bgzip not found; cannot compress pileup.vcf fallback")
    part = Path(str(dst_vcf_gz) + ".part")
    log_path = log_dir / "bgzip_pileup_fallback.log"
    with part.open("wb") as out, log_path.open("a", encoding="utf-8") as log_handle:
        log_handle.write("\n$ " + shell_join([bgzip_bin, "-c", str(src_vcf)]) + f" > {part}\n")
        log_handle.flush()
        result = subprocess.run([bgzip_bin, "-c", str(src_vcf)], stdout=out, stderr=log_handle)
    if result.returncode != 0:
        part.unlink(missing_ok=True)
        raise RuntimeError(f"bgzip pileup fallback failed ({result.returncode}); see log: {log_path}")
    part.replace(dst_vcf_gz)


def ensure_vcf_index(vcf: Path, config: Clair3LongReadConfig, log_dir: Path) -> None:
    if Path(str(vcf) + ".tbi").exists() or Path(str(vcf) + ".csi").exists():
        return
    tabix = shutil.which(config.tabix_bin) if "/" not in config.tabix_bin else config.tabix_bin
    if tabix:
        run_logged_command([config.tabix_bin, "-f", "-p", "vcf", str(vcf)], log_dir / "tabix_clair3_vcf.log")
        return
    # Fall back to samtools-compatible htslib through bcftools if users put it in PATH.
    bcftools = shutil.which("bcftools")
    if bcftools:
        run_logged_command([bcftools, "index", "-f", "-t", str(vcf)], log_dir / "bcftools_index_clair3_vcf.log")
        return
    raise FileNotFoundError(f"VCF index missing and neither tabix nor bcftools is available: {vcf}")


def export_sequence_vcf_if_requested(
    config: Clair3LongReadConfig,
    sample: str,
    source_vcf: Path,
    log_dir: Path,
) -> dict[str, Any]:
    """Export Clair3's long-read-derived VCF to the FASTA-consensus VCF slot.

    Downstream personal FASTA generation consumes ``<sample>.seq.vcf.gz`` from
    the per-sample VCF directory. For long-read samples, v2 policy uses the
    VCF called from the same long-read BAM/CRAM instead of the Illumina panel.
    This helper performs that path normalization and writes a sidecar metadata
    JSON so stage 40 can verify the source.
    """

    if not config.export_sequence_vcf:
        return {}
    if config.sequence_vcf_dir is None:
        raise ValueError("export_sequence_vcf=True requires sequence_vcf_dir")
    require_file(source_vcf, "Clair3 sequence-source VCF")
    ensure_vcf_index(source_vcf, config, log_dir)
    dest_dir = config.sequence_vcf_dir.resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{sample}.seq.vcf.gz"
    sidecar = dest_dir / f"{sample}.seq.vcf.source.json"
    mode = config.sequence_vcf_export_mode
    action = "created"
    source_resolved = str(source_vcf.resolve())

    if dest.exists() and has_vcf_index(dest) and not config.overwrite:
        if sidecar.exists():
            previous = json.loads(sidecar.read_text(encoding="utf-8"))
            if previous.get("source_vcf") == source_resolved and previous.get("source") == "clair3_longread":
                stats = inspect_vcf_phase_stats(dest)
                validate_sequence_vcf_phase(config, dest, stats)
                return {
                    "source": "clair3_longread",
                    "action": "existing_valid",
                    "dest_vcf": str(dest),
                    "source_vcf": source_resolved,
                    "export_mode": mode,
                    "phase_stats": stats,
                    "sidecar": str(sidecar),
                }
        raise RuntimeError(
            f"existing FASTA sequence VCF does not match the Clair3 source metadata: {dest}. "
            "Rerun with --overwrite or remove/relocate the stale VCF."
        )

    unlink_vcf_bundle(dest)
    sidecar.unlink(missing_ok=True)
    if mode == "symlink":
        dest.symlink_to(source_vcf.resolve())
        link_vcf_index(source_vcf, dest)
        action = "symlinked"
    else:
        part = Path(str(dest) + ".part")
        part.unlink(missing_ok=True)
        shutil.copy2(source_vcf, part)
        part.replace(dest)
        copy_vcf_index(source_vcf, dest)
        action = "copied"
    if not has_vcf_index(dest):
        ensure_vcf_index(dest, config, log_dir)
    stats = inspect_vcf_phase_stats(dest)
    validate_sequence_vcf_phase(config, dest, stats)
    payload = {
        "schema": "mfm_longread_sequence_vcf_export_v1",
        "sample": sample,
        "source": "clair3_longread",
        "source_vcf": source_resolved,
        "dest_vcf": str(dest),
        "export_mode": mode,
        "action": action,
        "reference_fasta_used_for_calling": str(config.reference_fasta.resolve()),
        "note": "This VCF is used by stage 40 bcftools consensus to generate HAP1/HAP2 FASTA for long-read samples.",
        "phase_stats": stats,
    }
    write_json(sidecar, payload)
    return {
        "source": "clair3_longread",
        "action": action,
        "dest_vcf": str(dest),
        "source_vcf": source_resolved,
        "export_mode": mode,
        "phase_stats": stats,
        "sidecar": str(sidecar),
    }


def has_vcf_index(vcf: Path) -> bool:
    return Path(str(vcf) + ".tbi").exists() or Path(str(vcf) + ".csi").exists()


def vcf_index_path(vcf: Path) -> Path | None:
    for suffix in (".tbi", ".csi"):
        candidate = Path(str(vcf) + suffix)
        if candidate.exists():
            return candidate
    return None


def unlink_vcf_bundle(vcf: Path) -> None:
    for path in (vcf, Path(str(vcf) + ".tbi"), Path(str(vcf) + ".csi"), Path(str(vcf) + ".part")):
        path.unlink(missing_ok=True)


def copy_vcf_index(source_vcf: Path, dest_vcf: Path) -> None:
    source_index = vcf_index_path(source_vcf)
    if source_index is None:
        return
    suffix = ".tbi" if str(source_index).endswith(".tbi") else ".csi"
    dest_index = Path(str(dest_vcf) + suffix)
    for stale_suffix in (".tbi", ".csi"):
        Path(str(dest_vcf) + stale_suffix).unlink(missing_ok=True)
    shutil.copy2(source_index, dest_index)


def link_vcf_index(source_vcf: Path, dest_vcf: Path) -> None:
    source_index = vcf_index_path(source_vcf)
    if source_index is None:
        return
    suffix = ".tbi" if str(source_index).endswith(".tbi") else ".csi"
    dest_index = Path(str(dest_vcf) + suffix)
    for stale_suffix in (".tbi", ".csi"):
        Path(str(dest_vcf) + stale_suffix).unlink(missing_ok=True)
    dest_index.symlink_to(source_index.resolve())


def inspect_vcf_phase_stats(vcf: Path) -> dict[str, Any]:
    """Count basic GT phasing statistics from a single-sample VCF/BCF-like text VCF.

    ``bcftools consensus -H 1/2`` needs meaningful haplotype-resolved GTs for
    heterozygous sites. Homozygous alternate sites do not need phasing, so the
    main QC value is ``phased_het_fraction``.
    """

    opener = gzip.open if str(vcf).endswith(".gz") else open
    records = gt_records = het_records = phased_het = phased_gt = nonref_gt = 0
    with opener(vcf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            records += 1
            fmt = parts[8].split(":")
            try:
                gt_index = fmt.index("GT")
            except ValueError:
                gt_index = 0
            sample_values = parts[9].split(":")
            if gt_index >= len(sample_values):
                continue
            gt = sample_values[gt_index]
            if gt in {".", "./.", ".|."}:
                continue
            gt_records += 1
            sep = "|" if "|" in gt else "/" if "/" in gt else ""
            alleles = gt.replace("|", "/").split("/") if sep else [gt]
            called = [a for a in alleles if a != "."]
            if any(a != "0" for a in called):
                nonref_gt += 1
            if "|" in gt:
                phased_gt += 1
            if len(set(called)) > 1:
                het_records += 1
                if "|" in gt:
                    phased_het += 1
    return {
        "records": records,
        "gt_records": gt_records,
        "nonref_gt_records": nonref_gt,
        "het_records": het_records,
        "phased_gt_records": phased_gt,
        "phased_het_records": phased_het,
        "phased_gt_fraction": round(phased_gt / gt_records, 6) if gt_records else 0.0,
        "phased_het_fraction": round(phased_het / het_records, 6) if het_records else 1.0,
    }


def validate_sequence_vcf_phase(config: Clair3LongReadConfig, vcf: Path, stats: dict[str, Any]) -> None:
    if config.allow_unphased_sequence_vcf:
        return
    if int(stats.get("het_records", 0)) > 0 and int(stats.get("phased_het_records", 0)) == 0:
        raise RuntimeError(
            f"sequence VCF has heterozygous records but no phased heterozygous GTs: {vcf}. "
            "This is unsafe for HAP1/HAP2 FASTA consensus. Use a phased Clair3 VCF or set "
            "--allow-unphased-sequence-vcf only for diagnostic runs."
        )


def inspect_optional_hp_ps(config: Clair3LongReadConfig, alignment: Path) -> dict[str, int]:
    if not alignment.exists():
        return {"inspected": 0, "hp_present": 0, "ps_present": 0}
    return inspect_hp_ps(config.samtools_bin, alignment, config.reference_fasta.resolve(), config.inspect_reads)


def cleanup_known_outputs(calling_dir: Path, hap_bam: Path, taglist: Path, qc_json: Path) -> None:
    # 삭제 범위는 이 wrapper가 생성하는 sample-level 산출물로 제한한다.
    # raw BAM/CRAM, reference, FASTA, VCF 입력은 건드리지 않는다.
    if calling_dir.exists():
        for path in calling_dir.iterdir():
            if path.is_file() or path.is_symlink():
                path.unlink()
            elif path.is_dir() and path.name in {"tmp", "log"}:
                shutil.rmtree(path)
    for path in [hap_bam, Path(str(hap_bam) + ".bai"), Path(str(hap_bam) + ".csi"), hap_bam.with_suffix(".bai"), hap_bam.with_suffix(hap_bam.suffix + ".part"), taglist, qc_json]:
        path.unlink(missing_ok=True)



def external_command_env(cmd: list[str]) -> dict[str, str]:
    """Build an execution environment for bundled bioinformatics tools.

    Clair3's Python pipeline receives explicit paths for python/pypy/samtools/
    whatshap/parallel, but some nested Clair3 utilities still call `bgzip` and
    `tabix` by name through the inherited PATH. If the parent process is running
    from `jj_mfm`, those binaries may not be visible even though they exist next
    to the Clair3 executables. We therefore prepend every absolute executable
    directory found in the command arguments to PATH. This does not modify the
    conda environment; it only makes the current subprocess tree resolve the
    intended tools.
    """

    env = os.environ.copy()
    path_dirs: list[str] = []
    for token in cmd:
        candidate = token
        if "=" in token:
            candidate = token.split("=", 1)[1]
        if not candidate.startswith("/"):
            continue
        candidate_path = Path(candidate)
        if candidate_path.exists():
            directory = str(candidate_path.parent)
            if directory not in path_dirs:
                path_dirs.append(directory)
    current_path = env.get("PATH", "")
    env["PATH"] = ":".join(path_dirs + ([current_path] if current_path else []))
    return env


def run_logged_command(cmd: list[str], log_path: Path) -> None:
    """Run an external command while streaming stdout/stderr to a log file.

    Clair3 and WhatsHap can print substantial progress logs. Capturing all output
    in memory is unnecessary and can become expensive during cohort-scale runs.
    """

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("\n$ " + shell_join(cmd) + "\n")
        handle.flush()
        result = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, text=True, env=external_command_env(cmd))
    if result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}); see log: {log_path}")


def finalize_part_output(output_bam: Path) -> None:
    """Promote WhatsHap .part output to the final BAM path when needed."""

    part = output_bam.with_suffix(output_bam.suffix + ".part")
    if part.exists():
        part.replace(output_bam)


def summary_row(
    config: Clair3LongReadConfig,
    sample: str,
    status: str,
    action: str,
    alignment: Path,
    calling_dir: Path,
    merged_vcf: Path,
    hap_bam: Path,
    taglist: Path,
    after: dict[str, int],
    elapsed_sec: float,
    message: str,
    sequence_export: dict[str, Any] | None = None,
) -> Clair3LongReadSummary:
    sequence_export = sequence_export or {}
    phase_stats = sequence_export.get("phase_stats", {}) if isinstance(sequence_export, dict) else {}
    return Clair3LongReadSummary(
        cohort=config.cohort,
        sample=sample,
        status=status,
        action=action,
        input_alignment=str(alignment),
        clair3_dir=str(calling_dir),
        merged_vcf=str(merged_vcf),
        sequence_vcf=str(sequence_export.get("dest_vcf", "")),
        sequence_vcf_source=str(sequence_export.get("source_vcf", "")),
        sequence_vcf_phased_het_fraction=float(phase_stats.get("phased_het_fraction", 0.0)),
        haplotagged_bam=str(hap_bam) if config.haplotag_mode == "external_whatshap" else "",
        taglist=str(taglist) if config.haplotag_mode == "external_whatshap" else "",
        hp_present_after=int(after.get("hp_present", 0)),
        ps_present_after=int(after.get("ps_present", 0)),
        inspected_after=int(after.get("inspected", 0)),
        elapsed_sec=round(elapsed_sec, 3),
        message=message,
    )


def build_qc_payload(
    config: Clair3LongReadConfig,
    sample: str,
    alignment: Path,
    calling_dir: Path,
    merged_vcf: Path,
    hap_bam: Path,
    taglist: Path,
    commands: list[list[str]],
    after: dict[str, int],
    status: str,
    sequence_export: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "schema": "mfm_longread_clair3_sample_qc_v1",
        "cohort": config.cohort,
        "sample": sample,
        "status": status,
        "policy": "Clair3 long-read self-calling followed by WhatsHap haplotagging; no Illumina VCF input for long-read haplotagging or default FASTA generation",
        "sequence_vcf_export": sequence_export or {},
        "input_alignment": str(alignment),
        "reference_fasta": str(config.reference_fasta.resolve()),
        "model_path": str(config.model_path.resolve()),
        "platform_clair3": normalize_clair3_platform(config.platform),
        "clair3_dir": str(calling_dir),
        "merged_vcf": str(merged_vcf),
        "haplotagged_bam": str(hap_bam) if config.haplotag_mode == "external_whatshap" else "",
        "taglist": str(taglist) if config.haplotag_mode == "external_whatshap" else "",
        "haplotag_mode": config.haplotag_mode,
        "tag_supplementary": bool(config.tag_supplementary),
        "ignore_read_groups": bool(config.ignore_read_groups),
        "contigs": config.contigs,
        "bed_fn": str(config.bed_fn.resolve()) if config.bed_fn else "",
        "hp_ps_after": after,
        "commands": [shell_join(cmd) for cmd in commands],
    }


def resolve_alignment(config: Clair3LongReadConfig, sample: str) -> Path:
    if config.alignment is not None:
        return config.alignment.resolve()
    assert config.raw_root is not None
    raw_root = config.raw_root.resolve()
    path = raw_root / config.alignment_pattern.format(sample=sample)
    if path.exists():
        return path
    candidates = [
        raw_root / sample / "wgs" / f"{sample}.hg38.cram",
        raw_root / sample / "wgs" / f"{sample}.ont.hg38.cram",
        raw_root / sample / "wgs" / f"{sample}.hifi.hg38.cram",
        raw_root / sample / "wgs" / f"{sample}.hg38.bam",
        raw_root / sample / "wgs" / f"{sample}.ont.hg38.bam",
        raw_root / sample / "wgs" / f"{sample}.hifi.hg38.bam",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return path.resolve()


def normalize_clair3_platform(platform: str) -> str:
    key = (platform or "").strip().lower().replace(" ", "_")
    return PLATFORM_ALIASES.get(key, key)


def normalize_contig_argument(contigs: str) -> str:
    # CLI의 다른 단계들은 chr1-chr22 같은 range 표현을 허용한다. Clair3는
    # ctg_name에 comma-separated contig string을 받을 수 있으므로 동일 표현을
    # chr1,chr2,...,chr22로 확장해 넘긴다.
    items: list[str] = []
    for token in contigs.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            left, right = token.split("-", 1)
            prefix = "chr" if left.startswith("chr") or right.startswith("chr") else ""
            a = int(left.replace("chr", ""))
            b = int(right.replace("chr", ""))
            step = 1 if a <= b else -1
            items.extend(f"{prefix}{idx}" for idx in range(a, b + step, step))
        else:
            items.append(token)
    return ",".join(items)


def require_file(path: Path, label: str) -> None:
    if not path.exists() or not path.stat().st_size:
        raise FileNotFoundError(f"{label} is missing: {path}")


def require_path_exists(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} is missing: {path}")


def require_alignment_index(alignment: Path) -> None:
    suffix = alignment.suffix.lower()
    candidates = [Path(str(alignment) + ".bai"), Path(str(alignment) + ".csi"), alignment.with_suffix(".bai")]
    if suffix == ".cram":
        candidates.extend([Path(str(alignment) + ".crai"), alignment.with_suffix(".crai")])
    if any(path.exists() for path in candidates):
        return
    raise FileNotFoundError(f"alignment index is missing for {alignment}")


def require_executable(command: str, label: str) -> None:
    if "/" in command:
        path = Path(command)
        if path.exists() and path.is_file():
            return
        raise FileNotFoundError(f"{label} executable is missing: {command}")
    if shutil.which(command):
        return
    raise FileNotFoundError(f"{label} executable is not in PATH: {command}")


def write_summary_tsv(path: Path, rows: list[Clair3LongReadSummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(Clair3LongReadSummary.__dataclass_fields__.keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def shell_join(command: list[str]) -> str:
    return " ".join(shlex.quote(str(item)) for item in command)
