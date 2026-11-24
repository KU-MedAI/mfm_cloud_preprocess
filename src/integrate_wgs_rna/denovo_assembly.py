from __future__ import annotations

import csv
import gzip
import json
import os
import re
import shlex
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .cohort_vcf import expand_chroms, load_sample_ids


DEFAULT_SAMTOOLS = "samtools"
DEFAULT_HIFIASM = "hifiasm"


@dataclass
class DenovoAssemblyConfig:
    """Run long-read de novo assembly benchmark/production-prep.

    This stage intentionally does not pretend that an arbitrary ONT assembly is
    already a chromosome-level HAP1/HAP2 FASTA. It extracts primary long reads,
    runs the selected assembler, converts GFA segments to FASTA for inspection,
    and records whether haplotype-resolved outputs are present. A later
    anchoring/scaffolding stage must decide whether these outputs can be promoted
    to ``individual_fasta/HAP1`` and ``individual_fasta/HAP2``.
    """

    cohort: str
    output_root: Path
    reference_fasta: Path
    platform: str
    sample_list: Path | None = None
    samples: list[str] | None = None
    raw_root: Path | None = None
    alignment: Path | None = None
    alignment_pattern: str = "{sample}/wgs/{sample}.hg38.cram"
    contigs: list[str] | None = None
    label: str = ""
    assembler: str = "hifiasm"
    hifiasm_bin: str = DEFAULT_HIFIASM
    samtools_bin: str = DEFAULT_SAMTOOLS
    threads: int = 48
    fastq_threads: int = 8
    exclude_flags: int = 0x904  # unmapped + secondary + supplementary
    min_contig_len: int = 100_000
    keep_fastq: bool = True
    keep_gfa: bool = True
    keep_bin: bool = False
    overwrite: bool = False
    dry_run: bool = False
    summary_tsv: Path | None = None
    metadata_json: Path | None = None


@dataclass
class DenovoAssemblyRow:
    cohort: str
    sample: str
    status: str
    assembler: str
    platform: str
    label: str
    output_dir: str
    fastq_path: str
    fastq_size_bytes: int
    gfa_count: int
    fasta_count: int
    hap1_fasta: str
    hap2_fasta: str
    primary_fasta: str
    total_contigs: int
    total_bases: int
    n50: int
    elapsed_sec: float
    message: str = ""


def run_denovo_assembly(config: DenovoAssemblyConfig) -> dict[str, Any]:
    started = time.time()
    validate_config(config)
    samples = load_sample_ids(config.sample_list, config.samples)
    if not samples:
        raise ValueError("sample list is empty")
    contigs = expand_chroms(config.contigs or ["chr1-chr22"])
    label = safe_label(config.label or f"{config.assembler}_{normalize_platform(config.platform)}_{scope_label(contigs)}")
    config.output_root.mkdir(parents=True, exist_ok=True)

    rows: list[DenovoAssemblyRow] = []
    for sample in samples:
        rows.append(process_sample(config, sample, contigs, label))

    summary_tsv = (config.summary_tsv or (config.output_root / f"{config.cohort}.denovo_assembly.summary.tsv")).resolve()
    metadata_json = (config.metadata_json or (config.output_root / f"{config.cohort}.denovo_assembly.metadata.json")).resolve()
    write_summary(summary_tsv, rows)
    metadata = {
        "schema": "mfm_denovo_assembly_v1",
        "cohort": config.cohort,
        "sample_count": len(samples),
        "chroms": contigs,
        "label": label,
        "assembler": config.assembler,
        "platform": config.platform,
        "reference_fasta": str(config.reference_fasta),
        "output_root": str(config.output_root),
        "hifiasm_bin": config.hifiasm_bin,
        "samtools_bin": config.samtools_bin,
        "threads": config.threads,
        "fastq_threads": config.fastq_threads,
        "exclude_flags": config.exclude_flags,
        "min_contig_len": config.min_contig_len,
        "keep_fastq": config.keep_fastq,
        "keep_gfa": config.keep_gfa,
        "keep_bin": config.keep_bin,
        "dry_run": config.dry_run,
        "started_at_unix": started,
        "completed_at_unix": time.time(),
        "elapsed_sec": round(time.time() - started, 3),
        "ok_count": sum(1 for row in rows if row.status == "OK"),
        "skip_count": sum(1 for row in rows if row.status == "SKIP"),
        "fail_count": sum(1 for row in rows if row.status == "FAIL"),
        "summary_tsv": str(summary_tsv),
    }
    write_json(metadata_json, metadata)
    failed = [row for row in rows if row.status == "FAIL"]
    if failed:
        raise RuntimeError(f"de novo assembly failed for {len(failed)} sample(s); see {summary_tsv}")
    return metadata


def process_sample(config: DenovoAssemblyConfig, sample: str, contigs: list[str], label: str) -> DenovoAssemblyRow:
    t0 = time.time()
    sample_root = config.output_root / sample
    out_dir = sample_root / "denovo_assembly" / label
    reads_dir = out_dir / "reads"
    asm_dir = out_dir / "assembly"
    fasta_dir = out_dir / "fastas"
    log_dir = out_dir / "logs"
    qc_json = out_dir / "denovo_assembly_qc.json"
    for path in (reads_dir, asm_dir, fasta_dir, log_dir):
        path.mkdir(parents=True, exist_ok=True)

    fastq = reads_dir / f"{sample}.{label}.fastq.gz"
    prefix = asm_dir / f"{sample}.{label}"
    try:
        alignment = resolve_alignment(config, sample)
        existing_gfas = sorted(asm_dir.glob(f"{prefix.name}*.gfa"))
        if qc_json.exists() and existing_gfas and not config.overwrite:
            row = row_from_existing(config, sample, label, out_dir, fastq, fasta_dir, qc_json, t0, "SKIP")
            return row

        if config.dry_run:
            plan = {
                "schema": "mfm_denovo_assembly_plan_v1",
                "sample": sample,
                "alignment": str(alignment),
                "contigs": contigs,
                "fastq": str(fastq),
                "assembler": config.assembler,
                "hifiasm_command": build_hifiasm_command(config, fastq, prefix),
            }
            write_json(qc_json, plan)
            return DenovoAssemblyRow(
                cohort=config.cohort,
                sample=sample,
                status="OK",
                assembler=config.assembler,
                platform=config.platform,
                label=label,
                output_dir=str(out_dir),
                fastq_path=str(fastq),
                fastq_size_bytes=0,
                gfa_count=0,
                fasta_count=0,
                hap1_fasta="",
                hap2_fasta="",
                primary_fasta="",
                total_contigs=0,
                total_bases=0,
                n50=0,
                elapsed_sec=round(time.time() - t0, 3),
                message="dry-run",
            )

        extract_fastq(config, alignment, contigs, fastq, log_dir / "01_extract_fastq.log")
        run_assembler(config, fastq, prefix, log_dir / "02_hifiasm.log")
        gfas = sorted(asm_dir.glob(f"{prefix.name}*.gfa"))
        if not gfas:
            raise RuntimeError(f"hifiasm finished but no GFA files were found under {asm_dir}")
        fasta_outputs = convert_gfas_to_fastas(gfas, fasta_dir, config.min_contig_len)
        role_fastas = classify_fastas(fasta_outputs)
        combined_stats = combined_fasta_stats(list(fasta_outputs.values()))
        qc = {
            "schema": "mfm_denovo_assembly_sample_qc_v1",
            "cohort": config.cohort,
            "sample": sample,
            "label": label,
            "assembler": config.assembler,
            "platform": config.platform,
            "alignment": str(alignment),
            "reference_fasta": str(config.reference_fasta),
            "contigs": contigs,
            "fastq": str(fastq),
            "fastq_size_bytes": file_size(fastq),
            "gfa_files": [str(path) for path in gfas],
            "fasta_outputs": {str(key): str(value) for key, value in fasta_outputs.items()},
            "role_fastas": role_fastas,
            "hap_resolved_ready": bool(role_fastas.get("hap1") and role_fastas.get("hap2")),
            "stats": combined_stats,
            "keep_fastq": config.keep_fastq,
            "keep_gfa": config.keep_gfa,
            "keep_bin": config.keep_bin,
            "elapsed_sec": round(time.time() - t0, 3),
            "status": "pass",
            "note": "HAP1/HAP2 promotion is allowed only when hap_resolved_ready=true and anchoring/scaffolding QC passes.",
        }
        write_json(qc_json, qc)
        maybe_cleanup(config, fastq, gfas, asm_dir, prefix.name)
        return DenovoAssemblyRow(
            cohort=config.cohort,
            sample=sample,
            status="OK",
            assembler=config.assembler,
            platform=config.platform,
            label=label,
            output_dir=str(out_dir),
            fastq_path=str(fastq),
            fastq_size_bytes=file_size(fastq),
            gfa_count=len(gfas),
            fasta_count=len(fasta_outputs),
            hap1_fasta=role_fastas.get("hap1", ""),
            hap2_fasta=role_fastas.get("hap2", ""),
            primary_fasta=role_fastas.get("primary", ""),
            total_contigs=int(combined_stats["contigs"]),
            total_bases=int(combined_stats["total_bases"]),
            n50=int(combined_stats["n50"]),
            elapsed_sec=round(time.time() - t0, 3),
            message="ok" if role_fastas.get("hap1") and role_fastas.get("hap2") else "assembly produced; hap1/hap2 not both detected",
        )
    except Exception as exc:  # noqa: BLE001 - persist row-level failures.
        return DenovoAssemblyRow(
            cohort=config.cohort,
            sample=sample,
            status="FAIL",
            assembler=config.assembler,
            platform=config.platform,
            label=label,
            output_dir=str(out_dir),
            fastq_path=str(fastq),
            fastq_size_bytes=file_size(fastq),
            gfa_count=len(list(asm_dir.glob(f"{prefix.name}*.gfa"))) if asm_dir.exists() else 0,
            fasta_count=len(list(fasta_dir.glob("*.fa"))) if fasta_dir.exists() else 0,
            hap1_fasta="",
            hap2_fasta="",
            primary_fasta="",
            total_contigs=0,
            total_bases=0,
            n50=0,
            elapsed_sec=round(time.time() - t0, 3),
            message=f"{type(exc).__name__}: {exc}",
        )


def validate_config(config: DenovoAssemblyConfig) -> None:
    if normalize_platform(config.platform) not in {"ont", "hifi", "pacbio_hifi", "pacbio_ccs", "ccs"}:
        raise ValueError("de novo assembly stage is intended for long-read platforms only")
    if config.assembler != "hifiasm":
        raise ValueError("currently only --assembler hifiasm is implemented")
    require_exe(config.samtools_bin, "samtools")
    require_exe(config.hifiasm_bin, "hifiasm")
    if not config.reference_fasta.exists():
        raise FileNotFoundError(f"reference FASTA not found: {config.reference_fasta}")
    if config.alignment is None and config.raw_root is None:
        raise ValueError("provide --alignment or --raw-root/--alignment-pattern")
    if config.alignment is not None:
        samples = load_sample_ids(config.sample_list, config.samples)
        if len(samples) != 1:
            raise ValueError("--alignment may be used only with exactly one sample")


def resolve_alignment(config: DenovoAssemblyConfig, sample: str) -> Path:
    if config.alignment is not None:
        path = config.alignment
    else:
        assert config.raw_root is not None
        path = config.raw_root / config.alignment_pattern.format(sample=sample)
    if not path.exists():
        raise FileNotFoundError(f"alignment not found: {path}")
    return path.resolve()


def extract_fastq(config: DenovoAssemblyConfig, alignment: Path, contigs: list[str], fastq: Path, log_path: Path) -> None:
    if fastq.exists() and fastq.stat().st_size > 0 and not config.overwrite:
        return
    fastq.unlink(missing_ok=True)
    regions = " ".join(shlex.quote(contig) for contig in contigs)
    cmd = (
        "set -euo pipefail\n"
        f"{shlex.quote(config.samtools_bin)} view -@ {max(1, config.fastq_threads)} -T {shlex.quote(str(config.reference_fasta))} "
        f"-u -F {int(config.exclude_flags)} {shlex.quote(str(alignment))} {regions} | "
        f"{shlex.quote(config.samtools_bin)} fastq -@ {max(1, config.fastq_threads)} --reference {shlex.quote(str(config.reference_fasta))} "
        f"-n -0 {shlex.quote(str(fastq))} -1 /dev/null -2 /dev/null -s /dev/null -\n"
    )
    run_bash(cmd, log_path)
    if not fastq.exists() or fastq.stat().st_size == 0:
        raise RuntimeError(f"FASTQ extraction produced empty output: {fastq}")


def run_assembler(config: DenovoAssemblyConfig, fastq: Path, prefix: Path, log_path: Path) -> None:
    expected = list(prefix.parent.glob(f"{prefix.name}*.gfa"))
    if expected and not config.overwrite:
        return
    for path in expected:
        path.unlink(missing_ok=True)
    cmd = build_hifiasm_command(config, fastq, prefix)
    run_exec(cmd, log_path)


def build_hifiasm_command(config: DenovoAssemblyConfig, fastq: Path, prefix: Path) -> list[str]:
    cmd = [config.hifiasm_bin, "-t", str(max(1, config.threads)), "-o", str(prefix)]
    if normalize_platform(config.platform) == "ont":
        cmd.append("--ont")
    cmd.append(str(fastq))
    return cmd


def convert_gfas_to_fastas(gfas: list[Path], fasta_dir: Path, min_contig_len: int) -> dict[Path, Path]:
    outputs: dict[Path, Path] = {}
    for gfa in gfas:
        fasta = fasta_dir / f"{gfa.stem}.min{min_contig_len}.fa"
        stats = gfa_to_fasta(gfa, fasta, min_contig_len)
        if stats["contigs"] > 0:
            outputs[gfa] = fasta
    return outputs


def gfa_to_fasta(gfa: Path, fasta: Path, min_len: int) -> dict[str, int]:
    contigs = 0
    bases = 0
    with gfa.open("rt", encoding="utf-8", errors="replace") as src, fasta.open("wt", encoding="utf-8") as out:
        for line in src:
            if not line.startswith("S\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name, seq = parts[1], parts[2]
            if seq == "*" or len(seq) < min_len:
                continue
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i : i + 80] + "\n")
            contigs += 1
            bases += len(seq)
    if contigs == 0:
        fasta.unlink(missing_ok=True)
    return {"contigs": contigs, "bases": bases}


def classify_fastas(fasta_outputs: dict[Path, Path]) -> dict[str, str]:
    roles: dict[str, str] = {}
    for gfa, fasta in fasta_outputs.items():
        name = gfa.name.lower()
        if "hap1" in name and "hap1" not in roles:
            roles["hap1"] = str(fasta)
        elif "hap2" in name and "hap2" not in roles:
            roles["hap2"] = str(fasta)
        elif ("p_ctg" in name or "primary" in name) and "primary" not in roles:
            roles["primary"] = str(fasta)
        elif ("a_ctg" in name or "alternate" in name or "alt" in name) and "alternate" not in roles:
            roles["alternate"] = str(fasta)
    if "primary" not in roles and fasta_outputs:
        first = sorted(fasta_outputs.values())[0]
        roles["primary"] = str(first)
    return roles


def combined_fasta_stats(fastas: list[Path]) -> dict[str, int]:
    lengths: list[int] = []
    for fasta in fastas:
        lengths.extend(read_fasta_lengths(fasta))
    total = sum(lengths)
    return {"contigs": len(lengths), "total_bases": total, "n50": n50(lengths), "max_contig": max(lengths) if lengths else 0}


def read_fasta_lengths(fasta: Path) -> list[int]:
    lengths: list[int] = []
    current = 0
    opener = gzip.open if fasta.suffix == ".gz" else open
    with opener(fasta, "rt") as fh:  # type: ignore[arg-type]
        for line in fh:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    return lengths


def n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    half = sum(lengths) / 2
    acc = 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc >= half:
            return int(length)
    return 0


def row_from_existing(config: DenovoAssemblyConfig, sample: str, label: str, out_dir: Path, fastq: Path, fasta_dir: Path, qc_json: Path, t0: float, status: str) -> DenovoAssemblyRow:
    qc = json.loads(qc_json.read_text(encoding="utf-8"))
    stats = qc.get("stats", {})
    roles = qc.get("role_fastas", {})
    return DenovoAssemblyRow(
        cohort=config.cohort,
        sample=sample,
        status=status,
        assembler=config.assembler,
        platform=config.platform,
        label=label,
        output_dir=str(out_dir),
        fastq_path=str(fastq),
        fastq_size_bytes=file_size(fastq),
        gfa_count=len(qc.get("gfa_files", [])),
        fasta_count=len(list(fasta_dir.glob("*.fa"))),
        hap1_fasta=roles.get("hap1", ""),
        hap2_fasta=roles.get("hap2", ""),
        primary_fasta=roles.get("primary", ""),
        total_contigs=int(stats.get("contigs", 0)),
        total_bases=int(stats.get("total_bases", 0)),
        n50=int(stats.get("n50", 0)),
        elapsed_sec=round(time.time() - t0, 3),
        message="existing valid assembly QC",
    )


def maybe_cleanup(config: DenovoAssemblyConfig, fastq: Path, gfas: list[Path], asm_dir: Path, prefix_name: str) -> None:
    if not config.keep_fastq:
        fastq.unlink(missing_ok=True)
    if not config.keep_gfa:
        for path in gfas:
            path.unlink(missing_ok=True)
    if not config.keep_bin:
        # hifiasm writes large restart/intermediate binary files. The original
        # CRAM/BAM and final GFA/FASTA/QC are enough for this preprocessing
        # benchmark, so remove these by default to avoid avoidable storage growth.
        for pattern in (f"{prefix_name}*.bin",):
            for path in asm_dir.glob(pattern):
                path.unlink(missing_ok=True)


def run_bash(script: str, log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("wt", encoding="utf-8") as log:
        log.write("$ bash -lc <script>\n")
        log.write(script + "\n")
        log.flush()
        proc = subprocess.run(["bash", "-lc", script], stdout=log, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed with exit code {proc.returncode}; see {log_path}")


def run_exec(cmd: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    exe_parent = str(Path(cmd[0]).resolve().parent) if "/" in cmd[0] else ""
    if exe_parent:
        env["PATH"] = exe_parent + os.pathsep + env.get("PATH", "")
    with log_path.open("wt", encoding="utf-8") as log:
        log.write("$ " + " ".join(shlex.quote(item) for item in cmd) + "\n")
        log.flush()
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed with exit code {proc.returncode}; see {log_path}")


def write_summary(path: Path, rows: list[DenovoAssemblyRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(asdict(rows[0]).keys()) if rows else [field for field in DenovoAssemblyRow.__dataclass_fields__]
    with path.open("wt", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def normalize_platform(value: str) -> str:
    value = (value or "").strip().lower().replace("-", "_").replace(" ", "_")
    if value in {"nanopore", "oxfordnanopore", "oxford_nanopore"} or "ont" in value or "nanopore" in value:
        return "ont"
    if "hifi" in value or "ccs" in value or "pacbio" in value:
        return "hifi"
    return value


def scope_label(contigs: list[str]) -> str:
    if not contigs:
        return "all"
    if len(contigs) == 1:
        return contigs[0]
    return f"{contigs[0]}_{contigs[-1]}"


def safe_label(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "denovo"


def require_exe(value: str, label: str) -> None:
    path = Path(value)
    if path.exists() and os.access(path, os.X_OK):
        return
    from shutil import which

    if which(value):
        return
    raise FileNotFoundError(f"{label} executable not found: {value}")


def file_size(path: Path) -> int:
    return path.stat().st_size if path.exists() else 0
