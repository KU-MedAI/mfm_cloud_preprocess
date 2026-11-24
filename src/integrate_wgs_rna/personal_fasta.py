from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from .cpu_executor import cpu_results, WorkerBudget
from pathlib import Path

import pysam

from .cohort_vcf import DEFAULT_BCFTOOLS, expand_chroms, has_vcf_index, load_sample_ids, run_capture, run_command
from .sample_vcf import normalize_sequence_variant_policy, policy_vtypes


DEFAULT_SAMTOOLS = "samtools"
CONSENSUS_ALPHABET = frozenset("ACGTN")


@dataclass
class PersonalFastaConfig:
    """Configuration for process_code-style personal FASTA generation.

    This command corresponds to process_code/03_personal_fasta.sh.
    It consumes the per-sample sequence VCF produced by step 02 and writes
    HAP1/HAP2 FASTA files to the final preprocessed-data layout.
    """

    cohort: str
    seq_vcf_dir: Path
    output_root: Path
    reference_fasta: Path
    chroms: list[str]
    sequence_variant_policy: str
    sample_list: Path | None = None
    samples: list[str] | None = None
    bcftools_bin: str = DEFAULT_BCFTOOLS
    samtools_bin: str = DEFAULT_SAMTOOLS
    jobs: int = 8
    overwrite: bool = False
    sanitize_invalid_bases: bool = False
    summary_tsv: Path | None = None
    metadata_json: Path | None = None
    parallel_contigs: bool = False


@dataclass
class PersonalFastaSummary:
    cohort: str
    sample: str
    sequence_variant_policy: str
    status: str
    individual_fasta_dir: str
    seq_vcf: str
    refok_vcf: str
    chain_required: bool
    expected_nonref_allele_count: int
    applied_variant_count_created: int
    dropped_ref_mismatch_count: int
    created_contig_count: int
    hap1_same_length_count: int
    hap1_diff_length_count: int
    hap2_same_length_count: int
    hap2_diff_length_count: int
    invalid_base_count: int
    sanitized_invalid_base_count: int
    elapsed_sec: float
    message: str = ""


def run_personal_fasta_generation(config: PersonalFastaConfig) -> dict[str, object]:
    """Create HAP1/HAP2 personal FASTA and optional chain files."""

    seq_vcf_dir = config.seq_vcf_dir.resolve()
    output_root = config.output_root.resolve()
    reference = config.reference_fasta.resolve()
    require_reference(reference)
    chroms = expand_chroms(config.chroms)
    samples = load_sample_ids(config.sample_list, config.samples)
    if not samples:
        raise ValueError("sample list is empty")
    config = replace(config,parallel_contigs=len(samples)==1)
    policy = normalize_sequence_variant_policy(config.sequence_variant_policy)
    chain_required = policy != "snp_only"
    output_root.mkdir(parents=True, exist_ok=True)

    started_at = time.time()
    summaries: list[PersonalFastaSummary] = []
    with ThreadPoolExecutor(max_workers=max(1, config.jobs)) as executor:
        futures = {
            executor.submit(process_one_sample, config, sample, seq_vcf_dir, output_root, reference, chroms, policy, chain_required): sample
            for sample in samples
        }
        for future in as_completed(futures):
            summaries.append(future.result())

    summaries.sort(key=lambda row: row.sample)
    summary_tsv = (config.summary_tsv or (output_root / f"{config.cohort}.personal_fasta.summary.tsv")).resolve()
    metadata_json = (config.metadata_json or (output_root / f"{config.cohort}.personal_fasta.metadata.json")).resolve()
    write_summary_tsv(summary_tsv, summaries)
    failed = [row for row in summaries if row.status == "FAIL"]
    completed_at = time.time()
    metadata = {
        "schema": "mfm_personal_fasta_v1",
        "cohort": config.cohort,
        "sample_count": len(samples),
        "chroms": chroms,
        "seq_vcf_dir": str(seq_vcf_dir),
        "output_root": str(output_root),
        "reference_fasta": str(reference),
        "reference_fingerprint": reference_fingerprint(reference),
        "sequence_variant_policy": policy,
        "sequence_vcf_types": policy_vtypes(policy),
        "chain_required": chain_required,
        "invalid_base_policy": "sanitize_to_N" if config.sanitize_invalid_bases else "fail",
        "bcftools_bin": config.bcftools_bin,
        "samtools_bin": config.samtools_bin,
        "jobs": config.jobs,
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
        raise RuntimeError(f"personal FASTA generation failed: {failed_labels}")
    return metadata


def process_one_sample(
    config: PersonalFastaConfig,
    sample: str,
    seq_vcf_dir: Path,
    output_root: Path,
    reference: Path,
    chroms: list[str],
    policy: str,
    chain_required: bool,
) -> PersonalFastaSummary:
    t0 = time.time()
    individual_dir = output_root / sample / "individual_fasta"
    seq_vcf = seq_vcf_dir / f"{sample}.seq.vcf.gz"
    refok_vcf = seq_vcf_dir / f"{sample}.seq.refok.vcf.gz"
    sample_qc_path = individual_dir / "personal_fasta_qc.json"
    try:
        if config.parallel_contigs and os.environ.get('MFM_FASTA_CONTIG_JOBS','auto')=='auto' and os.environ.get('MFM_CPU_POOL_LEASE_FILE'):
            # 직렬 VCF 검사 때는 부모만 유지하고 실제 contig 계산에서 다시 확보한다.
            WorkerBudget(1,1,True).capacity(1,grow=False)
        require_vcf(seq_vcf, "sequence VCF")
        raw_count = vcf_count(config.bcftools_bin, seq_vcf)
        refok_vcf = prepare_refok_vcf(config, seq_vcf, refok_vcf, reference)
        refok_count = vcf_count(config.bcftools_bin, refok_vcf)
        dropped_ref_mismatch = raw_count - refok_count
        expected_nonref_alleles = count_nonref_alleles(config.bcftools_bin, refok_vcf)
        if expected_nonref_alleles <= 0:
            raise RuntimeError("적용 대상 non-reference allele이 0개입니다")

        hap_results: dict[str, list[dict[str, object]]] = {"HAP1": [], "HAP2": []}
        created_contigs = 0
        applied_created = 0
        sanitized_total = 0
        invalid_total = 0
        overlap_skip_by_haplotype = {"HAP1": 0, "HAP2": 0}

        tasks=[]
        for hap_label,hap_number in (("HAP1","1"),("HAP2","2")):
            hap_dir=individual_dir/hap_label;hap_dir.mkdir(parents=True,exist_ok=True)
            for contig in chroms:
                tasks.append((dict(config=config,reference=reference,refok_vcf=refok_vcf,
                    hap_dir=hap_dir,hap_label=hap_label,hap_number=hap_number,
                    contig=contig,chain_required=chain_required),))
        request=os.environ.get('MFM_FASTA_CONTIG_JOBS','auto')
        if request!='auto' and (not request.isdigit() or int(request)<1):
            raise ValueError('MFM_FASTA_CONTIG_JOBS는 auto 또는 양수여야 함')
        budget=int(os.environ.get('MFM_SAMPLE_CPU_BUDGET','1'))
        parallel=config.parallel_contigs and request!='1'
        workers=max(1,budget//2) if parallel else 1
        if request!='auto':workers=min(workers,int(request))
        results={}
        import multiprocessing
        for index,result,error in cpu_results(_personal_contig_task,tasks,workers,
                elastic=bool(parallel and request=='auto' and os.environ.get('MFM_CPU_POOL_LEASE_FILE')),
                max_workers=len(tasks) if request=='auto' and parallel else workers,
                mp_context=multiprocessing.get_context('spawn'),memory_per_worker_bytes=4*1024**3,
                cpus_per_worker=2,wait_memory_bytes=512*1024**2):
            if error is not None:raise error
            results[index]=result
        for index in range(len(tasks)):
            hap_label,result=results[index]
            hap_results[hap_label].append(result)
            created_contigs+=int(bool(result['created']))
            applied_created+=int(result['applied_variant_count'])
            sanitized_total+=int(result['sanitized_invalid_bases'])
            invalid_total+=int(result['invalid_base_count'])
            overlap_skip_by_haplotype[hap_label]+=int(result['overlap_skip_count'])

        length_qc = validate_length_policy(policy, hap_results)
        if created_contigs > 0 and applied_created <= 0:
            raise RuntimeError("새 FASTA를 생성했지만 bcftools consensus 적용 변이 수가 0입니다")

        sample_qc = {
            "schema": "mfm_personal_fasta_sample_qc_v1",
            "cohort": config.cohort,
            "sample": sample,
            "sequence_variant_policy": policy,
            "reference_fasta": str(reference),
            "reference_fingerprint": reference_fingerprint(reference),
            "seq_vcf": str(seq_vcf),
            "refok_vcf": str(refok_vcf),
            "raw_variant_count": raw_count,
            "refok_variant_count": refok_count,
            "dropped_ref_mismatch_count": dropped_ref_mismatch,
            "expected_nonref_allele_count": expected_nonref_alleles,
            "chain_required": chain_required,
            "length_qc": length_qc,
            "hap_results": hap_results,
            "overlap_skip_count_by_haplotype": overlap_skip_by_haplotype,
            "overlap_skip_count_total": sum(overlap_skip_by_haplotype.values()),
            "invalid_base_policy": "sanitize_to_N" if config.sanitize_invalid_bases else "fail",
            "status": "pass",
        }
        write_json(sample_qc_path, sample_qc)
        status = "SKIP" if created_contigs == 0 else "OK"
        return PersonalFastaSummary(
            cohort=config.cohort,
            sample=sample,
            sequence_variant_policy=policy,
            status=status,
            individual_fasta_dir=str(individual_dir),
            seq_vcf=str(seq_vcf),
            refok_vcf=str(refok_vcf),
            chain_required=chain_required,
            expected_nonref_allele_count=expected_nonref_alleles,
            applied_variant_count_created=applied_created,
            dropped_ref_mismatch_count=dropped_ref_mismatch,
            created_contig_count=created_contigs,
            hap1_same_length_count=length_qc["HAP1"]["same_length"],
            hap1_diff_length_count=length_qc["HAP1"]["different_length"],
            hap2_same_length_count=length_qc["HAP2"]["same_length"],
            hap2_diff_length_count=length_qc["HAP2"]["different_length"],
            invalid_base_count=invalid_total,
            sanitized_invalid_base_count=sanitized_total,
            elapsed_sec=round(time.time() - t0, 3),
            message="created" if created_contigs else "existing valid FASTA",
        )
    except Exception as exc:  # noqa: BLE001 - record per-sample failures in summary TSV.
        return PersonalFastaSummary(
            cohort=config.cohort,
            sample=sample,
            sequence_variant_policy=policy,
            status="FAIL",
            individual_fasta_dir=str(individual_dir),
            seq_vcf=str(seq_vcf),
            refok_vcf=str(refok_vcf),
            chain_required=chain_required,
            expected_nonref_allele_count=0,
            applied_variant_count_created=0,
            dropped_ref_mismatch_count=0,
            created_contig_count=0,
            hap1_same_length_count=0,
            hap1_diff_length_count=0,
            hap2_same_length_count=0,
            hap2_diff_length_count=0,
            invalid_base_count=0,
            sanitized_invalid_base_count=0,
            elapsed_sec=round(time.time() - t0, 3),
            message=str(exc),
        )


def _personal_contig_task(arguments):
    """서로다른HAP/chrom의최종FASTA/chain만독립적으로작성함."""
    result=create_or_validate_contig(**arguments)
    return arguments['hap_label'],result


def prepare_refok_vcf(config: PersonalFastaConfig, seq_vcf: Path, refok_vcf: Path, reference: Path) -> Path:
    if refok_vcf.exists() and has_vcf_index(refok_vcf) and not config.overwrite:
        return refok_vcf
    unlink_vcf_bundle(refok_vcf)
    part = Path(str(refok_vcf) + ".part")
    cmd = [
        config.bcftools_bin,
        "norm",
        "--check-ref",
        "x",
        "-f",
        str(reference),
        "-Oz",
        "-o",
        str(part),
        str(seq_vcf),
    ]
    run_command(cmd)
    part.replace(refok_vcf)
    run_command([config.bcftools_bin, "index", "-t", str(refok_vcf)])
    return refok_vcf


def create_or_validate_contig(
    *,
    config: PersonalFastaConfig,
    reference: Path,
    refok_vcf: Path,
    hap_dir: Path,
    hap_label: str,
    hap_number: str,
    contig: str,
    chain_required: bool,
) -> dict[str, object]:
    fasta = hap_dir / f"{contig}.fa"
    chain = hap_dir / f"{contig}.chain" if chain_required else None
    fasta_index = Path(str(fasta) + ".fai")
    stderr_path = hap_dir / f"{contig}.{hap_label}.bcftools_consensus.stderr"
    existing_ok = fasta.exists() and fasta_index.exists() and (not chain_required or (chain is not None and chain.exists() and chain.stat().st_size))
    if existing_ok and not config.overwrite:
        result = validate_contig_outputs(
            config=config,
            reference=reference,
            fasta=fasta,
            chain=chain,
            contig=contig,
            chain_required=chain_required,
        )
        result["created"] = False
        result["applied_variant_count"] = 0
        result["stderr"] = str(stderr_path) if stderr_path.exists() else ""
        stderr_text = stderr_path.read_text(encoding="utf-8", errors="replace") if stderr_path.exists() else ""
        result["overlap_skip_count"] = parse_overlap_skip_count(stderr_text)
        return result

    unlink_fasta_bundle(fasta)
    if chain is not None:
        chain.unlink(missing_ok=True)
    fasta_part = Path(str(fasta) + ".part")
    chain_part = Path(str(chain) + ".part") if chain is not None else None
    applied_count = run_consensus(
        config=config,
        reference=reference,
        refok_vcf=refok_vcf,
        contig=contig,
        hap_number=hap_number,
        fasta_part=fasta_part,
        chain_part=chain_part,
        stderr_path=stderr_path,
    )
    if not fasta_part.exists() or not fasta_part.stat().st_size:
        raise RuntimeError(f"bcftools consensus가 빈 FASTA를 생성했습니다: {fasta_part}")
    os.replace(fasta_part, fasta)
    if chain_part is not None and chain is not None:
        if not chain_part.exists() or not chain_part.stat().st_size:
            raise RuntimeError(f"bcftools consensus가 빈 chain을 생성했습니다: {chain_part}")
        os.replace(chain_part, chain)
    run_command([config.samtools_bin, "faidx", str(fasta)])
    result = validate_contig_outputs(
        config=config,
        reference=reference,
        fasta=fasta,
        chain=chain,
        contig=contig,
        chain_required=chain_required,
    )
    result["created"] = True
    result["applied_variant_count"] = applied_count
    result["stderr"] = str(stderr_path)
    result["overlap_skip_count"] = parse_overlap_skip_count(
        stderr_path.read_text(encoding="utf-8", errors="replace")
    )
    return result


def run_consensus(
    *,
    config: PersonalFastaConfig,
    reference: Path,
    refok_vcf: Path,
    contig: str,
    hap_number: str,
    fasta_part: Path,
    chain_part: Path | None,
    stderr_path: Path,
) -> int:
    cmd_ref = [config.samtools_bin, "faidx", str(reference), contig]
    cmd_consensus = [config.bcftools_bin, "consensus", "-H", hap_number]
    if chain_part is not None:
        cmd_consensus.extend(["--chain", str(chain_part)])
    cmd_consensus.append(str(refok_vcf))
    with fasta_part.open("wb") as fasta_handle, stderr_path.open("wb") as stderr_handle:
        ref_proc = subprocess.Popen(cmd_ref, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert ref_proc.stdout is not None
        consensus_proc = subprocess.Popen(cmd_consensus, stdin=ref_proc.stdout, stdout=fasta_handle, stderr=stderr_handle)
        ref_proc.stdout.close()
        _, ref_stderr = ref_proc.communicate()
        consensus_return = consensus_proc.wait()
    stderr_text = stderr_path.read_text(encoding="utf-8", errors="replace") if stderr_path.exists() else ""
    if ref_proc.returncode != 0 or consensus_return != 0:
        raise RuntimeError(
            f"consensus pipeline failed: samtools={ref_proc.returncode}, bcftools={consensus_return}, "
            f"samtools_stderr={ref_stderr.decode(errors='replace')}, bcftools_stderr={stderr_text}"
        )
    return parse_applied_variant_count(stderr_text)


def validate_contig_outputs(
    *,
    config: PersonalFastaConfig,
    reference: Path,
    fasta: Path,
    chain: Path | None,
    contig: str,
    chain_required: bool,
) -> dict[str, object]:
    reference_length = fasta_length(reference, contig)
    consensus_length = fasta_length(fasta, contig)
    invalid_counts = count_invalid_fasta_bases(fasta)
    sanitized = 0
    if invalid_counts:
        if config.sanitize_invalid_bases:
            sanitized = sanitize_fasta_to_acgtn(fasta)
            run_command([config.samtools_bin, "faidx", str(fasta)])
            invalid_counts = count_invalid_fasta_bases(fasta)
        if invalid_counts:
            raise RuntimeError(f"FASTA에 ACGTN 외 문자/소문자가 있습니다: {fasta}, invalid={invalid_counts}")
    chain_qc = validate_chain(chain, contig, reference_length, consensus_length) if chain_required else None
    return {
        "contig": contig,
        "fasta": str(fasta),
        "chain": str(chain) if chain else "",
        "reference_length": reference_length,
        "consensus_length": consensus_length,
        "same_length": reference_length == consensus_length,
        "invalid_base_count": sum(invalid_counts.values()),
        "invalid_bases": invalid_counts,
        "sanitized_invalid_bases": sanitized,
        "chain_qc": chain_qc,
    }


def validate_length_policy(policy: str, hap_results: dict[str, list[dict[str, object]]]) -> dict[str, dict[str, int]]:
    qc: dict[str, dict[str, int]] = {}
    for hap, rows in hap_results.items():
        same = sum(1 for row in rows if bool(row["same_length"]))
        diff = len(rows) - same
        qc[hap] = {"contigs": len(rows), "same_length": same, "different_length": diff}
        if policy == "snp_only" and diff:
            raise RuntimeError(f"snp_only 정책인데 {hap}에서 reference와 길이가 다른 contig가 있습니다: {diff}")
        if policy == "snv_indel" and diff == 0:
            raise RuntimeError(f"snv_indel 정책인데 {hap}의 모든 contig 길이가 reference와 같습니다")
    return qc


def validate_chain(chain: Path | None, contig: str, reference_length: int, consensus_length: int) -> dict[str, object]:
    if chain is None or not chain.exists() or not chain.stat().st_size:
        raise RuntimeError(f"chain 파일이 없거나 비었습니다: {chain}")
    with chain.open("r", encoding="ascii", errors="replace") as handle:
        header = handle.readline().strip().split()
    if len(header) < 13 or header[0] != "chain":
        raise RuntimeError(f"chain header가 올바르지 않습니다: {chain}")
    qc = {
        "chain": str(chain),
        "target_name": header[2],
        "target_size": int(header[3]),
        "target_strand": header[4],
        "query_name": header[7],
        "query_size": int(header[8]),
        "query_strand": header[9],
    }
    if qc["target_size"] != reference_length:
        raise RuntimeError(f"chain tSize가 reference 길이와 다릅니다: {qc}")
    if qc["query_size"] != consensus_length:
        raise RuntimeError(f"chain qSize가 개인 FASTA 길이와 다릅니다: {qc}")
    if qc["target_strand"] != "+" or qc["query_strand"] != "+":
        raise RuntimeError(f"chain strand가 plus/plus가 아닙니다: {qc}")
    if qc["target_name"] != contig or qc["query_name"] != contig:
        raise RuntimeError(f"chain contig 이름이 예상과 다릅니다: {qc}, expected={contig}")
    qc["status"] = "pass"
    return qc


def require_reference(reference: Path) -> None:
    if not reference.exists() or not reference.stat().st_size:
        raise FileNotFoundError(f"reference FASTA가 없습니다: {reference}")
    if not Path(str(reference) + ".fai").exists():
        raise FileNotFoundError(f"reference FASTA index가 없습니다: {reference}.fai")


def require_vcf(vcf: Path, label: str) -> None:
    if not vcf.exists() or not vcf.stat().st_size:
        raise FileNotFoundError(f"{label}가 없습니다: {vcf}")
    if not has_vcf_index(vcf):
        raise FileNotFoundError(f"{label} index가 없습니다: {vcf}.tbi 또는 .csi")


def vcf_count(bcftools_bin: str, vcf: Path) -> int:
    return int(run_capture([bcftools_bin, "index", "-n", str(vcf)]).strip())


def count_nonref_alleles(bcftools_bin: str, vcf: Path) -> int:
    output = run_capture([bcftools_bin, "query", "-f", "[%GT]\n", str(vcf)])
    count = 0
    for gt in output.splitlines():
        for allele in re.split(r"[|/]", gt.strip()):
            if allele and allele not in {"0", "."}:
                count += 1
    return count


def parse_applied_variant_count(stderr_text: str) -> int:
    matches = re.findall(r"Applied\s+(\d+)", stderr_text)
    return int(matches[-1]) if matches else 0



def parse_overlap_skip_count(stderr_text: str) -> int:
    return len(re.findall(r"overlaps with another variant", stderr_text, flags=re.IGNORECASE))


def reference_fingerprint(reference: Path) -> dict[str, object]:
    """Cheap immutable-enough contract: resolved FASTA stat plus SHA-256 of its FAI."""
    reference = reference.resolve()
    fai = Path(str(reference) + ".fai")
    digest = hashlib.sha256()
    with fai.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    stat = reference.stat()
    return {
        "path": str(reference),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "fai_path": str(fai.resolve()),
        "fai_size_bytes": fai.stat().st_size,
        "fai_sha256": digest.hexdigest(),
    }

def fasta_length(fasta: Path, contig: str) -> int:
    with pysam.FastaFile(str(fasta)) as handle:
        if contig not in handle.references:
            raise RuntimeError(f"{contig}가 FASTA에 없습니다: {fasta}")
        return int(handle.get_reference_length(contig))


def count_invalid_fasta_bases(fasta: Path) -> dict[str, int]:
    invalid: dict[str, int] = {}
    with fasta.open("r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            for base in line.rstrip("\r\n"):
                if base not in CONSENSUS_ALPHABET:
                    invalid[base] = invalid.get(base, 0) + 1
    return invalid


def sanitize_fasta_to_acgtn(fasta: Path) -> int:
    partial = Path(str(fasta) + ".sanitize.part")
    count = 0
    with fasta.open("r", encoding="ascii", errors="replace") as src, partial.open("w", encoding="ascii") as dst:
        for line in src:
            if line.startswith(">"):
                dst.write(line)
                continue
            newline = "\n" if line.endswith("\n") else ""
            out = []
            for base in line.rstrip("\r\n").upper():
                if base in CONSENSUS_ALPHABET:
                    out.append(base)
                else:
                    out.append("N")
                    count += 1
            dst.write("".join(out) + newline)
    os.replace(partial, fasta)
    Path(str(fasta) + ".fai").unlink(missing_ok=True)
    return count


def unlink_vcf_bundle(vcf: Path) -> None:
    for path in (vcf, Path(str(vcf) + ".tbi"), Path(str(vcf) + ".csi"), Path(str(vcf) + ".part")):
        path.unlink(missing_ok=True)


def unlink_fasta_bundle(fasta: Path) -> None:
    for path in (fasta, Path(str(fasta) + ".fai"), Path(str(fasta) + ".part")):
        path.unlink(missing_ok=True)


def write_summary_tsv(path: Path, rows: list[PersonalFastaSummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(PersonalFastaSummary.__dataclass_fields__.keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

