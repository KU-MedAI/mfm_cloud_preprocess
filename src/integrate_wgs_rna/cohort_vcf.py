from __future__ import annotations

import csv
import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import pysam


DEFAULT_PANEL_BASE = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
)
DEFAULT_PANEL_FILE_TEMPLATE = (
    "1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)
DEFAULT_BCFTOOLS = "bcftools"


@dataclass
class CohortVcfConfig:
    """Configuration for process_code-style cohort common VCF extraction.

    This step corresponds to process_code/01_extract_panel_vcf.sh.  It creates
    one VCF per chromosome containing all samples in a cohort, instead of
    repeatedly streaming the remote panel once per sample.
    """

    cohort: str
    output_dir: Path
    chroms: list[str]
    sample_list: Path | None = None
    samples: list[str] | None = None
    panel_base: str = DEFAULT_PANEL_BASE
    panel_file_template: str = DEFAULT_PANEL_FILE_TEMPLATE
    bcftools_bin: str = DEFAULT_BCFTOOLS
    filter_status: str = "PASS,."
    min_ac: int = 1
    jobs: int = 6
    threads_per_chrom: int = 2
    overwrite: bool = False
    summary_tsv: Path | None = None
    metadata_json: Path | None = None
    idx_cache_dir: Path | None = None


@dataclass
class ChromVcfSummary:
    cohort: str
    contig: str
    source: str
    output_vcf: str
    status: str
    sample_count: int
    variant_count: int
    all_ref_ac0_count: int
    first_record_phased: bool
    size_bytes: int
    elapsed_sec: float
    message: str = ""


def run_cohort_vcf_extract(config: CohortVcfConfig) -> dict[str, object]:
    """Extract cohort common VCFs and return a run summary.

    The function intentionally accepts a panel URL/path template instead of
    hard-coding one data source.  The current default matches process_code's
    2022 1kGP Illumina phased panel, but future ONT-specific VCF sources can be
    selected by changing only CLI arguments.
    """

    output_dir = config.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    idx_cache_dir = (config.idx_cache_dir or (output_dir / ".idxcache")).resolve()
    idx_cache_dir.mkdir(parents=True, exist_ok=True)
    sample_ids = load_sample_ids(config.sample_list, config.samples)
    if not sample_ids:
        raise ValueError("sample list is empty")
    sample_list_path = materialize_sample_list(output_dir, config.cohort, sample_ids, config.sample_list)
    chroms = expand_chroms(config.chroms)
    summary_tsv = (config.summary_tsv or (output_dir / f"{config.cohort}.common_vcf.summary.tsv")).resolve()
    metadata_json = (config.metadata_json or (output_dir / f"{config.cohort}.common_vcf.metadata.json")).resolve()

    started_at = time.time()
    summaries: list[ChromVcfSummary] = []
    with ThreadPoolExecutor(max_workers=max(1, config.jobs)) as executor:
        futures = {
            executor.submit(
                extract_one_chrom,
                config,
                contig,
                sample_list_path,
                len(sample_ids),
                output_dir,
                idx_cache_dir,
            ): contig
            for contig in chroms
        }
        for future in as_completed(futures):
            summaries.append(future.result())

    summaries.sort(key=lambda row: chrom_sort_key(row.contig))
    write_summary_tsv(summary_tsv, summaries)
    failed = [row for row in summaries if row.status == "FAIL"]
    completed_at = time.time()
    metadata = {
        "schema": "mfm_cohort_common_vcf_v1",
        "cohort": config.cohort,
        "sample_count": len(sample_ids),
        "chroms": chroms,
        "output_dir": str(output_dir),
        "sample_list": str(sample_list_path),
        "panel_base": config.panel_base,
        "panel_file_template": config.panel_file_template,
        "bcftools_bin": config.bcftools_bin,
        "filter_status": config.filter_status,
        "min_ac": config.min_ac,
        "jobs": config.jobs,
        "threads_per_chrom": config.threads_per_chrom,
        "overwrite": config.overwrite,
        "summary_tsv": str(summary_tsv),
        "started_at_unix": started_at,
        "completed_at_unix": completed_at,
        "elapsed_sec": round(completed_at - started_at, 3),
        "ok_count": sum(1 for row in summaries if row.status == "OK"),
        "skip_count": sum(1 for row in summaries if row.status == "SKIP"),
        "fail_count": len(failed),
        "total_size_bytes": sum(row.size_bytes for row in summaries if row.status != "FAIL"),
    }
    write_json(metadata_json, metadata)
    if failed:
        failed_labels = ", ".join(f"{row.contig}:{row.message}" for row in failed[:5])
        raise RuntimeError(f"cohort common VCF extraction failed: {failed_labels}")
    return metadata


def extract_one_chrom(
    config: CohortVcfConfig,
    contig: str,
    sample_list_path: Path,
    expected_sample_count: int,
    output_dir: Path,
    idx_cache_dir: Path,
) -> ChromVcfSummary:
    t0 = time.time()
    source = panel_source(config.panel_base, config.panel_file_template, contig)
    output_vcf = output_dir / f"{config.cohort}.{contig}.vcf.gz"
    try:
        if output_vcf.exists() and has_vcf_index(output_vcf) and not config.overwrite:
            summary = validate_output_vcf(config, output_vcf, expected_sample_count)
            return ChromVcfSummary(
                cohort=config.cohort,
                contig=contig,
                source=source,
                output_vcf=str(output_vcf),
                status="SKIP",
                sample_count=summary["sample_count"],
                variant_count=summary["variant_count"],
                all_ref_ac0_count=summary["all_ref_ac0_count"],
                first_record_phased=summary["first_record_phased"],
                size_bytes=output_vcf.stat().st_size,
                elapsed_sec=round(time.time() - t0, 3),
                message="existing valid VCF",
            )

        part = output_vcf.with_suffix(output_vcf.suffix + ".part")
        part.unlink(missing_ok=True)
        cmd = [
            config.bcftools_bin,
            "view",
            "-r",
            contig,
            "-S",
            str(sample_list_path),
            "-f",
            config.filter_status,
            "-c",
            str(config.min_ac),
            "--threads",
            str(config.threads_per_chrom),
            "-Oz",
            "-o",
            str(part),
            source,
        ]
        run_command(cmd, cwd=idx_cache_dir)
        part.replace(output_vcf)
        run_command([config.bcftools_bin, "index", "-t", str(output_vcf)], cwd=idx_cache_dir)
        summary = validate_output_vcf(config, output_vcf, expected_sample_count)
        return ChromVcfSummary(
            cohort=config.cohort,
            contig=contig,
            source=source,
            output_vcf=str(output_vcf),
            status="OK",
            sample_count=summary["sample_count"],
            variant_count=summary["variant_count"],
            all_ref_ac0_count=summary["all_ref_ac0_count"],
            first_record_phased=summary["first_record_phased"],
            size_bytes=output_vcf.stat().st_size,
            elapsed_sec=round(time.time() - t0, 3),
            message="created",
        )
    except Exception as exc:  # noqa: BLE001 - record per-chromosome failure in summary TSV.
        return ChromVcfSummary(
            cohort=config.cohort,
            contig=contig,
            source=source,
            output_vcf=str(output_vcf),
            status="FAIL",
            sample_count=0,
            variant_count=0,
            all_ref_ac0_count=0,
            first_record_phased=False,
            size_bytes=output_vcf.stat().st_size if output_vcf.exists() else 0,
            elapsed_sec=round(time.time() - t0, 3),
            message=str(exc),
        )


def validate_output_vcf(config: CohortVcfConfig, vcf_path: Path, expected_sample_count: int) -> dict[str, int | bool]:
    if not vcf_path.exists() or not vcf_path.stat().st_size:
        raise FileNotFoundError(f"VCF is empty or missing: {vcf_path}")
    if not has_vcf_index(vcf_path):
        raise FileNotFoundError(f"VCF index missing: {vcf_path}.tbi or .csi")

    variant_count = int(run_capture([config.bcftools_bin, "index", "-n", str(vcf_path)]).strip())
    if variant_count <= 0:
        raise RuntimeError(f"variant count is zero: {vcf_path}")

    sample_count = len(run_capture([config.bcftools_bin, "query", "-l", str(vcf_path)]).splitlines())
    if sample_count != expected_sample_count:
        raise RuntimeError(f"sample count mismatch: {sample_count} != {expected_sample_count}")

    all_ref_ac0_count = count_all_ref_ac0(config.bcftools_bin, vcf_path)
    if all_ref_ac0_count != 0:
        raise RuntimeError(f"AC=0 records remain after -c {config.min_ac}: {all_ref_ac0_count}")

    first_record_phased = first_record_is_phased(vcf_path)
    if not first_record_phased:
        raise RuntimeError("first VCF record is not phased; expected GT with '|'")

    return {
        "sample_count": sample_count,
        "variant_count": variant_count,
        "all_ref_ac0_count": all_ref_ac0_count,
        "first_record_phased": first_record_phased,
    }


def load_sample_ids(sample_list: Path | None, samples: Iterable[str] | None) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    if sample_list is not None:
        with sample_list.open("r", encoding="utf-8") as handle:
            for line in handle:
                value = line.strip()
                if not value or value.startswith("#"):
                    continue
                sample = value.split()[0]
                if sample not in seen:
                    seen.add(sample)
                    ordered.append(sample)
    for sample in samples or []:
        value = sample.strip()
        if value and value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


def materialize_sample_list(output_dir: Path, cohort: str, sample_ids: list[str], source_sample_list: Path | None) -> Path:
    if source_sample_list is not None:
        return source_sample_list.resolve()
    path = output_dir / f"{cohort}.samples.txt"
    path.write_text("\n".join(sample_ids) + "\n", encoding="utf-8")
    return path.resolve()


def expand_chroms(values: list[str]) -> list[str]:
    if not values:
        values = ["chr1-chr22"]
    chroms: list[str] = []
    for raw in values:
        for item in raw.split(","):
            token = item.strip()
            if not token:
                continue
            chroms.extend(expand_one_chrom_token(token))
    seen: set[str] = set()
    result: list[str] = []
    for chrom in chroms:
        if chrom not in seen:
            seen.add(chrom)
            result.append(chrom)
    return result


def expand_one_chrom_token(token: str) -> list[str]:
    if "-" not in token:
        return [normalize_contig(token)]
    left, right = token.split("-", 1)
    left_n = chrom_number(normalize_contig(left))
    right_n = chrom_number(normalize_contig(right))
    if left_n is None or right_n is None:
        raise ValueError(f"only numeric chromosome ranges are supported: {token}")
    step = 1 if left_n <= right_n else -1
    return [f"chr{i}" for i in range(left_n, right_n + step, step)]


def normalize_contig(value: str) -> str:
    value = value.strip()
    if value.lower().startswith("chr"):
        return "chr" + value[3:]
    return "chr" + value


def chrom_number(contig: str) -> int | None:
    label = contig.removeprefix("chr")
    return int(label) if label.isdigit() else None


def chrom_sort_key(contig: str) -> tuple[int, str]:
    number = chrom_number(contig)
    if number is not None:
        return (number, "")
    return (10_000, contig)


def panel_source(panel_base: str, template: str, contig: str) -> str:
    chrom = contig.removeprefix("chr")
    if "{" in template:
        filename = template.format(chrom=chrom, contig=contig, chromosome=contig, chrom_num=chrom)
    elif "%s" in template:
        filename = template % chrom
    else:
        filename = template
    if not panel_base:
        return filename
    base = panel_base.rstrip("/")
    if "://" in base:
        return f"{base}/{filename}"
    return str(Path(base) / filename)


def has_vcf_index(vcf_path: Path) -> bool:
    return Path(str(vcf_path) + ".tbi").exists() or Path(str(vcf_path) + ".csi").exists()


def count_all_ref_ac0(bcftools_bin: str, vcf_path: Path) -> int:
    cmd = [bcftools_bin, "query", "-f", "%AC\n", str(vcf_path)]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    count = 0
    for line in process.stdout:
        if line.strip() == "0":
            count += 1
    _, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(f"bcftools query AC failed: {stderr.strip()}")
    return count


def first_record_is_phased(vcf_path: Path) -> bool:
    with pysam.VariantFile(str(vcf_path)) as handle:
        if not handle.header.samples:
            raise RuntimeError(f"VCF has no samples: {vcf_path}")
        sample = next(iter(handle.header.samples))
        for record in handle:
            gt = record.samples[sample].get("GT")
            phased = bool(record.samples[sample].phased)
            return bool(gt) and phased
    raise RuntimeError(f"VCF has no records: {vcf_path}")


def run_command(cmd: list[str], cwd: Path | None = None) -> None:
    result = subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )


def run_capture(cmd: list[str]) -> str:
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result.stdout


def write_summary_tsv(path: Path, rows: list[ChromVcfSummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()) if rows else list(ChromVcfSummary.__dataclass_fields__.keys()), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

