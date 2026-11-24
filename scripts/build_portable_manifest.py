#!/usr/bin/env python3
"""동결 registry를 외부 서버용 전처리 manifest와 sample list로 변환한다."""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from integrate_wgs_rna.library_provenance import (  # noqa: E402
    PROVENANCE_COLUMNS,
    project_manifest_fields,
)

BASE_COLUMNS = [
    "sample_id", "platform", "read_type", "dataset_role", "project_id", "cohort",
    "alignment_path", "raw_root", "rna_path", "final_root", "work_root",
    "common_vcf_dir", "persample_vcf_dir", "rna_bigwigs",
]
EXTRA_COLUMNS = [
    "data_group", "source_sample_id", *PROVENANCE_COLUMNS,
    "read_files_manifest", "source_decode_reference_fasta",
    "analysis_reference_fasta", "vcf_source", "public_vcf_path",
    "public_vcf_sample_id", "vcf_reference_fasta", "vcf_callset", "phase_source",
    "wgs_url", "wgs_index_url", "rna_url", "rna_index_url",
    "wgs_input_bytes", "wgs_input_gib", "wgs_index_bytes",
    "rna_input_bytes", "rna_input_gib", "rna_index_bytes",
    "vcf_input_bytes", "reference_input_bytes",
    "known_input_total_bytes", "known_input_total_gib",
    "size_completeness", "size_basis", "reference_id",
]
FILE_COLUMNS = [
    "sample_id", "run_id", "file_role", "format", "alignment_state", "path",
    "source_url", "size_bytes", "md5", "source_reference_fasta", "read_group_id",
]

# BAM @RG SM의 공개 source alias는 output sample ID와 구분한다. 이 목록은 evidence를
# 확인한 경우만 추가하며, VCF의 public_vcf_sample_id에는 사용하지 않는다.
SOURCE_SAMPLE_ID_OVERRIDES = {"NA12878": "HG001"}

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)

def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    temporary.write_text(text)
    os.replace(temporary, path)

def basename(url: str) -> str:
    return Path(urlparse(url).path).name

def as_int(value: str | None) -> int | None:
    value = (value or "").strip()
    return int(value) if value.isdigit() else None

def gib(value: int) -> str:
    return f"{value / 1073741824:.3f}"

def safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)

def is_remote_url(value: str) -> bool:
    """실제 다운로드 URL만 총량 집계에 포함하고 registry 설명값은 제외한다."""
    return urlparse(value).scheme.lower() in {"http", "https", "ftp", "s3"}

def local_path(kind: str, project: str, sample: str, url: str) -> str:
    return f"{kind}/{safe(project)}/{sample}/{basename(url)}"

def platform_fields(selected: str) -> tuple[str, str]:
    key = selected.lower()
    if "illumina" in key:
        return "illumina", "short_read"
    if "clr" in key:
        return "clr", "long_read"
    return "hifi", "long_read"

def file_format(row: dict[str, str]) -> tuple[str, str, str]:
    key = row.get("file_type", "").lower()
    state = row.get("alignment_state", "").lower()
    if "cram" in key:
        return "cram", "aligned", "alignment"
    if state == "aligned" and "bam" in key:
        return "bam", "aligned", "alignment"
    if "bam" in key:
        return "ubam", "unaligned", "alignment"
    return ("fastq_gz" if key.endswith("gz") or "fastq" in key else key), "unaligned", "long_reads"

def main() -> int:
    here = Path(__file__).resolve().parents[1]
    share = here / "share_package"
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", type=Path, default=share / "samples_3212.tsv")
    parser.add_argument("--pacbio-files", type=Path, default=share / "manifests/pacbio_file_manifest.tsv")
    parser.add_argument("--vcf-files", type=Path, default=share / "manifests/public_vcf_file_manifest.tsv")
    parser.add_argument("--reference-files", type=Path, default=share / "manifests/reference_file_manifest.tsv")
    parser.add_argument("--out-dir", type=Path, default=here)
    args = parser.parse_args()

    registry = read_tsv(args.registry)
    pacbio_rows = read_tsv(args.pacbio_files)
    vcf_rows = read_tsv(args.vcf_files)
    reference_rows = read_tsv(args.reference_files)
    if len(registry) != len({row["sample_id"] for row in registry}):
        raise SystemExit("registry sample_id가 중복됨")

    pacbio_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in pacbio_rows:
        pacbio_by_sample[row["sample_id"]].append(row)
    vcf_by_project_sample: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in vcf_rows:
        vcf_by_project_sample[(row["project"], row.get("sample_id", ""))].append(row)
    refs_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in reference_rows:
        refs_by_id[row["reference_id"]].append(row)

    manifest_rows: list[dict[str, object]] = []
    read_file_rows: list[dict[str, object]] = []
    lists: dict[str, list[str]] = {"all_3212": [], "wgs_only_2763": [], "wgs_rna_pair_449": []}
    project_lists: dict[str, list[str]] = defaultdict(list)
    unique_resource_urls: dict[str, int] = {}

    for source in sorted(registry, key=lambda row: row["sample_id"]):
        sample = source["sample_id"]
        project = source["selected_project"]
        platform, read_type = platform_fields(source["selected_platform"])
        is_pair = source["wgs_rna_pair_status"] == "RNA_BAM_pair"
        group = "wgs_rna_pair" if is_pair else "wgs"
        dataset_role = "post_training" if is_pair else "pre_training"
        wgs_url = source["wgs_primary_url"]
        wgs_index_url = source["wgs_index_url"]
        rna_url = source["rna_primary_url"] if is_pair else ""
        rna_index_url = source["rna_index_url"] if is_pair else ""
        wgs_bytes = as_int(source["wgs_primary_bytes"])
        rna_bytes = as_int(source["rna_primary_bytes"]) if is_pair else 0
        if wgs_bytes is None or rna_bytes is None:
            raise SystemExit(f"primary size 누락: {sample}")

        reference_id = source.get("vcf_reference", "") if platform != "illumina" else "ILLUMINA_GRCH38_DECODE"
        if platform == "illumina":
            source_decode_reference_path = "references/illumina_decode.fa"
            analysis_reference_path = source_decode_reference_path
            alignment_path = local_path("wgs", project, sample, wgs_url)
            read_files_manifest = ""
            vcf_source = "illumina_panel"
            public_vcf_paths: list[str] = []
            vcf_callset = source.get("vcf_callset", "1KGP_2022")
            common_vcf_dir = "vcf/1KGP_ILLUMINA"
            wgs_index_bytes: int | None = None
        else:
            # CRAM/BAM decode에는 원 source reference를, public VCF/HAP consensus와
            # raw-read canonical alignment에는 allele-compatible 22-autosome FASTA를 쓴다.
            source_decode_reference_path = "references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
            analysis_reference_path = "references/hg38_autosomes.fa"
            common_vcf_dir = ""
            reads = pacbio_by_sample.get(sample, [])
            if not reads:
                raise SystemExit(f"PacBio file manifest 누락: {sample}")
            aligned = [row for row in reads if row.get("alignment_state") == "aligned" and "bam" in row.get("file_type", "").lower()]
            if len(reads) == 1 and len(aligned) == 1:
                alignment_path = local_path("wgs", project, sample, reads[0]["file_url"])
                read_files_manifest = ""
            else:
                alignment_path = ""
                read_files_manifest = "read_files.tsv"
            wgs_index_values = [as_int(row.get("index_size_bytes")) for row in reads if row.get("index_url")]
            if not any(row.get("index_url") for row in reads):
                wgs_index_bytes = 0
            else:
                wgs_index_bytes = sum(value for value in wgs_index_values if value is not None) if wgs_index_values and all(value is not None for value in wgs_index_values) else None
            for index, row in enumerate(reads, start=1):
                fmt, state, role = file_format(row)
                run_id = row.get("run_accession") or row.get("source_record_id") or f"{sample}_{index}"
                read_file_rows.append({
                    "sample_id": sample, "run_id": safe(run_id), "file_role": role,
                    "format": fmt, "alignment_state": state,
                    "path": local_path("wgs", project, sample, row["file_url"]),
                    "source_url": row["file_url"], "size_bytes": row["size_bytes"],
                    "md5": row.get("checksum_value", ""),
                    "source_reference_fasta": source_decode_reference_path if fmt == "cram" else "",
                    "read_group_id": safe(run_id),
                })
                unique_resource_urls[row["file_url"]] = int(row["size_bytes"])
                if row.get("index_url") and as_int(row.get("index_size_bytes")) is not None:
                    unique_resource_urls[row["index_url"]] = int(row["index_size_bytes"])

            candidates = vcf_by_project_sample.get((project, sample), []) + vcf_by_project_sample.get((project, ""), [])
            public_vcf_rows = [row for row in candidates if row.get("file_type") == "vcf"]
            if not public_vcf_rows:
                raise SystemExit(f"public VCF 누락: {sample} project={project}")
            public_vcf_paths = [f"public_vcf/{safe(project)}/{basename(row['file_url'])}" for row in public_vcf_rows]
            vcf_source = "public"
            vcf_callset = ";".join(sorted({row.get("caller", "unknown") for row in public_vcf_rows}))

        selected_vcf_rows = [] if platform == "illumina" else (vcf_by_project_sample.get((project, sample), []) + vcf_by_project_sample.get((project, ""), []))
        vcf_bytes = sum(as_int(row.get("size_bytes")) or 0 for row in selected_vcf_rows)
        for row in selected_vcf_rows:
            if as_int(row.get("size_bytes")) is not None:
                unique_resource_urls[row["file_url"]] = int(row["size_bytes"])
        ref_rows = refs_by_id.get(source.get("vcf_reference", ""), []) if platform != "illumina" else []
        reference_bytes = sum(as_int(row.get("size_bytes")) or 0 for row in ref_rows)
        for row in ref_rows:
            if as_int(row.get("size_bytes")) is not None:
                unique_resource_urls[row["file_url"]] = int(row["size_bytes"])

        rna_index_bytes = None if is_pair else 0
        known_total = wgs_bytes + rna_bytes + vcf_bytes + reference_bytes
        if wgs_index_bytes is not None:
            known_total += wgs_index_bytes
        if rna_index_bytes is not None:
            known_total += rna_index_bytes
        complete = "full_known" if wgs_index_bytes is not None and rna_index_bytes is not None else "primary_only"
        row: dict[str, object] = {column: "" for column in BASE_COLUMNS + EXTRA_COLUMNS}
        row.update({
            "sample_id": sample, "platform": platform, "read_type": read_type,
            "dataset_role": dataset_role, "project_id": project, "cohort": project,
            "alignment_path": alignment_path, "raw_root": "",
            "rna_path": local_path("rna", source.get("rna_project", "GEUVADIS"), sample, rna_url) if is_pair else "",
            "final_root": "", "work_root": "", "common_vcf_dir": common_vcf_dir,
            "persample_vcf_dir": "", "rna_bigwigs": "", "data_group": group,
            "read_files_manifest": read_files_manifest,
            "source_sample_id": SOURCE_SAMPLE_ID_OVERRIDES.get(sample, sample),
            "source_decode_reference_fasta": source_decode_reference_path,
            "analysis_reference_fasta": analysis_reference_path,
            "vcf_source": vcf_source, "public_vcf_path": ";".join(public_vcf_paths),
            "public_vcf_sample_id": sample if platform != "illumina" else "",
            "vcf_reference_fasta": analysis_reference_path if platform != "illumina" else "",
            "vcf_callset": vcf_callset, "phase_source": "provided",
            "wgs_url": wgs_url, "wgs_index_url": wgs_index_url,
            "rna_url": rna_url, "rna_index_url": rna_index_url,
            "wgs_input_bytes": wgs_bytes, "wgs_input_gib": gib(wgs_bytes),
            "wgs_index_bytes": "" if wgs_index_bytes is None else wgs_index_bytes,
            "rna_input_bytes": rna_bytes, "rna_input_gib": gib(rna_bytes),
            "rna_index_bytes": "" if rna_index_bytes is None else rna_index_bytes,
            "vcf_input_bytes": vcf_bytes, "reference_input_bytes": reference_bytes,
            "known_input_total_bytes": known_total, "known_input_total_gib": gib(known_total),
            "size_completeness": complete, "size_basis": source.get("wgs_size_basis", "source_manifest_bytes"),
            "reference_id": reference_id,
        })
        row.update(project_manifest_fields(project, is_pair))
        manifest_rows.append(row)
        lists["all_3212"].append(sample)
        lists["wgs_rna_pair_449" if is_pair else "wgs_only_2763"].append(sample)
        project_lists[project].append(sample)
        # PacBio multi-file의 설명값은 실제 file URL loop에서 이미 집계했다.
        if is_remote_url(wgs_url):
            unique_resource_urls[wgs_url] = wgs_bytes
        if rna_url and is_remote_url(rna_url):
            unique_resource_urls[rna_url] = rna_bytes

    out = args.out_dir.resolve()
    write_tsv(out / "manifests/preprocess_manifest.tsv", BASE_COLUMNS + EXTRA_COLUMNS, manifest_rows)
    write_tsv(out / "manifests/read_files.tsv", FILE_COLUMNS, read_file_rows)
    lists["all_pair_first_3212"] = lists["wgs_rna_pair_449"] + lists["wgs_only_2763"]
    for name, samples in lists.items():
        path = out / "sample_lists" / f"{name}.txt"
        write_text(path, "\n".join(samples) + "\n")
    for project, samples in sorted(project_lists.items()):
        path = out / "sample_lists/project" / f"{safe(project)}.txt"
        write_text(path, "\n".join(samples) + "\n")
    pilot = [sample for sample in ("HG00281", "HG01204") if sample in set(lists["all_3212"])]
    write_text(out / "sample_lists/pilot_2.txt", "\n".join(pilot) + "\n")

    wgs_primary_total = sum(int(row["wgs_input_bytes"]) for row in manifest_rows)
    rna_primary_total = sum(int(row["rna_input_bytes"]) for row in manifest_rows)
    primary_total = wgs_primary_total + rna_primary_total
    summary = {
        "schema": "mfm_portable_manifest_v3",
        "samples": len(manifest_rows),
        "wgs_only": len(lists["wgs_only_2763"]),
        "wgs_rna_pair": len(lists["wgs_rna_pair_449"]),
        "wgs_primary_bytes": wgs_primary_total,
        "rna_primary_bytes": rna_primary_total,
        "primary_input_bytes": primary_total,
        "unique_known_resource_bytes": sum(unique_resource_urls.values()),
        "resource_count_deduplicated_by_url": len(unique_resource_urls),
        "resource_total_semantics": "known remote resources deduplicated by exact URL; unknown indexes and unfrozen panel files excluded",
        "unknown_wgs_index_samples": sum(row["wgs_index_bytes"] == "" for row in manifest_rows),
        "unknown_rna_index_samples": sum(row["data_group"] == "wgs_rna_pair" and row["rna_index_bytes"] == "" for row in manifest_rows),
        "wgs_library_amplification_counts": dict(sorted(Counter(
            str(row["wgs_library_amplification"]) for row in manifest_rows
        ).items())),
        "rna_library_amplification_counts": dict(sorted(Counter(
            str(row["rna_library_amplification"]) for row in manifest_rows
        ).items())),
        "manifest": "manifests/preprocess_manifest.tsv",
        "read_files_manifest": "manifests/read_files.tsv",
    }
    write_text(out / "manifests/portable_summary.json", json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
