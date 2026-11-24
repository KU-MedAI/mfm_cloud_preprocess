"""다중 형식 WGS 입력 manifest 계약 검증기.

이 모듈은 Stage 04 정렬을 수행하지 않는다. 원본 입력, 분석용 정렬 산출물,
VCF와 세 reference 역할이 모순 없이 기술됐는지만 검사한다.
"""
from __future__ import annotations

import csv
import json
import os
import tempfile
from pathlib import Path
from typing import Any

from .library_provenance import LibraryProvenanceError, validate_library_provenance


SAMPLE_REQUIRED_COLUMNS = (
    "sample_id", "platform", "read_type", "dataset_role", "project_id", "cohort", "alignment_path",
)
FILE_COLUMNS = (
    "sample_id", "run_id", "file_role", "format", "alignment_state", "path", "source_url",
    "size_bytes", "md5", "source_reference_fasta", "read_group_id",
)
PLATFORM_READ_TYPE = {
    "illumina": "short_read", "ont": "long_read", "hifi": "long_read", "clr": "long_read",
}
FORMAT_STATES = {
    "cram": "aligned", "bam": "aligned", "ubam": "unaligned",
    "fastq": "unaligned", "fastq_gz": "unaligned",
}


class InputManifestError(ValueError):
    """입력 manifest 계약 위반."""


def _present(value: str | None) -> bool:
    return bool(value and value.strip() and value.strip() != "-")


def _normalise_platform(value: str) -> str:
    key = value.lower().replace("-", "_").replace(" ", "")
    aliases = {
        "nanopore": "ont", "oxfordnanopore": "ont", "oxford_nanopore": "ont",
        "pacbio_hifi": "hifi", "pacbiohifi": "hifi", "pacbio_ccs": "hifi",
        "pacbioccs": "hifi", "ccs": "hifi",
        "pacbio_clr": "clr", "pacbioclr": "clr", "pb_clr": "clr",
        "ilmn": "illumina",
    }
    return aliases.get(key, key)


def _resolve(value: str | None, base: Path, *, use_input_root: bool = True) -> str | None:
    if not _present(value):
        return None
    path = Path(str(value))
    if path.is_absolute():
        return str(path)
    input_root = os.environ.get("INPUT_ROOT") if use_input_root else None
    return str((Path(input_root) / path).resolve() if input_root else (base / path).resolve())


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise InputManifestError(f"manifest가 없습니다: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise InputManifestError(f"TSV header가 없습니다: {path}")
        rows = list(reader)
    return list(reader.fieldnames), rows


def _require_file(path: str | None, label: str, allow_missing: bool) -> None:
    if not path:
        raise InputManifestError(f"{label} 경로가 비어 있습니다")
    if not allow_missing and not Path(path).is_file():
        raise InputManifestError(f"{label} 파일이 없습니다: {path}")


def _sample_row(manifest: Path, sample: str) -> tuple[dict[str, str], Path]:
    columns, rows = _read_tsv(manifest)
    if tuple(columns[:7]) != SAMPLE_REQUIRED_COLUMNS:
        raise InputManifestError(
            "sample manifest 첫 7열은 다음 순서여야 합니다: " + ",".join(SAMPLE_REQUIRED_COLUMNS)
        )
    matched = [row for row in rows if row.get("sample_id") == sample]
    if len(matched) != 1:
        raise InputManifestError(f"sample 행은 정확히 하나여야 합니다: sample={sample} count={len(matched)}")
    row = matched[0]
    platform = _normalise_platform(row.get("platform", ""))
    expected = PLATFORM_READ_TYPE.get(platform)
    if expected is None or row.get("read_type") != expected:
        raise InputManifestError(
            f"platform/read_type 불일치: platform={platform} read_type={row.get('read_type')}"
        )
    if row.get("dataset_role") not in {"pre_training", "post_training"}:
        raise InputManifestError(f"dataset_role 오류: {row.get('dataset_role')}")
    return row, manifest.parent.resolve()


def _validate_file_rows(
    file_manifest: Path, sample: str, platform: str, allow_missing: bool,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    columns, rows = _read_tsv(file_manifest)
    if tuple(columns) != FILE_COLUMNS:
        raise InputManifestError("read file manifest 열이 계약과 다릅니다: " + ",".join(FILE_COLUMNS))
    rows = [row for row in rows if row.get("sample_id") == sample]
    if not rows:
        raise InputManifestError(f"read file manifest에 sample 행이 없습니다: {sample}")

    states: set[str] = set()
    formats: set[str] = set()
    seen_paths: set[str] = set()
    run_roles: dict[str, list[str]] = {}
    resolved_rows: list[dict[str, str]] = []
    for row in rows:
        fmt = row.get("format", "").lower()
        state = row.get("alignment_state", "").lower()
        role = row.get("file_role", "").lower()
        expected_state = FORMAT_STATES.get(fmt)
        if expected_state is None:
            raise InputManifestError(f"지원하지 않는 format: {fmt}")
        if state != expected_state:
            raise InputManifestError(f"format/alignment_state 불일치: format={fmt} state={state}")
        allowed_roles = {"alignment", "long_reads"} if platform in {"ont", "hifi", "clr"} else {
            "alignment", "read1", "read2", "singletons"
        }
        if role not in allowed_roles:
            raise InputManifestError(f"platform/file_role 불일치: platform={platform} role={role}")
        if fmt in {"cram", "bam", "ubam"} and role != "alignment":
            raise InputManifestError(f"{fmt} file_role은 alignment여야 합니다")
        if fmt in {"fastq", "fastq_gz"} and role == "alignment":
            raise InputManifestError(f"{fmt} file_role은 read role이어야 합니다")
        if fmt == "cram" and not _present(row.get("source_reference_fasta")):
            raise InputManifestError("CRAM에는 source_reference_fasta가 필요합니다")
        resolved_path = _resolve(row.get("path"), file_manifest.parent.resolve())
        _require_file(resolved_path, "read input", allow_missing)
        if resolved_path in seen_paths:
            raise InputManifestError(f"read input path 중복: {resolved_path}")
        seen_paths.add(str(resolved_path))
        states.add(state)
        formats.add(fmt)
        run_id = row.get("run_id", "").strip()
        if not run_id:
            raise InputManifestError("run_id가 비어 있습니다")
        run_roles.setdefault(run_id, []).append(role)
        copied = dict(row)
        copied["path"] = str(resolved_path)
        copied["source_reference_fasta"] = _resolve(
            row.get("source_reference_fasta"), file_manifest.parent.resolve()
        ) or ""
        if copied["source_reference_fasta"]:
            _require_file(copied["source_reference_fasta"], "source reference", allow_missing)
        resolved_rows.append(copied)

    if len(states) != 1:
        raise InputManifestError(f"한 sample에 aligned/unaligned 입력을 섞을 수 없습니다: {sorted(states)}")
    if platform == "illumina" and states == {"unaligned"}:
        for run_id, roles in run_roles.items():
            n1, n2 = roles.count("read1"), roles.count("read2")
            if (n1, n2) not in {(1, 1), (0, 0)}:
                raise InputManifestError(
                    f"paired FASTQ는 run마다 read1/read2가 하나씩 필요합니다: run={run_id} roles={roles}"
                )
            if n1 == n2 == 0 and roles.count("singletons") != 1 and roles != ["alignment"]:
                raise InputManifestError(
                    f"short-read run은 paired, singleton 또는 uBAM alignment여야 합니다: run={run_id} roles={roles}"
                )
    elif platform in {"hifi", "ont", "clr"} and states == {"unaligned"}:
        for run_id, roles in run_roles.items():
            if len(roles) != 1 or roles[0] not in {"alignment", "long_reads"}:
                raise InputManifestError(
                    "raw long-read run은 uBAM/FASTQ 파일 하나여야 합니다: "
                    f"run={run_id} roles={roles}"
                )
    return resolved_rows, {
        "file_count": len(resolved_rows), "formats": sorted(formats), "alignment_states": sorted(states),
    }


def validate_wgs_input_manifest(
    manifest: Path, sample: str, *, read_files_manifest: Path | None = None,
    allow_missing_files: bool = False,
) -> dict[str, Any]:
    """한 sample의 다중 입력 계약을 검사하고 정규화된 계약을 반환한다."""
    manifest = manifest.resolve()
    row, base = _sample_row(manifest, sample)
    platform = _normalise_platform(row["platform"])
    # sample_id는 최종 산출물/모델용 canonical ID이고, source_sample_id는 BAM @RG SM에
    # 기록된 공개 원본 ID다(예: NA12878 <- HG001). VCF donor ID와는 별도 필드다.
    source_sample_id = (row.get("source_sample_id") or sample).strip()
    if not source_sample_id:
        raise InputManifestError("source_sample_id가 비어 있습니다")
    has_rna = _present(row.get("rna_path")) or (row.get("data_group") or "").strip() == "wgs_rna_pair"
    try:
        library_provenance = validate_library_provenance(row, has_rna=has_rna)
    except LibraryProvenanceError as exc:
        raise InputManifestError(str(exc)) from exc
    legacy_alignment = _resolve(row.get("alignment_path"), base)
    configured_file_manifest = _resolve(row.get("read_files_manifest"), base, use_input_root=False)
    if read_files_manifest is not None:
        explicit = str(read_files_manifest.resolve())
        if configured_file_manifest and explicit != configured_file_manifest:
            raise InputManifestError("CLI와 sample manifest의 read_files_manifest가 서로 다릅니다")
        configured_file_manifest = explicit
    if bool(legacy_alignment) == bool(configured_file_manifest):
        raise InputManifestError("alignment_path와 read_files_manifest 중 정확히 하나만 지정해야 합니다")

    analysis_reference = _resolve(row.get("analysis_reference_fasta"), base)
    source_decode_reference = _resolve(row.get("source_decode_reference_fasta"), base)
    vcf_reference = _resolve(row.get("vcf_reference_fasta"), base)
    prepared = _resolve(row.get("prepared_alignment_path"), base)
    if prepared and Path(prepared).suffix.lower() != ".bam":
        raise InputManifestError(
            "prepared_alignment_path는 Stage 04 canonical BAM이어야 합니다: " + prepared
        )
    file_rows: list[dict[str, str]] = []
    if legacy_alignment:
        _require_file(legacy_alignment, "legacy alignment", allow_missing_files)
        file_summary: dict[str, Any] = {
            "file_count": 1,
            "formats": [Path(legacy_alignment).suffix.lower().lstrip(".") or "unknown"],
            "alignment_states": ["aligned"],
        }
        input_mode = "legacy_alignment"
        if Path(legacy_alignment).suffix.lower() == ".cram" and not source_decode_reference:
            raise InputManifestError(
                "legacy CRAM 입력에는 source_decode_reference_fasta가 필요합니다: "
                + legacy_alignment
            )
    else:
        assert configured_file_manifest is not None
        _require_file(configured_file_manifest, "read_files_manifest", allow_missing_files)
        file_rows, file_summary = _validate_file_rows(
            Path(configured_file_manifest), sample, platform, allow_missing_files
        )
        if not analysis_reference:
            raise InputManifestError("raw/file manifest 입력에는 analysis_reference_fasta가 필요합니다")
        input_mode = "read_files_manifest"

    vcf_source = (row.get("vcf_source") or "").strip().lower()
    if platform in {"hifi", "ont", "clr"} and vcf_source not in {"public", "clair3"}:
        raise InputManifestError("long-read에는 vcf_source=public이 필요합니다(Clair3 재현은 명시적 clair3)")
    if platform == "illumina" and vcf_source not in {"", "illumina_panel"}:
        raise InputManifestError("short-read vcf_source는 illumina_panel이어야 합니다")
    public_vcfs = [
        _resolve(item.strip(), base) for item in (row.get("public_vcf_path") or "").split(";") if _present(item)
    ]
    if vcf_source == "public":
        if not public_vcfs:
            raise InputManifestError("vcf_source=public에는 public_vcf_path가 필요합니다")
        if not vcf_reference:
            raise InputManifestError("vcf_source=public에는 vcf_reference_fasta가 필요합니다")
        if not _present(row.get("public_vcf_sample_id")):
            raise InputManifestError("vcf_source=public에는 public_vcf_sample_id가 필요합니다")
        if not _present(row.get("vcf_callset")):
            raise InputManifestError("vcf_source=public에는 vcf_callset이 필요합니다")
        if (row.get("phase_source") or "").strip().lower() != "provided":
            raise InputManifestError("vcf_source=public에는 phase_source=provided가 필요합니다")
        for value in public_vcfs:
            _require_file(value, "public VCF", allow_missing_files)

    for value, label in (
        (analysis_reference, "analysis reference"),
        (source_decode_reference, "source decode reference"),
        (vcf_reference, "VCF reference"),
    ):
        if value:
            _require_file(value, label, allow_missing_files)

    return {
        "schema": "mfm_wgs_input_contract_v2", "sample_id": sample,
        "source_sample_id": source_sample_id, "platform": platform,
        "read_type": row["read_type"], "input_mode": input_mode,
        "library_provenance": library_provenance,
        "source": {
            "alignment": legacy_alignment, "read_files_manifest": configured_file_manifest,
            "files": file_rows, **file_summary,
        },
        "references": {
            "source_decode": source_decode_reference, "analysis": analysis_reference, "vcf": vcf_reference,
        },
        "prepared_alignment": prepared,
        "vcf": {
            "source": vcf_source or None, "paths": public_vcfs,
            "callset": row.get("vcf_callset") or None, "phase_source": row.get("phase_source") or None,
        },
    }


def write_contract_json(contract: dict[str, Any], output: Path) -> None:
    """부분 JSON이 남지 않도록 같은 디렉터리에서 원자적으로 교체한다."""
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{output.name}.", dir=output.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(contract, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        os.replace(temporary, output)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
